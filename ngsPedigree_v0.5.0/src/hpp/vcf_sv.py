"""
Stdlib SV-VCF reader for Delly2 / Manta DEL calls.

Reads uncompressed or .gz VCF. Filters to PASS DEL records. Returns
one record per called variant with per-sample genotype, plus a few
INFO fields commonly used for caller-merging (END, SVLEN, FILTER,
QUAL).

No pysam / no cyvcf2 dependency on purpose: the goal is to run this
on hatchery laptops without an HPC environment.
"""

from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


SV_TYPE_DEL = "DEL"


@dataclass(frozen=True)
class SvRecord:
    chrom: str
    pos: int                      # 1-based VCF POS (left breakpoint)
    end: int                      # 1-based END (right breakpoint, from INFO)
    sv_type: str                  # "DEL" | "DUP" | "INV" | "INS" | "BND" | ...
    svlen: Optional[int]
    qual: Optional[float]
    filter_pass: bool
    caller: str                   # "delly" | "manta" | other
    genotypes: Dict[str, str] = field(default_factory=dict)  # sample_id → "0/0" | "0/1" | "1/1" | "./."
    raw_id: str = ""

    @property
    def marker_key(self) -> Tuple[str, int, int]:
        """Internal key for hashing — useful for de-duping per-caller."""
        return (self.chrom, self.pos, self.end)


# ----------------------------------------------------------------------
# Opener (uncompressed or gzip transparently).
# ----------------------------------------------------------------------


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")


# ----------------------------------------------------------------------
# Genotype normalisation.
# ----------------------------------------------------------------------


def normalise_gt(raw: Optional[str]) -> str:
    if raw is None:
        return "./."
    g = raw.strip()
    if not g or g in (".", "./.", ".|."):
        return "./."
    # GT can be "0/1", "0|1", "1", etc. Collapse phased separators.
    g = g.split(":")[0]   # take only the GT subfield even if FORMAT mismatch
    g = g.replace("|", "/")
    # Defensive: only accept canonical biallelic genotypes; anything else
    # → missing. This lines up with downstream Mendelian code which
    # treats non-standard GTs as missing.
    if g in ("0/0", "0/1", "1/0", "1/1"):
        return "0/1" if g == "1/0" else g
    return "./."


# ----------------------------------------------------------------------
# INFO parser.
# ----------------------------------------------------------------------


def _parse_info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info_field in (".", ""):
        return out
    for token in info_field.split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = "1"
    return out


# ----------------------------------------------------------------------
# Caller heuristic — pulled from VCF headers (##source) or filename.
# ----------------------------------------------------------------------


def _detect_caller(path: Path, headers: List[str]) -> str:
    name = path.name.lower()
    for h in headers:
        hl = h.lower()
        if "manta" in hl:
            return "manta"
        if "delly" in hl:
            return "delly"
    if "manta" in name:
        return "manta"
    if "delly" in name:
        return "delly"
    return "unknown"


# ----------------------------------------------------------------------
# Iterator: yield one SvRecord per DEL call in the file.
# ----------------------------------------------------------------------


def iter_del_calls(
    vcf_path: str | Path,
    *,
    require_pass: bool = True,
    caller_override: Optional[str] = None,
) -> Iterator[SvRecord]:
    """Stream DEL records out of a Delly or Manta VCF.

    Skips non-DEL SV types and (by default) anything not PASS.
    """
    path = Path(vcf_path)
    samples: List[str] = []
    headers: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if not samples and not line.startswith("#"):
                raise ValueError(f"{path}: data line before #CHROM header")
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            raw_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual_raw = fields[5]
            filt = fields[6]
            info = _parse_info(fields[7])

            sv_type = info.get("SVTYPE", "")
            if sv_type != SV_TYPE_DEL:
                # Some VCFs have ALT="<DEL>" without explicit SVTYPE.
                if not (alt.upper() == "<DEL>"):
                    continue
                sv_type = SV_TYPE_DEL

            pass_ok = (filt == "PASS" or filt == "." or filt == "")
            if require_pass and not pass_ok:
                continue

            try:
                end = int(info.get("END", "0")) or pos + 1
            except ValueError:
                end = pos + 1
            try:
                svlen = abs(int(info.get("SVLEN", "0"))) or (end - pos)
            except ValueError:
                svlen = end - pos
            try:
                qual = float(qual_raw) if qual_raw not in ("", ".") else None
            except ValueError:
                qual = None

            # Parse FORMAT + per-sample columns for GT only.
            if len(fields) < 9:
                continue
            fmt = fields[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            gts: Dict[str, str] = {}
            for sample, col in zip(samples, fields[9:]):
                sub = col.split(":")
                raw_gt = sub[gt_idx] if gt_idx < len(sub) else ""
                gts[sample] = normalise_gt(raw_gt)

            yield SvRecord(
                chrom=chrom,
                pos=pos,
                end=end,
                sv_type=sv_type,
                svlen=svlen,
                qual=qual,
                filter_pass=pass_ok,
                caller=caller_override or _detect_caller(path, headers),
                genotypes=gts,
                raw_id=raw_id,
            )


def read_del_calls(
    vcf_path: str | Path,
    *,
    require_pass: bool = True,
    caller_override: Optional[str] = None,
) -> List[SvRecord]:
    return list(iter_del_calls(
        vcf_path,
        require_pass=require_pass,
        caller_override=caller_override,
    ))


def read_sample_ids(vcf_path: str | Path) -> List[str]:
    path = Path(vcf_path)
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []
