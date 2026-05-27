"""
Minimal stdlib VCF reader for HPP MVP 1 synthetic fixtures.

This is *not* a general VCF parser. It reads uncompressed, biallelic,
sample-genotype-only VCF files of the form produced by the test
fixtures. The real implementation in MVP 2+ will use cyvcf2 / pysam
against the joint VCF.

Scope:
  - parses #CHROM header to learn sample order;
  - yields VariantRecord(chrom, pos, ref, alt, variant_id, genotypes)
    where genotypes is {sample_id -> "0/0" | "0/1" | "1/1" | "./."};
  - normalises any phased separator ("|") to "/" for genotype strings,
    since HPP MVP 1 treats parental phase as carried by the sidecar
    parent_phase table, not the VCF.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    genotypes: Dict[str, str]

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


def _normalise_gt(raw: str) -> str:
    if raw in (".", "./.", ".|."):
        return "./."
    return raw.replace("|", "/")


def read_vcf(path: str | Path) -> Iterator[VariantRecord]:
    samples: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if not samples:
                raise ValueError(f"VCF {path}: data line before #CHROM header")
            fields = line.split("\t")
            chrom, pos, _id, ref, alt = fields[:5]
            fmt = fields[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError as exc:
                raise ValueError(
                    f"VCF {path}: FORMAT field has no GT at {chrom}:{pos}"
                ) from exc
            if "," in alt:
                raise ValueError(
                    f"VCF {path}: multiallelic site at {chrom}:{pos} — split upstream"
                )
            genotypes: Dict[str, str] = {}
            for sample, sample_field in zip(samples, fields[9:]):
                raw_gt = sample_field.split(":")[gt_idx]
                genotypes[sample] = _normalise_gt(raw_gt)
            yield VariantRecord(
                chrom=chrom,
                pos=int(pos),
                ref=ref,
                alt=alt,
                genotypes=genotypes,
            )
