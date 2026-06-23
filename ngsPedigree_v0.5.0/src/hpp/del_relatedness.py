"""
Pairwise relatedness from DEL markers — the broke-grad-student
ngsRelate-replacement.

Given a (marker × sample) DEL-genotype matrix, computes per-pair:
  - theta    (KING-robust kinship)
  - IBS0     (opposite-homozygote rate)
  - n_informative

These are the same two coefficients STEP_PED_01_annotate_relationships
uses to classify edges. So this module's output table can be fed
directly into Stage 1 (or into the stdlib shadow classifier) without
any further conversion.

KING-robust estimator (Manichaikul 2010):

    theta = (N_AaAa - 2 * N_AaBB - 2 * N_AAbb)
            / (N_Aa_x + N_x_Aa)

where:
    N_AaAa     = both samples het
    N_AaBB     = sample a het, sample b hom-alt
    N_AAbb     = sample a hom-ref, sample b hom-alt
                  (opposite-homozygote = IBS0)
    N_Aa_x     = sample a het (b anything called)
    N_x_Aa     = sample b het (a anything called)

For the catfish cohort, the marker set is DELs from Delly + Manta;
the estimator is identical to its SNP form.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class PairCoeffs:
    sample_a: str
    sample_b: str
    n_markers: int
    n_informative: int          # both samples called (non-missing)
    theta: Optional[float]
    IBS0: Optional[float]
    n_aa_aa: int                # both het
    n_opposite_hom: int          # IBS0 markers
    n_a_het: int                # sample a het (b anything called)
    n_b_het: int                # sample b het (a anything called)

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


def _is_called(gt: str) -> bool:
    return gt in ("0/0", "0/1", "1/1")


def pair_coeffs(
    sample_a: str,
    sample_b: str,
    genotype_matrix: Dict[str, Dict[str, str]],
) -> PairCoeffs:
    """Compute KING-robust (theta, IBS0) for one pair."""
    n_markers = len(genotype_matrix)
    n_inf = 0
    n_aa_aa = 0          # both het
    n_opp_hom = 0
    n_a_het = 0
    n_b_het = 0
    for calls in genotype_matrix.values():
        ga = calls.get(sample_a, "./.")
        gb = calls.get(sample_b, "./.")
        if not (_is_called(ga) and _is_called(gb)):
            continue
        n_inf += 1
        ha = (ga == "0/1")
        hb = (gb == "0/1")
        if ha and hb:
            n_aa_aa += 1
        if (ga == "0/0" and gb == "1/1") or (ga == "1/1" and gb == "0/0"):
            n_opp_hom += 1
        if ha:
            n_a_het += 1
        if hb:
            n_b_het += 1
    if n_inf == 0 or (n_a_het + n_b_het) == 0:
        return PairCoeffs(
            sample_a=sample_a, sample_b=sample_b,
            n_markers=n_markers, n_informative=n_inf,
            theta=None, IBS0=None,
            n_aa_aa=n_aa_aa, n_opposite_hom=n_opp_hom,
            n_a_het=n_a_het, n_b_het=n_b_het,
        )
    theta = (n_aa_aa - 2 * n_opp_hom) / (n_a_het + n_b_het)
    ibs0 = n_opp_hom / n_inf if n_inf else 0.0
    return PairCoeffs(
        sample_a=sample_a, sample_b=sample_b,
        n_markers=n_markers, n_informative=n_inf,
        theta=theta, IBS0=ibs0,
        n_aa_aa=n_aa_aa, n_opposite_hom=n_opp_hom,
        n_a_het=n_a_het, n_b_het=n_b_het,
    )


def all_pairs(
    samples: Sequence[str],
    genotype_matrix: Dict[str, Dict[str, str]],
    *,
    min_informative: int = 50,
) -> List[PairCoeffs]:
    out: List[PairCoeffs] = []
    samples = list(samples)
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            pc = pair_coeffs(a, b, genotype_matrix)
            if pc.n_informative >= min_informative:
                out.append(pc)
    return out


# ----------------------------------------------------------------------
# Write Stage-1-compatible pairwise TSV.
# ----------------------------------------------------------------------


def write_pairwise_tsv(
    path,
    pairs: Sequence[PairCoeffs],
    samples: Sequence[str],
) -> None:
    """Emit a pairwise table in the same shape Stage 1 produces.

    Columns: sample_a, sample_b, a, b, nSites, theta, IBS0, KING.

    The classifier downstream uses theta, IBS0 (and optionally KING).
    We set KING = theta as a placeholder — for KING-robust on biallelic
    markers, theta IS the KING-robust kinship.
    """
    from pathlib import Path
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    sample_index = {s: i for i, s in enumerate(samples)}
    with open(p, "w") as fh:
        fh.write("sample_a\tsample_b\ta\tb\tnSites\ttheta\tIBS0\tKING\n")
        for pc in pairs:
            theta_s = "" if pc.theta is None else f"{pc.theta:.6f}"
            ibs0_s = "" if pc.IBS0 is None else f"{pc.IBS0:.6f}"
            king_s = theta_s   # KING-robust ≡ theta on biallelic data
            fh.write(
                f"{pc.sample_a}\t{pc.sample_b}\t"
                f"{sample_index.get(pc.sample_a, '')}\t"
                f"{sample_index.get(pc.sample_b, '')}\t"
                f"{pc.n_informative}\t{theta_s}\t{ibs0_s}\t{king_s}\n"
            )
