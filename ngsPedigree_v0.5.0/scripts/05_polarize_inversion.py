#!/usr/bin/env python3
"""
05_polarize_inversion.py — read a karyotype-calls IN JSON, run inversion
polarization + transmission calling + the Mendelian drive test, and
write the polarized-transmissions OUT JSON for ngsTracts.

Usage:
  python 05_polarize_inversion.py --in PATH --out PATH

The OUT JSON is the hand-off contract to ngsTracts; ngsTracts reads
transmitted_arrangement per (parent, offspring, inversion) and uses it
as the phase polarity for marker-level CO / NCO / DCO scanning over the
inversion interval. ngsPedigree stops at the end of polarization;
recombination tract calling is ngsTracts' job.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.inversion_polarization import (  # noqa: E402
    call_transmissions,
    drive_test,
    polarize,
)
from hpp.mtdna_check import (  # noqa: E402
    HAMMING_THRESHOLD_DEFAULT,
    build_not_supplied_block,
    build_validation_block,
    check_pedigree,
    filter_by_mtdna,
    load_mtdna_haplotypes,
)
from hpp.ngstracts_io import (  # noqa: E402
    load_karyotype_calls,
    write_polarized_transmissions,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="input_path", required=True,
                    help="IN JSON (ngspedigree_karyotype_calls_in_v1)")
    ap.add_argument("--out", dest="output_path", required=True,
                    help="OUT JSON for ngsTracts (ngspedigree_polarized_transmissions_v1)")
    ap.add_argument("--mtdna", dest="mtdna_path", default=None,
                    help=("optional mtDNA haplotype TSV "
                          "(sample_id, mtdna_haplotype, [mtdna_sequence, "
                          "mtdna_n_sites]). Maternal-line pedigree pre-flight; "
                          "incompatible (mother, offspring) pairs are excluded "
                          "from polarization."))
    ap.add_argument("--mtdna-hamming-threshold", type=int,
                    default=HAMMING_THRESHOLD_DEFAULT,
                    help=(f"Hamming-distance tolerance when both samples carry "
                          f"mtDNA sequences (default {HAMMING_THRESHOLD_DEFAULT})."))
    args = ap.parse_args()

    inversion_id, polarity_hint, calls, dyads, triads = load_karyotype_calls(args.input_path)
    print(f"[polarize] inversion={inversion_id}  hint={polarity_hint}  "
          f"n_calls={len(calls)}  n_dyads={len(dyads)}  n_triads={len(triads)}")

    # --- mtDNA pre-flight (layer 2) -------------------------------------
    if args.mtdna_path:
        mt_records = load_mtdna_haplotypes(args.mtdna_path)
        n_triads_before = len(triads)
        n_maternal_dyads_before = sum(1 for d in dyads if d.parent_sex == "female")
        mt_checks = check_pedigree(
            mt_records, dyads, triads,
            hamming_threshold=args.mtdna_hamming_threshold,
        )
        dyads, triads = filter_by_mtdna(dyads, triads, mt_checks)
        n_triads_excluded = n_triads_before - len(triads)
        n_maternal_dyads_after = sum(1 for d in dyads if d.parent_sex == "female")
        n_dyads_excluded = n_maternal_dyads_before - n_maternal_dyads_after
        mtdna_validation = build_validation_block(
            mt_checks,
            n_triads_excluded=n_triads_excluded,
            n_dyads_excluded=n_dyads_excluded,
            hamming_threshold=args.mtdna_hamming_threshold,
        )
        n_inc = mtdna_validation["n_incompatible"]
        print(f"[polarize] mtDNA: checked={len(mt_checks)}  "
              f"compatible={mtdna_validation['n_compatible']}  "
              f"incompatible={n_inc}  "
              f"ambiguous={mtdna_validation['n_ambiguous']}  "
              f"triads_excluded={n_triads_excluded}  "
              f"maternal_dyads_excluded={n_dyads_excluded}")
    else:
        mtdna_validation = build_not_supplied_block()

    res = polarize(
        inversion_id=inversion_id,
        karyotype_calls=calls,
        dyads=dyads, triads=triads,
        polarity_hint=polarity_hint,
    )
    print(f"[polarize] chosen={res.chosen_polarity}  "
          f"dyad_incompatible={res.dyad_compat[res.chosen_polarity].n_incompatible}  "
          f"triad_incompatible={res.triad_compat[res.chosen_polarity].n_incompatible}  "
          f"symmetric={res.polarities_symmetric}")

    transmissions = call_transmissions(
        inversion_id=inversion_id,
        karyotype_calls=calls,
        dyads=dyads, triads=triads,
        chosen_polarity=res.chosen_polarity,
    )
    stats = drive_test(transmissions)
    print(f"[polarize] transmissions={len(transmissions)}  "
          f"informative={stats.n_informative_transmissions}  "
          f"INV_rate={stats.INV_transmission_rate}  "
          f"binomial_p={stats.binomial_pvalue}")

    out = write_polarized_transmissions(
        args.output_path,
        result=res,
        transmissions=transmissions,
        drive_stats=stats,
        mtdna_validation=mtdna_validation,
    )
    print(f"[polarize] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
