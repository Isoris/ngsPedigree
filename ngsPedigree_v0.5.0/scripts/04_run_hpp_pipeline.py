#!/usr/bin/env python3
"""
04_run_hpp_pipeline.py — end-to-end HPP smoke run on the synthetic fixtures.

Wires the four placeholder adapters + MVP 1-4 into one command and emits
Tables A, B, and (C or D) to TSV.

Usage:
  # dyad
  python 04_run_hpp_pipeline.py \\
      --relationship-type dyad \\
      --inheritance-map  fixtures/synthetic_dyad/inheritance_map_dyad.tsv \\
      --parent-phase     fixtures/synthetic_dyad/parent_phase.tsv \\
      --joint-vcf        fixtures/synthetic_dyad/joint.vcf \\
      --variant-master   fixtures/synthetic_dyad/variant_master.tsv \\
      --damaging-tier    T1 \\
      --outdir           /tmp/hpp_out_dyad

  # triad
  python 04_run_hpp_pipeline.py --relationship-type triad ...

The script is the placeholder pipeline — a real-cohort run swaps the
adapter sources but the script body does not change.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp import io  # noqa: E402
from hpp.gene_status import TABLE_B_COLUMNS, classify_gene_status  # noqa: E402
from hpp.parental_haps import build_for_parents, build_parental_hap_variants  # noqa: E402
from hpp.project import (  # noqa: E402
    TABLE_A_COLUMNS,
    project_dyad_to_offspring,
    project_triad_to_offspring,
)
from hpp.stage3_placeholder import (  # noqa: E402
    load_dyad_map,
    load_parent_phase,
    load_triad_map,
)
from hpp.summary import (  # noqa: E402
    TABLE_C_COLUMNS,
    TABLE_D_COLUMNS,
    summarise_dyad,
    summarise_triad,
)
from hpp.variant_master import PlaceholderVariantMaster, VariantAnnotation  # noqa: E402
from hpp.vcf_lite import read_vcf  # noqa: E402


def load_variant_master_tsv(path: Path) -> PlaceholderVariantMaster:
    """Tiny TSV loader for the placeholder variant_master.

    Real adapter (MVP 2b) will be a pandas/polars-backed reader against
    MODULE_CONSERVATION STEP 16 output.
    """
    vm = PlaceholderVariantMaster()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
            vm.add(VariantAnnotation(
                variant_id=row["variant_id"],
                gene_id=row.get("gene_id") or None,
                transcript_id=row.get("transcript_id") or None,
                consequence=row.get("consequence") or None,
                impact=row.get("impact") or None,
                sift_class=row.get("sift_class") or None,
                vesm_llr=float(row["vesm_llr"]) if row.get("vesm_llr") else None,
                splice_subclass=row.get("splice_subclass") or None,
            ))
    return vm


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--relationship-type", choices=["dyad", "triad"], required=True)
    ap.add_argument("--inheritance-map", required=True)
    ap.add_argument("--parent-phase", required=True)
    ap.add_argument("--joint-vcf", required=True)
    ap.add_argument("--variant-master", required=True)
    ap.add_argument("--damaging-tier", choices=["T1", "T2", "T3"], default="T1")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[hpp] tier={args.damaging_tier}  relationship={args.relationship_type}")

    # --- Load Stage 3 placeholder + joint VCF + variant_master ---
    phase = load_parent_phase(args.parent_phase)
    variants = list(read_vcf(args.joint_vcf))
    vm = load_variant_master_tsv(Path(args.variant_master))

    if args.relationship_type == "dyad":
        segments = load_dyad_map(args.inheritance_map)
        parent_id = segments[0].parent_sample_id
        parent_haps = build_parental_hap_variants(parent_id, variants, phase)
        print(f"[hpp] loaded {len(segments)} dyad segments, "
              f"{len(variants)} variants, "
              f"{len(vm._records)} variant_master annotations")

        # Table A — projection
        table_a = project_dyad_to_offspring(
            segments=segments, parent_haps=parent_haps,
            all_variants=variants, variant_master=vm,
            damaging_tier=args.damaging_tier,
        )
        # Table B — gene status
        table_b = classify_gene_status(
            table_a_rows=table_a, segments=segments,
            variant_master=vm, damaging_tier=args.damaging_tier,
        )
        # Table C — dyad summary
        table_c = summarise_dyad(
            segments=segments, table_a_rows=table_a, table_b_rows=table_b,
            all_variants=variants, variant_master=vm,
            damaging_tier=args.damaging_tier,
        )

        io.write_tsv(outdir / "table_A_offspring_haplotype_variants.tsv",
                     TABLE_A_COLUMNS, (r.to_dict() for r in table_a))
        io.write_tsv(outdir / "table_B_offspring_gene_status.tsv",
                     TABLE_B_COLUMNS, (r.to_dict() for r in table_b))
        io.write_tsv(outdir / "table_C_dyad_transmission_summary.tsv",
                     TABLE_C_COLUMNS, [table_c.to_dict()])

        print(f"[hpp] table A rows: {len(table_a)}")
        print(f"[hpp] table B rows: {len(table_b)}")
        print(f"[hpp] table C status: {table_c.mendelian_consistency_status} "
              f"({table_c.mendelian_inconsistent_sites} inconsistent sites)")

    else:  # triad
        segments = load_triad_map(args.inheritance_map)
        pat_id = segments[0].paternal_sample_id
        mat_id = segments[0].maternal_sample_id
        per_parent = build_for_parents([pat_id, mat_id], variants, phase)
        print(f"[hpp] loaded {len(segments)} triad segments, "
              f"{len(variants)} variants, "
              f"{len(vm._records)} variant_master annotations")

        table_a = project_triad_to_offspring(
            segments=segments,
            paternal_haps=per_parent[pat_id],
            maternal_haps=per_parent[mat_id],
            all_variants=variants, variant_master=vm,
            damaging_tier=args.damaging_tier,
        )
        table_b = classify_gene_status(
            table_a_rows=table_a, segments=segments,
            variant_master=vm, damaging_tier=args.damaging_tier,
        )
        table_d = summarise_triad(
            segments=segments, table_a_rows=table_a, table_b_rows=table_b,
            all_variants=variants, variant_master=vm,
            damaging_tier=args.damaging_tier,
        )

        io.write_tsv(outdir / "table_A_offspring_haplotype_variants.tsv",
                     TABLE_A_COLUMNS, (r.to_dict() for r in table_a))
        io.write_tsv(outdir / "table_B_offspring_gene_status.tsv",
                     TABLE_B_COLUMNS, (r.to_dict() for r in table_b))
        io.write_tsv(outdir / "table_D_triad_transmission_summary.tsv",
                     TABLE_D_COLUMNS, [table_d.to_dict()])

        print(f"[hpp] table A rows: {len(table_a)}")
        print(f"[hpp] table B rows: {len(table_b)}")
        print(f"[hpp] table D status: {table_d.mendelian_consistency_status} "
              f"({table_d.mendelian_inconsistent_sites} inconsistent / "
              f"{table_d.n_de_novo_candidates} de novo candidates)")

    print(f"[hpp] outputs in {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
