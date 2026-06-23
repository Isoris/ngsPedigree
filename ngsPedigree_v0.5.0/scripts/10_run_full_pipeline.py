#!/usr/bin/env python3
"""
10_run_full_pipeline.py — one-command end-to-end runner.

Chains the broke-grad-student stack:

  Delly + Manta VCFs                                                  [bloc 12]
       ↓ catalogue merge
  merged DEL catalogue
       ↓ KING-robust relatedness on DEL markers
  pairwise table + candidate PO pairs + triads (forced-offspring rule)
       ↓ chromosome inheritance map per (parent, offspring)            [bloc 12]
  inheritance_segments.tsv
       ↓ (optional) de novo LRR discovery from DEL correlations        [NEW bloc 16]
  candidate LRR list  ← --list-of-lrr overrides this when supplied
       ↓ family-based OR enrichment against matched background          [bloc 13]
  lrr_enrichment.tsv
       ↓ (when a karyotype catalogue is supplied)
       ↓ per-DEL × per-LRR arrangement-linkage classifier               [NEW bloc 17]
  del_arrangement_linkage.tsv   ← surfaces the "situation 1" cases:
                                  hemizygous-only DEL where the host
                                  homokaryotype is absent or depleted

Usage:
  python 10_run_full_pipeline.py \
      --delly delly.vcf[.gz] \
      --manta manta.vcf[.gz] \
      [--list-of-lrr lrr_list.tsv | --discover-lrrs] \
      [--karyotype-catalogue catalogue.json] \
      --outdir results/

Notes:
  - When neither --list-of-lrr nor --discover-lrrs is supplied, the LRR-side
    of the pipeline (blocs 13, 16, 17) is skipped — the family graph and
    inheritance map are still emitted.
  - When --karyotype-catalogue is omitted, bloc 17 (arrangement-linkage)
    is skipped because per-sample HOM_REF/HET/HOM_INV labels are needed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.catalogue_merge import (  # noqa: E402
    cohort_samples, merge_two_callers, to_genotype_matrix,
)
from hpp.del_arrangement_linkage import classify_all as classify_linkage_all  # noqa: E402
from hpp.del_arrangement_linkage import write_linkage_tsv  # noqa: E402
from hpp.del_inheritance import (  # noqa: E402
    DelMarkerLocus,
    build_inheritance_map_for_dyad,
    build_inheritance_map_for_triad,
    exclude_as_parent,
)
from hpp.del_relatedness import all_pairs  # noqa: E402
from hpp.karyotype_catalogue import load_catalogue  # noqa: E402
from hpp.lrr_discovery import (  # noqa: E402
    discover_candidate_lrrs, write_candidate_lrr_tsv,
)
from hpp.lrr_enrichment import (  # noqa: E402
    LRRInterval, Relationship,
    compute_enrichment, load_lrr_list, write_enrichment_tsv,
)
from hpp.relatedness_sim import classify_edge_stdlib  # noqa: E402
from hpp.vcf_sv import read_del_calls  # noqa: E402


def _candidates_to_intervals(cands):
    return [LRRInterval(c.lrr_id, c.chrom, c.start, c.end) for c in cands]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delly", required=True)
    ap.add_argument("--manta", required=True)
    ap.add_argument("--outdir", required=True)

    # LRR side
    lrr_group = ap.add_mutually_exclusive_group()
    lrr_group.add_argument("--list-of-lrr", "--list_of_LRR", dest="lrr_list",
                            help="curated LRR TSV (lrr_id, chrom, start, end)")
    lrr_group.add_argument("--discover-lrrs", action="store_true",
                            help=("de novo LRR discovery from DEL correlation "
                                  "patterns; emits candidate LRR list and "
                                  "feeds it into the enrichment step"))

    # Discovery knobs
    ap.add_argument("--window-size", type=int, default=1_000_000)
    ap.add_argument("--min-markers-per-window", type=int, default=4)
    ap.add_argument("--correlation-threshold", type=float, default=0.50)

    # Arrangement-linkage side
    ap.add_argument("--karyotype-catalogue", default=None,
                    help=("per-(sample, LRR) karyotype catalogue JSON "
                          "(ngspedigree_karyotype_catalogue_v1). When "
                          "supplied, runs the per-DEL × per-LRR "
                          "arrangement-linkage classifier (bloc 17)."))

    # Common knobs
    ap.add_argument("--bp-tolerance", type=int, default=500)
    ap.add_argument("--reciprocal-overlap", type=float, default=0.5)
    ap.add_argument("--min-informative-pair", type=int, default=50)
    ap.add_argument("--min-excluding", type=int, default=1)
    ap.add_argument("--n-background", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- (1) read & merge SV VCFs ---
    print(f"[full-pipeline] reading {args.delly}")
    delly = read_del_calls(args.delly, caller_override="delly")
    print(f"[full-pipeline] reading {args.manta}")
    manta = read_del_calls(args.manta, caller_override="manta")
    merged = merge_two_callers(
        delly, manta,
        bp_tolerance=args.bp_tolerance,
        reciprocal_overlap=args.reciprocal_overlap,
    )
    print(f"[full-pipeline] merged catalogue: {len(merged)} markers, "
          f"{sum(1 for m in merged if m.n_callers == 2)} two-caller, "
          f"{sum(1 for m in merged if m.n_callers == 1)} single-caller")
    with open(outdir / "merged_del_catalogue.json", "w") as fh:
        json.dump(
            {"schema": "ngspedigree_merged_del_catalogue_v1",
             "markers": [m.to_dict() for m in merged]},
            fh, indent=2,
        )

    gmatrix = to_genotype_matrix(merged)
    samples = cohort_samples(merged)
    loci = [DelMarkerLocus(m.marker_id, m.chrom, m.start, m.end) for m in merged]
    print(f"[full-pipeline] cohort: {len(samples)} samples")

    # --- (2) KING-robust pairwise + edge classification ---
    pairs = all_pairs(samples, gmatrix, min_informative=args.min_informative_pair)
    sample_index = {s: i for i, s in enumerate(samples)}
    with open(outdir / "pairwise_relationship_classification.tsv", "w") as fh:
        fh.write("sample_a\tsample_b\ta\tb\tnSites\ttheta\tIBS0\tKING\t"
                 "edge_class\tconfidence\treasons\n")
        po_pairs = []
        for p in pairs:
            v = classify_edge_stdlib({"theta": p.theta, "IBS0": p.IBS0,
                                       "KING": p.theta, "nSites": p.n_informative})
            theta_s = "" if p.theta is None else f"{p.theta:.6f}"
            ibs0_s = "" if p.IBS0 is None else f"{p.IBS0:.6f}"
            fh.write(f"{p.sample_a}\t{p.sample_b}\t"
                     f"{sample_index[p.sample_a]}\t{sample_index[p.sample_b]}\t"
                     f"{p.n_informative}\t{theta_s}\t{ibs0_s}\t{theta_s}\t"
                     f"{v['edge_class']}\t{v['confidence']}\t{v['reasons']}\n")
            if v["edge_class"] == "parent_offspring":
                po_pairs.append((p.sample_a, p.sample_b))

    # --- (3) Mendelian exclusion + triad detection ---
    edge_class_by_pair = {}
    for p in pairs:
        v = classify_edge_stdlib({"theta": p.theta, "IBS0": p.IBS0,
                                    "KING": p.theta, "nSites": p.n_informative})
        edge_class_by_pair[(p.sample_a, p.sample_b)] = v["edge_class"]
        edge_class_by_pair[(p.sample_b, p.sample_a)] = v["edge_class"]
    po_partners = {}
    for a, b in po_pairs:
        po_partners.setdefault(a, set()).add(b)
        po_partners.setdefault(b, set()).add(a)
    triads_detected = []
    for offspring, partners in po_partners.items():
        plist = sorted(partners)
        for i, p1 in enumerate(plist):
            for p2 in plist[i + 1:]:
                if edge_class_by_pair.get((p1, p2)) in {"unrelated",
                                                          "third_degree"}:
                    triads_detected.append((p1, p2, offspring))

    print(f"[full-pipeline] candidate PO edges: {len(po_pairs)}  "
          f"detected triads: {len(triads_detected)}")
    with open(outdir / "candidate_PO_pairs.tsv", "w") as fh:
        fh.write("sample_a\tsample_b\tedge_class\n")
        for a, b in po_pairs:
            fh.write(f"{a}\t{b}\tparent_offspring\n")
    with open(outdir / "detected_triads.tsv", "w") as fh:
        fh.write("triad_id\tpaternal_sample_id\tmaternal_sample_id\toffspring_sample_id\n")
        for i, (p1, p2, off) in enumerate(triads_detected, start=1):
            fh.write(f"triad_{i:03d}\t{p1}\t{p2}\t{off}\n")

    # --- (4) Chromosome inheritance map ---
    inh_rows = []
    for p1, p2, off in triads_detected:
        m = build_inheritance_map_for_triad(
            triad_id=f"triad_{p1}_{p2}_{off}",
            paternal_sample_id=p1, maternal_sample_id=p2,
            offspring_sample_id=off,
            loci=loci, genotype_matrix=gmatrix,
        )
        for side, segs in m.items():
            for s in segs:
                inh_rows.append((s, f"triad_{side}"))
    with open(outdir / "inheritance_segments.tsv", "w") as fh:
        fh.write("relationship_id\trelationship_type\tparent_sample_id\t"
                 "offspring_sample_id\tchrom\tseg_start\tseg_end\t"
                 "transmitted_allele\tn_informative_markers\tconfidence\t"
                 "recomb_event_left\trecomb_event_right\t"
                 "parental_hap_inherited\n")
        for s, rt in inh_rows:
            fh.write(f"{s.dyad_or_triad_id}\t{rt}\t{s.parent_sample_id}\t"
                     f"{s.offspring_sample_id}\t{s.chrom}\t{s.seg_start}\t"
                     f"{s.seg_end}\t{s.transmitted_allele}\t"
                     f"{s.n_informative_markers}\t{s.confidence}\t"
                     f"{str(s.recomb_event_left).lower()}\t"
                     f"{str(s.recomb_event_right).lower()}\t"
                     f"{s.parental_hap_inherited}\n")

    # --- (5) LRR set: curated list or de novo discovery ---
    lrrs = []
    if args.lrr_list:
        lrrs = load_lrr_list(args.lrr_list)
        print(f"[full-pipeline] LRRs from --list-of-lrr: {len(lrrs)}")
    elif args.discover_lrrs:
        cands = discover_candidate_lrrs(
            loci=loci, genotype_matrix=gmatrix, samples=samples,
            window_size=args.window_size,
            min_markers_per_window=args.min_markers_per_window,
            correlation_threshold=args.correlation_threshold,
        )
        write_candidate_lrr_tsv(outdir / "discovered_lrrs.tsv", cands)
        lrrs = _candidates_to_intervals(cands)
        print(f"[full-pipeline] de novo discovered LRRs: {len(lrrs)}")
    else:
        print("[full-pipeline] no --list-of-lrr or --discover-lrrs; "
              "skipping LRR-side analyses")

    # --- (6) LRR odds-ratio enrichment ---
    if lrrs and triads_detected:
        rels = [Relationship(
            relationship_id=f"triad_{p1}_{p2}_{off}",
            relationship_type="triad",
            paternal_sample_id=p1, maternal_sample_id=p2,
            parent_sample_id=None, offspring_sample_id=off,
        ) for p1, p2, off in triads_detected]
        enrichments, _ = compute_enrichment(
            lrrs=lrrs, relationships=rels,
            loci=loci, genotype_matrix=gmatrix,
            n_background_per_lrr=args.n_background,
            seed=args.seed,
        )
        write_enrichment_tsv(outdir / "lrr_enrichment.tsv", enrichments)
        print(f"[full-pipeline] wrote lrr_enrichment.tsv ({len(enrichments)} LRRs)")

    # --- (7) Per-DEL × per-LRR arrangement linkage (if karyotype catalogue) ---
    if lrrs and args.karyotype_catalogue:
        cat = load_catalogue(args.karyotype_catalogue)
        sample_arr_by_lrr = {}
        for lrr in lrrs:
            arr = {}
            for cr in cat.filter_to_inversion(lrr.lrr_id):
                band_to_arr = {0: "HOM_REF", 1: "HET", 2: "HOM_INV"}
                arr[cr.sample_id] = band_to_arr[cr.band]
            sample_arr_by_lrr[lrr.lrr_id] = arr
        records = classify_linkage_all(
            lrrs=lrrs, loci=loci, genotype_matrix=gmatrix,
            sample_arrangement_by_lrr=sample_arr_by_lrr,
        )
        write_linkage_tsv(outdir / "del_arrangement_linkage.tsv", records)
        print(f"[full-pipeline] wrote del_arrangement_linkage.tsv "
              f"({len(records)} (DEL, LRR) records)")

    print(f"[full-pipeline] done; outputs in {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
