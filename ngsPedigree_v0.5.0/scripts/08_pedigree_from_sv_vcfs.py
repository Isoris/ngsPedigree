#!/usr/bin/env python3
"""
08_pedigree_from_sv_vcfs.py — broke-grad-student pipeline.

Read Delly + Manta VCF catalogues (no ngsRelate, no BEAGLE, no HPC),
merge into a unified DEL marker catalogue, compute KING-robust theta
and IBS0 on the DEL genotypes, identify candidate first-degree edges
and parent–offspring pairs by Mendelian exclusion, and emit a
chromosome inheritance map per confirmed parent-offspring relationship.

Usage:
  python 08_pedigree_from_sv_vcfs.py \
      --delly  delly.vcf[.gz] \
      --manta  manta.vcf[.gz] \
      --outdir results/

Outputs in --outdir:
  pairwise_relationship_classification.tsv  (Stage 1 schema; feed into
                                              STEP_PED_01 or downstream
                                              06_build_polarization_input)
  pairwise_coefficients.tsv  (raw theta / IBS0 / counts per pair)
  candidate_PO_pairs.tsv     (pairs whose Mendelian exclusion test passes)
  inheritance_segments.tsv   (per-chrom transmission trace for each
                              candidate PO dyad)
  merged_del_catalogue.json  (the merged Delly + Manta catalogue)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.catalogue_merge import (  # noqa: E402
    cohort_samples,
    merge_two_callers,
    to_genotype_matrix,
)
from hpp.del_inheritance import (  # noqa: E402
    DelMarkerLocus,
    build_inheritance_map_for_dyad,
    build_inheritance_map_for_triad,
    exclude_as_parent,
)
from hpp.del_relatedness import all_pairs, write_pairwise_tsv  # noqa: E402
from hpp.relatedness_sim import classify_edge_stdlib  # noqa: E402
from hpp.vcf_sv import read_del_calls  # noqa: E402


def _loci_from_catalogue(catalogue) -> list:
    return [
        DelMarkerLocus(
            marker_id=m.marker_id, chrom=m.chrom,
            start=m.start, end=m.end,
        )
        for m in catalogue
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delly", required=True, help="Delly DEL VCF (.vcf or .vcf.gz)")
    ap.add_argument("--manta", required=True, help="Manta DEL VCF (.vcf or .vcf.gz)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bp-tolerance", type=int, default=500)
    ap.add_argument("--reciprocal-overlap", type=float, default=0.5)
    ap.add_argument("--min-informative-pair", type=int, default=50,
                    help="minimum informative DEL markers for a pair to be reported")
    ap.add_argument("--min-excluding", type=int, default=1,
                    help=("minimum opposite-homozygote markers to exclude a "
                          "candidate parent (default 1 = strict)"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Read both VCFs ---
    print(f"[sv-pedigree] reading {args.delly}")
    delly_calls = read_del_calls(args.delly, caller_override="delly")
    print(f"[sv-pedigree]   {len(delly_calls)} Delly DEL records")
    print(f"[sv-pedigree] reading {args.manta}")
    manta_calls = read_del_calls(args.manta, caller_override="manta")
    print(f"[sv-pedigree]   {len(manta_calls)} Manta DEL records")

    # --- Merge catalogues ---
    merged = merge_two_callers(
        delly_calls, manta_calls,
        bp_tolerance=args.bp_tolerance,
        reciprocal_overlap=args.reciprocal_overlap,
    )
    n_both = sum(1 for m in merged if m.n_callers == 2)
    n_delly_only = sum(1 for m in merged
                       if m.callers == ("delly",))
    n_manta_only = sum(1 for m in merged
                       if m.callers == ("manta",))
    print(f"[sv-pedigree] merged catalogue: {len(merged)} markers "
          f"(both={n_both}, delly-only={n_delly_only}, manta-only={n_manta_only})")

    with open(outdir / "merged_del_catalogue.json", "w") as fh:
        json.dump(
            {"schema": "ngspedigree_merged_del_catalogue_v1",
             "markers": [m.to_dict() for m in merged]},
            fh, indent=2,
        )

    gmatrix = to_genotype_matrix(merged)
    samples = cohort_samples(merged)
    print(f"[sv-pedigree] cohort: {len(samples)} samples")

    # --- Pairwise relatedness on DEL markers ---
    pairs = all_pairs(samples, gmatrix, min_informative=args.min_informative_pair)
    print(f"[sv-pedigree] computed coefficients for {len(pairs)} pairs")
    write_pairwise_tsv(outdir / "pairwise_coefficients.tsv", pairs, samples)

    # --- Classify edges (Stage 1 logic, stdlib shadow) ---
    pairwise_path = outdir / "pairwise_relationship_classification.tsv"
    sample_index = {s: i for i, s in enumerate(samples)}
    with open(pairwise_path, "w") as fh:
        fh.write("sample_a\tsample_b\ta\tb\tnSites\ttheta\tIBS0\tKING\t"
                 "edge_class\tconfidence\treasons\n")
        for p in pairs:
            row = {
                "theta": p.theta,
                "IBS0": p.IBS0,
                "KING": p.theta,
                "nSites": p.n_informative,
            }
            v = classify_edge_stdlib(row)
            fh.write(
                f"{p.sample_a}\t{p.sample_b}\t"
                f"{sample_index.get(p.sample_a, '')}\t"
                f"{sample_index.get(p.sample_b, '')}\t"
                f"{p.n_informative}\t"
                f"{'' if p.theta is None else p.theta:.6f}\t"
                if False else  # bypass — handle formatting below
                ""
            )
    # Re-write cleanly (the f-string above had a ternary issue).
    with open(pairwise_path, "w") as fh:
        fh.write("sample_a\tsample_b\ta\tb\tnSites\ttheta\tIBS0\tKING\t"
                 "edge_class\tconfidence\treasons\n")
        for p in pairs:
            row = {
                "theta": p.theta,
                "IBS0": p.IBS0,
                "KING": p.theta,
                "nSites": p.n_informative,
            }
            v = classify_edge_stdlib(row)
            theta_s = "" if p.theta is None else f"{p.theta:.6f}"
            ibs0_s = "" if p.IBS0 is None else f"{p.IBS0:.6f}"
            king_s = theta_s
            fh.write(
                f"{p.sample_a}\t{p.sample_b}\t"
                f"{sample_index.get(p.sample_a, '')}\t"
                f"{sample_index.get(p.sample_b, '')}\t"
                f"{p.n_informative}\t{theta_s}\t{ibs0_s}\t{king_s}\t"
                f"{v['edge_class']}\t{v['confidence']}\t{v['reasons']}\n"
            )
    n_po = sum(1 for p in pairs
               if classify_edge_stdlib({"theta": p.theta, "IBS0": p.IBS0,
                                          "KING": p.theta,
                                          "nSites": p.n_informative})["edge_class"]
               == "parent_offspring")
    print(f"[sv-pedigree] candidate parent-offspring edges: {n_po}")

    # --- Mendelian-exclusion test on every candidate PO pair ---
    po_pairs = []
    candidate_lines = ["parent_sample_id\toffspring_sample_id\tdirection_score_A\tdirection_score_B\tn_informative\tn_excluding_A\tn_excluding_B\tassigned_parent"]
    for p in pairs:
        v = classify_edge_stdlib({"theta": p.theta, "IBS0": p.IBS0,
                                    "KING": p.theta,
                                    "nSites": p.n_informative})
        if v["edge_class"] != "parent_offspring":
            continue
        a, b = p.sample_a, p.sample_b
        # Try both directions (a→b and b→a) and pick the one with zero
        # exclusions; if both pass, leave direction ambiguous (po_dyad_only).
        ab = exclude_as_parent(a, b, gmatrix, min_excluding=args.min_excluding)
        ba = exclude_as_parent(b, a, gmatrix, min_excluding=args.min_excluding)
        assigned = "ambiguous"
        if ab.can_be_parent and not ba.can_be_parent:
            assigned = a
        elif ba.can_be_parent and not ab.can_be_parent:
            assigned = b
        elif ab.can_be_parent and ba.can_be_parent:
            assigned = "ambiguous"  # both directions compatible (dyad-only)
        else:
            assigned = "neither_compatible"
        candidate_lines.append(
            f"{a}\t{b}\t{ab.n_excluding_markers}\t{ba.n_excluding_markers}\t"
            f"{ab.n_informative}\t{ab.n_excluding_markers}\t"
            f"{ba.n_excluding_markers}\t{assigned}"
        )
        po_pairs.append((a, b, assigned))
    with open(outdir / "candidate_PO_pairs.tsv", "w") as fh:
        fh.write("\n".join(candidate_lines) + "\n")

    # --- Detect triads from candidate PO pairs ---
    # A sample with 2+ PO partners who are themselves unrelated is the
    # offspring; the two PO partners are its parents. This is the same
    # forced-offspring rule as Stage 1, expressed against the candidate
    # PO set we just built.
    po_partners: dict = {}
    for a, b, _ in po_pairs:
        po_partners.setdefault(a, set()).add(b)
        po_partners.setdefault(b, set()).add(a)
    edge_class_by_pair: dict = {}
    for p in pairs:
        v = classify_edge_stdlib({"theta": p.theta, "IBS0": p.IBS0,
                                    "KING": p.theta,
                                    "nSites": p.n_informative})
        edge_class_by_pair[(p.sample_a, p.sample_b)] = v["edge_class"]
        edge_class_by_pair[(p.sample_b, p.sample_a)] = v["edge_class"]

    triads_detected = []
    for offspring, partners in po_partners.items():
        if len(partners) < 2:
            continue
        partner_list = sorted(partners)
        for i, p1 in enumerate(partner_list):
            for p2 in partner_list[i + 1:]:
                if edge_class_by_pair.get((p1, p2)) in {"unrelated",
                                                         "third_degree"}:
                    triads_detected.append((p1, p2, offspring))

    # --- Build inheritance maps for dyads where direction was resolved ---
    loci = _loci_from_catalogue(merged)
    inh_lines = ["relationship_id\trelationship_type\tparent_sample_id\t"
                 "offspring_sample_id\tchrom\tseg_start\tseg_end\t"
                 "transmitted_allele\tn_informative_markers\tconfidence\t"
                 "recomb_event_left\trecomb_event_right\t"
                 "parental_hap_inherited\tnotes"]
    n_dyads_resolved = 0
    n_triads_resolved = 0

    # First emit triads (highest direction confidence).
    for p1, p2, offspring in triads_detected:
        # Without external sex info, treat p1 as paternal and p2 as
        # maternal by convention; downstream consumers can flip.
        triad_id = f"triad_{p1}_{p2}_{offspring}"
        maps = build_inheritance_map_for_triad(
            triad_id=triad_id,
            paternal_sample_id=p1, maternal_sample_id=p2,
            offspring_sample_id=offspring,
            loci=loci, genotype_matrix=gmatrix,
        )
        n_triads_resolved += 1
        for side, segs in maps.items():
            for s in segs:
                inh_lines.append(
                    f"{s.dyad_or_triad_id}\ttriad_{side}\t"
                    f"{s.parent_sample_id}\t{s.offspring_sample_id}\t"
                    f"{s.chrom}\t{s.seg_start}\t{s.seg_end}\t"
                    f"{s.transmitted_allele}\t{s.n_informative_markers}\t"
                    f"{s.confidence}\t{str(s.recomb_event_left).lower()}\t"
                    f"{str(s.recomb_event_right).lower()}\t"
                    f"{s.parental_hap_inherited}\t{s.notes}"
                )

    # Mark samples that participated in a triad so we don't double-emit.
    in_triad = set()
    for p1, p2, off in triads_detected:
        in_triad.add((p1, off))
        in_triad.add((p2, off))
    for a, b, assigned in po_pairs:
        if assigned in ("ambiguous", "neither_compatible"):
            continue
        parent_id = assigned
        offspring_id = b if parent_id == a else a
        # Skip if already covered by a triad emission above.
        if (parent_id, offspring_id) in in_triad:
            continue
        dyad_id = f"dyad_{parent_id}_{offspring_id}"
        segs = build_inheritance_map_for_dyad(
            dyad_id=dyad_id,
            parent_sample_id=parent_id,
            offspring_sample_id=offspring_id,
            loci=loci,
            genotype_matrix=gmatrix,
        )
        n_dyads_resolved += 1
        for s in segs:
            inh_lines.append(
                f"{s.dyad_or_triad_id}\tdyad\t"
                f"{s.parent_sample_id}\t{s.offspring_sample_id}\t"
                f"{s.chrom}\t{s.seg_start}\t{s.seg_end}\t"
                f"{s.transmitted_allele}\t{s.n_informative_markers}\t"
                f"{s.confidence}\t{str(s.recomb_event_left).lower()}\t"
                f"{str(s.recomb_event_right).lower()}\t"
                f"{s.parental_hap_inherited}\t{s.notes}"
            )
    with open(outdir / "inheritance_segments.tsv", "w") as fh:
        fh.write("\n".join(inh_lines) + "\n")
    print(f"[sv-pedigree] inheritance map: "
          f"{n_triads_resolved} triads + {n_dyads_resolved} resolved dyads")

    print(f"[sv-pedigree] outputs in {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
