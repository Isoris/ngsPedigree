#!/usr/bin/env python3
"""
09_lrr_enrichment.py — family-based odds-ratio enrichment for candidate
LRR intervals against matched background.

Tests whether candidate LRRs (low-recombination regions / putative
inversion arrangements) are enriched for transmission-compatible
block inheritance across dyads/triads, compared with random matched
background regions on the same chromosome.

The output is NOT a proof of the inversion breakpoint. It is evidence
that the LRR behaves as an inherited recombination-suppressed haplotype
block.

Usage:
  python 09_lrr_enrichment.py \
      --catalogue   merged_del_catalogue.json \
      --list-of-lrr lrr_list.tsv \
      --triads      triads.tsv \    # optional, lines "triad_id<TAB>paternal<TAB>maternal<TAB>offspring"
      --dyads       dyads.tsv \     # optional, lines "dyad_id<TAB>parent<TAB>offspring"
      --n-background 10 \
      --seed 42 \
      --out lrr_enrichment.tsv

Outputs:
  lrr_enrichment.tsv         per-LRR ORs (combined / triad-only / dyad-only)
  lrr_classifications.tsv    raw (relationship × region) block-compatibility table
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.del_inheritance import DelMarkerLocus  # noqa: E402
from hpp.lrr_enrichment import (  # noqa: E402
    Relationship,
    compute_enrichment,
    load_lrr_list,
    write_enrichment_tsv,
)


def _load_catalogue(path):
    with open(path) as fh:
        doc = json.load(fh)
    markers = doc["markers"]
    loci = [
        DelMarkerLocus(m["marker_id"], m["chrom"], m["start"], m["end"])
        for m in markers
    ]
    gmatrix = {m["marker_id"]: dict(m["genotypes"]) for m in markers}
    return loci, gmatrix


def _load_triads(path):
    rels = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            d = dict(zip(header, parts))
            rels.append(Relationship(
                relationship_id=d.get("triad_id",
                                        f"{d['paternal_sample_id']}+{d['maternal_sample_id']}+{d['offspring_sample_id']}"),
                relationship_type="triad",
                paternal_sample_id=d["paternal_sample_id"],
                maternal_sample_id=d["maternal_sample_id"],
                parent_sample_id=None,
                offspring_sample_id=d["offspring_sample_id"],
            ))
    return rels


def _load_dyads(path):
    rels = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            d = dict(zip(header, parts))
            rels.append(Relationship(
                relationship_id=d.get("dyad_id",
                                        f"{d['parent_sample_id']}->{d['offspring_sample_id']}"),
                relationship_type="dyad",
                paternal_sample_id=None,
                maternal_sample_id=None,
                parent_sample_id=d["parent_sample_id"],
                offspring_sample_id=d["offspring_sample_id"],
            ))
    return rels


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", required=True,
                    help="merged_del_catalogue.json from script 08")
    ap.add_argument("--list-of-lrr", "--list_of_LRR", dest="lrr_list", required=True,
                    help="TSV: lrr_id, chrom, start, end (header required)")
    ap.add_argument("--triads", default=None,
                    help="TSV: triad_id, paternal_sample_id, maternal_sample_id, offspring_sample_id")
    ap.add_argument("--dyads", default=None,
                    help="TSV: dyad_id, parent_sample_id, offspring_sample_id")
    ap.add_argument("--n-background", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-markers", type=int, default=3)
    ap.add_argument("--dominance-threshold", type=float, default=0.8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    loci, gmatrix = _load_catalogue(args.catalogue)
    lrrs = load_lrr_list(args.lrr_list)
    rels = []
    if args.triads:
        rels.extend(_load_triads(args.triads))
    if args.dyads:
        rels.extend(_load_dyads(args.dyads))
    if not rels:
        print("[lrr-enrich] error: no triads or dyads supplied",
              file=sys.stderr)
        return 2

    print(f"[lrr-enrich] catalogue: {len(loci)} markers")
    print(f"[lrr-enrich] LRRs: {len(lrrs)}")
    print(f"[lrr-enrich] relationships: "
          f"{sum(1 for r in rels if r.relationship_type == 'triad')} triads, "
          f"{sum(1 for r in rels if r.relationship_type == 'dyad')} dyads")

    enrichments, classifications = compute_enrichment(
        lrrs=lrrs, relationships=rels,
        loci=loci, genotype_matrix=gmatrix,
        n_background_per_lrr=args.n_background,
        seed=args.seed,
        min_markers=args.min_markers,
        dominance_threshold=args.dominance_threshold,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_enrichment_tsv(out, enrichments)
    print(f"[lrr-enrich] wrote {out}")

    classif_path = out.with_name(out.stem + "_classifications.tsv")
    with open(classif_path, "w") as fh:
        cols = ["relationship_id", "relationship_type", "parent_sample_id",
                "offspring_sample_id", "region_id", "region_kind",
                "chrom", "start", "end",
                "n_informative", "n_dominant_allele", "dominance",
                "mendelian_errors", "block_compatible"]
        fh.write("\t".join(cols) + "\n")
        for c in classifications:
            row = c.to_dict()
            fh.write("\t".join(
                f"{row[k]:.4f}" if isinstance(row[k], float) else str(row[k])
                for k in cols
            ) + "\n")
    print(f"[lrr-enrich] wrote {classif_path}")

    # Console summary.
    print()
    for e in enrichments:
        cor = e.combined.odds_ratio
        cor_s = f"{cor:.3f}" if cor is not None else "NA"
        ci = (f"[{e.combined.ci_low:.2f}, {e.combined.ci_high:.2f}]"
              if e.combined.ci_low is not None else "[NA]")
        tcor = e.triad_only.odds_ratio if e.triad_only else None
        tcor_s = f"{tcor:.3f}" if tcor is not None else "NA"
        print(f"  {e.lrr_id} ({e.chrom}:{e.start}-{e.end}): "
              f"combined OR = {cor_s} {ci}  "
              f"triad OR = {tcor_s}  "
              f"haldane={e.combined.haldane_corrected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
