# ngsPedigree v0.5.0 — Stage 4: HPP (Haplotype Projection from Pedigree)

**Status:** SPEC + MVP 1 only. **Not ready for production.**
Awaiting (a) Stage 3 inheritance-map schema lock and (b) audit of `SPEC_HPP.md`.

## What this is

Stage 4 of ngsPedigree. Consumes Stage 3's per-dyad / per-triad inheritance
maps and projects parental consequence-annotated variants onto offspring
haplotypes, producing per-offspring × per-segment haplotype-resolved
deleterious-variant assignments.

See `SPEC_HPP.md` for the full specification, `../docs/hpp/HANDOFF.md` for
the implementation hand-off, and `../docs/hpp/UMBRELLA_README.md` for how
HPP sits in the four-sibling deleterious-burden stack (KBC / HPP / HAPS).

## What's built (MVP 1)

- `src/hpp/stage3_placeholder.py` — placeholder loader for the §3.2
  inheritance-map schema (swappable when Stage 3 ships).
- `src/hpp/parental_haps.py` — `build_parental_hap_variants()`:
  emits the per-parent hap-1 / hap-2 variant lists used by the
  projection step.
- `src/hpp/vcf_lite.py` — stdlib-only VCF reader for synthetic fixtures
  (no cyvcf2/pysam dependency at MVP 1).
- `schemas/inheritance_map_dyad.placeholder.schema.json`,
  `schemas/inheritance_map_triad.placeholder.schema.json`,
  `schemas/parent_phase.placeholder.schema.json` — placeholder schemas.
- `scripts/02_ingest_stage3.py`, `scripts/03_build_parental_haps.py` —
  thin CLI wrappers.
- `tests/` — synthetic dyad + synthetic triad fixtures with known truth.

## What's NOT built yet

MVP 2 (projection → table A), MVP 3 (per-gene status → table B),
MVP 4 (triad Mendelian → tables C, D), MVP 5 (KBC cross-check → table E),
MVP 6 (swap placeholder → real Stage 3). See `SPEC_HPP.md` §9.

## Run tests

```
python ngsPedigree_v0.5.0/tests/run_tests.py
```

No external deps — stdlib only.

## Hard rules (from `docs/hpp/HANDOFF.md`)

1. Three-cohort rule absolute (226-sample pure *C. gariepinus* hatchery).
2. No family-hub aggregation — dyad / triad only.
3. No cohort-level claim — that's KBC.
4. No fitness / balancing-selection / overdominance claim.
5. Confidence labels mandatory on every assignment.
6. Stage 3 placeholder schema is a contract — do not change casually.
7. No HWE-as-evidence.
