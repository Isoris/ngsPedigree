# Status, MVP roadmap, open questions

## MVP roadmap

| MVP | Scope | Built in this session? | Blocker |
|---|---|---|---|
| 1 | Stage 3 placeholder loader + parental haplotype builder + synthetic fixture | **yes** | — |
| 2 | Dyad + triad projection → Table A | **yes** (synthetic-fixture path; real joint-VCF adapter still MVP 2b) | — |
| 3 | Per-gene status → Table B | no | — |
| 4 | Triad Mendelian + transmission summaries → Tables C, D | partial (triad projection done; Mendelian + rollup pending) | — |
| 5 | KBC cross-check → Table E | no | KBC must ship in `catfish-variant-analysis` |
| 6 | Swap `stage3_placeholder` → `stage3_real` | no | ngsPedigree Stage 3 (v0.4.0 roadmap) |

## Built in this session

- `src/hpp/stage3_placeholder.py` — strict-validated loader
- `src/hpp/parental_haps.py` — `build_parental_hap_variants`
- `src/hpp/vcf_lite.py` — stdlib synthetic-fixture VCF reader
- `src/hpp/io.py` — adapter `Protocol`s + default wiring + `write_tsv`
- `src/hpp/stage3_real.py` — swap-target stub
- `src/hpp/variant_master.py` — placeholder + T1/T2/T3 classifier + real-loader stub
- `src/hpp/kbc_adapter.py` — placeholder + real-loader stub
- 8 JSON Schemas (3 inputs + 5 outputs A–E)
- `config/hpp_config.sh` — project paths + runtime knobs
- `scripts/02_ingest_stage3.py`, `scripts/03_build_parental_haps.py` — CLI smoke wrappers
- Synthetic dyad + triad fixtures with EXPECTED.md
- 43 unit tests, all passing
- 4 atlas-core JSONL blocks + tarball + smoke test
- 9 documentation pages

## Open questions (carried over from HANDOFF.md)

1. **Parent-het phase representation in real Stage 3 output.** MVP 1
   parks this in a sidecar `parent_phase.tsv`. When Stage 3 ships,
   either (a) parent_phase folds into the inheritance-map record per
   segment, or (b) it stays as a parallel file. Affects
   `stage3_real.py` signature.
2. **Bronze-segment default policy.** `include` (current) vs `exclude`
   from headline counts. Settle at audit.
3. **KBC cross-check timing.** Cross-check is optional, so HPP can
   ship MVP 4 without waiting for KBC. MVP 5 wires the optional
   terminal bloc.
4. **`compound_het_trans` deeper inspection.** Should HPP also produce
   a frame-aware coding consequence (apply variants to CDS and
   translate)? Current answer: no — that's HAPS's job
   (sibling spec, post-manuscript).
5. **Dyad partial-Mendelian thresholds.** `pass` / `warn` / `fail`
   cutoffs for the parent-hom subset. Wants empirical calibration.

## Three-cohort rule check

The 226-sample pure *C. gariepinus* hatchery cohort is the only input
scope for HPP. F1 hybrid and *C. macrocephalus* wild cohorts are
explicitly out. Adapters filter at the boundary.

## What HPP does NOT do — re-confirmed

- No cohort-level statistical claim (KBC's job)
- No family-hub aggregation (dyad / triad only)
- No SnpEff / SIFT4G / VESM re-running (HAPS's job)
- No EGO / GERP / phastCons / phyloP / orthology / CAFE inputs
- No statistical phasing of the offspring's own genotypes
- No SV burden
- No claim of fitness, balancing selection, or overdominance
- No HWE-deviation as evidence
- No cohort conflation
