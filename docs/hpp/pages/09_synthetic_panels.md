# Bloc 09 — Synthetic-panel recovery tests

| | |
|---|---|
| **modules** | `src/hpp/relatedness_sim.py`, `src/hpp/synthetic_panels.py`, `src/hpp/recovery_harness.py` |
| **CLI** | `scripts/07_run_synthetic_panels.py` |
| **suite** | `tests/test_synthetic_panels.py` (33 assertions) |
| **purpose** | regression-proof the entire ngsRelate → ngsTracts backbone against the topology mixes the real cohort will hit |

## What this is

A deterministic, stdlib-only simulator that generates ground-truth
pedigrees, ngsRelate-style coefficients, PCAngsd karyotype bands and
mtDNA haplotypes — then runs the full pipeline (Stage 1 shadow
classifier → roster construction → pedigree extract → mtDNA check →
polarization → transmission calling) and asserts recovery metrics
against the truth. Every panel is seeded so the test suite catches
regressions exactly.

## The six panels

| Panel | Families | Topology | mtDNA coverage |
|---|---|---|---|
| A — many_small_mixed | 50 | 40% dyads / 30% single-offspring triads / 30% small sibships (2–3 offspring), 1–5 catfish per family | 50% |
| B — few_large_triads | 10 | full triads with 6–10 offspring each (8–12 total) | 100% |
| C — dyad_only | 30 | isolated parent–offspring pairs | 70% |
| D — triad_only | 20 | single-offspring triads | 100% |
| E — many_small_dyads | 60 | isolated PO dyads with deliberate mtDNA swaps (8%) | 80% |
| F — eighty_twenty | 70 | 80% dyads / 20% triads (1–3 offspring) | 60% |

## Coefficient simulator

Per pair, `simulate_pair_coeffs(rel, rng, cfg)` draws a noisy
`(theta, IBS0, KING)` triple from a truncated Gaussian around the
true pedigree-relationship expectation:

| Relationship | theta | IBS0 |
|---|---|---|
| duplicate / MZ | 0.50 | 0.0000 |
| parent–offspring | 0.25 | 0.0002 |
| full sibling | 0.25 | 0.0125 |
| half sibling | 0.125 | 0.04 |
| first cousin | 0.0625 | 0.08 |
| unrelated | 0.0 | 0.12 |

Default noise (`SimConfig`): `sigma_theta = 0.010`, `sigma_ibs0 = 0.0008`.
These match genome-wide ngsRelate at ~100k sites, where PO and FS are
well-separated by IBS0.

## Shadow classifier

`classify_edge_stdlib` mirrors `STEP_PED_01_annotate_relationships.classify_edge`
exactly (same decision order, same thresholds), without the pandas
dependency. This lets the recovery tests run in any environment.

## Karyotype + mtDNA simulators

- `simulate_karyotype` walks the pedigree topologically, draws founder
  alleles from a Bernoulli at the inversion frequency, and inherits one
  allele per parent at random — emitting Mendelian-correct bands
  0 / 1 / 2 per individual.
- `simulate_mtdna` assigns each founder female a unique haplotype and
  propagates strictly matrilineally; optionally accepts a list of
  offspring whose mtDNA is swapped to inject incompatible
  mother-offspring pairs.

## Recovery metrics (per panel)

`PipelineReport` carries:

- **edge-classification accuracy** and per-class recall (PO, FS,
  half-sib, unrelated, …)
- **ground-truth vs recovered** dyad and triad counts
- **mtDNA swap injection vs detection** — the recall of the maternal
  check when swaps were introduced
- **polarization output** — dyad/triad Mendelian incompatibility counts
  and transmission-row totals

## Documented limitation surfaced by panel C / E

For families that are isolated PO pairs (`po_dyad_only` hub), Stage 1
emits `ambiguous_first_degree_PO` because **direction is genuinely
unrecoverable from a single PO edge**. So `extract_dyads` skips them
and recovered_dyads ≈ 0. Panels C and E are designed to exhibit this
behaviour. The mtDNA check on such hubs also yields 0 detected swaps
because parent_sex is unknown (blind mode), so the maternal-dyad path
does not fire. **Both are limitations, not bugs.** The simulator
documents them explicitly.

## Run the report

```
python ngsPedigree_v0.5.0/scripts/07_run_synthetic_panels.py [--panel A_many_small_mixed]
```

Current pass: all six panels at 100% edge accuracy and 100% PO recall
under the default noise; triad recovery is 100% on every panel where
triads exist.

## How this insures the PhD

Every change to the pipeline runs through these six panels in the test
suite. A regression that drops edge accuracy on any panel below the
configured floor fails the suite. The floors are set tight enough that
the catfish cohort (which sits between Panel A and Panel F in topology)
should land in safe territory if every panel passes.
