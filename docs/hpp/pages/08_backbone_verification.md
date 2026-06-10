# Bloc 08 — Backbone: ngsRelate `.res` → ngsTracts (verified)

End-to-end pipeline from per-chromosome ngsRelate output to the
ngsTracts OUT JSON. **Every stage below is implemented in this
repository and verified by the test suite running on synthetic
fixtures.**

## Pipeline

```
ngsRelate genome-wide .res ─┐
ngsRelate per-chrom   .res ─┤
                            ▼
                    STAGE 1 (scripts/STEP_PED_01_annotate_relationships.py)
                    │  pairwise_relationship_classification.tsv
                    │  family_hub_roster.tsv
                    │  ngspedigree_run_envelope.json
                            │
                            ▼
                    STAGE 2 (scripts/STEP_PED_02_per_chromosome_qc.py)
                    │  pairwise_relationship_classification.tsv  ← extended w/ pair_review_flag
                    │  per_chromosome_qc_flags.tsv
                            │
                            ▼
PCAngsd K=3 karyotype calls (TSV)  ── ─┤
mtDNA haplotypes (optional TSV)    ── ─┤
                            ▼
                    STAGE 4a (scripts/06_build_polarization_input.py)
                    │  karyotype_calls.in.json  (ngspedigree_karyotype_calls_in_v1)
                            │
                            ▼
                    STAGE 4b (scripts/05_polarize_inversion.py [--mtdna ...])
                    │  polarized_transmissions.out.json
                            │           (ngspedigree_polarized_transmissions_v1)
                            ▼
                    NGSTRACTS  (separate repo — layer 4: CO/NCO/DCO)
```

## Stage-by-stage status

| Stage | Built | Source |
|---|---|---|
| Stage 1 — edge classifier + hub topology + roles | ✓ | `scripts/STEP_PED_01_annotate_relationships.py` |
| Stage 2 — per-chrom QC + `pair_review_flag` | ✓ | `scripts/STEP_PED_02_per_chromosome_qc.py` |
| Stage 3 (inheritance map) | not needed for this path | — |
| Pedigree extraction (Stage 1/2 → dyads/triads + karyotype JSON) | ✓ | `scripts/06_build_polarization_input.py`, `src/hpp/pedigree_extract.py` |
| mtDNA maternal pre-flight (`--mtdna`) | ✓ | `src/hpp/mtdna_check.py` |
| Polarization + transmission calling + drive test | ✓ | `src/hpp/inversion_polarization.py` |
| ngsTracts OUT JSON | ✓ | `src/hpp/ngstracts_io.py` |

## Two-command run

```bash
# (1) assemble polarization IN JSON from Stage 1/2 outputs + karyotype TSV
python ngsPedigree_v0.5.0/scripts/06_build_polarization_input.py \
    --stage1-edges  ngspedigree_stage1/pairwise_relationship_classification.tsv \
    --stage1-roster ngspedigree_stage1/family_hub_roster.tsv \
    --stage2-edges  ngspedigree_stage2/pairwise_relationship_classification.tsv \
    --karyotype     pcangsd_k3_calls.tsv \
    --inversion-id  inv_LG01_pod \
    --polarity-hint band_0_is_REF \
    --out           polarization_in.json

# (2) polarize + transmissions + drive + (optional) mtDNA pre-flight
python ngsPedigree_v0.5.0/scripts/05_polarize_inversion.py \
    --in    polarization_in.json \
    --mtdna mtdna_haplotypes.tsv \
    --out   ngstracts_input.json
```

`ngstracts_input.json` is the contract ngsTracts consumes.

## What Stage 2's QC flag does to the input set

When `--stage2-edges` is supplied:

- Edges with `pair_review_flag != "OK"` are dropped from the PO
  neighbour set used for dyad/triad construction.
- A dyad whose underlying edge was dropped is omitted.
- A triad loses its parent-offspring support if either parent edge
  was dropped, and so is omitted.

On the synthetic_pedigree_pipeline fixture this drops triad 3
(S007) because Stage 2 flagged `S004↔S007`.

## What mtDNA does to the input set

When `--mtdna` is supplied:

- A maternal (triad-mother or female-flagged dyad) edge whose mtDNA
  haplotype check returns `incompatible` is dropped before
  polarization runs.
- Paternal edges are never affected.

In **blind** mode (no sex info), the `extract_dyads` step sets
`parent_sex = None` for `parent_a/parent_b` roles, so mtDNA can only
fire through the **triad** lineage (where the triad assembler
labels one parent as `maternal_sample_id`). The synthetic_inversion
fixture exercises this path.

## Verification on synthetic fixtures

| Fixture | Tests | What it verifies |
|---|---|---|
| `synthetic_pedigree_pipeline/` | 13 (`test_pedigree_extract.py`) | Stage 1/2 → IN-JSON converter (loaders, dyad/triad extraction, blind-mode warning, Stage 2 QC filter), plus 3 end-to-end CLI integration tests chaining `06 → 05` with and without `--stage2-edges` / `--mtdna` |
| `synthetic_inversion/` | (covered via end-to-end) | mtDNA TSV + IN JSON |
| `synthetic_dyad/`, `synthetic_triad/`, `synthetic_triad_violation/` | (HPP MVP 1-4) | the coding-burden side of the chain |

132 assertions in the v0.5.0 suite; all pass.

## What's still external / out of scope

- Computing `theta` / `IBS0` / `KING` from BAM/genotype data
  (ngsRelate / ANGSD).
- Computing PCAngsd K=3 karyotype bands per inversion candidate
  (PCAngsd + Hungarian arrangement-label matching).
- Computing mtDNA haplotypes from reads (an upstream consensus +
  haplogrouping step).
- Calling CO / NCO / DCO recombination tracts (ngsTracts — layer 4).
