# ngsPedigree

A graph annotator + per-chromosome QC tool for low-coverage population
genomics. Takes the filtered first-degree relatedness graph (from ngsRelate)
and annotates every edge with a relationship class and every node with its
possible role within its family hub. Operates in **blind mode by default** —
no metadata required. Stage 2 adds per-chromosome QC of those classifications.

Designed for the MS_Inversions_North_african_catfish manuscript (226-sample
pure *C. gariepinus* hatchery cohort, LANTA HPC).

## Stages

```
ngsRelate (whole-genome .res)
        │
        ▼
Stage 1 — STEP_PED_01_annotate_relationships.py
        │   classify edges, infer hubs, assign blind roles
        │
        ▼
ngsRelate (per-chromosome .res files)
        │
        ▼
Stage 2 — STEP_PED_02_per_chromosome_qc.py
        │   add per-chromosome columns; flag suspicious chromosomes
        │
        ▼
BEAGLE GLs (whole-genome .beagle.gz)
        │
        ▼
Stage 3 — STEP_PED_03_inheritance_map.py            [TODO next session]
            per-window microhaplotype state for confirmed PO pairs
        │
        ▼
Stage 4 — HPP: Haplotype Projection from Pedigree   [SPEC + MVP 1; ngsPedigree_v0.5.0/]
            per-offspring, per-segment, haplotype-resolved
            consequence-annotated variant projection
        │
        ▼
ngsTracts (NCO/DCO calling)                         [separate repo, downstream]
```

## Repository layout

```
ngsPedigree/
  scripts/
    STEP_PED_01_annotate_relationships.py    # Stage 1: graph annotator
    STEP_PED_02_per_chromosome_qc.py         # Stage 2: per-chromosome QC
    query_pedigree.py                         # forward/reverse query CLI
    adapt_catfish226_json.py                  # adapter for theta-only JSON inputs
  lanta/
    STEP_A07b_relatedness_per_chrom.sh       # NEW: per-chrom ngsRelate dispatch
    SLURM_A07b_relatedness_per_chrom.sh      # NEW: array job for per-chrom .res
    _REFERENCE_STEP_A03_panels.sh            # already-done; included for context
    _REFERENCE_STEP_A04_beagles.sh           # already-done; included for context
    _REFERENCE_STEP_A05_merge_beagles.sh     # already-done; included for context
    _REFERENCE_STEP_A07_relatedness.sh       # already-done; included for context
  tests/
    run_all_tests.py                          # master test runner (108 assertions)
    synthetic_12sample_fixture/               # Stage 1 unit tests, 57 assertions
    synthetic_226_realistic/                  # Stage 1 cohort-scale, 32 assertions
    synthetic_per_chrom/                      # Stage 2 per-chrom, 19 assertions
  pages/
    Relatedness_atlas.html                    # Family/Individual Evidence Hub UI
    Relatedness_atlas.js                      # page logic, Mendelian segregation tester
  docs/
    SCHEMA.md                                 # output column documentation
    SLURM_TEMPLATES.md                        # how to run on LANTA
  README.md                                   # this file
```

## What goes where

### `scripts/` — analysis logic
Pure Python. Runs anywhere with pandas. The two stage scripts and the
helpers. Stage 2 imports Stage 1's classifier directly so the per-chromosome
logic is provably the same.

### `lanta/` — LANTA pipeline pieces
Bash + SLURM. Two new files (`STEP_A07b_*`) and four reference copies
of existing pipeline scripts. The reference copies are included so this
repo is self-documenting about the prep chain that produces Stage 2's
inputs (per-RF BEAGLEs at thin-500). DO NOT RUN the reference copies —
their outputs already exist on LANTA.

### `tests/` — synthetic-only test suite
108 assertions across 3 fixtures, 0 dependencies on real cohort data.
Run with `python tests/run_all_tests.py`.

### `pages/` — browser UI for the Family/Individual Evidence Hub
Standalone HTML+JS page (`Relatedness_atlas.html` + `Relatedness_atlas.js`)
that consumes ngsPedigree outputs (Stage 1 + 2) plus the inversion
candidate TSV from the Inversion Atlas. Visual style cloned 1:1 from
Population_atlas.html (same CSS variable system, same dark/light/academic
themes, same monospace numerics, same atlas-dropdown chrome). Atlas brand
color is violet/indigo so it slots into the four-atlas palette as a
distinct fifth color. The Mendelian sub-tab is the headline analytical
feature — it tests Mendelian segregation for selected dyads or triads at
every inversion candidate, using a real binomial-exact P-value for dyads
(parent-het transmission test) and chi-square 1:2:1 for het×het triads,
plus Stouffer combination across cohorts. Demo data baked in so the page
renders fully on first open. Open by double-clicking the HTML or serving
the `pages/` directory; the four sibling atlas links assume the file
lives at `pages/` relative to the other atlas HTML files.

## Workflow on LANTA

```
1. Pre-existing pipeline (already done — produces inputs)
   ─────────────────────────────────────────────────────
   STEP_A03_panels.sh          → site panels (thin-500 etc.)
   STEP_A04_beagles.sh         → per-RF BEAGLE files at thin-500
   STEP_A05_merge_beagles.sh   → whole-genome BEAGLE
   STEP_A07_relatedness.sh     → whole-genome catfish_226_relatedness.res

2. NEW: per-chromosome ngsRelate (run once)
   ────────────────────────────────────────
   bash lanta/STEP_A07b_relatedness_per_chrom.sh
   sbatch --array=1-${N_RF}%8 lanta/SLURM_A07b_relatedness_per_chrom.sh
   → produces ${THIN_DIR}/06_relatedness/per_chrom/<rf>.res

3. ngsPedigree Stage 1 (post-process whole-genome .res)
   ────────────────────────────────────────────────────
   python scripts/STEP_PED_01_annotate_relationships.py \
     --res ${THIN_DIR}/06_relatedness/catfish_226_relatedness.res \
     --samples ${SAMPLE_LIST} \
     --outdir ${THIN_DIR}/06_relatedness/ngspedigree_stage1/ \
     --run-id cohort_226_full_v1

4. ngsPedigree Stage 2 (per-chromosome QC)
   ────────────────────────────────────────
   python scripts/STEP_PED_02_per_chromosome_qc.py \
     --per-chrom-dir ${THIN_DIR}/06_relatedness/per_chrom/ \
     --stage1-outdir ${THIN_DIR}/06_relatedness/ngspedigree_stage1/ \
     --samples ${SAMPLE_LIST} \
     --outdir ${THIN_DIR}/06_relatedness/ngspedigree_stage2/
```

## Stage 1 outputs

Three files in `--outdir`:
- `pairwise_relationship_classification.tsv` — every pair, edge_class +
  confidence + reasons
- `family_hub_roster.tsv` — every sample's hub assignment + role +
  confidence + reason
- `ngspedigree_run_envelope.json` — self-describing wrapper

## Stage 2 outputs

Four files in `--outdir`:
- `pairwise_relationship_classification.tsv` — Stage 1's table extended
  with one `edge_class_<chrom>` column per chromosome, plus
  `n_chrom_compared`, `n_chrom_disagreements`, `frac_disagreement`,
  `pair_review_flag`
- `per_chromosome_qc_flags.tsv` — long-format flags for chromosome×pair
  cells where local class disagrees with genome-wide or n_sites is low
- `per_chromosome_summary.tsv` — per-chromosome counts and mean n_sites
- `ngspedigree_run_envelope.json` — Stage 2 envelope (embeds Stage 1's)

## Edge classification

Decision order (first match wins):

| Order | Class | Rule |
|-------|-------|------|
| 1 | `duplicate_or_clone` | `theta ≥ theta_dup_min` AND `IBS0 < 0.001` |
| 2 | `parent_offspring` | `theta_first ≤ theta < theta_dup_min` AND `IBS0 < ibs0_po_max` |
| 3 | `full_sibling` | `theta_first ≤ theta < theta_dup_min` AND `IBS0 ≥ ibs0_po_max` |
| 4 | `ambiguous_first_degree` | first-degree by theta but `IBS0` missing |
| 5 | `second_degree` | `theta_second ≤ theta < theta_first` |
| 6 | `third_degree` | `theta_third ≤ theta < theta_second` |
| 7 | `unrelated` | `theta < theta_third` |

Stage 2 uses the same classifier with one threshold tweak: `ibs0_po_max`
default is 0.008 instead of 0.005, because per-chromosome IBS0 has more
sampling noise than genome-wide.

## Hub topology rules (blind mode)

| Hub type | Condition |
|----------|-----------|
| `two_parents_with_sibship` | Two unrelated nodes each have PO edges to a shared FS-clique of ≥2 |
| `parent_with_sibship` | One node has PO edges to ≥2 mutually-FS others |
| `sibship_only` | All edges FS, no PO anchor |
| `po_dyad_only` | Exactly one PO edge, n=2 |
| `duplicate_pair` | One duplicate edge, n=2 |
| `mixed_or_complex` | Doesn't fit any pattern above |

Forced-parent rule: if a node has PO edges to ≥3 nodes that form an
FS-clique, that node *must* be the parent (mathematical necessity — a
single individual can't have 3+ parents). With PO to exactly 2 mutually-FS
nodes, the rule still applies but at lower confidence.

## Test suite

```bash
python tests/run_all_tests.py
```

Runs three suites, 108 assertions total, all on synthetic data:
- 13-sample basic fixture (every topology branch + sex-promotion)
- 226-sample realistic fixture (matches real cohort hub-size distribution)
- 30-chromosome QC fixture (Stage 2; deliberate disagreements on LG14;
  low-data on LG29, LG30)

## Versioning

- v0.1.0 (2026-05-10): Stage 1 graph annotator + sex promotion + query CLI.
- v0.2.0 (2026-05-10): Stage 2 per-chromosome QC + LANTA bash bundle.
- v0.3.0 (2026-05-10): Browser UI in `pages/` (Relatedness_atlas.html +
  Relatedness_atlas.js) — Family/Individual Evidence Hub with Mendelian
  segregation tester for dyads/triads vs inversion candidate karyotypes.
- v0.3.1 (2026-05-10): Inversion-table v2 — clickable rows with expanded
  family-level inheritance roster. Each candidate scored across all
  triads with PASS/WARN/FAIL counts, diagnostic-family count, and a
  support tier (strong/moderate/weak/conflict). Strong-tier rows get a
  gold shimmer animation. New Compatibility sub-tab — breeding-partner
  finder with sex-aware mode and per-chromosome scope. TSV export for
  both inversion rosters and compatibility results.
- v0.4.0 (2026-05-10): Inversion scoring v3 — four-stage hierarchy.
  (1) family validity from genome-wide trio QC (filters out suspect
  trios so their failures don't contaminate per-inversion verdicts).
  (2) local Mendelian compatibility per valid family. (3) X-of-Y
  aggregate. (4) transmission test (binomial against 50:50 on het-
  parent gametes; concordance across families flags drive). New
  categories: PASS, WARN_CALL, WARN_FAMILY, LOCAL_CONFLICT,
  TRANSMISSION_SKEW, DRIVE_CANDIDATE, NEEDS_CROSSES. Visual rewards:
  ruby aura for STRONG-tier inversions (red shimmer, intensity scales
  with pass fraction); diamond aura for DRIVE_CANDIDATE (prismatic
  spectrum animation — biologically extraordinary inversions get the
  most special visual); black ❓ for NEEDS_CROSSES so you immediately
  see which candidates need experimental work. New "Export
  experimental cross design" button generates a priority-ordered cross
  matrix TSV with cohort-availability counts, recommended n offspring,
  and sex-aware reciprocal crosses for distinguishing maternal vs
  paternal drive.
- v0.5.0 (TODO): Stage 3 chromosome inheritance map from Beagle dosages.
- v0.5.0/ (this session): HPP (Stage 4) SPEC + MVP 1 staged under
  `ngsPedigree_v0.5.0/` — placeholder Stage 3 loader + parental haplotype
  builder + synthetic dyad/triad fixtures. SPEC ONLY beyond MVP 1;
  awaiting Stage 3 schema lock and audit. See `docs/hpp/HANDOFF.md` and
  `ngsPedigree_v0.5.0/SPEC_HPP.md`. Note: version-number reconciliation
  between the planned Stage-3 (was v0.5.0) and HPP-Stage-4
  (HANDOFF proposes v0.5.0) is deferred to the project owner.

## Provenance

Built for the MS_Inversions_North_african_catfish manuscript (226-sample
pure *C. gariepinus* hatchery cohort, LANTA HPC, account lt200308).
Stage 1 reframes the dyad-table framework per the May 10 epistemic
correction: theta-only data cannot distinguish PO from FS, so Stage 1
emits `ambiguous_first_degree` when IBS0 is missing rather than guessing.
Stage 2 unblocks the FAIL_PO_ONLY tier in the family-segregation gate
spec by providing per-chromosome compatibility checks.

LANTA pipeline pieces in `lanta/` are adapted from
`STEP_A07_relatedness.sh` (Quentin Andres) — the per-chromosome variant
mirrors the whole-genome workflow exactly, omitting only NAToRA culling,
greedy pruning, and the 3-panel figure (which are meaningful only at
genome-wide scale).
