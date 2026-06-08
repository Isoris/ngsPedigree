# Bloc 06 — `inversion_polarization` (ngsPedigree → ngsTracts)

| | |
|---|---|
| **module** | `ngspedigree_hpp` (sibling of the coding-burden HPP chain) |
| **scope** | layers 2–3 of the four-layer inversion-inheritance stack |
| **produces** | `polarized_transmissions.out.json` (ngsTracts consumes) |
| **schemas** | `schemas/karyotype_calls.in.schema.json`, `schemas/polarized_transmissions.out.schema.json` |
| **MVP** | **built** |

## Where this bloc stops

The four-layer stack runs:

```
Layer 1: karyotype call  (sample × inversion → band 0/1/2)
Layer 2: pedigree compatibility (dyad/triad → compatible/incompatible)
Layer 3: transmission polarization (parent → transmitted REF/INV)
Layer 4: recombination tract calling (marker-level CO/NCO/DCO)
```

**ngsPedigree implements layers 2–3 and stops.** Layer 4 is
ngsTracts' job. The hand-off between them is the OUT JSON written by
this bloc.

## 1. Inputs (IN JSON contract)

`karyotype_calls.in.json` carries:

- `inversion_id`
- `polarity_hint ∈ {band_0_is_REF, band_0_is_INV}` — see §1.4 for why
- `karyotype_calls`: per-(sample) band assignment (`0 / 1 / 2` =
  lower / middle / upper PC1 cluster)
- `dyads`: confirmed parent → offspring pairs with optional `parent_sex`
- `triads`: confirmed P1 + P2 → O trios

Loader: `src/hpp/ngstracts_io.load_karyotype_calls()` — validates the
schema string, polarity hint, band enum, and required fields; raises
`NgsTractsAdapterError` on contract violation.

## 2. Dyad / triad compatibility (layer 2)

Polarity-applied arrangement labels (`HOM_REF`, `HET`, `HOM_INV`) feed
the Mendelian-compatibility predicate.

**Dyad rule.** A dyad is compatible iff the parent's allele set
intersects the offspring's allele set (`src/hpp/inversion_polarization.py`
`_dyad_compatible`). Concretely:

| Parent | Offspring | Verdict |
|---|---|---|
| HOM_REF | HOM_REF | ✓ (parent transmits REF) |
| HOM_REF | HET | ✓ |
| HOM_REF | HOM_INV | ✗ |
| HOM_INV | HOM_REF | ✗ |
| HOM_INV | HET | ✓ |
| HOM_INV | HOM_INV | ✓ |
| HET | any | ✓ (HET parent can transmit either allele) |

**Triad rule.** From the Punnett enumeration
(`src/hpp/inversion_polarization.py` `_TRIAD_ALLOWED_OFFSPRING`):

| P1 × P2 | Allowed offspring |
|---|---|
| HOM_REF × HOM_REF | {HOM_REF} |
| HOM_REF × HET | {HOM_REF, HET} |
| HOM_REF × HOM_INV | {HET} |
| HET × HET | {HOM_REF, HET, HOM_INV} |
| HET × HOM_INV | {HET, HOM_INV} |
| HOM_INV × HOM_INV | {HOM_INV} |

## 3. Transmission calling (layer 3)

`call_transmissions()` emits one row per (parent, offspring)
relationship. **Triads override the bare dyads they contain** — when a
(parent, offspring) pair appears in both a triad and a dyad, the
triad's inference is used and the bare dyad row is suppressed.

Transmission inference uses the obvious rules:

- HOM_REF parent ⇒ transmitted REF;
- HOM_INV parent ⇒ transmitted INV;
- HET parent + HOM_REF offspring (or HOM_REF co-parent) ⇒ transmitted REF;
- HET parent + HOM_INV offspring (or HOM_INV co-parent) ⇒ transmitted INV;
- HET × HET parents producing HET offspring ⇒ **ambiguous** (cannot
  decide which parent gave which allele);
- triad contradictions tag both parents as `contradiction`.

A transmission is **informative for the Mendelian drive test** only when
the parent is HET and the call resolves to REF or INV (not `ambiguous`
or `contradiction`).

## 4. The Mendelian drive test

`binomial_two_sided_pvalue(k, n, 0.5)` is the exact two-sided binomial
against H0: `P(INV transmitted from HET parents) = 0.5`, computed in
stdlib (no scipy dependency). The aggregator
(`drive_test`) reports the overall test plus per-parent-sex
stratifications (paternal and maternal) when `parent_sex` was supplied.

## 5. Why the polarity hint is required: a symmetry caveat

Under a global relabel of band-0 ↔ band-2, the Mendelian-compatibility
predicate is **invariant**: the same (parent-band, offspring-band)
combinations are incompatible under either orientation. So:

> **Dyad/triad compatibility from band data alone cannot choose which
> homozygote band corresponds to which arrangement.**

Polarization needs an external anchor. The standard anchor is the
**reference assembly's arrangement** at the inversion (whichever
homozygote band matches the reference is HOM_REF). The IN-JSON
`polarity_hint` carries that anchor. The bloc reports compatibility
counts under BOTH orientations in
`polarization.dyad_compat` / `polarization.triad_compat`; they should
be equal up to per-individual confidence weighting (the
`polarities_symmetric` flag confirms this). The chosen polarity comes
straight from the hint.

(`polarities_symmetric = false` is a red flag: either the input was
weighted asymmetrically or one of the compatibility primitives drifted
from spec. The synthetic_inversion fixture verifies the symmetric case.)

## 6. Outputs (OUT JSON contract — ngsTracts hand-off)

`polarized_transmissions.out.json` (schema
`ngspedigree_polarized_transmissions_v1`) carries:

- `chosen_polarity` (= hint)
- `polarization.dyad_compat` / `triad_compat` under both orientations
- `polarization.polarities_symmetric` (QC flag)
- `polarization.incompatible_dyads`, `incompatible_triads` (per-pair lists)
- `transmissions[]` — one row per (parent, offspring):
  - `transmitted_arrangement ∈ {REF, INV, ambiguous, contradiction}`
  - `informative_for_drive` boolean
  - `relationship_type ∈ {dyad, triad}`
  - `parent_sex` when known
- `drive_stats` — overall + paternal + maternal counts, rates, binomial p-values

ngsTracts reads `transmissions[*]` and uses
`transmitted_arrangement` as the phase polarity per (parent, offspring,
inversion) — at each marker inside the inversion interval it knows
which arrangement the offspring inherited from each parent at the
inversion call, and so it can detect switches (CO), short
leave-return tracts (NCO), and long leave-return blocks (DCO) over the
interval.

**The key novelty:** karyotype-derived polarization provides the
phase anchor that classical three-generation pedigrees obtain from
grandparental haplotypes. That is what makes a two-generation
recombination map possible for inversion intervals.

## 7. Running the bloc

```
python ngsPedigree_v0.5.0/scripts/05_polarize_inversion.py \
    --in  ngsPedigree_v0.5.0/tests/fixtures/synthetic_inversion/karyotype_calls.json \
    --out /tmp/inv_polarization_out/polarized_transmissions.out.json
```

On the synthetic_inversion fixture: 8 samples / 8 dyads / 4 triads /
1 dyad contradiction (S005→S008) / 1 triad contradiction (S008 trio) /
8 transmission rows (all triad-resolved) / 2 informative HET-parent
transmissions / binomial p = 1.0.

## 8. What this bloc does NOT do

- It does **not** call recombination tracts (CO / NCO / DCO). That is
  ngsTracts' layer 4.
- It does **not** alter or override the karyotype calls. Per-individual
  band labels are taken as given; samples whose Mendelian status is
  contradictory under a chosen polarity are flagged (in
  `incompatible_dyads` / `incompatible_triads` and through
  `transmitted_arrangement = "contradiction"` rows) rather than
  silently re-classified.
- It does **not** integrate Stage 3 inheritance maps. Inversion
  polarization is a separate, parallel structural-mutational signal
  to the per-variant projection chain (blocs 01–04); the two
  meet downstream, not here.
- It does **not** choose the polarity from pedigree data. That decision
  is external (reference assembly or sentinel markers), supplied via
  `polarity_hint`.

## 9. Open questions

- **Cross-inversion polarity coherence.** If multiple inversions share a
  hatchery cohort, the polarities ought to be consistent in the sense
  that a sample's REF/INV genotype at each inversion is correctly
  resolved against the same reference assembly. This bloc treats each
  inversion independently; a wrapper that loops over inversions and
  flags cross-inversion polarity inconsistencies is a candidate next
  refinement.
- **Confidence-weighted compatibility.** The current compatibility
  predicate is unweighted (a contradiction is a contradiction). A
  confidence-weighted version would let low-confidence karyotype calls
  be downweighted rather than treated as hard contradictions; this is
  also where genuine polarity asymmetry between the two orientations
  could emerge.
