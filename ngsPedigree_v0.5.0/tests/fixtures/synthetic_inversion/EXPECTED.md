# synthetic_inversion — known truth

One candidate inversion `inv_LG01_pod`. 8 samples in 4 triads:

| triad | father | mother | offspring | parental bands | offspring band |
|---|---|---|---|---|---|
| 1 | S001 (b0) | S002 (b2) | S003 (b1) | HOM_REF × HOM_INV | HET ✓ |
| 2 | S004 (b1) | S005 (b2) | S006 (b1) | HET × HOM_INV | HET ✓ |
| 3 | S004 (b1) | S005 (b2) | S007 (b2) | HET × HOM_INV | HOM_INV ✓ |
| 4 | S004 (b1) | S005 (b2) | **S008 (b0)** | HET × HOM_INV | **HOM_REF — IMPOSSIBLE** |

Polarity hint = `band_0_is_REF`. Under that:

## Compatibility counts (under both orientations — symmetric)

- Dyads: 8 tested, 7 compatible, 1 incompatible (S005→S008: HOM_INV→HOM_REF)
- Triads: 4 tested, 3 compatible, 1 incompatible (triad 4)

Under `band_0_is_INV` the orientation flips but the same observations
are incompatible (the contradiction set is invariant under polarity
flip — see module docstring).

## Transmission calls (`chosen_polarity = band_0_is_REF`)

Triads override the bare dyads they contain.

| triad | parent | transmitted | informative for drive? |
|---|---|---|---|
| 1 | S001 (HOM_REF) | REF | no (HOM parent) |
| 1 | S002 (HOM_INV) | INV | no (HOM parent) |
| 2 | S004 (HET) | **REF** | **yes** |
| 2 | S005 (HOM_INV) | INV | no |
| 3 | S004 (HET) | **INV** | **yes** |
| 3 | S005 (HOM_INV) | INV | no |
| 4 | S004 (HET) | contradiction | no |
| 4 | S005 (HOM_INV) | contradiction | no |

Total transmissions emitted: 8 (all triad-resolved; no bare-dyad
fallback fires).

## Drive test against H0: P(INV) = 0.5

Informative HET transmissions: 2 (1 REF, 1 INV). Trivial small-n,
binomial two-sided p-value = 1.0. Paternal n = 2, INV rate 0.5;
maternal n = 0.

## mtDNA pre-flight (when `--mtdna mtdna_haplotypes.tsv` is supplied)

mtDNA labels:

| sample | role | mtDNA |
|---|---|---|
| S001 | father in triad 1 | `M_F1` (irrelevant) |
| S002 | mother in triad 1 | `M_A` |
| S003 | offspring triad 1 | `M_A` |
| S004 | father in triads 2–4 | `M_F2` (irrelevant) |
| S005 | mother in triads 2–4 | `M_B` |
| S006 | offspring triad 2 | `M_B` |
| S007 | offspring triad 3 | **`M_C`** — does not match S005 |
| S008 | offspring triad 4 | `M_B` |

Unique (mother, offspring) checks: 4.

| pair | mother | offspring | status |
|---|---|---|---|
| (S002, S003) | M_A | M_A | compatible |
| (S005, S006) | M_B | M_B | compatible |
| (S005, S007) | M_B | M_C | **incompatible** |
| (S005, S008) | M_B | M_B | compatible |

Filtering by mtDNA removes triad 3 and the maternal-dyad (S005, S007).
After filtering: 3 triads, 7 dyads. Triad 4 remains for the nuclear
contradiction test.

Two independent failure modes are now visible in the OUT JSON:
- **mtDNA pedigree failure**: S007 (excluded from polarization, listed
  in `mtdna_validation.checks` with `status="incompatible"`).
- **nuclear Mendelian failure**: S008 (still in polarization input,
  flagged in `polarization.incompatible_triads`).

After mtDNA filtering the transmission set becomes 7 rows: 6
triad-resolved + 1 bare paternal dyad S004→S007 (paternal-side
mtDNA-irrelevant; HET parent + HOM_INV offspring → transmitted INV,
informative). Drive-informative count stays at 2 (one REF from triad
2 via S004, one INV from the S004→S007 dyad).

## What ngsTracts receives

`/tmp/inv_polarization_out/polarized_transmissions.out.json` carries:

- `chosen_polarity = band_0_is_REF` (from hint)
- `polarities_symmetric = true` (the two orientations agree on
  contradiction counts)
- 8 transmission rows; 6 with `transmitted_arrangement ∈ {REF, INV}`,
  2 with `transmitted_arrangement = contradiction`
- 2 rows have `informative_for_drive = true` — these are the rows
  ngsTracts will use to seed its marker-level scan with a known
  phase polarity per (parent, offspring, inversion).
