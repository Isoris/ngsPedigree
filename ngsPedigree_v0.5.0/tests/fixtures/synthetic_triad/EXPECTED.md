# synthetic_triad — known truth

One triad `triad_T1`: paternal `P_pat` + maternal `P_mat` → offspring
`O_T1`. Two segments on `LG01` — a Gold segment with both parental haps
resolved, and a Bronze segment where the maternal hap is ambiguous.

## Parent P_pat — variants

| variant_id | GT | parent_phase | expected bucket |
|---|---|---|---|
| `LG01:100:A:T` | 0/1 | hap 1 | hap1 only |
| `LG01:300:G:C` | 1/1 | (n/a) | hap1 AND hap2 |
| `LG01:500:C:G` | 0/1 | hap 2 | hap2 only |
| `LG01:800:T:A` | 0/0 | (n/a) | skipped |

Expected: hap1 = {100, 300}, hap2 = {300, 500}, unphased = {}.

## Parent P_mat — variants

| variant_id | GT | parent_phase | expected bucket |
|---|---|---|---|
| `LG01:100:A:T` | 0/1 | hap 2 | hap2 only |
| `LG01:300:G:C` | 0/0 | (n/a) | skipped |
| `LG01:500:C:G` | 1/1 | (n/a) | hap1 AND hap2 |
| `LG01:800:T:A` | 0/1 | hap 1 | hap1 only |

Expected: hap1 = {500, 800}, hap2 = {100, 500}, unphased = {}.

## Expected Stage 3 ingestion

- `triad_segments` = 2
- `parent_phase` rows = 4
- `triad_ids` = ["triad_T1"]
- segment 1: Gold, paternal=1, maternal=2
- segment 2: Bronze, paternal=2, maternal=ambiguous
