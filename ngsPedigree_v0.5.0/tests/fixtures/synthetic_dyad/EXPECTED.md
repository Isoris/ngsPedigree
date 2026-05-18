# synthetic_dyad — known truth

One dyad `dyad_PA_OA`: parent `P_A` → offspring `O_A`. One chromosome
`LG01` split into two segments at 1000 bp: Gold (hap 1 inherited),
Silver (hap 2 inherited).

## Parent P_A — variants

| variant_id | parent GT | parent_phase | expected bucket |
|---|---|---|---|
| `LG01:100:A:T` | 0/1 | hap 1   | hap1 only |
| `LG01:200:G:C` | 0/1 | hap 2   | hap2 only |
| `LG01:300:C:A` | 1/1 | (n/a)   | hap1 AND hap2 |
| `LG01:400:T:G` | 0/1 | hap 1   | hap1 only |
| `LG01:500:A:C` | 0/0 | (n/a)   | skipped (no alt) |
| `LG01:600:G:T` | 0/1 | (none)  | unphased |

## Expected `ParentHapVariants(P_A)`

- `hap1` = {`LG01:100:A:T`, `LG01:300:C:A`, `LG01:400:T:G`} — 3 variants
- `hap2` = {`LG01:200:G:C`, `LG01:300:C:A`} — 2 variants
- `unphased` = {`LG01:600:G:T`} — 1 variant
- `n_total` = 6 (the hom call is counted twice — once per hap)

## Expected Stage 3 ingestion

- `dyad_segments` = 2
- `parent_phase` rows = 3
- `dyad_ids` = ["dyad_PA_OA"]
- segment 1: `Gold`, `parental_hap_inherited="1"`, `n_informative_markers=42`
- segment 2: `Silver`, `parental_hap_inherited="2"`, `n_informative_markers=19`
