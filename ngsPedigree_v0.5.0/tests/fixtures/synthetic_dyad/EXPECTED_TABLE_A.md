# Expected Table A — synthetic_dyad at damaging_tier=T1

Dyad `dyad_PA_OA`: P_A → O_A. Both inheritance-map segments resolve
P_A's hap 1 (Gold, seg1) and hap 2 (Silver, seg2). All six VCF
variants fall in seg1 (positions 100–600 < 1000), so seg2 produces
no rows.

## Filter (damaging_tier=T1 + synonymous)

| variant_id | consequence | kept? | reason |
|---|---|---|---|
| LG01:100:A:T | stop_gained | yes | T1 LoF |
| LG01:200:G:C | missense_variant | **no** | not T1 LoF, not synonymous |
| LG01:300:C:A | synonymous_variant | yes | synonymous control |
| LG01:400:T:G | frameshift_variant | yes | T1 LoF |
| LG01:500:A:C | (not in variant_master) | **no** | unannotated — and parent is 0/0 anyway |
| LG01:600:G:T | synonymous_variant | yes | synonymous control |

## Expected rows (4)

For seg1 (parental_hap_inherited=1, Gold):

| variant_id | hap_copy | allele_state | projection_source | confidence |
|---|---|---|---|---|
| LG01:100:A:T | hap_from_P1 | alt | parent_heterozygous_phased | high |
| LG01:300:C:A | hap_from_P1 | alt | parent_homozygous | high |
| LG01:400:T:G | hap_from_P1 | alt | parent_heterozygous_phased | high |
| LG01:600:G:T | unassigned | unknown | parent_heterozygous_unphased | unresolved |

Note row 4: LG01:600 is a het in P_A with no parent_phase entry, so
even though the segment is Gold, this variant projects as
`parent_heterozygous_unphased` → `unassigned` → `unresolved`.

Note: no row for LG01:200 because the variant is not in the T1
damaging set (it's a Tier-3 missense). At `damaging_tier=T3` it would
be kept and emit a `hap_from_P1` row with `allele_state=ref` (since
parent_phase puts it on hap 2, not the inherited hap 1).
