# Expected Table A — synthetic_triad at damaging_tier=T1

Triad `triad_T1`: P_pat + P_mat → O_T1. Two segments on LG01:
seg1 (0–1500, Gold, paternal=1, maternal=2) and seg2 (1500–3000,
Bronze, paternal=2, maternal=ambiguous). All four VCF variants
(positions 100, 300, 500, 800) are in seg1; seg2 produces no rows.

## Filter (damaging_tier=T1 + synonymous)

| variant_id | consequence | kept? |
|---|---|---|
| LG01:100:A:T | stop_gained | yes |
| LG01:300:G:C | missense_variant | **no** (T3 only) |
| LG01:500:C:G | frameshift_variant | yes |
| LG01:800:T:A | synonymous_variant | yes |

## Expected rows (6 — 3 paternal + 3 maternal)

For seg1, paternal=1, maternal=2:

| variant_id | parent | hap_copy | allele_state | projection_source | confidence |
|---|---|---|---|---|---|
| LG01:100:A:T | P_pat (0/1 → hap 1) | hap_from_P1 | alt | parent_heterozygous_phased | high |
| LG01:500:C:G | P_pat (0/1 → hap 2) | hap_from_P1 | ref | parent_heterozygous_phased | high |
| LG01:800:T:A | P_pat (0/0) | — | — | (skipped) | — |
| LG01:100:A:T | P_mat (0/1 → hap 2) | hap_from_P2 | alt | parent_heterozygous_phased | high |
| LG01:500:C:G | P_mat (1/1) | hap_from_P2 | alt | parent_homozygous | high |
| LG01:800:T:A | P_mat (0/1 → hap 1) | hap_from_P2 | ref | parent_heterozygous_phased | high |

Paternal contributes 2 rows (LG01:100 hap 1, LG01:500 hap 2 → ref).
Maternal contributes 3 rows. Total 5 rows in Table A under T1.

Note: LG01:500 P_pat is ref on the inherited hap because P_pat is
phased to hap 2 but offspring inherited paternal hap 1 — the alt
allele lives on the *other* paternal hap.
