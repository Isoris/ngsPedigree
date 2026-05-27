# synthetic_triad_violation — deliberate Mendelian-violation fixture

Triad `triad_violator`: P_v_pat + P_v_mat → O_violator. Designed to
exercise the warn/fail Mendelian path and the de-novo candidate flag.

One LG02 segment (0–2000, Gold, paternal=1, maternal=2).

## Variants and expected Mendelian verdict

| variant | consequence | P_pat | P_mat | O_violator | expected_set | inconsistent? | de novo? |
|---|---|---|---|---|---|---|---|
| LG02:100:A:T | stop_gained (T1) | 0/0 | 0/0 | 0/1 | {0/0} | **yes** | **yes** (novel "1") |
| LG02:200:C:G | frameshift (T1) | 0/0 | 0/1 | 1/1 | {0/0, 0/1} | **yes** | no (1 is in P_mat) |
| LG02:300:A:T | synonymous | 0/1 | 0/1 | 0/0 | {0/0, 0/1, 1/1} | no | — |
| LG02:400:G:C | missense (T3 only) | 1/1 | 0/0 | 0/0 | {0/1} | **yes** | no (0 is in P_mat) |

## Expected counts at damaging_tier=T1

- `testable`: 4 (all both-parents-called)
- `mendelian_inconsistent_sites`: 3 (LG02:100, LG02:200, LG02:400)
- `mendelian_inconsistent_damaging_sites`: 2 (LG02:100, LG02:200 — LG02:400 is missense, not T1 damaging)
- `n_de_novo_candidates`: 1 (LG02:100)
- `mendelian_consistency_status`: **fail** (3 ≥ 3 threshold)

## Expected counts at damaging_tier=T3

- `inconsistent_damaging_sites` rises to 3 (LG02:400 is now T3 damaging)
- Everything else unchanged
