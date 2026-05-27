# Bloc 04 — `hpp_kbc_arrangement_crosscheck`

| | |
|---|---|
| **analysis_id** | `hpp_kbc_arrangement_crosscheck` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_kbc_crosscheck` (Table E) |
| **schema** | `ngsPedigree_v0.5.0/schemas/E_hpp_kbc_crosscheck.schema.json` |
| **MVP** | 5 (depends on KBC having shipped) |

## Goal

For damaging variants inside inversion intervals, cross-check the
per-offspring **HPP haplotype assignment** against the **KBC
karyotype-derived arrangement-background assignment**. Concordance is
the expected case; discordances are diagnostic.

## Inputs

| Dimension | Adapter |
|---|---|
| `hpp_offspring_haplotype_variants` | Table A (Bloc 01) |
| `kbc_variant_arrangement_assignments` | `kbc_adapter.PlaceholderKbcAdapter` → real loader at MVP 5 |
| `inversion_karyotype_calls` | PCAngsd K=3 + Hungarian arrangement-label matching |

## Algorithm

For each (offspring, variant) inside an inversion interval:

1. Read HPP's `hap_copy` (`hap_from_P1` / `hap_from_P2` / `unassigned`).
2. Read offspring's karyotype class (AA / AB / BB / NA).
3. Translate `hap_copy` → arrangement (`A` / `B`) under the karyotype.
4. Compare against KBC's `arrangement_background` for that variant.
5. Emit `concordant = (hpp_arr == kbc_arr)`.

## What discordances mean

| Where | Likely cause |
|---|---|
| inversion edges | KBC arrangement-assignment boundary effects |
| Bronze HPP segments | Stage 3 inheritance-map noise |
| Class M segments of double-CO PODs | recombinant haplotypes — expected disagreement |

Discordances are **reported as a diagnostic**, never as a failure of
either method. The two methods answer different questions.

## What's NOT built yet (MVP 5)

Everything. KBC must ship first (`catfish-variant-analysis` repo).
Then the cross-check algorithm lands at `src/hpp/crosscheck_kbc.py`.

## Open questions

- Should discordance above a threshold rate automatically downgrade
  the relevant Stage 3 segments to Bronze in a second pass, or just
  be reported? Current draft: report only — Stage 3 confidence is
  Stage 3's responsibility, not HPP's to retroactively edit.
