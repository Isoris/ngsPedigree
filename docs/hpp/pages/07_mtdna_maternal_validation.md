# Bloc 07 — `mtdna_maternal_validation` (pedigree pre-flight)

| | |
|---|---|
| **module** | `src/hpp/mtdna_check.py` |
| **layer** | 2 of the four-layer inversion-inheritance stack (pedigree validation) |
| **CLI** | `scripts/05_polarize_inversion.py --mtdna PATH` |
| **IN schema** | `schemas/mtdna_haplotypes.in.schema.json` |
| **OUT contract** | additional `mtdna_validation` block in `polarized_transmissions.out.json` |
| **MVP** | **built** |

## What it does

Validates **mother → offspring** pedigree edges against mitochondrial
haplotypes. Because mtDNA is maternally inherited with no paternal
contribution and no recombination, it is asymmetrically informative:

| Observation | Verdict |
|---|---|
| different mtDNA haplotype between candidate mother and offspring | **strong rejection** of maternity |
| same haplotype | compatible (NOT unique proof — multiple females in the cohort may share a haplotype) |
| paternal-side edges | **never checked** — mtDNA carries no paternity signal |
| nuclear inversion REF / INV polarity | **never inferred** — mtDNA carries no nuclear-inversion signal |

This bloc runs as a **pre-flight check** before layer 3 (polarization
+ transmission calling). Maternal (dyad, triad) relationships whose
mtDNA check returns `incompatible` are **filtered out** of the
polarization input set so the Mendelian-compatibility counts and the
drive test are not contaminated by likely-wrong-mother edges. Paternal
dyads pass through untouched.

## Input contract (`--mtdna PATH`)

TSV. Required columns: `sample_id`, `mtdna_haplotype`. Optional:
`mtdna_sequence`, `mtdna_n_sites`.

```
sample_id   mtdna_haplotype   [mtdna_sequence]   [mtdna_n_sites]
S001        M_F1
S002        M_A
...
```

Schema: `schemas/mtdna_haplotypes.in.schema.json`. The loader is
strict — missing required columns or contradictory duplicate rows
raise `MtdnaContractError`.

## Decision rule (`check_pair`)

```
1. either record missing                           → ambiguous
2. haplotype labels equal                          → compatible (distance 0)
3. labels differ, sequences both present + equal length,
   Hamming distance d ≤ threshold (default 2)      → compatible (close)
4. labels differ, sequences present, Hamming > threshold → incompatible
5. labels differ, no usable sequence pair          → incompatible
```

The Hamming tolerance is for the boundary case where two named
haplotype clusters carry near-identical sequences and the label
disagreement is essentially a clustering edge; the default threshold
(2) is conservative.

## Filter rule

Only `incompatible` pairs filter their (triad, maternal-dyad)
relationships out of the polarization input set. `compatible` and
`ambiguous` keep the relationship. Ambiguity is not evidence of an
error — absence of mtDNA is absence of evidence, not evidence of
absence.

## OUT JSON addition

`polarized_transmissions.out.json` gains an `mtdna_validation` block:

```json
"mtdna_validation": {
  "supplied": true,
  "n_relationships_checked": 4,
  "n_compatible": 3,
  "n_incompatible": 1,
  "n_ambiguous": 0,
  "n_triads_excluded": 1,
  "n_dyads_excluded": 1,
  "hamming_threshold": 2,
  "checks": [
    {
      "mother_sample_id": "S005",
      "offspring_sample_id": "S007",
      "relationship_type": "triad",
      "mother_haplotype": "M_B",
      "offspring_haplotype": "M_C",
      "distance": null,
      "status": "incompatible",
      "note": "haplotype labels differ; no sequence available for distance"
    },
    ...
  ]
}
```

When `--mtdna` is not supplied, the block is just
`{"supplied": false}`.

## Two independent failure paths surface cleanly

On the synthetic_inversion fixture with mtDNA enabled, the OUT JSON
distinguishes mtDNA-pedigree failure from nuclear-Mendelian failure:

| sample | failure mode | where reported |
|---|---|---|
| S007 | mtDNA-incompatible (mother M_B → offspring M_C) | `mtdna_validation.checks` (status `incompatible`); triad 3 dropped from polarization |
| S008 | nuclear-Mendelian-incompatible (HOM_INV × HET cannot yield HOM_REF) | `polarization.incompatible_triads`, `incompatible_dyads` |

S007 passes the nuclear band check (HET × HOM_INV → HOM_INV is OK), so
nuclear alone would not flag it. S008 passes the mtDNA check (mother
and offspring both M_B), so mtDNA alone would not flag it. The two
layers catch different errors.

## Method-paragraph language

> Mitochondrial haplotypes were used as a maternal-lineage control
> for pedigree reconstruction. Offspring were expected to share the
> mitochondrial haplotype of their assigned mother; mother–offspring
> dyads and the maternal side of triads with incompatible mitochondrial
> haplotypes were flagged as pedigree conflicts and excluded before
> nuclear Mendelian compatibility testing and inversion-karyotype
> transmission calling. Mitochondrial information was not used to
> polarize nuclear inversion arrangements (mtDNA is uninformative for
> nuclear orientation); rather, the maternal-side validation improved
> the separation of paternal versus maternal transmission for
> downstream tests of Mendelian segregation and sex-specific
> transmission distortion.

## What this bloc does NOT do

- It does **not** confirm fatherhood (mtDNA cannot).
- It does **not** uniquely identify the correct mother when multiple
  cohort females share a haplotype (mtDNA can reject more strongly
  than it can prove).
- It does **not** alter polarity (REF vs INV) — that is a nuclear
  question, anchored by the polarity hint to bloc 06.
- It does **not** modify or re-classify karyotype calls. Samples with
  mtDNA-incompatible maternal labels keep their band assignments;
  only the *relationship* is excluded from polarization.

## Run

```
python ngsPedigree_v0.5.0/scripts/05_polarize_inversion.py \
    --in    ngsPedigree_v0.5.0/tests/fixtures/synthetic_inversion/karyotype_calls.json \
    --mtdna ngsPedigree_v0.5.0/tests/fixtures/synthetic_inversion/mtdna_haplotypes.tsv \
    --out   /tmp/inv_polarization_out/polarized_with_mtdna.out.json
```

The CLI also accepts `--mtdna-hamming-threshold N` to override the
default tolerance.
