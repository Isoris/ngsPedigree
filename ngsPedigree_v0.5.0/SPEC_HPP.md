# SPEC — HPP: Haplotype Projection from Pedigree

**Status:** SPEC ONLY — awaiting audit before implementation.
**Version:** v0.5 (split from IHC v0.4 — sibling of KBC)
**Sibling spec:** `KBC_SPEC.md` (Karyotype Burden Contrast)
**Scope:** `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, ~9× WGS, joint VCF, with ngsPedigree
Stage 3 inheritance maps once those are available.
**Position in stack:** downstream of `MODULE_CONSERVATION`, `ngsPedigree`
Stages 1 / 2 / 3, and (optionally) KBC.
**Repo home (proposed):** `catfish-variant-analysis/Modules/NN+1_hpp/`.
**Manuscript role:** **mechanism / case studies / supplementary** — adds
per-individual haplotype-resolved deleterious-variant assignments where
the inheritance map supports them.
**Blocking dependency:** ngsPedigree Stage 3 inheritance maps. HPP can be
implemented now using the placeholder schema in §3.2 and switched to
real Stage 3 inputs without changing downstream contracts.

---

## 0. One-line goal

> For each ngsPedigree-confirmed dyad (parent → offspring) or triad
> (two parents → offspring), project the parents' consequence-annotated
> variants onto the offspring's chromosome copies using the Stage 3
> inheritance map, producing **per-individual, per-segment,
> haplotype-resolved deleterious-variant assignments** for that
> offspring.

This is per-individual genotype reasoning, not a cohort-level test.
HPP does not generalise to inversions as a class; KBC does that. HPP
generalises to specific offspring whose inheritance map is high
quality.

---

## 1. Conceptual overview

### 1.1 What HPP is

A per-dyad / per-triad inference procedure. The input is the
inheritance map for one offspring-parent pair: a segment-by-segment
record of which parental haplotype the offspring inherited. The output
is, for that offspring, a haplotype-resolved variant list with each
variant tagged as residing on hap-1 (from parent 1) or hap-2 (from
parent 2), and the gene-level consequence labelled per-copy.

### 1.2 Why this is different from KBC

KBC is **cohort-level**: it tests whether AA / AB / BB karyotype classes
differ in burden, marginalising over individual identity. It does not
need to know who is related to whom; it corrects for relatedness as a
confound.

HPP is **individual-level**: for one specific offspring, it states what
that fish carries on which specific chromosome copy, at the segments
where the inheritance map is reliable. It needs explicit dyads or triads
from ngsPedigree, and it works dyad-by-dyad — never aggregating across
family hubs.

The two methods are siblings, not tiers. Stage 3 inheritance maps do not
"upgrade" KBC's confidence; they enable HPP, which answers a different
question.

### 1.3 Why family hubs are not the unit

ngsRelate first-degree edge clustering produces hubs that often contain
multiple unrelated parental pairs that contributed offspring to the same
hatchery batch. A hub is therefore not a pedigree family in the
single-cross sense. HPP works at the dyad / triad level only — the
explicit, directionally known parent-offspring relationships from
ngsPedigree Stage 2 / 3. There is no `family_transmission_summary` table
in HPP; the summary unit is the dyad (or the triad).

### 1.4 What HPP does and does not claim

**Does claim** (for an offspring O with high-quality Stage 3 map):

- "Offspring O's hap-1, in segment Seg, carries this set of damaging
  variants inherited from parent P1."
- "Therefore, at gene G in segment Seg, offspring O is het-masked
  (one copy intact from P2) or hom-exposed (both parents transmitted
  a damaging copy) or compound-het-trans (different damaging variants
  on each transmitted copy)."

**Does not claim**:

- balancing selection, overdominance, fitness;
- any cohort-level inference — that's KBC;
- anything at segments where the inheritance map is below confidence
  threshold;
- anything for samples not in a confirmed dyad or triad.

### 1.5 Three-cohort rule

HPP operates exclusively on the 226-sample pure *C. gariepinus*
hatchery cohort.

---

## 2. The biology of inheritance projection

For a dyad (parent P → offspring O):
- at each segment of the inheritance map, O's hap-from-P is identified
  as either P's hap-1 or P's hap-2;
- P's hap-1 and hap-2 each have a set of called variants (from the
  joint VCF + P's genotypes);
- in regions where P is homozygous, the inherited hap carries the
  homozygous allele unambiguously;
- in regions where P is heterozygous, the inheritance map's
  "which parental hap" call is what makes the projection possible — if
  the map says O inherited P's hap-1, then at a P-het site, O carries
  whichever allele was assigned to P's hap-1 in the upstream phasing
  used by Stage 3.

For a triad (P1 + P2 → O):
- the same projection runs from both sides;
- O's two haps are now both resolved (one from P1, one from P2);
- Mendelian-consistency checks become available at every site, and
  inconsistencies flag either a Stage 3 error or a de novo / genotyping
  error.

For a dyad only (one parent known, the other inferred):
- O's hap-from-P is projected from P;
- O's hap-from-other-parent is **left unassigned** at variants where P
  is homozygous (since the inherited allele could come from either
  parental hap); at sites where P is heterozygous, the other-parent
  allele can be deduced if O's genotype is known;
- HPP marks the unassigned side explicitly rather than guessing.

### 2.1 Confidence depends on the inheritance map

Stage 3 emits per-segment confidence (Gold / Silver / Bronze in the
project's vocabulary — see ngsPedigree Stage 3 docs). HPP inherits
that confidence at the segment level. Within a Gold segment, the
projection is treated as exact; within a Bronze segment, the
projection is treated as a guess flagged with `confidence = low`.

---

## 3. Inputs

### 3.1 Required

| File | Source | Notes |
|---|---|---|
| `variant_master_scored.tsv` | `MODULE_CONSERVATION` STEP 16 | canonical CSQ + scoring |
| Joint multisample VCF/BCF | `MODULE_CONSERVATION` STEP 03 | for parental genotypes |
| ngsPedigree Stage 2 dyad / triad table | ngsPedigree | confirmed relationships only |
| ngsPedigree Stage 3 inheritance maps | ngsPedigree (placeholder for now) | per-dyad, per-segment, per-chromosome |
| Reference FASTA | `fClaHyb_Gar_LG.fa` | |
| Sample metadata | project sample sheet | |

### 3.2 Stage 3 inheritance-map placeholder schema

Because ngsPedigree Stage 3 is not yet ready, HPP is implemented
against this placeholder schema. When Stage 3 ships, the placeholder
loader is swapped for the real reader; downstream contracts do not
change.

`inheritance_map_dyad.tsv` (one file per dyad):

| col | type | req | description |
|---|---|---|---|
| `dyad_id` | str | yes | unique, e.g. `dyad_P-S001_O-S045` |
| `parent_sample_id` | str | yes | |
| `offspring_sample_id` | str | yes | |
| `chrom` | str | yes | |
| `seg_start` | int | yes | 0-based |
| `seg_end` | int | yes | half-open |
| `parental_hap_inherited` | enum | yes | `1` / `2` / `ambiguous` (which of the parent's haps O inherited at this segment) |
| `segment_confidence` | enum | yes | `Gold` / `Silver` / `Bronze` |
| `recomb_event_left` | bool | yes | true if a recombination event borders this segment on the left |
| `recomb_event_right` | bool | yes | true if recomb on the right |
| `n_informative_markers` | int | yes | informative markers backing this segment |
| `notes` | str | opt | |

For triads, the equivalent `inheritance_map_triad.tsv` carries two
columns — `paternal_hap_inherited` and `maternal_hap_inherited` — and
otherwise the same schema.

### 3.3 Optional

| File | Source | Adds |
|---|---|---|
| KBC table B (`kbc_variant_arrangement_assignments.tsv`) | KBC output | lets HPP cross-check inversion-interval variants against the karyotype-derived arrangement assignment, flagging discordances |
| Inversion karyotype calls | PCAngsd K=3 + Hungarian | sanity-check the projection inside inversion intervals |
| ROH BED per sample | `MODULE_3` | parents' ROH is informative for the homozygous regions where projection is unambiguous |
| Splice subclass | splice module | annotation only |

HPP, like KBC, reads only the consequence / impact / SIFT class / VESM
LLR / splice subclass fields from `variant_master_scored.tsv`. No EGO
dependency.

---

## 4. Outputs

All outputs are TSV with sidecar JSON schemas. Paths under
`${BASE}/results/catfish-variant-analysis/NN+1_hpp/`.

### A. `hpp_offspring_haplotype_variants.tsv`

For each offspring × variant × hap-copy combination.

| col | type | req | description |
|---|---|---|---|
| `offspring_sample_id` | str | yes | |
| `dyad_id` or `triad_id` | str | yes | source relationship |
| `variant_id` | str | yes | chrom:pos:ref:alt |
| `chrom` / `pos` / `ref` / `alt` | yes | |
| `gene_id` / `transcript_id` | opt | |
| `consequence` / `impact` / `sift_class` / `vesm_llr` / `splice_subclass` | yes | from variant_master |
| `hap_copy` | enum | yes | `hap_from_P1` / `hap_from_P2` / `unassigned` |
| `allele_state` | enum | yes | `ref` / `alt` / `unknown` |
| `projection_source` | enum | yes | `parent_homozygous` (unambiguous) / `parent_heterozygous_phased` (Stage 3 map resolved) / `parent_heterozygous_unphased` (cannot resolve from this dyad alone) |
| `segment_confidence` | enum | yes | `Gold` / `Silver` / `Bronze` from Stage 3 |
| `confidence` | enum | yes | composite per §5 |

### B. `hpp_offspring_gene_status.tsv`

Per offspring × gene × segment.

| col | type | req | description |
|---|---|---|---|
| `offspring_sample_id` | str | yes | |
| `gene_id` / `transcript_id` | yes | |
| `chrom` / `seg_start` / `seg_end` | yes | inherits Stage 3 segment |
| `hap_from_P1_damaging_variants` | str | yes | semicolon-separated variant_ids |
| `hap_from_P2_damaging_variants` | str | yes | semicolon-separated variant_ids |
| `predicted_gene_status` | enum | yes | `reference_like` / `het_masked` / `hom_exposed_same_variant` / `compound_het_trans` (different damaging variants on each copy → both copies hit) / `compound_het_cis` (multiple damaging variants on the same copy, other copy intact) / `partially_resolved` (one copy resolved, other unassigned) / `unresolved` |
| `inside_inversion` | str | opt | inversion_id if applicable |
| `inversion_karyotype_class` | enum | opt | AA / AB / BB / NA (for cross-check with KBC) |
| `confidence` | enum | yes | from §5 |

### C. `hpp_dyad_transmission_summary.tsv` (replaces the old family_transmission table)

Per dyad — strictly per dyad, not aggregated by hub.

| col | type | req | description |
|---|---|---|---|
| `dyad_id` | str | yes | |
| `parent_sample_id` | str | yes | |
| `offspring_sample_id` | str | yes | |
| `n_segments_total` / `n_segments_Gold` / `n_segments_Silver` / `n_segments_Bronze` | int | yes | from Stage 3 |
| `n_damaging_variants_in_parent` | int | yes | total damaging variants in P, genome-wide |
| `n_damaging_variants_transmitted` | int | yes | projected onto O's hap-from-P |
| `n_damaging_variants_resolved` | int | yes | projection_source ∈ {parent_homozygous, parent_heterozygous_phased} |
| `n_damaging_variants_unresolved` | int | yes | parent_heterozygous_unphased |
| `n_genes_het_masked` / `n_genes_hom_exposed_same_variant` / `n_genes_compound_het_trans` / `n_genes_compound_het_cis` / `n_genes_partially_resolved` | int | yes | from table B |
| `mendelian_consistency_status` | enum | yes | `pass` / `warn` / `fail` / `untestable` (only triads are fully testable) |
| `mendelian_inconsistent_sites` | int | yes | |

### D. `hpp_triad_transmission_summary.tsv`

Per triad — same as C but with both parents resolved and Mendelian
checks running on every site.

| col | type | req | description |
|---|---|---|---|
| `triad_id` / `paternal_sample_id` / `maternal_sample_id` / `offspring_sample_id` | yes | |
| same segment / variant / gene counts as C | yes | |
| `mendelian_consistency_status` | enum | yes | full test available |
| `mendelian_inconsistent_sites` / `mendelian_inconsistent_damaging_sites` | int | yes | the latter is the headline check |
| `n_de_novo_candidates` | int | yes | sites where O carries an allele neither parent carries (flag, do not over-claim — could be Stage 3 error or genotyping noise at 9×) |

### E. `hpp_kbc_crosscheck.tsv` (when KBC is available)

For damaging variants inside inversion intervals.

| col | type | req | description |
|---|---|---|---|
| `variant_id` | yes | |
| `inversion_id` | yes | |
| `offspring_sample_id` | yes | |
| `hpp_hap_assignment` | enum | yes | `A` / `B` / `unassigned` — translated from `hap_copy` using the offspring's karyotype |
| `kbc_arrangement_background` | enum | yes | `A_private` / `B_private` / `shared` / `unassigned` from KBC table B |
| `concordant` | bool | yes | true if HPP and KBC agree on the arrangement background |
| `confidence` | enum | yes | min of HPP segment confidence and KBC assignment confidence |

Discordances are diagnostic — they flag either Stage 3 errors at the
inversion interval (Bronze segments) or arrangement-assignment
boundary effects in KBC (variants near inversion edges).

---

## 5. Confidence rules

Per-variant projection confidence (`hpp_offspring_haplotype_variants.confidence`):

```
if segment_confidence == 'Bronze':
    confidence = 'low'
else:
    if projection_source == 'parent_homozygous':
        confidence = 'high'                   # unambiguous from genotype
    elif projection_source == 'parent_heterozygous_phased':
        if segment_confidence == 'Gold':
            confidence = 'high'
        else:                                  # Silver
            confidence = 'medium'
    else:                                       # parent_heterozygous_unphased
        confidence = 'unresolved'
```

Per-gene status confidence (`hpp_offspring_gene_status.confidence`):
the minimum of the per-variant confidences contributing to either copy.

---

## 6. Core algorithm

### Step 1 — Ingest dyad / triad list and inheritance maps

- Validate every relationship has a Stage 3 file (or placeholder).
- Validate every offspring's chromosomes are fully covered by segments.
- Validate segment_confidence values are in the allowed enum.

### Step 2 — Build parental haplotype variant lists

For each parent in any dyad / triad:
- pull the parent's genotype at every variant from the joint VCF;
- if homozygous, both haps carry the homozygous allele;
- if heterozygous, the assignment to hap-1 vs hap-2 comes from whatever
  upstream phasing Stage 3 used (placeholder: assume a `parent_phase`
  field is available in the inheritance-map record at hetero sites; the
  real Stage 3 output will provide this);
- emit `parent_hap_variants[(parent_id, hap_no)]` keyed structures.

### Step 3 — Project onto offspring

For each (dyad, segment, variant):
- read `parental_hap_inherited` from the Stage 3 record;
- if `ambiguous`: emit `hap_copy = unassigned`,
  `projection_source = parent_heterozygous_unphased`;
- else look up the variant in
  `parent_hap_variants[(parent, parental_hap_inherited)]`;
- emit row in table A.

For triads, run the projection for both parents and tag rows by which
hap the variant came from.

### Step 4 — Per-gene status

For each (offspring, gene, segment):
- collect damaging variants assigned to hap-from-P1 and hap-from-P2;
- apply the classification rules in §7;
- emit row in table B.

### Step 5 — Mendelian consistency

For triads, at every site where both parents have called genotypes:
- compute the Mendelian-expected genotype set for the offspring;
- if the offspring's observed genotype is outside that set:
  - if at a damaging variant: flag for the headline
    `mendelian_inconsistent_damaging_sites` count;
  - if the offspring carries an allele neither parent carries: add to
    `n_de_novo_candidates` (flag, do not claim).

For dyads, the Mendelian test runs only at sites where the parent is
homozygous (the inherited allele is fixed) — full Mendelian check is
only available for triads.

### Step 6 — Per-dyad and per-triad summaries

Aggregate to tables C and D.

### Step 7 — Optional KBC cross-check

If KBC's table B is available, for each damaging variant inside an
inversion interval, compare HPP's hap assignment against KBC's
arrangement-background assignment. Emit table E.

---

## 7. Per-gene status classification rules

For each (offspring O, gene G, segment Seg):

```
d1 = damaging variants on hap_from_P1 in G ∩ Seg
d2 = damaging variants on hap_from_P2 in G ∩ Seg
unresolved = damaging variants in G ∩ Seg with hap_copy = unassigned

case 0: any unresolved
        if d1 or d2:
            → partially_resolved
        else:
            → unresolved

case 1: not d1 and not d2
        → reference_like

case 2: d1 and not d2
        if len(d1) == 1:
            → het_masked
        else:
            → compound_het_cis

case 3: d2 and not d1
        symmetric of case 2

case 4: d1 and d2
        if d1 == d2 and len(d1) == 1:
            → hom_exposed_same_variant
        else:
            → compound_het_trans
```

Notes:
- `hom_exposed_same_variant` means both copies carry the same homozygous
  damaging variant.
- `compound_het_trans` means at least one damaging variant on each copy;
  the gene has no intact copy.
- `compound_het_cis` means multiple damaging variants on a single copy,
  the other copy intact at this gene.

This is the same set of categories KBC uses for its per-gene status,
but here the assignments are made with **per-individual exactness** at
the limits of the Stage 3 confidence.

---

## 8. Implementation design

### 8.1 Directory layout

```
catfish-variant-analysis/Modules/NN+1_hpp/
  README.md
  SPEC_HPP.md
  config/
    hpp_config.sh
    stage3_placeholder/         ← swappable when real Stage 3 lands
  schemas/
    A_hpp_offspring_haplotype_variants.schema.json
    B_hpp_offspring_gene_status.schema.json
    C_hpp_dyad_transmission_summary.schema.json
    D_hpp_triad_transmission_summary.schema.json
    E_hpp_kbc_crosscheck.schema.json
    inheritance_map_dyad.placeholder.schema.json
    inheritance_map_triad.placeholder.schema.json
  scripts/
    01_validate_inputs.sh
    02_ingest_stage3.py
    03_build_parental_haps.py
    04_project_to_offspring.py
    05_per_gene_status.py
    06_mendelian_check.py
    07_dyad_summary.py
    08_triad_summary.py
    09_kbc_crosscheck.py
  src/
    hpp/
      __init__.py
      io.py
      stage3_placeholder.py     ← swap target
      stage3_real.py            ← created when Stage 3 ships
      project.py
      gene_status.py
      mendelian.py
      summary.py
      crosscheck_kbc.py
  tests/
    fixtures/                   ← synthetic dyad + triad with known truth
    test_project_homozygous.py
    test_project_heterozygous.py
    test_unresolved.py
    test_per_gene_status.py
    test_mendelian_pass.py
    test_mendelian_fail.py
    test_kbc_crosscheck.py
  outputs/
```

### 8.2 Function names

```
hpp_load_config()
hpp_validate_inputs()
hpp_load_stage3_dyad()
hpp_load_stage3_triad()
hpp_build_parental_hap_variants()
hpp_project_dyad_to_offspring()
hpp_project_triad_to_offspring()
hpp_classify_per_gene_status()
hpp_mendelian_check_dyad()
hpp_mendelian_check_triad()
hpp_summarise_dyad()
hpp_summarise_triad()
hpp_crosscheck_kbc()
hpp_export_tables()
hpp_write_report()
```

### 8.3 Language + dependencies

Same as KBC: Python 3.11, pysam / cyvcf2, pandas / polars, scipy.
Conda env `assembly` plus `envs/hpp.yaml` for whatever's missing.

---

## 9. Staged MVP plan

### MVP 1 — Stage 3 placeholder loader + parental haplotype builder

- Implement the placeholder schema in §3.2.
- Build `parent_hap_variants` structures.
- Smoke test on a synthetic dyad with known truth.

### MVP 2 — Dyad projection + table A

- Project parental variants onto offspring at every segment.
- Emit table A with confidence labels.

### MVP 3 — Per-gene status + table B

- Apply §7 classification.
- Cross-check that all damaging genes inside ROH-overlapping segments
  with parent_homozygous source come out hom_exposed_same_variant.

### MVP 4 — Triad projection + Mendelian check + tables C, D

### MVP 5 — KBC cross-check + table E

- Once KBC has shipped, run the cross-check; report discordance rates.

### MVP 6 — Switch placeholder → real Stage 3

- When ngsPedigree Stage 3 lands, swap `stage3_placeholder.py` for
  `stage3_real.py`. No downstream contracts change.

---

## 10. Pseudocode

### 10.1 Build parental hap variants

```python
def hpp_build_parental_hap_variants(parent_id, joint_vcf, parent_phase_map):
    """
    parent_phase_map[(parent_id, variant_id)] -> 1 | 2 | None
    For homozygous sites, both haps carry the allele (None is fine).
    For heterozygous sites, the upstream phasing in Stage 3 must
    provide the hap assignment.
    """
    out = {1: [], 2: []}
    for record in joint_vcf:
        gt = record.genotype(parent_id)
        if gt == '0/0':
            continue
        if gt == '1/1':
            out[1].append(record)
            out[2].append(record)
        elif gt == '0/1':
            hap = parent_phase_map.get((parent_id, record.variant_id))
            if hap == 1:   out[1].append(record)
            elif hap == 2: out[2].append(record)
            else:          pass   # unphased — emitted with projection_source = parent_heterozygous_unphased downstream
    return out
```

### 10.2 Project dyad to offspring

```python
def hpp_project_dyad_to_offspring(dyad, inheritance_map, parent_haps,
                                  joint_vcf, variant_master):
    rows = []
    for segment in inheritance_map[dyad]:
        inherited_hap = segment.parental_hap_inherited
        seg_conf      = segment.segment_confidence
        for v in variants_in_segment(joint_vcf, segment):
            if not is_damaging_or_synonymous(v, variant_master):
                continue
            parent_gt = parent_genotype(v, dyad.parent)
            if parent_gt == '1/1':
                rows.append(emit(v, hap_copy='hap_from_P1',
                                 allele_state='alt',
                                 source='parent_homozygous',
                                 seg_conf=seg_conf))
            elif parent_gt == '0/0':
                continue                # nothing to project
            elif parent_gt == '0/1':
                if inherited_hap == 'ambiguous':
                    rows.append(emit(v, hap_copy='unassigned',
                                     allele_state='unknown',
                                     source='parent_heterozygous_unphased',
                                     seg_conf=seg_conf))
                else:
                    on_inherited = v in parent_haps[inherited_hap]
                    rows.append(emit(v,
                                     hap_copy='hap_from_P1',
                                     allele_state='alt' if on_inherited else 'ref',
                                     source='parent_heterozygous_phased',
                                     seg_conf=seg_conf))
    return rows
```

### 10.3 Per-gene status

```python
def hpp_classify_per_gene_status(offspring, gene, segment, table_A):
    rows = filter_table_A(table_A, offspring, gene, segment)
    d1 = [r for r in rows if r.hap_copy == 'hap_from_P1'
                              and r.allele_state == 'alt'
                              and is_damaging(r)]
    d2 = [r for r in rows if r.hap_copy == 'hap_from_P2'
                              and r.allele_state == 'alt'
                              and is_damaging(r)]
    unresolved = [r for r in rows if r.hap_copy == 'unassigned'
                                       and is_damaging(r)]

    if unresolved:
        return 'partially_resolved' if d1 or d2 else 'unresolved'
    if not d1 and not d2:
        return 'reference_like'
    if d1 and not d2:
        return 'het_masked' if len(d1) == 1 else 'compound_het_cis'
    if d2 and not d1:
        return 'het_masked' if len(d2) == 1 else 'compound_het_cis'
    # d1 and d2
    if len(d1) == 1 and len(d2) == 1 and d1[0].variant_id == d2[0].variant_id:
        return 'hom_exposed_same_variant'
    return 'compound_het_trans'
```

### 10.4 Mendelian check (triad)

```python
def hpp_mendelian_check_triad(triad, joint_vcf, damaging_only=False):
    inconsistent = 0
    inconsistent_damaging = 0
    de_novo = 0
    for v in joint_vcf:
        gp1 = parent_genotype(v, triad.parent1)
        gp2 = parent_genotype(v, triad.parent2)
        go  = offspring_genotype(v, triad.offspring)
        if any(g is None for g in (gp1, gp2, go)):
            continue
        expected = mendelian_expected_set(gp1, gp2)
        if go not in expected:
            inconsistent += 1
            if is_damaging(v):
                inconsistent_damaging += 1
            if offspring_carries_novel_allele(go, gp1, gp2):
                de_novo += 1
    return inconsistent, inconsistent_damaging, de_novo
```

### 10.5 KBC cross-check

```python
def hpp_crosscheck_kbc(table_A_hpp, table_B_kbc, karyotypes):
    rows = []
    for r in table_A_hpp:
        if not r.variant_id in table_B_kbc:
            continue
        kbc_row = table_B_kbc[r.variant_id]
        offspring_karyotype = karyotypes[(r.offspring_sample_id,
                                          kbc_row.inversion_id)]
        hpp_arr = translate_hap_to_arrangement(r.hap_copy,
                                               offspring_karyotype,
                                               r.allele_state)
        concord = (hpp_arr == kbc_row.arrangement_background)
        rows.append(emit_crosscheck_row(r, kbc_row, hpp_arr, concord))
    return rows
```

---

## 11. QC and caveats

- **Stage 3 is the binding bottleneck.** HPP cannot run before Stage 3
  is ready; the placeholder loader exists only for early development
  and unit testing.
- **Inheritance-map quality propagates directly.** Bronze segments emit
  low-confidence projections. The headline numbers in C and D should
  be reported with and without Bronze segments included; the audit
  chat should decide which is the default.
- **Parent heterozygous sites without upstream phase** are the largest
  source of `unresolved` projections. The Stage 3 documentation must
  clarify whether Stage 3's BEAGLE phasing is exported per-parent-het
  in a form HPP can consume. If not, HPP's coverage drops sharply at
  heterozygous parental sites — a known limitation, not a bug.
- **Low coverage (~9×):** propagates to parental genotype calls; sites
  with low GQ / DP in the parent are dropped before projection.
- **Mendelian inconsistencies** can be Stage 3 errors, genotyping
  errors at 9× coverage, or genuine de novo events. HPP reports the
  count but does not classify. The `n_de_novo_candidates` field is a
  flag, not a claim.
- **Family hubs are not used.** Aggregation is strictly dyad and
  triad. The hatchery's mixed-family hub structure is a confound for
  cohort-level work (handled by KBC), not a unit of analysis for HPP.
- **KBC cross-check discordances are informative, not failures.**
  Discordance at inversion edges or in M segments of double-CO PODs is
  expected and useful diagnostic information.
- **No claim of balancing selection, fitness, or overdominance.** HPP
  reports per-individual structural-mutational state; biological
  interpretation stays at the manuscript level.
- **HWE deviation is not used as evidence** for anything.
- **Three-cohort rule:** 226-sample pure *C. gariepinus* hatchery only.

---

## 12. Out of scope

- No cohort-level burden contrast. That's KBC.
- No family-hub aggregation. Dyad / triad only.
- No SnpEff / SIFT4G / VESM re-running.
- No EGO / GERP / phastCons / phyloP / orthology / CAFE inputs.
- No statistical phasing of the offspring directly — HPP's phasing
  comes from inheritance, not from Beagle / SHAPEIT on the offspring's
  own genotypes.
- No SV burden.
- No claim of fitness or balancing selection.
- No use of HWE deviation as evidence.
- No conflation of cohorts.

---

## 13. Manuscript-facing language (supplementary)

> For a subset of offspring with high-confidence ngsPedigree Stage 3
> inheritance maps, we projected parental haplotype-resolved
> consequence-annotated variants onto the offspring genome
> segment-by-segment, producing per-individual diploid
> deleterious-variant assignments at the limits of the inheritance
> map's confidence (Gold / Silver / Bronze). For each offspring × gene,
> we classified the gene as reference-like, het-masked, hom-exposed,
> compound-heterozygous-trans, compound-heterozygous-cis, partially
> resolved, or unresolved based on the variant set inherited on each
> copy. In triads where both parents had been genotyped, Mendelian
> consistency was tested at every site, and the rate of inconsistent
> damaging-variant calls is reported per triad in Supplementary Table
> SX. Where these per-individual assignments fell inside candidate
> inversion intervals, hap-to-arrangement translations were
> cross-checked against the karyotype-derived A_private / B_private /
> shared assignment from the cohort-level KBC analysis; discordance
> rates are reported as a diagnostic, not a primary result. HPP
> provides individual-level mechanism behind the cohort-level
> POD-compatible pattern reported by KBC; it does not, and cannot,
> demonstrate balancing selection or fitness.

---

## 14. Open questions for the audit chat

1. The Stage 3 placeholder schema in §3.2 — does it match what
   ngsPedigree Stage 3 will actually emit? In particular, the
   `parent_phase` for heterozygous parental sites — where will that
   live in the real output?
2. Default Bronze-segment policy — include or exclude from headline
   counts? Current draft includes them but tags confidence; the audit
   chat may prefer the opposite default.
3. The `compound_het_trans` vs `compound_het_cis` distinction — should
   HPP also produce a Tier-4-style frame-aware coding consequence
   (apply the variants to the CDS and translate) for the trans cases,
   or is the categorical label enough for the manuscript role this
   spec has?
4. KBC cross-check — should discordances above a threshold rate
   automatically downgrade the relevant Stage 3 segments to Bronze in
   a second pass, or just be reported?
5. Triad de novo flag — what's the minimum confidence configuration
   for HPP to *report* a count of de novo candidates without that
   getting misread as a claim? Current draft says "flag, do not claim";
   want explicit text from the audit chat.
6. Repo placement — sibling module to KBC under
   `catfish-variant-analysis`, or separate repo? Current draft says
   sibling.

---

*End of HPP SPEC v0.5. SPEC ONLY — not implemented. Awaiting audit and
Stage 3 readiness.*
