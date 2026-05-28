# Pedigree-graph annotation and haplotype-resolved deleterious-burden projection in a low-coverage hatchery cohort: methods and a worked example

## Provenance note (read first)

This document describes the mathematics and decision rules of the analysis
chain implemented in this repository, from the relatedness encoding through to
per-individual haplotype-resolved coding-consequence calls. Three provenance
classes appear:

- **In-repo, executed:** the pedigree-graph edge classifier and hub-topology
  solver (Stage 1), the per-chromosome consistency pass (Stage 2), and the
  haplotype-projection chain (Stage 4: parental haplotype construction →
  projection → per-gene zygosity classification → Mendelian/transmission
  rollup). The Stage 4 chain has been executed end-to-end on synthetic
  fixtures with result artifacts on disk; Stage 1/2 are implemented and
  unit-tested but were **not** run on the real cohort within this repository.
- **In-repo, not yet built (contract only):** the chromosome inheritance map
  (Stage 3). Stage 4 consumes it through a placeholder loader whose schema is
  a frozen contract; the real loader raises until Stage 3 ships.
- **External engines (must be confirmed against their own sources):** the
  pairwise relatedness coefficients `theta`, `IBS0`, `KING` (ANGSD/ngsRelate),
  the genotype-likelihood panels (BEAGLE), and the per-variant consequence /
  deleteriousness annotation table consumed by Stage 4 (produced by an upstream
  conservation-annotation module: SnpEff / SIFT4G / VESM / a splice module).
  The estimators those engines use are stated where Stage-1/Stage-4 logic
  depends on them and are flagged as needing confirmation against the engine.

No stage has been run on the real 226-sample cohort inside this repository; the
Results section is therefore a synthetic-fixture validation with real on-disk
numbers, and is labelled as such.

---

## Abstract

We analyse a 226-sample, ~9× whole-genome cohort of pure *Clarias gariepinus*
(reference `fClaHyb_Gar_LG`) drawn from a mixed-family hatchery. The analytical
problem is twofold: (i) reconstruct enough of the first-degree relatedness
structure to identify directional parent→offspring relationships without
metadata, and (ii) use those relationships to state, for a given offspring,
which deleterious coding variants sit on which inherited chromosome copy, and
hence whether each gene is functionally exposed or masked. The chain encodes
pairwise relatedness as (`theta`, `IBS0`) pairs and applies a first-match-wins
threshold classifier to label each edge (parent–offspring, full-sibling,
duplicate, 2nd/3rd-degree, unrelated, or — when `IBS0` is absent —
ambiguous-first-degree); a union-find pass groups edges into hubs and a
forced-parent rule (a node with parent–offspring edges to ≥2 mutually
full-sib nodes must be the parent) assigns directional roles; a per-chromosome
replay of the same classifier flags pairs whose local class disagrees with the
genome-wide call. Given a (placeholder) per-segment inheritance map, parental
genotypes are split into haplotype-resolved variant sets, projected onto the
offspring, classified per (offspring, gene, segment) into a seven-state
zygosity-exposure enum, and rolled up with a Mendelian-consistency test. On the
synthetic triad fixture the chain produces the headline qualitative result it is
designed to detect: a stop-gained allele transmitted from both parents to the
offspring at gene `G1` is called `hom_exposed_same_variant`, while a maternal-only
frameshift at `G5` is called `het_masked`, with full Mendelian consistency
(`pass`, 0 inconsistent sites).

---

## 1. Introduction

In a hatchery cohort, individuals are related through a tangle of overlapping
crosses, and low coverage (~9×) makes genotypes uncertain. Two facts make the
deleterious-burden question hard. First, relatedness must be recovered from
summary coefficients (`theta`, `IBS0`) rather than from a known pedigree, and
the dominant first-degree ambiguity — parent–offspring versus full-sibling —
cannot be resolved from `theta` alone, because both have kinship ≈ ¼; only the
zero-IBS-sharing rate `IBS0` separates them (a parent and child share at least
one allele at essentially every site, so `IBS0 ≈ 0`, whereas full sibs are
homozygous-discordant at a non-trivial fraction of sites). Second, a deleterious
allele's functional impact depends on whether the *other* chromosome copy is
intact, which is a per-individual, per-haplotype fact that a genotype table
cannot express: a heterozygous genotype is masked if the trans copy is wild-type
and exposed if a second hit lands in trans. Resolving exposure therefore
requires knowing which parental haplotype each offspring inherited at each locus.

The remaining sections follow the data flow. §2 fixes the input encodings
(relatedness coefficients; genotype strings; the per-variant annotation
contract). §3 defines the core detection statistic — the edge classifier — and
its thresholds. §4 covers graph construction and the hub-topology / forced-parent
role assignment. §5 is the per-chromosome consistency pass. §6 states the
inheritance-map contract (the not-yet-built Stage 3) and the parental-haplotype
construction that consumes it. §7 is the projection onto offspring haplotypes.
§8 is the per-(offspring, gene, segment) zygosity-exposure classifier. §9 is the
Mendelian-consistency test and per-relationship transmission rollup. §10 covers
serialization and reproducibility. The Results section (§11) gives a worked
example on synthetic fixtures with on-disk numbers and an explicit
not-produced subsection.

---

## 2. Input encoding

### 2.1 Relatedness coefficients (external; ngsRelate)

The Stage-1 input is the standard 23-column ngsRelate table; the classifier
requires only `a, b, theta, IBS0` and uses `KING, R0, R1, J9, nSites` when
present (`scripts/STEP_PED_01_annotate_relationships.py:88`,
`scripts/STEP_PED_01_annotate_relationships.py:89`). `theta` is the
ngsRelate kinship coefficient (≈ relatedness/2), `IBS0` the fraction of sites
at which the pair shares zero alleles identical-by-state, and `KING` the
KING-robust kinship estimate. **These estimators are computed by ngsRelate, not
in this repo, and their definitions must be confirmed against the ngsRelate /
ANGSD source.** Column-name case is normalised on load
(`scripts/STEP_PED_01_annotate_relationships.py:496`).

A theta-only adapter exists for inputs that carry no IBS sharing rates: it
writes `IBS0` as `NA` and a constant `nSites=100000`, deliberately forcing the
first-degree ambiguity branch downstream rather than guessing direction
(`scripts/adapt_catfish226_json.py:52`).

### 2.2 Genotype encoding (Stage 4)

Stage 4 reads a multisample VCF. Genotype strings are normalised by collapsing
the phased separator (`|` → `/`) so that phase is *not* taken from the VCF
(`ngsPedigree_v0.5.0/src/hpp/vcf_lite.py:41`); haplotype phase enters only
through the inheritance map and a parent-phase sidecar (§6). A variant is keyed
by `chrom:pos:ref:alt` (`ngsPedigree_v0.5.0/src/hpp/vcf_lite.py:34`).
Multiallelic records are rejected and must be split upstream
(`ngsPedigree_v0.5.0/src/hpp/vcf_lite.py:69`). The fixture reader here is a
stdlib parser; a production run swaps in a cyvcf2/pysam-backed adapter exposing
the same record interface.

### 2.3 Per-variant annotation contract (external)

Stage 4 reads only the consequence/deleteriousness fields of an upstream
variant table: `consequence`, `impact`, `sift_class`, `vesm_llr`,
`splice_subclass` (`ngsPedigree_v0.5.0/src/hpp/variant_master.py:48`). The
damaging/neutral call is computed *in-repo* from these fields under a tiered
rule (§8.1), but the fields themselves are produced by an external
conservation-annotation engine (SnpEff/SIFT4G/VESM/splice) and **must be
confirmed against that engine.**

---

## 3. Core detection statistic: the relatedness-edge classifier

The atom of the chain is a pure function mapping one coefficient row to a
relationship class, a confidence, and diagnostic reasons
(`scripts/STEP_PED_01_annotate_relationships.py:109`). The decision order is
first-match-wins, so the order of the branches is itself part of the estimator:

```
1. duplicate_or_clone   : theta ≥ theta_dup_min  AND IBS0 < 0.001
2. parent_offspring     : theta ≥ theta_first AND IBS0 < ibs0_po_max
3. full_sibling         : theta ≥ theta_first AND IBS0 ≥ ibs0_po_max
4. ambiguous_first_degree: theta ≥ theta_first AND IBS0 missing
5. second_degree        : theta_second ≤ theta < theta_first
6. third_degree         : theta_third  ≤ theta < theta_second
7. unrelated            : theta < theta_third
```

Branch 1 (duplicate) is tested first because a clone also has first-degree-level
`theta` and would otherwise be captured by branch 2/3; the near-zero `IBS0`
distinguishes it (`scripts/STEP_PED_01_annotate_relationships.py:144`). Inside
the first-degree band, the parent–offspring versus full-sibling split is made
entirely by `IBS0`: below `ibs0_po_max` is parent–offspring, at-or-above is
full-sibling (`scripts/STEP_PED_01_annotate_relationships.py:164`,
`scripts/STEP_PED_01_annotate_relationships.py:174`). If `IBS0` is missing the
function refuses to guess and emits `ambiguous_first_degree` at low confidence
(`scripts/STEP_PED_01_annotate_relationships.py:159`) — this is the branch the
theta-only adapter (§2.1) deliberately triggers.

Default thresholds (the Manichaikul first-degree cut and its successive halvings)
are: `theta_first = 0.177`
(`scripts/STEP_PED_01_annotate_relationships.py:542`),
`theta_second = 0.0884` (`:544`), `theta_third = 0.0442` (`:546`),
`ibs0_po_max = 0.005` (`:548`), `theta_dup_min = 0.45` (`:550`).

Confidence is high by default and is downgraded by two rules: a site-count floor
(`nSites < 1000` → low confidence,
`scripts/STEP_PED_01_annotate_relationships.py:139`) and a KING cross-check
(a parent–offspring call with `KING < 0.15` is downgraded,
`scripts/STEP_PED_01_annotate_relationships.py:166`). An unexpectedly large
`IBS0` accompanying a duplicate-level `theta` is also flagged
(`scripts/STEP_PED_01_annotate_relationships.py:148`).

---

## 4. Graph construction and directional role assignment

### 4.1 First-degree subgraph and components

Edges whose class is in {parent_offspring, full_sibling, duplicate_or_clone,
ambiguous_first_degree} are treated as in-hub
(`scripts/STEP_PED_01_annotate_relationships.py:615`) and grouped into connected
components by a path-compressing union-find
(`scripts/STEP_PED_01_annotate_relationships.py:220`).

### 4.2 Hub topology and the forced-parent rule

Within each component the solver assigns roles
(`scripts/STEP_PED_01_annotate_relationships.py:243`). The decisive estimator is
the **forced-parent rule**: a node with parent–offspring edges to ≥2 other
nodes that themselves form a full-sibling clique must be the parent, because a
set of mutual full sibs cannot all be parents of one another's shared neighbour
(`scripts/STEP_PED_01_annotate_relationships.py:316`–`:331`). The clique test
is an explicit all-pairs check over the candidate's PO-neighbours
(`scripts/STEP_PED_01_annotate_relationships.py:323`–`:329`).

Topology outcomes:

- **two_parents_with_sibship** — two forced-parent nodes that are themselves
  unrelated (their edge is unrelated/2nd/3rd degree), sharing an offspring set
  (`scripts/STEP_PED_01_annotate_relationships.py:342`–`:365`).
- **parent_with_sibship** — one forced parent; the parent role is reported as
  `forced_parent` at high confidence only when it subtends ≥3 sibs, else
  `likely_parent` at medium (`scripts/STEP_PED_01_annotate_relationships.py:374`).
- **po_dyad_only** — a single parent–offspring edge; direction is *not*
  decidable and both nodes are left ambiguous
  (`scripts/STEP_PED_01_annotate_relationships.py:275`).
- **sibship_only**, **duplicate_pair**, **mixed_or_complex** for the remaining
  patterns (`scripts/STEP_PED_01_annotate_relationships.py:394`,
  `:268`, `:405`).

### 4.3 Sex-based label promotion (optional, metadata-only)

If a sex map is supplied, parent roles are promoted to mother/father only where
topology has already fixed the parental role and the sexes are unambiguous; sex
never changes an edge class or a direction
(`scripts/STEP_PED_01_annotate_relationships.py:422`).

---

## 5. Per-chromosome consistency pass

Stage 2 replays the *same* classifier on each chromosome's coefficients
(imported, not re-implemented:
`scripts/STEP_PED_02_per_chromosome_qc.py:86`), with one loosened parameter:
the per-chromosome parent–offspring `IBS0` ceiling is `0.008` rather than
`0.005`, because per-chromosome `IBS0` is noisier
(`scripts/STEP_PED_02_per_chromosome_qc.py:160`). Chromosomes with
`nSites < 1000` are marked `low_data` and excluded from comparison
(`scripts/STEP_PED_02_per_chromosome_qc.py:137`,
`scripts/STEP_PED_02_per_chromosome_qc.py:159`).

For each pair, the genome-wide class is compared against every non-low-data
chromosome class, yielding

```
frac_disagreement = n_chrom_disagreements / n_chrom_compared
```

(`scripts/STEP_PED_02_per_chromosome_qc.py:259`–`:261`). A pair is flagged for
review when `n_chrom_disagreements ≥ disagreement_threshold` (default 3)
(`scripts/STEP_PED_02_per_chromosome_qc.py:300`,
`scripts/STEP_PED_02_per_chromosome_qc.py:161`). Stage 2 adds QC columns only;
it does not re-solve hub topology.

---

## 6. Inheritance map (contract) and parental-haplotype construction

### 6.1 The inheritance-map contract (Stage 3 — not built in-repo)

Stage 4 consumes a per-segment inheritance map: for each (relationship, chrom,
segment), which parental haplotype the offspring inherited
(`1` / `2` / `ambiguous`) and a segment confidence (`Gold` / `Silver` /
`Bronze`). The placeholder loader validates these enums and the half-open
segment geometry (`seg_end > seg_start`)
(`ngsPedigree_v0.5.0/src/hpp/stage3_placeholder.py:23`,
`ngsPedigree_v0.5.0/src/hpp/stage3_placeholder.py:25`,
`ngsPedigree_v0.5.0/src/hpp/stage3_placeholder.py:188`). The map itself is
produced by Stage 3, which is **not implemented**; the real adapter raises
until it ships. Parent-heterozygous-site phase is carried in a sidecar table
keyed by (parent, variant) → haplotype
(`ngsPedigree_v0.5.0/src/hpp/stage3_placeholder.py:286`); its final location in
the Stage-3 output is an open contract question.

### 6.2 Parental-haplotype construction

For each parent, called genotypes are split into two haplotype variant lists
(`ngsPedigree_v0.5.0/src/hpp/parental_haps.py:47`):

- homozygous-alt (`1/1`) → the variant is placed on **both** haplotypes
  (`ngsPedigree_v0.5.0/src/hpp/parental_haps.py:57`–`:59`);
- heterozygous (`0/1`) with a sidecar phase of 1 or 2 → placed on that
  haplotype (`ngsPedigree_v0.5.0/src/hpp/parental_haps.py:61`–`:66`);
- heterozygous with no phase → placed in an `unphased` bucket
  (`ngsPedigree_v0.5.0/src/hpp/parental_haps.py:68`);
- homozygous-ref / missing → dropped.

---

## 7. Projection onto offspring haplotypes

For each segment, each candidate variant is first kept only if it is damaging
(under the active tier, §8.1) or synonymous — the synonymous class is retained
as a neutral comparator (`ngsPedigree_v0.5.0/src/hpp/project.py:87`). Variant–
segment membership uses the half-open interval on 0-based coordinates
(`pos−1 ∈ [seg_start, seg_end)`,
`ngsPedigree_v0.5.0/src/hpp/project.py:133`). The per-parent projection rule
(`ngsPedigree_v0.5.0/src/hpp/project.py:143`) is:

- parent `1/1` → offspring carries ALT on the inherited copy;
  `projection_source = parent_homozygous`
  (`ngsPedigree_v0.5.0/src/hpp/project.py:168`);
- parent `0/1`, segment phase known → ALT if the variant sits on the inherited
  haplotype, else REF; `projection_source = parent_heterozygous_phased`
  (`ngsPedigree_v0.5.0/src/hpp/project.py:182`–`:185`);
- parent `0/1`, segment `ambiguous` or variant unphased → emitted as
  `hap_copy = unassigned`, `allele_state = unknown`,
  `projection_source = parent_heterozygous_unphased`.

The per-assignment confidence is a deterministic function of projection source
and segment confidence (`ngsPedigree_v0.5.0/src/hpp/project.py:67`):

```
Bronze segment                         -> low
parent_homozygous                      -> high
parent_heterozygous_phased  & Gold     -> high
parent_heterozygous_phased  & Silver   -> medium
parent_heterozygous_unphased           -> unresolved
```

Dyad projection emits the known-parent copy only; triad projection runs both
parents, tagging paternal contributions as `hap_from_P1` and maternal as
`hap_from_P2` (`ngsPedigree_v0.5.0/src/hpp/project.py:222`,
`ngsPedigree_v0.5.0/src/hpp/project.py:269`).

---

## 8. Per-(offspring, gene, segment) zygosity-exposure classification

### 8.1 Damaging-variant tiers (in-repo decision over external annotations)

The damaging call is a nested three-tier rule
(`ngsPedigree_v0.5.0/src/hpp/variant_master.py:79`):

- **T1** (high-confidence loss-of-function): `consequence ∈ {stop_gained,
  frameshift_variant, start_lost, splice_donor_variant,
  splice_acceptor_variant}` (`ngsPedigree_v0.5.0/src/hpp/variant_master.py:23`);
- **T2** adds validated Class-A splice subclasses
  (`ngsPedigree_v0.5.0/src/hpp/variant_master.py:34`);
- **T3** adds model-dependent missense: SIFT-deleterious or
  `vesm_llr ≤ −7.0` (`ngsPedigree_v0.5.0/src/hpp/variant_master.py:42`,
  `ngsPedigree_v0.5.0/src/hpp/variant_master.py:43`).

The VESM/SIFT scores and the splice subclass are external inputs; the tier
thresholds and membership are in-repo.

### 8.2 The seven-state classifier

Group the kept Table-A rows by (offspring, gene, segment); let `d1`, `d2` be the
sets of damaging ALT variants on the two inherited copies and `unresolved` the
damaging variants with `hap_copy = unassigned`. The rule
(`ngsPedigree_v0.5.0/src/hpp/gene_status.py:104`) is:

```
any unresolved              -> partially_resolved  (if d1 or d2) else unresolved
not d1 and not d2           -> reference_like
d1 xor d2, |hit set| == 1   -> het_masked
d1 xor d2, |hit set| >  1   -> compound_het_cis
d1 and d2, same single var  -> hom_exposed_same_variant
d1 and d2, otherwise        -> compound_het_trans
```

(`ngsPedigree_v0.5.0/src/hpp/gene_status.py:107`–`:117`). A REF call on a copy
means that copy does not carry the damaging allele and so does not enter `d1`/`d2`;
this is what lets a heterozygous parent transmit either a hit or a wild-type copy.

Per-gene confidence is the **worst** (lowest-quality) of the contributing
per-assignment confidences, under the order high < medium < low < unresolved
(`ngsPedigree_v0.5.0/src/hpp/gene_status.py:59`,
`ngsPedigree_v0.5.0/src/hpp/gene_status.py:63`). Segment coordinates are
re-derived by locating each variant's position in the inheritance map
(`ngsPedigree_v0.5.0/src/hpp/gene_status.py:79`), since Table A does not carry
segment bounds.

---

## 9. Mendelian consistency and per-relationship transmission rollup

### 9.1 Mendelian primitives

The expected offspring-genotype set is the Punnett enumeration over one allele
drawn from each parent (`ngsPedigree_v0.5.0/src/hpp/mendelian.py:27`):

```
expected(gp1, gp2) = { sorted(a1,a2) : a1 ∈ alleles(gp1), a2 ∈ alleles(gp2) }
```

A triad site is inconsistent when the offspring genotype is outside this set
(uncalled sites are never inconsistencies)
(`ngsPedigree_v0.5.0/src/hpp/mendelian.py:41`). A de-novo candidate is flagged
when the offspring carries an allele present in neither parent — a flag, not a
claim (`ngsPedigree_v0.5.0/src/hpp/mendelian.py:47`). For dyads only the
parent-homozygous subset is testable: the offspring must carry the parent's
fixed allele (`ngsPedigree_v0.5.0/src/hpp/mendelian.py:54`).

The status map is (`ngsPedigree_v0.5.0/src/hpp/mendelian.py:76`):

```
0 testable sites      -> untestable
0 inconsistent        -> pass
1-2 inconsistent      -> warn
>= 3 inconsistent     -> fail
```

The 1/2/≥3 thresholds are provisional defaults pending empirical calibration.

### 9.2 Transmission rollup

The per-relationship summary aggregates segment-confidence counts, parental
damaging-variant counts (`ngsPedigree_v0.5.0/src/hpp/summary.py:111`),
projection counts — transmitted = ALT assignments, resolved =
homozygous-or-phased sources, unresolved = unphased
(`ngsPedigree_v0.5.0/src/hpp/summary.py:141`) — and the per-gene status tallies
(`ngsPedigree_v0.5.0/src/hpp/summary.py:127`). Dyads emit a 19-field record
with partial Mendelian (`ngsPedigree_v0.5.0/src/hpp/summary.py:167`); triads
emit a 23-field record with the full Mendelian test plus de-novo and
inconsistent-damaging counts (`ngsPedigree_v0.5.0/src/hpp/summary.py:232`).

---

## 10. Serialization and reproducibility

Stage 1 and Stage 2 each emit a self-describing JSON envelope recording the
exact thresholds, input paths, mode, and per-class / per-hub counts
(`scripts/STEP_PED_01_annotate_relationships.py:685`,
`scripts/STEP_PED_02_per_chromosome_qc.py:320`); Stage 2 nests the Stage 1
envelope so a run is fully traceable. Stage 4 tables are emitted as TSV whose
column order is fixed by the dataclass field order and checked against the
JSON-Schema column lists in tests. The end-to-end Stage-4 driver is a single
script (`ngsPedigree_v0.5.0/scripts/04_run_hpp_pipeline.py`); the unit suite
(114 assertions, stdlib-only) reproduces every number below.

---

## 11. Results (synthetic-fixture validation)

**Scope of this section.** No stage was run on the real 226-sample cohort in
this repository. The numbers below come from executing the Stage-4 driver
(§10) on the committed synthetic fixtures; the cited result files are the
TSVs it wrote. They demonstrate that the estimators produce the intended
qualitative calls on inputs with known truth.

### 11.1 Dyad worked example (`damaging_tier = T1`)

Inputs: one parent→offspring dyad (`P_A`→`O_A`), six genotyped sites in two
segments (Gold + Silver) on `LG01`. Result file:
`/tmp/hpp_dyad_out/table_A_offspring_haplotype_variants.tsv` (4 rows).

| variant | gene | consequence | hap_copy | allele_state | projection_source | confidence |
|---|---|---|---|---|---|---|
| LG01:100:A:T | G1 | stop_gained | hap_from_P1 | alt | parent_heterozygous_phased | high |
| LG01:300:C:A | G3 | synonymous_variant | hap_from_P1 | alt | parent_homozygous | high |
| LG01:400:T:G | G4 | frameshift_variant | hap_from_P1 | alt | parent_heterozygous_phased | high |
| LG01:600:G:T | G6 | synonymous_variant | unassigned | unknown | parent_heterozygous_unphased | unresolved |

The T3-only missense `LG01:200` is correctly excluded at T1, and the
parent-heterozygous site with no phase (`LG01:600`) correctly degrades to
`unassigned` / `unresolved`. Per-gene status
(`/tmp/hpp_dyad_out/table_B_offspring_gene_status.tsv`, 4 rows): `G1` and `G4`
→ `het_masked` (high); `G3` → `reference_like` (high, synonymous-only); `G6`
→ `reference_like` (unresolved confidence, inherited from the unphased site).
Dyad rollup (`/tmp/hpp_dyad_out/table_C_dyad_transmission_summary.tsv`):
`n_segments_total=2` (Gold 1, Silver 1), `n_damaging_variants_in_parent=2`,
`transmitted=2`, `resolved=2`, `unresolved=0`, `n_genes_het_masked=2`,
`mendelian_consistency_status=pass`, `mendelian_inconsistent_sites=0`.

### 11.2 Triad worked example (`damaging_tier = T1`) — the headline call

Inputs: one triad (`P_pat` + `P_mat` → `O_T1`), four genotyped sites, two
segments (Gold + Bronze) on `LG01`. Per-gene status
(`/tmp/hpp_triad_out/table_B_offspring_gene_status.tsv`, 3 rows):

| gene | hap_from_P1 damaging | hap_from_P2 damaging | predicted_gene_status | confidence |
|---|---|---|---|---|
| G1 | LG01:100:A:T | LG01:100:A:T | **hom_exposed_same_variant** | high |
| G5 | — | LG01:500:C:G | het_masked | high |
| G8 | — | — | reference_like | high |

This is the qualitative result the chain exists to detect: a stop-gained allele
transmitted on **both** inherited copies (one from each parent) at `G1` is
called `hom_exposed_same_variant` (no intact copy), while a maternal-only
frameshift at `G5` is `het_masked` (paternal copy intact). Triad rollup
(`/tmp/hpp_triad_out/table_D_triad_transmission_summary.tsv`):
`n_damaging_variants_in_paternal=2`, `…_in_maternal=2`, `transmitted=3`,
`resolved=4`, `n_genes_het_masked=1`, `n_genes_hom_exposed_same_variant=1`,
`mendelian_consistency_status=pass`, all of
`mendelian_inconsistent_sites` / `…_damaging_sites` / `n_de_novo_candidates` = 0.

### 11.3 Adversarial-Mendelian worked example (`damaging_tier = T1`)

A deliberately Mendelian-violating triad (`P_v_pat` + `P_v_mat` →
`O_violator`, four sites) exercises the failure path. Result file:
`/tmp/hpp_violation_out/table_D_triad_transmission_summary.tsv`:
`mendelian_consistency_status=fail`, `mendelian_inconsistent_sites=3`,
`mendelian_inconsistent_damaging_sites=2`, `n_de_novo_candidates=1`. The single
de-novo candidate is the site where both parents are `0/0` and the offspring
carries an ALT allele; the two damaging inconsistencies are the stop-gained and
frameshift sites, while the T3-only missense violation is correctly *not*
counted as damaging at T1.

### 11.4 What was NOT produced on this data

- **No real-cohort run exists in this repository.** Stages 1, 2, and 4 were
  not executed on the 226-sample cohort here; there are no cohort
  `pairwise_relationship_classification.tsv`, `family_hub_roster.tsv`, or
  Stage-4 tables on disk. The only `.res`/truth files present are synthetic
  test fixtures.
- **Stage 3 (the inheritance map) is not implemented.** All Stage-4 results
  above used the placeholder inheritance-map loader; the real loader raises by
  design.
- **No per-group population-genetic statistics (FST, dxy, π, Tajima's D) are
  computed in this repository.** Those, and the cohort-level karyotype burden
  contrast, live in sibling modules outside this repo and were not run here.
- **The external annotation inputs were synthetic.** The `consequence` /
  `sift_class` / `vesm_llr` / `splice_subclass` values in the worked examples
  are fixture values, not engine output; on real data they come from the
  upstream conservation module and must be confirmed against it.
- **The dyad unknown-parent copy is not reconstructed.** Dyad projection emits
  the known-parent copy only; the other copy is left unrepresented (a documented
  refinement, not a bug).

---

## Appendix A — Default parameters

| Parameter | Default | Source |
|---|---|---|
| `theta_first` (first-degree cut) | 0.177 | `scripts/STEP_PED_01_annotate_relationships.py:542` |
| `theta_second` | 0.0884 | `scripts/STEP_PED_01_annotate_relationships.py:544` |
| `theta_third` | 0.0442 | `scripts/STEP_PED_01_annotate_relationships.py:546` |
| `ibs0_po_max` (PO vs FS split) | 0.005 | `scripts/STEP_PED_01_annotate_relationships.py:548` |
| `theta_dup_min` (duplicate) | 0.45 | `scripts/STEP_PED_01_annotate_relationships.py:550` |
| duplicate `IBS0` ceiling | 0.001 | `scripts/STEP_PED_01_annotate_relationships.py:148` |
| `nSites` confidence floor | 1000 | `scripts/STEP_PED_01_annotate_relationships.py:139` |
| KING PO-downgrade threshold | 0.15 | `scripts/STEP_PED_01_annotate_relationships.py:166` |
| forced-parent: high-confidence min offspring | 3 | `scripts/STEP_PED_01_annotate_relationships.py:374` |
| per-chrom `ibs0_po_max` | 0.008 | `scripts/STEP_PED_02_per_chromosome_qc.py:160` |
| per-chrom low-data threshold | 1000 | `scripts/STEP_PED_02_per_chromosome_qc.py:159` |
| disagreement review threshold | 3 | `scripts/STEP_PED_02_per_chromosome_qc.py:161` |
| theta-only adapter `nSites` constant | 100000 | `scripts/adapt_catfish226_json.py:52` |
| Tier-1 LoF consequence set | (5 terms) | `ngsPedigree_v0.5.0/src/hpp/variant_master.py:23` |
| Tier-3 VESM LLR threshold | −7.0 | `ngsPedigree_v0.5.0/src/hpp/variant_master.py:42` |
| Tier-3 SIFT damaging classes | deleterious{,_low_confidence} | `ngsPedigree_v0.5.0/src/hpp/variant_master.py:43` |
| confidence order (worst-of) | high<medium<low<unresolved | `ngsPedigree_v0.5.0/src/hpp/gene_status.py:59` |
| Mendelian warn/fail cut | 1–2 / ≥3 | `ngsPedigree_v0.5.0/src/hpp/mendelian.py:94` |
| default damaging tier | T1 | `ngsPedigree_v0.5.0/src/hpp/variant_master.py:79` |

---

## Appendix B — Statistic provenance and run status

| Statistic / step | In-repo | External | Run on this data |
|---|---|---|---|
| `theta`, `IBS0`, `KING` coefficients | no | ngsRelate/ANGSD (confirm at source) | n/a (consumed, not run here) |
| Edge classifier (relationship class) | yes | — | not on cohort; unit-tested only |
| Hub topology / forced-parent role assignment | yes | — | not on cohort; unit-tested only |
| Per-chromosome consistency + disagreement frac | yes | — | not on cohort; unit-tested only |
| Genotype-likelihood panels (BEAGLE) | no | BEAGLE (confirm at source) | n/a |
| Inheritance map (per-segment hap inherited) | contract only | Stage 3 — **not built** | no (placeholder used) |
| Per-variant consequence / SIFT / VESM / splice | no | conservation module (confirm) | synthetic fixtures only |
| Damaging-tier (T1/T2/T3) decision | yes | — | yes (fixtures) |
| Parental-haplotype construction | yes | — | yes (fixtures) |
| Haplotype projection (Table A) | yes | — | yes (fixtures) |
| Per-gene zygosity-exposure (Table B) | yes | — | yes (fixtures) |
| Mendelian consistency + de-novo flag | yes | — | yes (fixtures) |
| Dyad/triad transmission rollup (Tables C/D) | yes | — | yes (fixtures) |
| Per-group popgen statistics (FST/dxy/π/Tajima) | no | sibling module (out of repo) | no |
| Cohort-level karyotype burden contrast | no | sibling module (out of repo) | no |

---

*End of draft.*
