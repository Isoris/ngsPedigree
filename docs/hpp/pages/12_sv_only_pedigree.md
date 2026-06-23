# Bloc 12 — SV-only pedigree + chromosome inheritance map (no ngsRelate, no BEAGLE, no HPC)

| | |
|---|---|
| **modules** | `vcf_sv.py`, `catalogue_merge.py`, `del_relatedness.py`, `del_inheritance.py` |
| **CLI** | `scripts/08_pedigree_from_sv_vcfs.py` |
| **inputs** | one Delly VCF + one Manta VCF (uncompressed or `.gz`) |
| **purpose** | resolve first-degree dyads / triads and the chromosome inheritance map using only the DEL marker catalogues already produced by Delly and Manta |

## Why this exists

The teacher pays for nothing. ngsRelate per-chromosome takes HPC time
we do not have. But Delly + Manta have already been run for the SV
manuscript. Their two VCFs are a free, segregating biallelic marker
catalogue that can drive the entire pedigree-side of the analysis.

## Pipeline (two commands all the way to ngsTracts)

```
delly.vcf + manta.vcf
   ↓
scripts/08_pedigree_from_sv_vcfs.py
   ↓
  pairwise_relationship_classification.tsv   (Stage 1 schema — drop-in)
  pairwise_coefficients.tsv                  (raw theta / IBS0 / counts)
  candidate_PO_pairs.tsv                     (exclusion-tested PO list)
  inheritance_segments.tsv                   (chromosome inheritance map)
  merged_del_catalogue.json
   ↓
(feed pairwise_relationship_classification.tsv into bloc 08's
 scripts/06_build_polarization_input.py to get the polarization IN JSON,
 then scripts/05_polarize_inversion.py [--mtdna ...] for ngsTracts.)
```

## Stage-by-stage

### 1. VCF SV reader (`vcf_sv.py`)

Stdlib parser for Delly2 and Manta VCFs. Reads `.vcf` or `.vcf.gz`,
filters to `SVTYPE=DEL` PASS records, extracts per-sample GT and
breakpoint coordinates. Detects caller from `##source` header or
filename. No pysam/cyvcf2 dependency.

### 2. Catalogue merger (`catalogue_merge.py`)

Cross-caller DEL unification. Two records merge if:
- same chromosome,
- breakpoint distance ≤ 500 bp on both sides (`--bp-tolerance`),
- reciprocal overlap ≥ 0.5 (`--reciprocal-overlap`).

Matched records get a consensus midpoint position and a unified
genotype per sample (agreement is kept; disagreement defaults to
`./.`). Unmatched records pass through with `callers=("delly",)` or
`("manta",)`.

### 3. KING-robust relatedness on DEL markers (`del_relatedness.py`)

For each pair, computes the standard KING-robust kinship and IBS0:

```
theta = (N_both_het - 2 * N_opposite_hom) / (N_a_het + N_b_het)
IBS0  = N_opposite_hom / N_informative
```

Identical to the SNP form (Manichaikul 2010). The same numbers Stage 1
expects, so the pairwise TSV is **drop-in compatible** with
`STEP_PED_01_annotate_relationships.py` and with the stdlib shadow
classifier (`classify_edge_stdlib`).

### 4. Mendelian exclusion-based parent ID (`del_inheritance.py`)

For each candidate PO pair, test both directions: at every marker
where (parent, offspring) genotypes are observed,

| Parent | Offspring | Verdict |
|---|---|---|
| 0/0 | 1/1 | **excluded** (cannot be parent) |
| 1/1 | 0/0 | **excluded** |
| else | else | compatible |

A single high-confidence opposite-homozygote excludes parentage
(`--min-excluding` to tune).

**Dyad limitation.** A symmetric PO genotype pattern leaves both
directions equally compatible. Without a third anchor (a co-parent,
sex info, or generation order), the direction is reported as
`ambiguous`. This is honest — direction from pure DEL genotypes alone
is not always recoverable, same as in the `po_dyad_only` hub from
bloc 09.

### 5. Triad detection (forced-offspring rule)

If a sample has two candidate PO partners that are themselves
**unrelated** (or 2nd/3rd degree), it is the offspring; the two
partners are the parents. This is the natural complement of Stage 1's
forced-parent rule and runs against the candidate PO set produced in
step 4.

Triads emit two parental traces per chromosome via
`build_inheritance_map_for_triad`, using the co-parent's homozygous
state to disambiguate parent-HET-offspring-HET markers.

### 6. Chromosome inheritance map (`del_inheritance.py`)

For every confirmed PO pair (dyad with known direction, or triad), walk
all DEL markers on each chromosome in coordinate order. At each
informative marker (parent het and offspring resolvable), record which
parental allele was transmitted (`REF` or `DEL`). Run-length encode
into segments and emit:

| Column | Meaning |
|---|---|
| `relationship_id` | `dyad_…` or `triad_…` |
| `relationship_type` | `dyad` / `triad_paternal` / `triad_maternal` |
| `chrom`, `seg_start`, `seg_end` | segment span |
| `transmitted_allele` | `REF` / `DEL` / `ambiguous` |
| `n_informative_markers` | parent-HET resolvable markers in this segment |
| `confidence` | `Gold` (≥8 markers) / `Silver` (3–7) / `Bronze` (1–2) |
| `recomb_event_left` / `_right` | true at segment boundaries |
| `parental_hap_inherited` | `1` (REF-transmitting hap) / `2` (DEL) / `ambiguous` |

A switch in `transmitted_allele` between adjacent segments is a
candidate recombination event — the same signal ngsTracts will refine.

A `Bronze` segment indicates few informative markers (typical of
chromosomes with low DEL density); these are emitted but flagged. A
small noise-smoothing rule (`--min-run-length`, default 2) collapses
isolated single-marker minority runs flanked by the same neighbour to
avoid spurious recombination calls from sequencing noise.

## Demo on the synthetic fixture

```
$ python 08_pedigree_from_sv_vcfs.py \
      --delly delly.vcf --manta manta.vcf --outdir out/

[sv-pedigree] reading delly.vcf
[sv-pedigree]   15 Delly DEL records
[sv-pedigree] reading manta.vcf
[sv-pedigree]   9 Manta DEL records
[sv-pedigree] merged catalogue: 18 markers (both=6, delly-only=9, manta-only=3)
[sv-pedigree] cohort: 4 samples
[sv-pedigree] computed coefficients for 6 pairs
[sv-pedigree] candidate parent-offspring edges: 2
[sv-pedigree] inheritance map: 1 triads + 0 resolved dyads
```

The triad (P_F + P_M → C) is detected automatically; the chromosome
inheritance map shows C inherited a REF→DEL switch on Chr1 from
P_F (a candidate recombination event between markers at ~40 kb and
~50 kb) and a coherent REF-then-DEL transmission on Chr2 from P_M.

## What this bloc does NOT do

- It does **not** call deletions from reads — Delly + Manta already
  did that. This bloc only consumes their VCFs.
- It does **not** orient ambiguous dyads. If both directions are
  Mendelian-compatible, direction is left as `ambiguous`; the user
  can resolve with a sex map (the existing `--sex` plumbing) or with
  mtDNA pre-flight (bloc 07).
- It does **not** call recombination tracts at marker resolution
  inside the inversion intervals — that is still ngsTracts' job. This
  bloc gives ngsTracts a usable phase polarity per (parent, offspring,
  chromosome segment) at DEL-marker resolution.
- It does **not** mask SNPs inside DEL intervals. The merged
  catalogue's `chrom/start/end` per marker is available for downstream
  consumers to do the masking themselves.

## Honesty pass

This pipeline replaces an HPC-scale ngsRelate per-chromosome + BEAGLE
inheritance map with a stdlib script driven by two cheap SV VCFs. The
trade-offs are:

- **Marker density:** SVs are sparser than SNPs (typically ≤ thousands
  of DELs vs millions of SNPs). The KING-robust estimator still
  separates first-degree from unrelated, but with bigger error bars.
  At ~50–100 informative DEL markers per pair the noise on theta is
  still well under the Manichaikul thresholds, so PO/FS/unrelated
  classification is robust.
- **Direction inference:** opposite-hom exclusion is strong evidence,
  but symmetric dyads (no exclusion in either direction) are left
  `ambiguous`. The triad forced-offspring rule (a sample with two
  PO partners who are unrelated to each other) recovers direction
  exactly. This matches the documented `po_dyad_only` limitation
  surfaced by panel C/E (bloc 09).
- **Inheritance-map resolution:** segments are at DEL-marker
  resolution (typically 10 kb–10 Mb), not SNP resolution. Good enough
  to seed ngsTracts, not a replacement for SNP-resolution inheritance
  mapping.

The PhD ships on this.
