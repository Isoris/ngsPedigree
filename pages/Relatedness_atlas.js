/* ===========================================================================
   Relatedness_atlas.js — page logic for Family / Individual Evidence Hub.

   Sibling deliverable to Population_atlas.html / Inversion_atlas.html etc.
   This is the consumer of ngsPedigree's outputs:
     • Stage 1 → pairwise_relationship_classification.tsv
                 family_hub_roster.tsv
     • Stage 2 → per-chromosome QC columns
     • Stage 3 → chromosome inheritance map  (consumed via the Mendelian tab)
   Plus an inversion karyotype TSV from the Inversion Atlas / scrubber.

   Structure:
     §1  Demo data (baked in so the page renders meaningfully on first open)
     §2  State container + small utilities
     §3  Atlas dropdown + theme + tab routing  (mirrors Population_atlas.html)
     §4  Population Browser tree (left column)
     §5  Population Summary (left column bottom)
     §6  Network panel (center, sub-tab #1)
     §7  Karyotypes table (center, sub-tab #2 + inline preview on Network)
     §8  Inversion candidates table (center, sub-tab #3 + inline preview)
     §9  Mendelian segregation tester (center, sub-tab #4) — the headline feature
     §10 Inspector / Stats (right column) — selection-driven
     §11 Loaded files list (right column)
     §12 Initial render boot
   =========================================================================== */

(function() {
'use strict';

/* ╔══════════════════════════════════════════════════════════════════════╗
   §1  DEMO DATA — baked-in synthetic that mirrors the screenshot exactly.
       Replaced wholesale by the data layer once real .res / beagle / inv
       files are wired in via the status pills.
   ╚══════════════════════════════════════════════════════════════════════╝ */

const DEMO = {
  cohort_meta: {
    n_samples: 226,
    species: 'Clarias gariepinus',
    cohort: 'broodstock-226',
  },

  // K=8 ancestry palette — matches the Inversion Atlas defaults.
  ancestry_palette: [
    '#3b82f6',  // blue
    '#f97316',  // orange
    '#10b981',  // emerald
    '#8b5cf6',  // violet
    '#ec4899',  // pink
    '#f59e0b',  // amber
    '#06b6d4',  // cyan
    '#94a3b8',  // slate
  ],

  // The screenshot shows a flat 12-individual demo with 3 families. We keep
  // the same scope so the layout reads identically. The real cohort would
  // wire in 226 samples via Stage 1's family_hub_roster.tsv.
  individuals: [
    'Ind_001','Ind_044','Ind_087','Ind_102','Ind_155','Ind_168','Ind_193',
    'Ind_204','Ind_221','Ind_240','Ind_295','Ind_322'
  ],

  // Family hub topology — matches the network diagram in the screenshot.
  // Hub Ind_001 is the central node with PO connections to 4 confirmed
  // (strong) offspring + 1 possible PO + 1 ambiguous + 1 Mendelian-conflict
  // (which is the kind of diagnostic anomaly Stage 2's per-chrom QC catches).
  families: [
    {
      family_id: 'Family 1',
      hub_individual: 'Ind_001',
      members: ['Ind_001','Ind_044','Ind_087','Ind_102','Ind_155'],
      n: 5,
    },
    {
      family_id: 'Family 2',
      hub_individual: null,
      members: ['Ind_168','Ind_193','Ind_204','Ind_221'],
      n: 4,
    },
    {
      family_id: 'Family 3',
      hub_individual: null,
      members: ['Ind_240','Ind_295','Ind_322'],
      n: 3,
    },
  ],
  ambiguous_clusters: [
    { id: 'cluster_A', members: ['Ind_168','Ind_193'] },
    { id: 'cluster_B', members: ['Ind_322'] },
  ],
  unassigned: [],

  // Sex assignments. In real cohorts this is sometimes provided as a sidecar
  // TSV; if absent everyone is 'unknown' and the sex-aware Compatibility
  // mode just says "no sex data — falling back to non-sex-aware".
  // For the demo we have 12 individuals with ~50/50 split; some unknown.
  sex: {
    'Ind_001':'F', 'Ind_044':'M', 'Ind_087':'F', 'Ind_102':'M', 'Ind_155':'F',
    'Ind_168':'M', 'Ind_193':'F', 'Ind_204':'M', 'Ind_221':'F',
    'Ind_240':'M', 'Ind_295':'F', 'Ind_322':'?',
  },

  // Triads — explicit (parent_a, parent_b, offspring) trios used by the
  // per-inversion family-roster and the Mendelian inheritance summary.
  // 13 triads — enough to make DRIVE_CANDIDATE fire on INV_056 with binomial
  // p < 0.01. In real cohorts these are derived from ngsPedigree Stage 1's
  // family_hub_roster.tsv (rows where role == forced_parent or
  // parent_a/parent_b drive the triad assembly).
  triads: [
    { id: 'T01', parent_a: 'Ind_001', parent_b: 'Ind_044', offspring: 'Ind_087' },
    { id: 'T02', parent_a: 'Ind_001', parent_b: 'Ind_044', offspring: 'Ind_102' },
    { id: 'T03', parent_a: 'Ind_001', parent_b: 'Ind_044', offspring: 'Ind_155' },
    { id: 'T04', parent_a: 'Ind_168', parent_b: 'Ind_193', offspring: 'Ind_204' },
    { id: 'T05', parent_a: 'Ind_168', parent_b: 'Ind_193', offspring: 'Ind_221' },
    { id: 'T06', parent_a: 'Ind_240', parent_b: 'Ind_295', offspring: 'Ind_322' },
    { id: 'T07', parent_a: 'Ind_001', parent_b: 'Ind_240', offspring: 'Ind_087' },
    { id: 'T08', parent_a: 'Ind_044', parent_b: 'Ind_295', offspring: 'Ind_155' },
    { id: 'T09', parent_a: 'Ind_168', parent_b: 'Ind_044', offspring: 'Ind_193' },
    { id: 'T10', parent_a: 'Ind_001', parent_b: 'Ind_193', offspring: 'Ind_322' },
    // 3 extra triads for stronger drive-test power on INV_056
    { id: 'T11', parent_a: 'Ind_001', parent_b: 'Ind_044', offspring: 'Ind_204' },
    { id: 'T12', parent_a: 'Ind_240', parent_b: 'Ind_295', offspring: 'Ind_221' },
    { id: 'T13', parent_a: 'Ind_168', parent_b: 'Ind_295', offspring: 'Ind_087' },
  ],

  // Chromosomes (28 autosomes — typical for catfish at ~28-29).
  chromosomes: (() => {
    const c = [];
    for (let i = 1; i <= 28; i++) {
      c.push('Chr' + String(i).padStart(2, '0'));
    }
    return c;
  })(),

  // PO edges shown in the network panel for hub Ind_001.
  // class: strong_po | possible_po | ambiguous | mendelian_conflict
  network_edges: [
    { a: 'Ind_001', b: 'Ind_044', class: 'strong_po' },
    { a: 'Ind_001', b: 'Ind_087', class: 'strong_po' },
    { a: 'Ind_001', b: 'Ind_102', class: 'strong_po' },
    { a: 'Ind_001', b: 'Ind_155', class: 'strong_po' },
    { a: 'Ind_001', b: 'Ind_168', class: 'possible_po' },
    { a: 'Ind_001', b: 'Ind_193', class: 'ambiguous' },
    // The Mendelian-conflict edge would be shown if any inversion test
    // produced a hard conflict for this dyad. Demo: none currently.
  ],

  // Per-individual ngsRelate-derived stats vs hub Ind_001. Used by Inspector.
  pairwise_stats: {
    'Ind_044': { k0: 0.012, k1: 0.955, k2: 0.033, kinship: 0.489, IBS0: 0.012,
                 PO_distance: 1, relationship_class: 'strong PO-compatible',
                 PASS: 24, WARN: 3, FAIL: 1, TOTAL: 28 },
    'Ind_087': { k0: 0.018, k1: 0.948, k2: 0.034, kinship: 0.483, IBS0: 0.014,
                 PO_distance: 1, relationship_class: 'strong PO-compatible',
                 PASS: 25, WARN: 2, FAIL: 1, TOTAL: 28 },
    'Ind_102': { k0: 0.024, k1: 0.940, k2: 0.036, kinship: 0.476, IBS0: 0.018,
                 PO_distance: 1, relationship_class: 'strong PO-compatible',
                 PASS: 23, WARN: 4, FAIL: 1, TOTAL: 28 },
    'Ind_155': { k0: 0.019, k1: 0.951, k2: 0.030, kinship: 0.485, IBS0: 0.016,
                 PO_distance: 1, relationship_class: 'strong PO-compatible',
                 PASS: 26, WARN: 1, FAIL: 1, TOTAL: 28 },
    'Ind_168': { k0: 0.083, k1: 0.812, k2: 0.105, kinship: 0.380, IBS0: 0.035,
                 PO_distance: 1.4, relationship_class: 'possible PO',
                 PASS: 18, WARN: 7, FAIL: 3, TOTAL: 28 },
    'Ind_193': { k0: 0.21, k1: 0.62, k2: 0.17, kinship: 0.250, IBS0: 0.062,
                 PO_distance: 2.1, relationship_class: 'ambiguous (FS/PO)',
                 PASS: 14, WARN: 9, FAIL: 5, TOTAL: 28 },
    // ↓ demo defaults for non-hub-related individuals so swapping selection
    //   in the Inspector still shows reasonable numbers.
    'Ind_001': { k0: 0,     k1: 0,     k2: 1,     kinship: 0.500, IBS0: 0,
                 PO_distance: 0, relationship_class: 'self',
                 PASS: 28, WARN: 0, FAIL: 0, TOTAL: 28 },
  },

  // Karyotype matrix: each row = an individual, each col = an inversion
  // candidate, value = '0/0' | '0/1' | '1/1' | 'NA'.
  // Some inversions have STD/INV named-allele convention; cell shows that.
  // Built deterministically below.
  inversion_candidates_full: [],   // filled procedurally
  karyotype_matrix: {},            // filled procedurally
  ancestry_q: {},                  // filled procedurally
};

// Generate inversion candidate set (248 to match the mockup's "248" badge).
(function genInversions() {
  const N_INV = 248;
  const chroms = DEMO.chromosomes;
  let inv_id = 1;
  for (let i = 0; i < N_INV; i++) {
    const chrom = chroms[i % chroms.length];
    // Distribute lengths between 0.4 Mb and 6 Mb plausibly.
    const len_mb = Math.round((0.4 + (i * 17 % 56) / 10) * 10) / 10;
    const start_mb = Math.round((((i * 31 + 13) % 100) + 5) * 10) / 10;
    const end_mb = Math.round((start_mb + len_mb) * 10) / 10;
    const freq = Math.round((0.04 + (i * 7 % 41) / 100) * 100) / 100;
    // Status: ~85% pass, 11% warn, 4% fail, mirroring the mockup's bars.
    const status_roll = (i * 41 + 17) % 100;
    let status = 'pass', notes = '';
    if (status_roll < 85)      { status = 'pass';
                                  notes = ['Strong PO-supported','High frequency','Stable across chrs',''][i%4]; }
    else if (status_roll < 96) { status = 'warn';
                                  notes = 'Mendelian warns in chr'+chrom.slice(-2); }
    else                       { status = 'fail';
                                  notes = 'Mendelian failures'; }
    if (freq < 0.10) notes = 'Low frequency';
    DEMO.inversion_candidates_full.push({
      candidate: 'INV_' + String(inv_id++).padStart(3, '0'),
      chromosome: chrom,
      start_mb, end_mb, length_mb: len_mb,
      frequency: freq,
      status,
      notes,
    });
  }
  // Override the first four to exactly match the mockup screenshot's table.
  Object.assign(DEMO.inversion_candidates_full[0], {
    candidate: 'INV_001', chromosome: 'Chr28',
    start_mb: 15.1, end_mb: 18.0, length_mb: 2.9, frequency: 0.23,
    status: 'pass', notes: 'Strong PO-supported',
  });
  Object.assign(DEMO.inversion_candidates_full[1], {
    candidate: 'INV_002', chromosome: 'Chr03',
    start_mb: 22.4, end_mb: 25.9, length_mb: 3.5, frequency: 0.41,
    status: 'warn', notes: 'Mendelian warns in chr17',
  });
  Object.assign(DEMO.inversion_candidates_full[2], {
    candidate: 'INV_003', chromosome: 'Chr17',
    start_mb: 7.2, end_mb: 9.8, length_mb: 2.6, frequency: 0.12,
    status: 'fail', notes: 'Mendelian failures',
  });
  Object.assign(DEMO.inversion_candidates_full[3], {
    candidate: 'INV_004', chromosome: 'Chr10',
    start_mb: 12.7, end_mb: 14.3, length_mb: 1.6, frequency: 0.08,
    status: 'pass', notes: 'Low frequency',
  });
})();

// Generate the karyotype matrix (genotypes at every candidate, for each
// individual). Deterministic per (sample, candidate).
(function genKaryotypes() {
  const k_options = ['0/0','0/1','1/1','NA'];
  // Distribution roughly: 35% 0/0, 40% 0/1, 20% 1/1, 5% NA.
  function pick(seed) {
    const r = Math.abs(Math.sin(seed * 12.9898)) * 43758.5453 % 1;
    if (r < 0.05) return 'NA';
    if (r < 0.40) return '0/0';
    if (r < 0.80) return '0/1';
    return '1/1';
  }
  for (const ind of DEMO.individuals) {
    const seed_base = hashStr(ind);
    DEMO.karyotype_matrix[ind] = {};
    DEMO.inversion_candidates_full.forEach((inv, j) => {
      DEMO.karyotype_matrix[ind][inv.candidate] = pick(seed_base * 0.001 + j);
    });
  }
  // ── Override karyotypes for INV_001 specifically. This is the showcase
  //    inversion in the user's example: across triads we want to see
  //    diagnostic 0/0 × 1/1 → 0/1 patterns (the "real result" the user
  //    wants to feel confident about), with 9 PASS, 1 WARN, 0 FAIL.
  // Karyotypes for INV_001 (Chr28 15.1-18.0 Mb):
  DEMO.karyotype_matrix['Ind_001']['INV_001'] = '0/0';   // Parent A type (Family 1)
  DEMO.karyotype_matrix['Ind_044']['INV_001'] = '1/1';   // Parent B type (Family 1)
  DEMO.karyotype_matrix['Ind_087']['INV_001'] = '0/1';   // Offspring (Family 1) — diagnostic ✓
  DEMO.karyotype_matrix['Ind_102']['INV_001'] = '0/1';   // Offspring (Family 1) — diagnostic ✓
  DEMO.karyotype_matrix['Ind_155']['INV_001'] = '0/1';   // Offspring (Family 1) — diagnostic ✓
  DEMO.karyotype_matrix['Ind_168']['INV_001'] = '0/0';   // Family 2 parent A
  DEMO.karyotype_matrix['Ind_193']['INV_001'] = '0/1';   // Family 2 parent B (het)
  DEMO.karyotype_matrix['Ind_204']['INV_001'] = '0/1';   // Family 2 offspring (informative)
  DEMO.karyotype_matrix['Ind_221']['INV_001'] = '0/1';   // Family 2 offspring (informative)
  DEMO.karyotype_matrix['Ind_240']['INV_001'] = '1/1';   // Family 3 parent
  DEMO.karyotype_matrix['Ind_295']['INV_001'] = '0/1';   // Family 3 parent (het)
  DEMO.karyotype_matrix['Ind_322']['INV_001'] = '0/1';   // Family 3 offspring + cross-family
  // T05 (Ind_168 × Ind_193 → Ind_221): 0/0 × 0/1 → expects {0/0, 0/1};
  //   observed 0/1 → PASS but borderline → mark as WARN via low-confidence flag.
  // T07 (cross-family Ind_001 × Ind_240 → Ind_087): 0/0 × 1/1 → must be 0/1;
  //   Ind_087 = 0/1 → PASS (diagnostic).
  // T09 (Ind_168 × Ind_044 → Ind_193): 0/0 × 1/1 → must be 0/1;
  //   Ind_193 = 0/1 → PASS (diagnostic).
  // T10 (Ind_001 × Ind_193 → Ind_322): 0/0 × 0/1 → expects {0/0, 0/1};
  //   Ind_322 = 0/1 → PASS.

  // Rig a "fail" inversion (INV_003 on Chr17) to demonstrate Conflict tier:
  // engineer one impossible inheritance pattern.
  DEMO.karyotype_matrix['Ind_001']['INV_003'] = '0/0';
  DEMO.karyotype_matrix['Ind_044']['INV_003'] = '0/0';
  DEMO.karyotype_matrix['Ind_087']['INV_003'] = '1/1';   // Hard FAIL: 0/0 × 0/0 → 1/1 impossible
  DEMO.karyotype_matrix['Ind_102']['INV_003'] = '0/0';
  DEMO.karyotype_matrix['Ind_155']['INV_003'] = '0/0';
})();

// Karyotype quality flags. 'low_confidence' karyotypes get downgraded
// from PASS to WARN even when consistent with Mendelian expectations.
// Demonstrates the "boundary uncertainty" / "ambiguous band" WARN case.
DEMO.karyotype_quality = {};
(function genQuality() {
  for (const ind of DEMO.individuals) {
    DEMO.karyotype_quality[ind] = {};
    DEMO.inversion_candidates_full.forEach(inv => {
      DEMO.karyotype_quality[ind][inv.candidate] = 'high';
    });
  }
  // Rig one low-confidence call on INV_001 for triad T05's offspring,
  // to produce the "WARN family" entry in the user's example (1/10 WARN).
  DEMO.karyotype_quality['Ind_221']['INV_001'] = 'low';
})();

// ──────────────────────────────────────────────────────────────────────
// Genome-wide trio QC. For each triad, this would normally come from a
// pipeline that scans every SNP across the genome and counts Mendelian
// errors (parent_a homozygous-ref + parent_b homozygous-ref + offspring
// homozygous-alt = 1 error). The error rate is the per-locus probability
// of a hard inheritance violation. Real pipelines compute this from the
// .res file's IBS0 column or from a vcftools --mendelian-pass run.
//
// The threshold below 0.5% means the trio is genome-wide consistent;
// above ~1% means the family annotation is suspect (sample swap, wrong
// parent, contamination, etc.). A flagged trio's failures at any
// individual inversion are likely just downstream consequences of the
// global problem, not local biology — so we DROP these trios from the
// per-inversion denominator before reporting "9 of 10 PASS".
//
// Demo: T08 is engineered as a family-suspect trio. Its failures across
// inversions should be attributed to family-annotation error, not real
// inversion-level conflicts.
// ──────────────────────────────────────────────────────────────────────
DEMO.trio_qc = {
  // sample_a-sample_b PO support is from the .res file; here we mock it.
  // genome_wide_mendelian_error: fraction of SNPs with hard violations
  // ancestry_midpoint_match: |Q_offspring - mean(Q_parents)|; lower = better
  'T01': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.003, anc_dist: 0.05, valid: true },
  'T02': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.004, anc_dist: 0.06, valid: true },
  'T03': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.002, anc_dist: 0.04, valid: true },
  'T04': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.005, anc_dist: 0.07, valid: true },
  'T05': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.006, anc_dist: 0.08, valid: true },
  'T06': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.003, anc_dist: 0.05, valid: true },
  'T07': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.004, anc_dist: 0.06, valid: true },
  'T08': { po_a: 'weak',   po_b: 'strong', gw_mend_error: 0.024, anc_dist: 0.31, valid: false },
                                          // ↑ family-suspect trio: high genome-wide error
                                          //   plus weak PO on parent A; failures at any
                                          //   inversion are downstream of THIS, not real.
  'T09': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.003, anc_dist: 0.05, valid: true },
  'T10': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.004, anc_dist: 0.06, valid: true },
  'T11': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.003, anc_dist: 0.05, valid: true },
  'T12': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.004, anc_dist: 0.06, valid: true },
  'T13': { po_a: 'strong', po_b: 'strong', gw_mend_error: 0.005, anc_dist: 0.07, valid: true },
};

// ──────────────────────────────────────────────────────────────────────
// Drive-candidate scaffolding. To showcase the TRANSMISSION_SKEW and
// DRIVE_CANDIDATE categories, we engineer one inversion (INV_056 on
// Chr28) where heterozygous parents preferentially transmit the INV
// allele to offspring across multiple families — i.e. simulated
// transmission distortion (which in real data could indicate meiotic
// drive, gametic competition, or genotype-dependent viability — all
// indistinguishable without controlled crosses).
//
// Mechanism: override the karyotype matrix so that for INV_056, every
// het × het cross produces het or hom-INV offspring (no hom-STD), and
// every het × hom-STD cross still gives a het offspring (the expected
// 50% — but we tip the gametic representation in favor of INV).
// ──────────────────────────────────────────────────────────────────────
(function rigDriveCandidate() {
  // INV_056 on Chr28. Override karyotypes so multiple heterozygous parents
  // over-transmit the INV (1) allele. We need ≥8 het-parent gametes total
  // and ≥3 concordant families for the DRIVE_CANDIDATE category.
  // Setup: heterozygous parents in 5 families all transmit "1" to offspring.
  //
  // T01 (Ind_001 × Ind_044 → Ind_087): 0/1 × 0/0 → expected 50:50; observed 0/1 (parent 1 transmitted "1")
  // T02 (Ind_001 × Ind_044 → Ind_102): same parents; observed 0/1 (parent 1 transmitted "1")
  // T03 (Ind_001 × Ind_044 → Ind_155): same parents; observed 0/1 again
  // T04 (Ind_168 × Ind_193 → Ind_204): 0/1 × 0/0 → observed 0/1
  // T05 (Ind_168 × Ind_193 → Ind_221): same; observed 0/1
  // T06 (Ind_240 × Ind_295 → Ind_322): 0/1 × 0/0 → observed 0/1
  // T07 (Ind_001 × Ind_240 → Ind_087): both parents 0/1 → 1/1 (informative
  //                                    only when offspring is hom)
  // T09 (Ind_168 × Ind_044 → Ind_193): 0/1 × 0/0 → observed 0/1
  // That's 7 families with INV transmitted (concordant), giving 7+ "1"
  // gametes. We add a couple of "0" gametes via T10 to make it not perfectly
  // 8/0 (which would be too clean — real drive is partial).
  DEMO.karyotype_matrix['Ind_001']['INV_056'] = '0/1';
  DEMO.karyotype_matrix['Ind_044']['INV_056'] = '0/0';
  DEMO.karyotype_matrix['Ind_087']['INV_056'] = '0/1';   // parent het transmitted "1"
  DEMO.karyotype_matrix['Ind_102']['INV_056'] = '0/1';   // parent het transmitted "1"
  DEMO.karyotype_matrix['Ind_155']['INV_056'] = '0/1';   // parent het transmitted "1"
  DEMO.karyotype_matrix['Ind_168']['INV_056'] = '0/1';
  DEMO.karyotype_matrix['Ind_193']['INV_056'] = '0/0';
  DEMO.karyotype_matrix['Ind_204']['INV_056'] = '0/1';   // transmitted "1"
  DEMO.karyotype_matrix['Ind_221']['INV_056'] = '0/1';   // transmitted "1"
  DEMO.karyotype_matrix['Ind_240']['INV_056'] = '0/1';
  DEMO.karyotype_matrix['Ind_295']['INV_056'] = '0/0';
  DEMO.karyotype_matrix['Ind_322']['INV_056'] = '0/1';   // transmitted "1"
  // T07: Ind_001 × Ind_240 → Ind_087  ; both het, offspring het → ambiguous
  // T09: Ind_168 × Ind_044 → Ind_193  ; 0/1 × 0/0 → Ind_193 = 0/0 (parent
  //   transmitted "0"; gives one "0" transmission to balance the picture
  //   slightly while keeping skew strong)
  DEMO.karyotype_matrix['Ind_193']['INV_056'] = '0/0';   // Ind_193 in T09 is now 0/0
  // Wait — Ind_193 is also in T04 as parent. Need to keep that = 0/0 for
  // T04 to be 0/1 × 0/0; offspring 0/1 means Ind_168 transmitted "1".
  // Since Ind_193 is offspring in T09 AND parent in T04, the same value
  // works: 0/0. (consistent).
  // T10: Ind_001 × Ind_193 → Ind_322. Ind_001 = 0/1, Ind_193 = 0/0.
  //   Ind_322 = 0/1 → Ind_001 transmitted "1". Already covered above.
  // Net result: many families with parent het 0/1 → offspring 0/1 (allele 1
  //   transmitted). Only T09 produces a "0" transmission. Heavily skewed.
})();

// Engineer INV_028 as NEEDS_CROSSES — too few diagnostic / informative
// families. Solution: most karyotype calls are NA so the inversion can't
// be tested.
(function rigNeedsCrosses() {
  // Set most individuals to NA at INV_028 (too few callable karyotypes
  // to evaluate any triad). Only Ind_001 is called; the rest unknown.
  DEMO.individuals.forEach(p => {
    if (p !== 'Ind_001') {
      DEMO.karyotype_matrix[p]['INV_028'] = 'NA';
    }
  });
  DEMO.karyotype_matrix['Ind_001']['INV_028'] = '0/0';
})();



// Generate ancestry-Q vectors (K=8 columns) for each individual.
// Each row sums to 1.0. Used by the Karyotypes ancestry stripe.
(function genAncestry() {
  const K = DEMO.ancestry_palette.length;
  for (const ind of DEMO.individuals) {
    const seed_base = hashStr(ind);
    const raw = [];
    let sum = 0;
    for (let k = 0; k < K; k++) {
      const r = Math.abs(Math.sin((seed_base + k) * 78.233)) * 43758.5 % 1;
      const v = Math.pow(r, 1.5);
      raw.push(v); sum += v;
    }
    // One dominant component
    const dominant = (seed_base | 0) % K;
    raw[dominant] += sum * 1.4;
    sum = raw.reduce((a, b) => a + b, 0);
    DEMO.ancestry_q[ind] = raw.map(v => v / sum);
  }
})();

// Hash a string to a small integer (deterministic, ≥ 0).
function hashStr(s) {
  let h = 5381;
  for (let i = 0; i < s.length; i++) h = ((h << 5) + h + s.charCodeAt(i)) | 0;
  return Math.abs(h);
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §2  STATE CONTAINER + small utilities
   ╚══════════════════════════════════════════════════════════════════════╝ */

const state = window._RA_state = {
  // What's selected in the population tree.
  selected_individual: 'Ind_001',
  selected_family: 'Family 1',
  selected_chromosome: null,             // null = "all chromosomes"
  // Pair currently shown in the Inspector / Mendelian-check panel.
  inspector_pair: { a: 'Ind_001', b: 'Ind_044' },
  // Inversion candidates table pagination.
  inv_page_inline: 1,
  inv_page_full:   1,
  inv_per_page:    4,                    // matches mockup's 4-row table
  // Mendelian-tester state (selections, results).
  mend: {
    mode: 'dyad',
    parent1: 'Ind_001',
    parent2: 'Ind_044',
    offspring: 'Ind_087',
    inv_subset: 'all',
    alpha: '0.05',
    last_result: null,
  },
  // Loaded files (mock — pretend the four files already loaded).
  loaded_files: [
    { kind: 'res',    path: '/data/project/population.res',     loaded: true },
    { kind: 'beagle', path: '/data/beagle/beagle.gz',            loaded: true },
    { kind: 'prune',  path: '/data/prune.in',                    loaded: true },
    { kind: 'inv',    path: '/data/inversion_candidates.tsv',    loaded: true },
  ],
  loaded_at: 'May 20, 2026  10:32:18',
};

// Expose for debugging / external scripts.
window.RA_state = state;

const $  = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));
const fmt = n => Number.isFinite(n) ? n.toFixed(3) : '—';

function el(tag, attrs = {}, ...children) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined) continue;
    if (k === 'class')        e.className = v;
    else if (k === 'html')    e.innerHTML = v;
    else if (k === 'text')    e.textContent = v;
    else if (k.startsWith('on') && typeof v === 'function')
                              e.addEventListener(k.slice(2), v);
    else if (k === 'style' && typeof v === 'object')
                              Object.assign(e.style, v);
    else                      e.setAttribute(k, v);
  }
  for (const c of children) {
    if (c == null) continue;
    if (typeof c === 'string')      e.appendChild(document.createTextNode(c));
    else                            e.appendChild(c);
  }
  return e;
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §3  Header chrome — atlas dropdown, theme toggle, tab routing
        (lifted directly from Population_atlas.html so the chrome behaves
         identically across the whole atlas suite.)
   ╚══════════════════════════════════════════════════════════════════════╝ */

// Sub-tab activation (Network / Karyotypes / Inversions / Mendelian)
$$('#subTabBar button[data-subtab]').forEach(btn => {
  btn.addEventListener('click', () => {
    const target = btn.getAttribute('data-subtab');
    $$('#subTabBar button[data-subtab]').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    $$('[data-subtab-page]').forEach(p => p.classList.remove('active'));
    const pg = $('#subpage-' + target);
    if (pg) pg.classList.add('active');
    try { localStorage.setItem('relatedness_atlas_v0.1.activeSubTab', target); } catch (_) {}
  });
});
try {
  const saved = localStorage.getItem('relatedness_atlas_v0.1.activeSubTab');
  if (saved) {
    const btn = document.querySelector(`#subTabBar button[data-subtab="${saved}"]`);
    if (btn) btn.click();
  }
} catch (_) {}

// Hash routing — `Relatedness_atlas.html#mendelian` lands on the Mendelian tab.
function _routeFromHash() {
  const h = (location.hash || '').replace(/^#/, '').trim();
  if (!h) return false;
  const btn = document.querySelector(`#subTabBar button[data-subtab="${h}"]`);
  if (btn) { btn.click(); return true; }
  return false;
}
window.addEventListener('hashchange', _routeFromHash);
_routeFromHash();

// Theme toggle — three-way cycle (dark → light → academic → dark).
const themeBtn = $('#themeToggleBtn');
function applyTheme(theme) {
  if (theme === 'dark') {
    document.documentElement.removeAttribute('data-theme');
    themeBtn.textContent = '☀ light';
  } else {
    document.documentElement.setAttribute('data-theme', theme);
    if (theme === 'light')    themeBtn.textContent = '📜 academic';
    if (theme === 'academic') themeBtn.textContent = '🌙 dark';
  }
  try { localStorage.setItem('relatedness_atlas_v0.1.theme', theme); } catch (_) {}
}
try {
  const saved = localStorage.getItem('relatedness_atlas_v0.1.theme');
  if (saved && ['dark','light','academic'].includes(saved)) {
    applyTheme(saved);
  } else {
    // Default to light mode since the mockup is light-themed.
    applyTheme('light');
  }
} catch (_) { applyTheme('light'); }
themeBtn.addEventListener('click', () => {
  const cur = document.documentElement.getAttribute('data-theme') || 'dark';
  const next = cur === 'dark' ? 'light' : (cur === 'light' ? 'academic' : 'dark');
  applyTheme(next);
});

// Section collapse/expand toggles (used in left + right columns).
$$('[data-toggle-section]').forEach(btn => {
  btn.addEventListener('click', (e) => {
    e.stopPropagation();
    const section = btn.closest('.col-section');
    if (!section) return;
    section.classList.toggle('collapsed');
    btn.textContent = section.classList.contains('collapsed') ? '⌃' : '⌄';
  });
});

// Pill click → would open a file picker in production. Demo: log only.
$$('.pill').forEach(p => {
  p.addEventListener('click', () => {
    const kind = p.getAttribute('data-pill');
    console.log('[RA] pill click:', kind, '— would open file picker for', kind);
  });
});

// Reload button (right column footer) — placeholder for re-run.
$('#reloadBtn').addEventListener('click', () => {
  state.loaded_at = new Date().toLocaleString();
  $('#loadedAt').textContent = 'Loaded: ' + state.loaded_at;
});


/* ╔══════════════════════════════════════════════════════════════════════╗
   §4  Population Browser tree (left column)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderPopTree() {
  const root = $('#popTree');
  root.innerHTML = '';

  const totalN = DEMO.cohort_meta.n_samples;

  // Population (root)
  const popNode = el('div', { class: 'tree-node', 'data-tree-key': 'population' },
    el('span', { class: 'tree-toggle expanded', text: '▶' }),
    el('span', { class: 'tree-icon population', text: '⚇' }),
    el('span', { class: 'tree-label', text: 'Population' }),
    el('span', { class: 'tree-count', text: 'n = ' + DEMO.individuals.length })
  );
  popNode.addEventListener('click', e => toggleTree(popNode, popKids));
  root.appendChild(popNode);

  const popKids = el('div', { class: 'tree-children' });
  root.appendChild(popKids);

  // Each family
  DEMO.families.forEach((fam, fi) => {
    const isFam1 = fam.family_id === 'Family 1';
    const famNode = el('div', { class: 'tree-node' },
      el('span', { class: 'tree-toggle' + (isFam1 ? ' expanded' : ''),
                   text: '▶' }),
      el('span', { class: 'tree-icon family', text: '⚇' }),
      el('span', { class: 'tree-label',
                   text: fam.family_id + (fam.hub_individual ? ' (hub)' : '') }),
      el('span', { class: 'tree-count', text: 'n = ' + fam.n })
    );
    popKids.appendChild(famNode);

    const famKids = el('div', { class: 'tree-children' + (isFam1 ? '' : ' collapsed') });
    popKids.appendChild(famKids);
    famNode.addEventListener('click', e => {
      e.stopPropagation();
      toggleTree(famNode, famKids);
      state.selected_family = fam.family_id;
      $('#bcFamily').textContent = fam.family_id;
    });

    // Hub (if any)
    if (fam.hub_individual) {
      const hubNode = el('div', { class: 'tree-node selected' },
        el('span', { class: 'tree-toggle expanded', text: '▶' }),
        el('span', { class: 'tree-icon hub', text: '◉' }),
        el('span', { class: 'tree-label', text: 'Hub: ' + fam.hub_individual })
      );
      famKids.appendChild(hubNode);
      const hubKids = el('div', { class: 'tree-children' });
      famKids.appendChild(hubKids);

      hubNode.addEventListener('click', e => {
        e.stopPropagation();
        toggleTree(hubNode, hubKids);
        selectIndividual(fam.hub_individual);
      });

      // Chromosome leaves under hub
      DEMO.chromosomes.slice(0, 4).forEach(ch => {
        const chNode = el('div', { class: 'tree-node dimmed' },
          el('span', { class: 'tree-toggle empty', text: '▶' }),
          el('span', { class: 'tree-icon chromosome', text: '⌬' }),
          el('span', { class: 'tree-label', text: ch })
        );
        chNode.addEventListener('click', e => {
          e.stopPropagation();
          state.selected_chromosome = ch;
          $('#chrFilter').value = ch;
          renderInversionTables();
        });
        hubKids.appendChild(chNode);
      });
      // Add a "..." indicator for the rest
      hubKids.appendChild(el('div', {
        class: 'tree-node dimmed',
        style: { fontStyle: 'italic', fontSize: '10px', color: 'var(--ink-dimmer)' }
      },
        el('span', { class: 'tree-toggle empty', text: '▶' }),
        el('span', { class: 'tree-icon chromosome', text: '⌬' }),
        el('span', { class: 'tree-label', text: 'Chr28' })
      ));
    }

    // Other members
    fam.members.filter(m => m !== fam.hub_individual).forEach(member => {
      const m_n = el('div', { class: 'tree-node dimmed' },
        el('span', { class: 'tree-toggle', text: '▶' }),
        el('span', { class: 'tree-icon individual', text: '○' }),
        el('span', { class: 'tree-label', text: member })
      );
      m_n.addEventListener('click', e => {
        e.stopPropagation();
        selectIndividual(member);
      });
      famKids.appendChild(m_n);
    });
  });

  // Ambiguous clusters
  const ambN = DEMO.ambiguous_clusters.reduce((s, c) => s + c.members.length, 0);
  const ambNode = el('div', { class: 'tree-node' },
    el('span', { class: 'tree-toggle', text: '▶' }),
    el('span', { class: 'tree-icon family', text: '◇' }),
    el('span', { class: 'tree-label', text: 'Ambiguous clusters' }),
    el('span', { class: 'tree-count', text: 'n = ' + ambN })
  );
  popKids.appendChild(ambNode);
  const ambKids = el('div', { class: 'tree-children collapsed' });
  popKids.appendChild(ambKids);
  ambNode.addEventListener('click', e => {
    e.stopPropagation(); toggleTree(ambNode, ambKids);
  });

  // Unassigned
  const unNode = el('div', { class: 'tree-node' },
    el('span', { class: 'tree-toggle empty', text: '▶' }),
    el('span', { class: 'tree-icon individual', text: '?' }),
    el('span', { class: 'tree-label', text: 'Unassigned' }),
    el('span', { class: 'tree-count', text: 'n = 1' })
  );
  popKids.appendChild(unNode);
}

function toggleTree(node, kidsContainer) {
  if (!kidsContainer) return;
  const collapsed = kidsContainer.classList.contains('collapsed');
  kidsContainer.classList.toggle('collapsed', !collapsed);
  const toggle = node.querySelector('.tree-toggle');
  if (toggle && !toggle.classList.contains('empty')) {
    toggle.classList.toggle('expanded', collapsed);
  }
}

function selectIndividual(ind) {
  state.selected_individual = ind;
  state.inspector_pair = { a: 'Ind_001', b: ind === 'Ind_001' ? 'Ind_044' : ind };
  // Re-render selected highlight
  $$('.tree-node').forEach(n => {
    const lbl = n.querySelector('.tree-label');
    if (!lbl) return;
    const matches = lbl.textContent.includes(ind);
    n.classList.toggle('selected', matches && !n.classList.contains('dimmed'));
  });
  // Update breadcrumb
  $('#bcLeaf').textContent = ind === DEMO.families[0].hub_individual
    ? 'Hub: ' + ind : ind;
  renderInspector();
  renderNetwork();
  renderInversionTables();   // status-pill counts may change with selection
}

// Tree search
$('#treeSearch').addEventListener('input', e => {
  const q = e.target.value.trim().toLowerCase();
  $$('.tree-node').forEach(n => {
    const lbl = n.querySelector('.tree-label');
    if (!lbl) return;
    if (!q) { n.style.display = ''; return; }
    n.style.display = lbl.textContent.toLowerCase().includes(q) ? '' : 'none';
  });
});


/* ╔══════════════════════════════════════════════════════════════════════╗
   §5  Population Summary (left column bottom)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderPopSummary() {
  const root = $('#popSummary');
  root.innerHTML = '';
  const rows = [
    ['Individuals',          12, null],
    ['PO edges (strong)',     6, '<span class="legend-line" style="border-top-width:2.5px;"></span>'],
    ['PO edges (possible)',   4, '<span class="legend-line possible" style="border-top-width:1.5px;"></span>'],
    ['Ambiguous edges',       3, '<span class="legend-line ambig"></span>'],
    ['Mendelian warnings',    5, '<span class="legend-icon warn">⚠</span>'],
    ['Mendelian failures',    2, '<span class="legend-icon fail">⚐</span>'],
  ];
  for (const [lbl, val, ico] of rows) {
    const lblWrap = el('div', { class: 'pop-summary-label' });
    if (ico) lblWrap.innerHTML = ico + ' ' + lbl;
    else     lblWrap.textContent = lbl;
    root.appendChild(lblWrap);
    root.appendChild(el('div', { class: 'pop-summary-value', text: String(val) }));
  }
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §6  Network panel (sub-tab #1)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderNetwork() {
  const svg = $('#networkSvg');
  while (svg.firstChild) svg.removeChild(svg.firstChild);

  const W = 800, H = 360;
  const cx = W / 2, cy = H / 2;
  const center = state.selected_individual || 'Ind_001';

  // Find edges incident to the focal node, plus edges in its hub.
  const incident = DEMO.network_edges.filter(e =>
    e.a === center || e.b === center);
  const neighbors = incident.map(e => e.a === center ? e.b : e.a);
  // Position center node + neighbors on a radial layout.
  const positions = { [center]: { x: cx, y: cy } };
  neighbors.forEach((n, i) => {
    const angle = (-Math.PI / 2) + (i * 2 * Math.PI / neighbors.length);
    const r = 130;
    positions[n] = { x: cx + r * Math.cos(angle), y: cy + r * Math.sin(angle) };
  });

  // SVG namespace helper
  const NS = 'http://www.w3.org/2000/svg';
  function svgEl(tag, attrs) {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) e.setAttribute(k, v);
    return e;
  }

  // Render edges
  for (const e of incident) {
    const pa = positions[e.a];
    const pb = positions[e.b];
    if (!pa || !pb) continue;
    const cls = 'net-edge ' + (
      e.class === 'strong_po'   ? 'strong-po' :
      e.class === 'possible_po' ? 'possible-po' :
      e.class === 'ambiguous'   ? 'ambig' : 'conflict'
    );
    const line = svgEl('line', {
      x1: pa.x, y1: pa.y, x2: pb.x, y2: pb.y,
      class: cls,
    });
    svg.appendChild(line);
  }

  // Render nodes
  for (const [name, p] of Object.entries(positions)) {
    const isCenter = name === center;
    const grp = svgEl('g', {
      class: 'net-node' + (isCenter ? ' selected' : ''),
      transform: `translate(${p.x},${p.y})`,
      'data-individual': name,
    });
    grp.appendChild(svgEl('circle', {
      r: isCenter ? 22 : 16,
    }));
    grp.appendChild(svgEl('text', {
      y: isCenter ? -32 : -22,
      class: 'net-label',
    })).textContent = name;
    grp.addEventListener('click', () => {
      selectIndividual(name);
    });
    svg.appendChild(grp);
  }
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §7  Karyotypes table (inline preview + full sub-page view)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderKaryotypeTable(slot, opts = { rows: null, columns: null }) {
  const target = typeof slot === 'string' ? $(slot) : slot;
  if (!target) return;
  target.innerHTML = '';

  const rows = opts.rows || DEMO.individuals.slice(0, 5);
  const cols = opts.columns || (() => {
    // Mockup shows Chr01..Chr06, then ..., then Chr17, Chr18, ..., Chr28.
    const sample = ['Chr01','Chr02','Chr03','Chr04','Chr05','Chr06',
                    null, 'Chr17','Chr18', null, 'Chr28'];
    return sample;
  })();

  // Pick the first inversion candidate per shown chromosome to populate a
  // single column per chromosome (genotype call for that representative
  // candidate).
  const chrCandidate = {};
  for (const inv of DEMO.inversion_candidates_full) {
    if (!chrCandidate[inv.chromosome]) chrCandidate[inv.chromosome] = inv.candidate;
  }

  const table = el('table', { class: 'data-table' });
  const thead = el('thead');
  const tr = el('tr');
  tr.appendChild(el('th', { text: 'Ancestry' }));
  cols.forEach(c => tr.appendChild(el('th',
    { text: c === null ? '…' : c, style: { textAlign: 'center' } })));
  thead.appendChild(tr);
  table.appendChild(thead);

  const tbody = el('tbody');
  for (const ind of rows) {
    const r = el('tr');
    // Ancestry cell: sample id + horizontal stacked-bar
    r.appendChild((function(){
      const td = el('td', null);
      const ancWrap = el('div', { style: { display: 'flex', alignItems: 'center', gap: '8px' } },
        el('span', { class: 'sample-id', text: ind }),
        renderAncestryStripe(ind),
      );
      td.appendChild(ancWrap);
      return td;
    })());
    cols.forEach(c => {
      const td = el('td', { style: { textAlign: 'center' } });
      if (c === null) {
        td.textContent = '…';
        td.style.color = 'var(--ink-dimmer)';
      } else {
        const cand = chrCandidate[c];
        const k = cand ? DEMO.karyotype_matrix[ind][cand] : null;
        let cls = 'kt-cell';
        if (k === '0/0') cls += ' kt-00';
        else if (k === '0/1') cls += ' kt-01';
        else if (k === '1/1') cls += ' kt-11';
        else cls += ' kt-na';
        // For Ind_155 row — use STD/INV named-allele convention to mirror
        // the mockup screenshot's last row.
        let display = k || 'NA';
        if (ind === 'Ind_155' && k && k !== 'NA') {
          if (k === '0/0') display = 'STD/STD';
          else if (k === '0/1') display = 'STD/INV';
          else if (k === '1/1') display = 'INV/INV';
        }
        td.appendChild(el('span', { class: cls, text: display }));
      }
      r.appendChild(td);
    });
    tbody.appendChild(r);
  }
  table.appendChild(tbody);
  target.appendChild(table);

  // Legend
  const legend = el('div', { class: 'kt-legend' },
    legendChip('kt-00', '0/0 or STD/STD'),
    legendChip('kt-01', '0/1 or STD/INV'),
    legendChip('kt-11', '1/1 or INV/INV'),
    legendChip('kt-na', 'No call / missing'),
  );
  target.appendChild(legend);
}

function legendChip(cls, label) {
  return el('div', { class: 'item' },
    el('span', { class: 'kt-cell ' + cls,
                 style: { width: '18px', height: '12px',
                          minWidth: '0', padding: '0' } }),
    el('span', { text: label })
  );
}

function renderAncestryStripe(ind) {
  const wrap = el('div', { class: 'anc-stripe' });
  const q = DEMO.ancestry_q[ind] || [];
  q.forEach((p, k) => {
    if (p < 0.005) return;
    const seg = el('span', {
      style: {
        width: (p * 100).toFixed(2) + '%',
        background: DEMO.ancestry_palette[k % DEMO.ancestry_palette.length],
      },
      title: 'K' + (k + 1) + ': ' + (p * 100).toFixed(1) + '%',
    });
    wrap.appendChild(seg);
  });
  return wrap;
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §8  Inversion candidates table (inline + full sub-page)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function getFilteredInversions() {
  let inv = DEMO.inversion_candidates_full;
  if (state.selected_chromosome && state.selected_chromosome !== 'all') {
    inv = inv.filter(i => i.chromosome === state.selected_chromosome);
  }
  return inv;
}

/* ─── Family-inheritance scoring for an inversion ────────────────────────
   Iterates every triad in DEMO.triads, evaluates Mendelian compatibility
   against the candidate's karyotypes, and aggregates into PASS/WARN/FAIL
   counts plus a support tier. This is the core of "is this inversion
   real?" — not one trio's evidence, but the whole-cohort consistency.

   A family is:
     PASS  — observed offspring karyotype ∈ Mendelian-expected set
             (all parents called, all karyotypes high-confidence)
     WARN  — observed in expected set but at least one karyotype call is
             low-confidence (boundary uncertainty)
     FAIL  — observed offspring karyotype ∉ Mendelian-expected set
             (hard impossibility under standard genetics)
     not informative — at least one parent or offspring is NA

   "Diagnostic" families are the subset where parent karyotypes uniquely
   determine the offspring karyotype (the most powerful kind of evidence):
     0/0 × 0/0 → 0/0 only
     1/1 × 1/1 → 1/1 only
     0/0 × 1/1 → 0/1 only

   Support tiers:
     strong:   PASS ≥ 90% AND FAIL = 0 AND informative ≥ 5
     moderate: PASS 70-90%, FAIL ≤ 1, informative ≥ 3
     weak:     PASS < 70% OR informative < 3 (no clear conclusion)
     conflict: FAIL ≥ 2 OR (FAIL = 1 AND informative ≤ 5)
*/

// Mendelian rule table: key = sorted "p1|p2", value = array of allowed
// offspring karyotypes. Sorted to canonicalize 0/0×0/1 == 0/1×0/0.
const MENDEL_RULES = (() => {
  const rules = {};
  function add(p1, p2, allowed) {
    const k = [p1, p2].sort().join('|');
    rules[k] = allowed;
  }
  add('0/0','0/0', ['0/0']);
  add('1/1','1/1', ['1/1']);
  add('0/0','1/1', ['0/1']);
  add('0/0','0/1', ['0/0','0/1']);
  add('1/1','0/1', ['0/1','1/1']);
  add('0/1','0/1', ['0/0','0/1','1/1']);
  return rules;
})();

function isDiagnosticParents(p1, p2) {
  // Diagnostic = parent karyotypes uniquely determine offspring
  // 0/0 × 0/0, 1/1 × 1/1, 0/0 × 1/1
  const k = [p1, p2].sort().join('|');
  return ['0/0|0/0', '1/1|1/1', '0/0|1/1'].includes(k);
}

function expectedOffspring(p1, p2) {
  if (!p1 || !p2 || p1 === 'NA' || p2 === 'NA') return null;
  return MENDEL_RULES[[p1, p2].sort().join('|')] || null;
}

function classifyFamilyForInversion(triad, invId) {
  const p1 = (DEMO.karyotype_matrix[triad.parent_a] || {})[invId];
  const p2 = (DEMO.karyotype_matrix[triad.parent_b] || {})[invId];
  const o  = (DEMO.karyotype_matrix[triad.offspring]  || {})[invId];
  // Quality flags
  const q1 = (DEMO.karyotype_quality[triad.parent_a] || {})[invId] || 'high';
  const q2 = (DEMO.karyotype_quality[triad.parent_b] || {})[invId] || 'high';
  const qo = (DEMO.karyotype_quality[triad.offspring]  || {})[invId] || 'high';
  const anyLowConf = (q1 === 'low' || q2 === 'low' || qo === 'low');

  // Family-level validity from genome-wide trio QC (Stage 1 of the
  // four-stage hierarchy). If gw_mend_error > 1% OR a PO edge is weak,
  // this family's local failures are downstream of an annotation problem,
  // not real biology — so we mark it but compute the local result anyway
  // for transparency.
  const qc = DEMO.trio_qc[triad.id] || { valid: true, gw_mend_error: 0, anc_dist: 0 };
  const family_valid = qc.valid !== false;

  if (!p1 || !p2 || !o || p1 === 'NA' || p2 === 'NA' || o === 'NA') {
    return {
      family_id: triad.id, parent_a: triad.parent_a, parent_b: triad.parent_b,
      offspring: triad.offspring,
      p1_kt: p1 || 'NA', p2_kt: p2 || 'NA', o_kt: o || 'NA',
      expected: null, status: 'not_informative',
      diagnostic: false, family_valid, qc,
    };
  }
  const expected = expectedOffspring(p1, p2);
  const diagnostic = isDiagnosticParents(p1, p2);
  const compatible = expected && expected.includes(o);
  let status;
  if (!compatible) {
    // Hard fail. But if the family is suspect, downgrade to family_warn
    // to prevent the inversion from being judged on bad annotation.
    status = family_valid ? 'fail' : 'family_warn';
  }
  else if (anyLowConf) status = 'warn';
  else status = 'pass';
  return {
    family_id: triad.id,
    parent_a: triad.parent_a, parent_b: triad.parent_b,
    offspring: triad.offspring,
    p1_kt: p1, p2_kt: p2, o_kt: o,
    expected, status, diagnostic, anyLowConf,
    family_valid, qc,
  };
}

function scoreInversion(invId) {
  // Four-stage hierarchical scoring (per the user's framework):
  //   Stage 1: family validation (genome-wide trio QC)
  //   Stage 2: local Mendelian compatibility per VALID family
  //   Stage 3: aggregate "X of Y families PASS" verdict
  //   Stage 4: transmission test (heterozygous parent gametes), only if
  //            stages 1-3 are clean. Repeated direction across families
  //            is the drive signal.
  //
  // The category emitted is one of:
  //   'PASS'              — Mendelian-compatible, no skew
  //   'WARN_CALL'         — Mendelian-compatible but low-confidence calls
  //   'WARN_FAMILY'       — failures concentrated in suspect trios
  //   'LOCAL_CONFLICT'    — valid families fail this specific inversion
  //   'TRANSMISSION_SKEW' — valid families compatible, but binomial test
  //                         on het-parent transmissions rejects p=0.5
  //   'DRIVE_CANDIDATE'   — TRANSMISSION_SKEW that repeats across ≥3
  //                         independent valid families with same direction
  //   'NEEDS_CROSSES'     — too few informative families; can't decide
  //                         (also: too few het parents to test drive)
  const families = DEMO.triads.map(t => classifyFamilyForInversion(t, invId));
  const informative_families = families.filter(f => f.status !== 'not_informative');
  const valid_inf = informative_families.filter(f => f.family_valid !== false);
  const suspect_inf = informative_families.filter(f => f.family_valid === false);
  const diagnostic_families = valid_inf.filter(f => f.diagnostic);

  // Stage 2/3: aggregate over VALID informative families only.
  const n_pass = valid_inf.filter(f => f.status === 'pass').length;
  const n_warn = valid_inf.filter(f => f.status === 'warn').length;
  const n_fail = valid_inf.filter(f => f.status === 'fail').length;
  const n_inf  = valid_inf.length;
  const pass_frac = n_inf > 0 ? n_pass / n_inf : 0;

  // Suspect-trio failures (downstream of family annotation errors)
  const n_suspect_fail = suspect_inf.filter(f =>
    f.status === 'fail' || f.status === 'family_warn').length;

  // Stage 4: transmission test. Only meaningful when n_inf is reasonable
  // and we have heterozygous parents.
  // For each VALID family with at least one het parent, count how the
  // offspring's allele was transmitted. parent het 0/1, offspring 0/0
  // → 0-allele transmitted; offspring 1/1 → 1-allele transmitted; offspring
  // 0/1 with hom-other partner → unambiguously the OTHER allele was
  // transmitted from the het parent. With both parents het, the
  // transmission is ambiguous (skip).
  let n_transmissions_0 = 0, n_transmissions_1 = 0;
  let n_het_parents_used = 0;
  // Track per-family transmission direction (for "repeated direction" rule)
  const family_directions = []; // 'over_1', 'over_0', 'balanced', 'na'
  for (const f of valid_inf) {
    let fam_t0 = 0, fam_t1 = 0;
    function tally(parent_kt, partner_kt, offspring_kt) {
      // Returns 0 (parent transmitted ref), 1 (parent transmitted alt),
      //         null (uninformative).
      if (parent_kt !== '0/1') return null;
      // Hom partner: offspring genotype reveals what the het parent gave.
      if (partner_kt === '0/0') {
        if (offspring_kt === '0/0') return 0;
        if (offspring_kt === '0/1') return 1;
        return null; // 1/1 impossible in clean cross — already failed Mendelian
      }
      if (partner_kt === '1/1') {
        if (offspring_kt === '0/1') return 0;
        if (offspring_kt === '1/1') return 1;
        return null;
      }
      // Both het: ambiguous unless offspring is hom (rare unambiguous case)
      if (partner_kt === '0/1') {
        if (offspring_kt === '0/0') return 0;
        if (offspring_kt === '1/1') return 1;
        return null;
      }
      return null;
    }
    const t_a = tally(f.p1_kt, f.p2_kt, f.o_kt);
    const t_b = tally(f.p2_kt, f.p1_kt, f.o_kt);
    if (t_a === 0) fam_t0++;
    if (t_a === 1) fam_t1++;
    if (t_b === 0) fam_t0++;
    if (t_b === 1) fam_t1++;
    n_transmissions_0 += fam_t0;
    n_transmissions_1 += fam_t1;
    const fam_total = fam_t0 + fam_t1;
    if (fam_total > 0) {
      n_het_parents_used += fam_total;
      const dir = fam_t1 > fam_t0 ? 'over_1' : (fam_t0 > fam_t1 ? 'over_0' : 'balanced');
      family_directions.push({ family_id: f.family_id, t0: fam_t0, t1: fam_t1, dir });
    } else {
      family_directions.push({ family_id: f.family_id, t0: 0, t1: 0, dir: 'na' });
    }
  }
  const n_total_t = n_transmissions_0 + n_transmissions_1;
  // Two-sided binomial-exact test against 50:50.
  let trans_p = NaN, trans_skew_dir = 'none';
  if (n_total_t >= 4) {
    const k = Math.min(n_transmissions_0, n_transmissions_1);
    trans_p = binomialPValueTwoSided(k, n_total_t, 0.5);
    if (n_transmissions_1 > n_transmissions_0)      trans_skew_dir = 'over_1';
    else if (n_transmissions_0 > n_transmissions_1) trans_skew_dir = 'over_0';
    else                                              trans_skew_dir = 'balanced';
  }
  // Repeatability across families: count how many families had the
  // same direction as the global skew (excluding 'balanced'/'na').
  const concordant_families = family_directions.filter(
    d => d.dir === trans_skew_dir).length;
  const family_dirs_observed = family_directions.filter(
    d => d.dir === 'over_0' || d.dir === 'over_1').length;

  // Decide category, in order
  let category, verdict;
  if (n_inf < 3) {
    category = 'NEEDS_CROSSES';
    verdict = 'Insufficient informative families — design crosses to validate';
  }
  else if (pass_frac < 0.70) {
    // Lots of failures in valid families
    if (n_suspect_fail > n_fail) {
      category = 'WARN_FAMILY';
      verdict = 'Failures concentrated in suspect trios — review annotation';
    } else {
      category = 'LOCAL_CONFLICT';
      verdict = 'Valid families show repeated Mendelian inconsistencies at this locus';
    }
  }
  else if (n_fail >= 2) {
    category = 'LOCAL_CONFLICT';
    verdict = 'Valid families show repeated Mendelian inconsistencies at this locus';
  }
  else if (Number.isFinite(trans_p) && trans_p < 0.02
           && n_total_t >= 8 && concordant_families >= 5) {
    category = 'DRIVE_CANDIDATE';
    verdict = `Transmission distortion repeated across ${concordant_families} `
            + `independent valid families (binomial p = ${trans_p.toExponential(2)})`;
  }
  else if (Number.isFinite(trans_p) && trans_p < 0.05 && n_total_t >= 6) {
    category = 'TRANSMISSION_SKEW';
    verdict = `Heterozygous parents transmit one allele at ${(Math.max(n_transmissions_0,n_transmissions_1)/n_total_t*100).toFixed(0)}% `
            + `(binomial p = ${trans_p.toFixed(3)}); needs more families before drive call`;
  }
  else if (n_warn > 0 && n_pass + n_warn === n_inf) {
    category = 'WARN_CALL';
    verdict = 'Mendelian-compatible across all valid families, but some calls are low-confidence';
  }
  else {
    category = 'PASS';
    verdict = 'Mendelian-compatible inversion candidate';
  }

  // Visual tier (drives row styling). Mapping from category → tier:
  //   PASS, WARN_CALL                  → 'strong' or 'moderate' depending on PASS frac
  //   DRIVE_CANDIDATE                  → 'drive'   (diamond aura)
  //   TRANSMISSION_SKEW                → 'skew'    (deeper amber)
  //   LOCAL_CONFLICT                   → 'conflict'
  //   WARN_FAMILY                      → 'warn_family'
  //   NEEDS_CROSSES                    → 'needs_crosses' (black ?)
  let tier;
  if      (category === 'DRIVE_CANDIDATE')   tier = 'drive';
  else if (category === 'TRANSMISSION_SKEW') tier = 'skew';
  else if (category === 'LOCAL_CONFLICT')    tier = 'conflict';
  else if (category === 'WARN_FAMILY')       tier = 'warn_family';
  else if (category === 'NEEDS_CROSSES')     tier = 'needs_crosses';
  else if (category === 'WARN_CALL')         tier = 'moderate';
  else if (n_inf >= 5 && pass_frac >= 0.90 && n_fail === 0)
                                              tier = 'strong';   // ruby aura
  else if (n_inf >= 3 && pass_frac >= 0.70)  tier = 'moderate';
  else                                        tier = 'weak';

  // Aura intensity (0..1) — bigger glow for stronger evidence.
  // Strong tier: scales with pass_frac and n_inf.
  // Drive tier: scales with skew magnitude and n_total_t.
  let aura_intensity = 0;
  if (tier === 'strong') {
    // 5/0/0 perfect → 1.0; 9/1/0 → ~0.85; 7/2/0 → ~0.65
    aura_intensity = Math.min(1.0, 0.5 + 0.5 * pass_frac
                                       + 0.05 * Math.min(n_inf - 5, 5));
    // Penalize for warns (small)
    aura_intensity -= 0.05 * n_warn;
    aura_intensity = Math.max(0.4, Math.min(1.0, aura_intensity));
  } else if (tier === 'drive') {
    const skew_mag = n_total_t > 0
      ? Math.abs(n_transmissions_1 - n_transmissions_0) / n_total_t : 0;
    aura_intensity = Math.min(1.0, 0.5 + 0.5 * skew_mag);
  }

  return {
    inv_id: invId,
    n_informative: n_inf,
    n_diagnostic: diagnostic_families.length,
    n_pass, n_warn, n_fail,
    n_suspect_inf: suspect_inf.length,
    n_suspect_fail,
    pass_frac, tier,
    category, verdict,
    families,
    valid_families: valid_inf,
    suspect_families: suspect_inf,
    // transmission-test outputs
    n_transmissions_0, n_transmissions_1, n_total_t,
    trans_p, trans_skew_dir,
    family_directions,
    concordant_families,
    aura_intensity,
  };
}

function renderInversionTables() {
  renderInversionTable('#invTableSlotInline', '#invPaginationInline', 'inv_page_inline');
  renderInversionTable('#invTableSlotFull',   '#invPaginationFull',   'inv_page_full');
  $('#pillInvValue').textContent = String(DEMO.inversion_candidates_full.length);
}

function renderInversionTable(slot, paginationSlot, pageStateKey) {
  const target = $(slot);
  if (!target) return;
  target.innerHTML = '';

  const all = getFilteredInversions();
  const total = all.length;
  const page = state[pageStateKey];
  const per = state.inv_per_page;
  const start = (page - 1) * per;
  const visible = all.slice(start, start + per);

  const table = el('table', { class: 'data-table' });
  const thead = el('thead');
  // Updated header — includes inheritance-summary columns
  // (Informative / PASS / WARN / FAIL / Support tier).
  const headers = ['Candidate','Chrom','Start','End','Length','Freq',
                   'Inform.','PASS','WARN','FAIL','Support','Notes'];
  const tr = el('tr');
  headers.forEach(h => tr.appendChild(el('th', { text: h })));
  thead.appendChild(tr);
  table.appendChild(thead);

  const tbody = el('tbody');
  for (const inv of visible) {
    // Score inversion against family roster
    const score = scoreInversion(inv.candidate);
    // Map tier name (snake_case) to CSS class (kebab-case)
    const tier_cls = score.tier.replace(/_/g, '-');
    const auraStyle = score.aura_intensity > 0
      ? '--aura-intensity:' + score.aura_intensity.toFixed(2) : '';
    const r = el('tr', { class: 'clickable tier-' + tier_cls,
                          'data-inv-id': inv.candidate,
                          style: auraStyle });
    r.appendChild(el('td', { class: 'sample-id', text: inv.candidate }));
    r.appendChild(el('td', { text: inv.chromosome }));
    r.appendChild(el('td', { class: 'num', text: inv.start_mb.toFixed(1) + ' Mb' }));
    r.appendChild(el('td', { class: 'num', text: inv.end_mb.toFixed(1) + ' Mb' }));
    r.appendChild(el('td', { class: 'num', text: inv.length_mb.toFixed(1) + ' Mb' }));
    r.appendChild(el('td', { class: 'num', text: inv.frequency.toFixed(2) }));
    r.appendChild(el('td', { class: 'num',
                              text: score.n_informative + (score.n_diagnostic > 0
                                ? ' (' + score.n_diagnostic + ' diag.)' : '') }));
    r.appendChild(el('td', { class: 'num',
                              text: String(score.n_pass),
                              style: { color: score.n_pass > 0 ? 'var(--good)' : 'var(--ink-dim)' } }));
    r.appendChild(el('td', { class: 'num',
                              text: String(score.n_warn),
                              style: { color: score.n_warn > 0 ? 'var(--warn)' : 'var(--ink-dim)' } }));
    r.appendChild(el('td', { class: 'num',
                              text: String(score.n_fail),
                              style: { color: score.n_fail > 0 ? 'var(--bad)' : 'var(--ink-dim)' } }));
    {
      const td = el('td');
      const pill_cls = score.tier.replace(/_/g, '-');
      td.appendChild(el('span', {
        class: 'tier-pill ' + pill_cls,
        text: score.category.replace(/_/g, ' ').toLowerCase(),
        title: score.verdict,
      }));
      r.appendChild(td);
    }
    r.appendChild(el('td', { text: inv.notes,
                              style: { color: 'var(--ink-dim)', fontSize: '9.5px' } }));
    // Make the row clickable: clicking expands the family-inheritance roster
    r.addEventListener('click', () => toggleInversionExpand(inv, score, r, tbody));
    tbody.appendChild(r);
  }
  table.appendChild(tbody);
  target.appendChild(table);

  if (paginationSlot) {
    const p = $(paginationSlot);
    if (p) renderPagination(p, total, page, per, n => {
      state[pageStateKey] = n;
      renderInversionTable(slot, paginationSlot, pageStateKey);
    });
  }
}

// Toggle the expanded family-inheritance row directly under the clicked row.
function toggleInversionExpand(inv, score, parentRow, tbody) {
  // If already expanded, collapse
  const next = parentRow.nextElementSibling;
  if (next && next.classList.contains('inv-expand-row')
        && next.dataset.parent === inv.candidate) {
    next.remove();
    parentRow.classList.remove('expanded');
    return;
  }
  // Remove any other expanded rows in this tbody
  Array.from(tbody.querySelectorAll('.inv-expand-row')).forEach(n => n.remove());
  Array.from(tbody.querySelectorAll('tr.expanded')).forEach(n => n.classList.remove('expanded'));
  // Add new expand row
  const expandTr = el('tr', { class: 'inv-expand-row',
                               'data-parent': inv.candidate });
  const expandTd = el('td', { colspan: 12 });
  expandTd.appendChild(buildInversionExpandPanel(inv, score));
  expandTr.appendChild(expandTd);
  parentRow.classList.add('expanded');
  parentRow.parentNode.insertBefore(expandTr, parentRow.nextSibling);
}

function buildInversionExpandPanel(inv, score) {
  const panel = el('div', { class: 'inv-expand-panel' });

  // Title block
  const titleText = inv.candidate + ' — '
                    + inv.chromosome + '  '
                    + inv.start_mb.toFixed(1) + '–' + inv.end_mb.toFixed(1) + ' Mb · '
                    + 'length ' + inv.length_mb.toFixed(1) + ' Mb · '
                    + 'freq ' + inv.frequency.toFixed(2);
  panel.appendChild(el('div', { class: 'ie-title', text: titleText }));
  panel.appendChild(el('div', { class: 'ie-subtitle',
    text: 'Four-stage scoring: (1) family validity from genome-wide trio QC, '
        + '(2) local Mendelian compatibility per valid family, '
        + '(3) aggregate "X of Y" verdict, '
        + '(4) transmission test for distortion/drive.' }));

  // ─── Stats grid ───
  const stats = el('div', { class: 'ie-stats' });
  stats.appendChild(statCell('Valid informative', score.n_informative));
  stats.appendChild(statCell('Diagnostic ★', score.n_diagnostic));
  stats.appendChild(statCell('PASS', score.n_pass, 'pass'));
  stats.appendChild(statCell('WARN', score.n_warn, 'warn'));
  stats.appendChild(statCell('FAIL', score.n_fail, 'fail'));
  stats.appendChild(statCell('Pass fraction', (score.pass_frac * 100).toFixed(0) + '%'));
  if (score.n_suspect_inf > 0) {
    stats.appendChild(statCell('Suspect trios ⚐', score.n_suspect_inf, 'warn'));
  }
  panel.appendChild(stats);

  // ─── Family-level roster ───
  panel.appendChild(el('div', {
    style: { fontSize: '11px', fontWeight: '600', color: 'var(--ink)',
             marginTop: '14px', marginBottom: '6px' },
    text: 'Stage 1+2: per-family roster (validity from genome-wide trio QC, '
        + 'local Mendelian inheritance at this inversion)',
  }));
  const rosterTable = el('table', { class: 'ie-roster-table' });
  const thead = el('thead');
  const tr = el('tr');
  ['Family','Trio QC','Parent A','kt','Parent B','kt','Offspring','kt','Expected','Verdict']
    .forEach(h => tr.appendChild(el('th', { text: h })));
  thead.appendChild(tr); rosterTable.appendChild(thead);
  const tbody = el('tbody');
  score.families.forEach(f => {
    const r = el('tr');
    if (f.family_valid === false) r.style.opacity = '0.55';
    r.appendChild(el('td', { text: f.family_id + (f.diagnostic ? ' ★' : '') }));
    {
      // Trio QC pill
      const td = el('td');
      const qc = f.qc || {};
      if (f.family_valid === false) {
        td.appendChild(el('span', {
          class: 'tier-pill warn-family', text: 'SUSPECT',
          title: `gw_mend_error=${(qc.gw_mend_error*100).toFixed(1)}% · `
               + `anc_dist=${qc.anc_dist?.toFixed(2) || '?'} · `
               + `PO support: ${qc.po_a || '?'}/${qc.po_b || '?'}`,
        }));
      } else {
        td.appendChild(el('span', {
          class: 'tier-pill', text: 'OK', title: 'genome-wide trio QC passed',
          style: { background: 'rgba(95,212,154,0.14)',
                    color: 'var(--good)', borderColor: 'rgba(95,212,154,0.40)' },
        }));
      }
      r.appendChild(td);
    }
    r.appendChild(el('td', { html: f.parent_a + ' ' + sexBadgeHtml(f.parent_a) }));
    r.appendChild(el('td', null,
      el('span', { class: 'karyo-cell ' + ktClass(f.p1_kt), text: f.p1_kt })));
    r.appendChild(el('td', { html: f.parent_b + ' ' + sexBadgeHtml(f.parent_b) }));
    r.appendChild(el('td', null,
      el('span', { class: 'karyo-cell ' + ktClass(f.p2_kt), text: f.p2_kt })));
    r.appendChild(el('td', { html: f.offspring + ' ' + sexBadgeHtml(f.offspring) }));
    r.appendChild(el('td', null,
      el('span', { class: 'karyo-cell ' + ktClass(f.o_kt), text: f.o_kt })));
    r.appendChild(el('td', { class: 'expected',
      text: f.expected ? f.expected.join(' / ') : '—' }));
    const tdV = el('td', { class: 'verdict-cell' });
    if (f.status === 'not_informative') {
      tdV.appendChild(el('span', { class: 'tier-pill weak', text: 'NA' }));
    } else if (f.status === 'family_warn') {
      tdV.appendChild(el('span', { class: 'tier-pill warn-family', text: 'FAM⚐' }));
    } else {
      tdV.appendChild(el('span', {
        class: 'status-pill-cell ' + f.status,
        text: f.status.toUpperCase(),
      }));
    }
    r.appendChild(tdV);
    tbody.appendChild(r);
  });
  rosterTable.appendChild(tbody);
  panel.appendChild(rosterTable);

  // Footnotes
  panel.appendChild(el('div', {
    style: { fontSize: '9px', color: 'var(--ink-dim)',
             marginTop: '6px', fontStyle: 'italic', lineHeight: '1.5' },
    text: '★ = diagnostic family (parent karyotypes uniquely determine the offspring '
        + 'outcome — most informative evidence). SUSPECT = trio failed genome-wide '
        + 'QC (likely sample-swap / annotation error); its local result is shown '
        + 'transparently but excluded from the X-of-Y aggregate.',
  }));

  // ─── Stage 4: Transmission test (only meaningful when relevant) ───
  if (score.n_total_t >= 4) {
    panel.appendChild(el('div', {
      style: { fontSize: '11px', fontWeight: '600', color: 'var(--ink)',
               marginTop: '14px', marginBottom: '6px' },
      text: 'Stage 4: transmission test — heterozygous-parent gametes '
          + '(binomial against 50:50)',
    }));
    const tStats = el('div', { class: 'ie-stats' });
    tStats.appendChild(statCell('Het-parent transmissions', score.n_total_t));
    tStats.appendChild(statCell('Reference allele transmitted', score.n_transmissions_0));
    tStats.appendChild(statCell('Inverted allele transmitted', score.n_transmissions_1));
    if (Number.isFinite(score.trans_p)) {
      const skew_pct = score.n_total_t > 0
        ? (Math.max(score.n_transmissions_0, score.n_transmissions_1) /
           score.n_total_t * 100).toFixed(0)
        : '50';
      const sev = score.trans_p < 0.01 ? 'fail'
                : score.trans_p < 0.05 ? 'warn' : 'pass';
      tStats.appendChild(statCell('Skew (favored allele %)', skew_pct + '%', sev));
      tStats.appendChild(statCell('Binomial p',
        score.trans_p < 0.001
          ? score.trans_p.toExponential(2)
          : score.trans_p.toFixed(3),
        sev));
      tStats.appendChild(statCell('Concordant families',
        score.concordant_families + '/' + score.family_directions.filter(
          d => d.dir === 'over_0' || d.dir === 'over_1').length));
    } else {
      tStats.appendChild(statCell('Status', 'not enough informative gametes'));
    }
    panel.appendChild(tStats);
  } else if (score.n_informative >= 3) {
    panel.appendChild(el('div', {
      style: { fontSize: '10px', color: 'var(--ink-dim)',
               marginTop: '14px', marginBottom: '6px',
               padding: '8px 12px', background: 'var(--panel)',
               border: '1px dashed var(--rule)', borderRadius: '3px' },
      text: 'Stage 4 (transmission test): not run — fewer than 4 informative '
          + 'het-parent gametes available. Drive cannot be tested at this '
          + 'inversion without more heterozygous parents.',
    }));
  }

  // ─── Conclusion paragraph ───
  const conclLabel = score.category.replace(/_/g, ' ');
  const concl = el('div', { class: 'ie-conclusion tier-' + score.tier.replace(/_/g, '-') });
  concl.appendChild(el('div', { class: 'verdict', text: conclLabel }));
  // Body text by category
  let body;
  switch (score.category) {
    case 'PASS':
      body = `Mendelian-compatible inversion candidate. ${score.n_pass} of ${score.n_informative} `
           + `valid informative families consistent with expected segregation `
           + `(${score.n_diagnostic} diagnostic), no hard conflicts. `
           + `Retain as high-confidence call.`; break;
    case 'WARN_CALL':
      body = `Inversion is Mendelian-compatible across all ${score.n_informative} valid `
           + `informative families, but ${score.n_warn} family/families had at least `
           + `one low-confidence karyotype call. Consider boundary refinement or `
           + `additional sequencing.`; break;
    case 'WARN_FAMILY':
      body = `Failures concentrated in ${score.n_suspect_inf} suspect trios with elevated `
           + `genome-wide Mendelian errors — likely an annotation error rather than a `
           + `real inversion problem. Re-examine those triads' parent assignments before `
           + `judging this inversion.`; break;
    case 'LOCAL_CONFLICT':
      body = `Valid families show ${score.n_fail} hard Mendelian conflict(s) at this `
           + `inversion specifically (genome-wide QC was OK for these trios). This is `
           + `the pattern of a real local problem: karyotyping artifact, mis-defined `
           + `boundary, or genuine biological complexity. Investigate before retaining.`; break;
    case 'TRANSMISSION_SKEW':
      body = `Valid families show no Mendelian incompatibilities, but heterozygous parents `
           + `transmit one allele at ${(Math.max(score.n_transmissions_0,score.n_transmissions_1)/score.n_total_t*100).toFixed(0)}% `
           + `(binomial p = ${score.trans_p.toFixed(3)}, n = ${score.n_total_t}). `
           + `Could indicate transmission distortion, but with only `
           + `${score.concordant_families} concordant families this is preliminary. `
           + `Add more triads before claiming drive.`; break;
    case 'DRIVE_CANDIDATE':
      body = `Strong, repeated transmission distortion across ${score.concordant_families} `
           + `independent valid families. Heterozygous parents transmit `
           + `${(Math.max(score.n_transmissions_0,score.n_transmissions_1)/score.n_total_t*100).toFixed(0)}% `
           + `of one allele (n = ${score.n_total_t}, binomial p = `
           + `${score.trans_p.toExponential(2)}). This is the signature of transmission `
           + `distortion — possibly meiotic drive, gametic competition, or genotype-`
           + `dependent viability. Controlled crosses with embryo/larval sampling are `
           + `required to distinguish these mechanisms.`; break;
    case 'NEEDS_CROSSES':
      body = `Insufficient informative families (only ${score.n_informative}). Cannot `
           + `judge Mendelian compatibility or transmission distortion from existing `
           + `data alone. Use the "Export cross design" button below to generate a `
           + `recommended experimental-cross matrix that would maximize informativeness.`; break;
    default:
      body = score.verdict;
  }
  concl.appendChild(el('div', null, body));
  panel.appendChild(concl);

  // ─── Action buttons ───
  const actions = el('div', { class: 'ie-actions' });
  actions.appendChild(el('button', {
    class: 'primary', text: '→ Run formal P-value test (Mendelian tab)',
    onclick: () => runMendelianFromInversion(inv.candidate),
  }));
  actions.appendChild(el('button', {
    text: '⤓ Export inheritance roster TSV',
    onclick: () => exportInversionRosterTsv(inv, score),
  }));
  actions.appendChild(el('button', {
    text: '⤓ Export experimental cross design',
    onclick: () => exportCrossDesignTsv(inv, score),
    title: 'Generate a recommended cross matrix to validate or extend this inversion',
  }));
  actions.appendChild(el('button', {
    text: '→ Compatibility planner',
    onclick: () => openCompatibilityForInversion(inv.candidate),
  }));
  panel.appendChild(actions);

  return panel;
}

function statCell(label, value, severity = '') {
  return el('div', { class: 'ie-stat' },
    el('div', { class: 'lbl', text: label }),
    el('div', { class: 'val' + (severity ? ' ' + severity : ''), text: String(value) })
  );
}

function ktClass(k) {
  if (k === '0/0') return 'k00';
  if (k === '0/1') return 'k01';
  if (k === '1/1') return 'k11';
  return 'kna';
}

function sexBadge(ind) {
  // returned as suffix-able HTML to attach next to an individual ID
  const sex = (DEMO.sex || {})[ind] || '?';
  return ' '; // sex pill rendered separately to avoid HTML-escape issues;
              // here we just return a space and let downstream renderers
              // add the badge if they want.
}

function exportInversionRosterTsv(inv, score) {
  const cols = ['family_id','parent_a','parent_a_kt','parent_b','parent_b_kt',
                'offspring','offspring_kt','expected','status','diagnostic'];
  const lines = [
    '# Inversion: ' + inv.candidate,
    '# Chromosome: ' + inv.chromosome + ' ' + inv.start_mb.toFixed(1)
      + '-' + inv.end_mb.toFixed(1) + ' Mb',
    '# Frequency: ' + inv.frequency.toFixed(2),
    '# Tier: ' + score.tier + ' — ' + score.verdict,
    '# Informative: ' + score.n_informative + ' (diagnostic: ' + score.n_diagnostic + ')',
    '# PASS=' + score.n_pass + '  WARN=' + score.n_warn + '  FAIL=' + score.n_fail,
    cols.join('\t'),
  ];
  score.families.forEach(f => {
    lines.push([f.family_id, f.parent_a, f.p1_kt, f.parent_b, f.p2_kt,
                f.offspring, f.o_kt,
                f.expected ? f.expected.join('|') : '',
                f.status,
                f.diagnostic ? 'yes' : 'no'].join('\t'));
  });
  downloadTsv('inversion_inheritance_' + inv.candidate + '.tsv',
              lines.join('\n') + '\n');
}

// ─────────────────────────────────────────────────────────────────────
// Experimental cross-design exporter.
//
// Looks at which karyotype combinations are present vs missing among the
// existing trios for this inversion, then recommends crosses that would
// fill the gap. The most diagnostic crosses are 0/0 × 1/1 (mandatory 0/1
// offspring, perfect Mendelian witness) and 0/1 × 0/0 / 0/1 × 1/1
// (gametic transmission test from het parents — needed for drive call).
//
// Output is a TSV with columns: cross_type, parent_a_karyotype,
// parent_b_karyotype, expected_offspring, n_offspring_recommended,
// purpose, priority. Plus a header block describing the design rationale.
// ─────────────────────────────────────────────────────────────────────
function exportCrossDesignTsv(inv, score) {
  // What parent-karyotype crosses already exist?
  const observedCrosses = new Set();
  score.valid_families.forEach(f => {
    observedCrosses.add([f.p1_kt, f.p2_kt].sort().join('×'));
  });

  // Frequency of each karyotype in the cohort, to know which parents are
  // available. (For demo: scan DEMO.karyotype_matrix.)
  const freq = { '0/0': 0, '0/1': 0, '1/1': 0, 'NA': 0 };
  DEMO.individuals.forEach(ind => {
    const k = (DEMO.karyotype_matrix[ind] || {})[inv.candidate];
    if (k && freq[k] !== undefined) freq[k]++;
  });
  const total_called = freq['0/0'] + freq['0/1'] + freq['1/1'];

  // The four cross types worth recommending, ordered by priority for
  // resolving the unknowns about this inversion.
  const recommendations = [];

  // 1. Most diagnostic: 0/0 × 1/1 → forced 0/1 (single-class, Mendelian witness)
  recommendations.push({
    cross_type: '0/0 × 1/1',
    parent_a_kt: '0/0',
    parent_b_kt: '1/1',
    expected_offspring: '100% 0/1',
    n_offspring_rec: 5,
    purpose: 'Diagnostic Mendelian witness — every offspring must be 0/1; '
           + 'any other genotype is a hard inheritance error',
    priority: 1,
    available_parents: `${freq['0/0']} hom-STD × ${freq['1/1']} hom-INV available`,
    already_done: observedCrosses.has(['0/0','1/1'].sort().join('×'))
                  ? 'YES (existing in cohort)' : 'NO',
  });

  // 2. Drive test (paternal): 0/1 male × 0/0 female → 50% 0/0 / 50% 0/1
  if (score.n_total_t < 12 || score.category === 'TRANSMISSION_SKEW'
                            || score.category === 'NEEDS_CROSSES') {
    recommendations.push({
      cross_type: '0/1 ♂ × 0/0 ♀',
      parent_a_kt: '0/1',
      parent_b_kt: '0/0',
      expected_offspring: '50% 0/0, 50% 0/1',
      n_offspring_rec: 12,
      purpose: 'Test paternal transmission ratio. Significant deviation from '
             + '50:50 implies meiotic drive, gametic competition, or post-'
             + 'fertilization viability selection.',
      priority: 2,
      available_parents: `${freq['0/1']} het × ${freq['0/0']} hom-STD available; `
                       + `pick one of each opposite sex if known`,
      already_done: 'PARTIAL (current trios provide '
                  + score.n_total_t + ' het-parent gametes; '
                  + 'aim for ≥30 to detect modest drive)',
    });
    // 3. Drive test (maternal) — symmetric
    recommendations.push({
      cross_type: '0/1 ♀ × 0/0 ♂',
      parent_a_kt: '0/1',
      parent_b_kt: '0/0',
      expected_offspring: '50% 0/0, 50% 0/1',
      n_offspring_rec: 12,
      purpose: 'Test maternal transmission ratio. Compare with paternal cross '
             + 'above to determine if drive is sex-of-origin specific.',
      priority: 3,
      available_parents: `Reciprocal of above; ensure female is the 0/1 parent`,
      already_done: 'see paternal note',
    });
  }

  // 4. Both-het cross (segregation distortion + viability)
  recommendations.push({
    cross_type: '0/1 × 0/1',
    parent_a_kt: '0/1',
    parent_b_kt: '0/1',
    expected_offspring: '25% 0/0, 50% 0/1, 25% 1/1',
    n_offspring_rec: 16,
    purpose: 'Test full 1:2:1 segregation. Deviation may indicate viability '
           + 'selection on a homozygous class (heterozygote advantage) or '
           + 'segregation distortion if both parents drive same direction.',
    priority: 4,
    available_parents: `${freq['0/1']} het individuals available`,
    already_done: observedCrosses.has(['0/1','0/1'].sort().join('×'))
                  ? 'YES (existing in cohort)' : 'NO',
  });

  // 5. If the inversion frequency is high enough to find rare 1/1 ×
  // 1/1 crosses, recommend that too — important for hom-INV viability.
  if (freq['1/1'] >= 2) {
    recommendations.push({
      cross_type: '1/1 × 1/1',
      parent_a_kt: '1/1',
      parent_b_kt: '1/1',
      expected_offspring: '100% 1/1',
      n_offspring_rec: 5,
      purpose: 'Test fertility / viability of homozygous-INV genotype. If '
             + 'INV/INV pairs produce reduced clutch sizes or skewed survival, '
             + 'this is evidence the inversion carries recessive deleterious '
             + 'alleles.',
      priority: 5,
      available_parents: `${freq['1/1']} hom-INV individuals available`,
      already_done: observedCrosses.has(['1/1','1/1'].sort().join('×'))
                    ? 'YES' : 'NO',
    });
  }

  // Build TSV
  const cols = ['cross_type','parent_a_kt','parent_b_kt','expected_offspring',
                'n_offspring_recommended','purpose','priority',
                'available_parents','already_done'];
  const lines = [
    '# Experimental cross design',
    '# Inversion: ' + inv.candidate + ' (' + inv.chromosome + ' '
      + inv.start_mb.toFixed(1) + '–' + inv.end_mb.toFixed(1) + ' Mb)',
    '# Frequency: ' + inv.frequency.toFixed(2),
    '# Current category: ' + score.category + ' (' + score.tier + ')',
    '# Current evidence: ' + score.n_pass + ' PASS / ' + score.n_warn + ' WARN / '
      + score.n_fail + ' FAIL across ' + score.n_informative + ' valid families',
    '# Het-parent gametes counted: ' + score.n_total_t,
    '# Cohort karyotype availability:',
    '#   0/0: ' + freq['0/0'] + ' individuals',
    '#   0/1: ' + freq['0/1'] + ' individuals',
    '#   1/1: ' + freq['1/1'] + ' individuals',
    '#   NA:  ' + freq['NA']  + ' individuals',
    '#   total called: ' + total_called,
    '#',
    '# Design rationale:',
    '#   This table lists experimental crosses ordered by priority that would',
    '#   either confirm Mendelian inheritance (priority 1: 0/0 × 1/1 forced',
    '#   0/1) or quantify transmission distortion (priority 2-3: 0/1 × 0/0',
    '#   reciprocal crosses to test maternal vs paternal transmission ratios).',
    '#   For DRIVE_CANDIDATE-tier inversions, controlled crosses with',
    '#   embryo/larval sampling are the gold standard for distinguishing',
    '#   meiotic drive from genotype-dependent viability.',
    '#',
    cols.join('\t'),
  ];
  recommendations.forEach(r => {
    lines.push([r.cross_type, r.parent_a_kt, r.parent_b_kt, r.expected_offspring,
                r.n_offspring_rec,
                '"' + r.purpose.replace(/"/g, '""') + '"',
                r.priority, r.available_parents, r.already_done].join('\t'));
  });

  downloadTsv('cross_design_' + inv.candidate + '.tsv',
              lines.join('\n') + '\n');
}

function downloadTsv(filename, content) {
  const blob = new Blob([content], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1500);
}

// Wire-up button: jump to the Mendelian tab pre-populated with this inversion
function runMendelianFromInversion(invId) {
  // Switch to Mendelian sub-tab
  const btn = document.querySelector('#subTabBar button[data-subtab="mendelian"]');
  if (btn) btn.click();
  // Set the inversion-subset filter to "current chromosome only"
  // and select an example family from the triads list.
  const t0 = DEMO.triads[0];
  $('#mendMode').value = 'triad';
  $('#mendMode').dispatchEvent(new Event('change'));
  $('#mendParent1').value = t0.parent_a;
  $('#mendParent2').value = t0.parent_b;
  $('#mendOffspring').value = t0.offspring;
  state.mend.parent1 = t0.parent_a;
  state.mend.parent2 = t0.parent_b;
  state.mend.offspring = t0.offspring;
  // Scroll into view
  const slot = $('#subpage-mendelian');
  if (slot) slot.scrollTop = 0;
  $('#mendRunBtn').click();
}

function openCompatibilityForInversion(invId) {
  // Switch to Compatibility sub-tab and pre-select this inversion
  const btn = document.querySelector('#subTabBar button[data-subtab="compatibility"]');
  if (btn) btn.click();
  $('#compatInvScope').value = 'single';
  $('#compatInvScope').dispatchEvent(new Event('change'));
  $('#compatInvSingle').value = invId;
  state.compat.scope = 'single';
  state.compat.inv_single = invId;
}



function renderPagination(target, total, page, per, onPage) {
  target.innerHTML = '';
  const nPages = Math.max(1, Math.ceil(total / per));
  const start = (page - 1) * per + 1;
  const end = Math.min(page * per, total);
  target.appendChild(el('div', { text: `Showing ${start}–${end} of ${total}` }));
  const pages = el('div', { class: 'pages' });
  pages.appendChild(makePageBtn('‹', () => onPage(Math.max(1, page - 1)), page === 1));
  // Page buttons: 1 ... (n-1) (n) (n+1) ... last  — collapsed for many pages
  const wantPages = new Set();
  wantPages.add(1);
  wantPages.add(2);
  wantPages.add(3);
  wantPages.add(nPages);
  if (page > 3) wantPages.add(page);
  let lastShown = 0;
  Array.from(wantPages).filter(n => n <= nPages).sort((a,b)=>a-b).forEach(n => {
    if (n - lastShown > 1) pages.appendChild(el('button', { text: '…', disabled: true,
        style: { background: 'transparent', border: 'none', cursor: 'default' } }));
    pages.appendChild(makePageBtn(String(n), () => onPage(n), false, n === page));
    lastShown = n;
  });
  pages.appendChild(makePageBtn('›', () => onPage(Math.min(nPages, page + 1)), page === nPages));
  target.appendChild(pages);
}

function makePageBtn(label, onClick, disabled = false, active = false) {
  return el('button', {
    text: label,
    onclick: disabled ? null : onClick,
    class: active ? 'active' : '',
    disabled: disabled ? true : null,
  });
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §9  Mendelian segregation tester (sub-tab #4)
       — the headline analytical feature.

       Methods:
         • DYAD — parent × offspring, no other parent. We test whether the
           offspring's allele draws from the parent's two alleles is
           consistent with the binomial(n, 0.5) expectation. Sample size is
           the number of inversions where the parent is heterozygous (the
           informative subset). Statistic: number of times the offspring
           inherits each allele; binomial-exact P-value, two-sided.

         • TRIAD — parent₁ × parent₂ → offspring. The parental karyotype
           pair determines the expected offspring distribution exactly via
           Mendel's first law:
              0/0 × 0/0 → 1.00 0/0
              0/0 × 0/1 → 0.50 0/0, 0.50 0/1
              0/0 × 1/1 → 1.00 0/1
              0/1 × 0/1 → 0.25 0/0, 0.50 0/1, 0.25 1/1
              0/1 × 1/1 → 0.50 0/1, 0.50 1/1
              1/1 × 1/1 → 1.00 1/1
           For each candidate where both parents are called, we derive the
           expected category and compare to observed. Cohort-level test:
           collect a binary "expected vs observed match" per candidate;
           binomial test against the no-error expectation (1.0). Or under
           a per-class chi-square if at least one parent is heterozygous.

         • ALL_DYADS / ALL_TRIADS — same logic, applied across every dyad
           or triad in the current hub. Output: per-pair P, plus a Fisher
           or Stouffer-like combined statistic.

       This is the in-page UI test only — for production, ngsPedigree's
       Stage 1 + 2 outputs are the source of truth and the manuscript
       methods refer to those tables, not this widget.
   ╚══════════════════════════════════════════════════════════════════════╝ */

// Fill the Mendelian-tab selectors from the current hub members.
function populateMendSelectors() {
  const fam = DEMO.families.find(f => f.family_id === state.selected_family)
              || DEMO.families[0];
  const opts = fam.members.map(m => `<option value="${m}">${m}</option>`).join('');
  $('#mendParent1').innerHTML = opts;
  $('#mendParent2').innerHTML = opts;
  $('#mendOffspring').innerHTML = opts;
  // Reasonable defaults — pick hub as parent1, first non-hub as offspring.
  if (fam.hub_individual) {
    $('#mendParent1').value = fam.hub_individual;
    state.mend.parent1 = fam.hub_individual;
  }
  const others = fam.members.filter(m => m !== fam.hub_individual);
  if (others.length) {
    $('#mendParent2').value = others[0];
    state.mend.parent2 = others[0];
  }
  if (others.length > 1) {
    $('#mendOffspring').value = others[1];
    state.mend.offspring = others[1];
  }
}

// Mode change toggles which selectors are shown.
function updateMendModeUI() {
  const mode = $('#mendMode').value;
  state.mend.mode = mode;
  const showParent1   = (mode === 'dyad' || mode === 'triad');
  const showParent2   = (mode === 'triad');
  const showOffspring = (mode === 'dyad' || mode === 'triad');
  $('#rowParent1').style.display   = showParent1   ? '' : 'none';
  $('#rowParent2').style.display   = showParent2   ? '' : 'none';
  $('#rowOffspring').style.display = showOffspring ? '' : 'none';
  $('#lblParent1').textContent = (mode === 'dyad') ? 'Parent' : 'Parent 1';
}

// Statistical primitives ────────────────────────────────────────────────

// Log of n choose k — uses the lgamma-ish recurrence good for our scale.
function logChoose(n, k) {
  if (k < 0 || k > n) return -Infinity;
  if (k === 0 || k === n) return 0;
  // Use Stirling-ish accumulation
  let s = 0;
  for (let i = 1; i <= k; i++) s += Math.log(n - k + i) - Math.log(i);
  return s;
}

// Two-sided binomial-exact P(X = k | n, p), summing over outcomes at most
// as likely as the observed.
function binomialPValueTwoSided(k, n, p) {
  if (n === 0) return 1.0;
  // P(X = j)
  function pmf(j) {
    return Math.exp(logChoose(n, j) + j * Math.log(p) + (n - j) * Math.log(1 - p));
  }
  const p_obs = pmf(k);
  let pval = 0;
  for (let j = 0; j <= n; j++) {
    const pj = pmf(j);
    if (pj <= p_obs + 1e-12) pval += pj;
  }
  return Math.min(1.0, pval);
}

// Chi-square P-value (Wilson–Hilferty approximation, fine for df ≥ 1 and
// the modest sample sizes we test here). Returns p = P(X² ≥ x | df).
function chiSquarePValue(x, df) {
  if (x <= 0) return 1.0;
  if (df <= 0) return NaN;
  // Wilson–Hilferty
  const t = Math.cbrt(x / df) - (1 - 2 / (9 * df));
  const z = t / Math.sqrt(2 / (9 * df));
  // 1 - Phi(z) for upper tail; use Abramowitz approximation.
  const a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741;
  const a4 = -1.453152027, a5 =  1.061405429, p_ = 0.3275911;
  const sign = z < 0 ? -1 : 1;
  const xz = Math.abs(z) / Math.SQRT2;
  const t_ = 1 / (1 + p_ * xz);
  const y = 1 - (((((a5 * t_ + a4) * t_) + a3) * t_ + a2) * t_ + a1) * t_ * Math.exp(-xz * xz);
  const phi = 0.5 * (1 + sign * y);
  return Math.max(0, Math.min(1, 1 - phi));
}

// For triads: given parent karyotypes, what's the expected offspring
// karyotype probability vector [P(0/0), P(0/1), P(1/1)]?
// Returns null if either parent karyotype is missing or 'NA'.
function expectedOffspringPrior(p1, p2) {
  if (!p1 || !p2 || p1 === 'NA' || p2 === 'NA') return null;
  const a = p1.split('/').map(Number);
  const b = p2.split('/').map(Number);
  // Each parent transmits a uniformly random allele from its two.
  const probs = [0, 0, 0];   // [0/0, 0/1, 1/1]
  for (const a1 of a) for (const b1 of b) {
    const sum = a1 + b1;
    probs[sum] += 0.25;      // 4 combinations, each prob 0.25
  }
  return probs;
}

// Run a dyad test for one (parent, offspring) pair across selected inversions.
// At each candidate where the parent is heterozygous (informative), we look
// at whether the offspring carries the parent's "0" allele or "1" allele.
//   parent het '0/1', offspring '0/0' → offspring inherited '0' (1 success counted as "0")
//   parent het '0/1', offspring '0/1' → offspring inherited either; we
//                                       treat this as half-success for the
//                                       binomial (or, more honestly, we
//                                       skip — we choose to skip since the
//                                       site is uninformative for dyad-only).
// We report:
//   n_informative = parent-het sites where offspring is homozygous
//   n_zero        = offspring 0/0 at those sites
//   n_one         = offspring 1/1 at those sites
//   binomial_p    = P-value for n_zero vs binomial(n_inf, 0.5)
function runDyadTest(parentId, offspringId, candidateList) {
  const pK = DEMO.karyotype_matrix[parentId]    || {};
  const oK = DEMO.karyotype_matrix[offspringId] || {};
  let n_inf = 0, n_zero = 0, n_one = 0;
  let n_consistent = 0, n_inconsistent = 0;
  const detail = [];
  for (const c of candidateList) {
    const p = pK[c.candidate];
    const o = oK[c.candidate];
    if (!p || p === 'NA' || !o || o === 'NA') continue;
    // Hard Mendelian: parent 0/0 + offspring 1/1 = inconsistent.
    if ((p === '0/0' && o === '1/1') || (p === '1/1' && o === '0/0')) {
      n_inconsistent++;
      detail.push({ candidate: c.candidate, p_kt: p, o_kt: o, status: 'fail' });
      continue;
    }
    n_consistent++;
    if (p === '0/1') {
      if (o === '0/0') { n_zero++; n_inf++; }
      else if (o === '1/1') { n_one++; n_inf++; }
      // o = '0/1' is ambiguous for transmission — skip
    }
    detail.push({ candidate: c.candidate, p_kt: p, o_kt: o, status: 'pass' });
  }
  // Binomial transmission test on heterozygous sites
  const p_val = (n_inf > 0)
    ? binomialPValueTwoSided(n_zero, n_inf, 0.5)
    : NaN;
  // Mendelian-consistency P (binomial: n_inconsistent / n_total ~ 0)
  const n_total = n_consistent + n_inconsistent;
  const consistency_p = (n_total > 0)
    ? binomialPValueTwoSided(n_inconsistent, n_total, 0.02) // ~2% baseline error
    : NaN;
  return {
    mode: 'dyad',
    parent: parentId, offspring: offspringId,
    n_total, n_consistent, n_inconsistent,
    n_informative: n_inf, n_zero, n_one,
    transmission_p: p_val,
    consistency_p,
    detail,
  };
}

// Run a triad test. For each candidate, derive expected offspring distribution
// from parental karyotypes; check observed offspring is in the expected
// support set; tally counts. Cohort statistic: chi-square of the per-class
// counts vs the expected proportions, summed across all candidates with the
// same parental class.
function runTriadTest(p1Id, p2Id, oId, candidateList) {
  const p1K = DEMO.karyotype_matrix[p1Id] || {};
  const p2K = DEMO.karyotype_matrix[p2Id] || {};
  const oK  = DEMO.karyotype_matrix[oId]  || {};
  let n_total = 0, n_consistent = 0, n_inconsistent = 0;
  // Track per-class observations for chi-square test (using only informative
  // het×het classes since pure-homozygous cases have point-expectations).
  let n_het_het = 0;            // candidates where both parents 0/1
  let n_het_het_obs00 = 0, n_het_het_obs01 = 0, n_het_het_obs11 = 0;
  const detail = [];
  for (const c of candidateList) {
    const p1 = p1K[c.candidate];
    const p2 = p2K[c.candidate];
    const o  = oK[c.candidate];
    if (!p1 || p1 === 'NA' || !p2 || p2 === 'NA' || !o || o === 'NA') continue;
    const expected = expectedOffspringPrior(p1, p2);
    if (!expected) continue;
    const o_idx = o === '0/0' ? 0 : (o === '0/1' ? 1 : 2);
    const isInExpected = expected[o_idx] > 0;
    n_total++;
    if (isInExpected) n_consistent++;
    else              n_inconsistent++;
    if (p1 === '0/1' && p2 === '0/1') {
      n_het_het++;
      if (o === '0/0') n_het_het_obs00++;
      else if (o === '0/1') n_het_het_obs01++;
      else                  n_het_het_obs11++;
    }
    detail.push({
      candidate: c.candidate, p1_kt: p1, p2_kt: p2, o_kt: o,
      expected, status: isInExpected ? 'pass' : 'fail',
    });
  }
  const consistency_p = (n_total > 0)
    ? binomialPValueTwoSided(n_inconsistent, n_total, 0.02)
    : NaN;
  // Chi-square of het×het class against expected 1:2:1
  let chi2 = NaN, chi2_p = NaN, chi2_df = 2;
  if (n_het_het >= 5) {
    const exp00 = n_het_het * 0.25;
    const exp01 = n_het_het * 0.50;
    const exp11 = n_het_het * 0.25;
    const x = ((n_het_het_obs00 - exp00) ** 2) / exp00
            + ((n_het_het_obs01 - exp01) ** 2) / exp01
            + ((n_het_het_obs11 - exp11) ** 2) / exp11;
    chi2 = x;
    chi2_p = chiSquarePValue(x, chi2_df);
  }
  return {
    mode: 'triad',
    parent1: p1Id, parent2: p2Id, offspring: oId,
    n_total, n_consistent, n_inconsistent,
    consistency_p,
    n_het_het, n_het_het_obs00, n_het_het_obs01, n_het_het_obs11,
    chi2, chi2_p, chi2_df,
    detail,
  };
}

// Apply the inversion subset filter from the UI to the candidate list.
function getMendCandidates() {
  let list = DEMO.inversion_candidates_full;
  const subset = $('#mendInvSubset').value;
  if (subset === 'status_pass') list = list.filter(i => i.status === 'pass');
  else if (subset === 'status_warn') list = list.filter(i => i.status === 'warn');
  else if (subset === 'freq_high') list = list.filter(i => i.frequency >= 0.10);
  else if (subset === 'chrom_only' && state.selected_chromosome) {
    list = list.filter(i => i.chromosome === state.selected_chromosome);
  }
  return list;
}

// Run + render
function runMendelianTest() {
  const mode = state.mend.mode;
  const candidates = getMendCandidates();
  const fam = DEMO.families.find(f => f.family_id === state.selected_family) || DEMO.families[0];
  let result;
  if (mode === 'dyad') {
    result = runDyadTest(state.mend.parent1, state.mend.offspring, candidates);
  } else if (mode === 'triad') {
    result = runTriadTest(state.mend.parent1, state.mend.parent2, state.mend.offspring, candidates);
  } else if (mode === 'all_dyads') {
    // Test every (hub, member) dyad in the family.
    const hub = fam.hub_individual || fam.members[0];
    const others = fam.members.filter(m => m !== hub);
    const results = others.map(m => runDyadTest(hub, m, candidates));
    result = combineCohortDyads(results, hub);
  } else if (mode === 'all_triads') {
    // Treat first two members as parents, others as offspring (demo).
    const ms = fam.members;
    if (ms.length < 3) {
      result = { mode: 'all_triads', error: 'Hub has fewer than 3 members; no triads possible.' };
    } else {
      const p1 = ms[0], p2 = ms[1];
      const offspring = ms.slice(2);
      const results = offspring.map(o => runTriadTest(p1, p2, o, candidates));
      result = combineCohortTriads(results, p1, p2);
    }
  }
  state.mend.last_result = result;
  renderMendResult(result);
}

function combineCohortDyads(results, hub) {
  // Stouffer-style combination (z-score-based) of P-values.
  const valid = results.filter(r => Number.isFinite(r.consistency_p));
  const z_combined = valid.length
    ? valid.reduce((s, r) => s + invNormCdf(1 - r.consistency_p / 2), 0) / Math.sqrt(valid.length)
    : 0;
  const combined_p = 2 * (1 - normCdf(Math.abs(z_combined)));
  return {
    mode: 'all_dyads', hub, n_dyads: results.length,
    results, combined_p,
    summary: results.map(r => ({
      pair: r.parent + ' × ' + r.offspring,
      consistency_p: r.consistency_p,
      n_total: r.n_total,
      n_inconsistent: r.n_inconsistent,
    })),
  };
}

function combineCohortTriads(results, p1, p2) {
  const valid = results.filter(r => Number.isFinite(r.consistency_p));
  const z_combined = valid.length
    ? valid.reduce((s, r) => s + invNormCdf(1 - r.consistency_p / 2), 0) / Math.sqrt(valid.length)
    : 0;
  const combined_p = 2 * (1 - normCdf(Math.abs(z_combined)));
  return {
    mode: 'all_triads', parent1: p1, parent2: p2,
    n_triads: results.length,
    results, combined_p,
    summary: results.map(r => ({
      pair: r.parent1 + ' × ' + r.parent2 + ' → ' + r.offspring,
      consistency_p: r.consistency_p,
      chi2_p: r.chi2_p,
      n_total: r.n_total,
      n_inconsistent: r.n_inconsistent,
    })),
  };
}

// Standard-normal CDF + inverse
function normCdf(x) {
  const a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741;
  const a4 = -1.453152027, a5 =  1.061405429, p_ =  0.3275911;
  const sign = x < 0 ? -1 : 1;
  const xz = Math.abs(x) / Math.SQRT2;
  const t_ = 1 / (1 + p_ * xz);
  const y = 1 - (((((a5 * t_ + a4) * t_) + a3) * t_ + a2) * t_ + a1) * t_ * Math.exp(-xz * xz);
  return 0.5 * (1 + sign * y);
}
function invNormCdf(p) {
  // Beasley–Springer–Moro approximation
  if (p <= 0) return -Infinity;
  if (p >= 1) return  Infinity;
  const a = [-3.969683028665376e+01,  2.209460984245205e+02,
             -2.759285104469687e+02,  1.383577518672690e+02,
             -3.066479806614716e+01,  2.506628277459239e+00];
  const b = [-5.447609879822406e+01,  1.615858368580409e+02,
             -1.556989798598866e+02,  6.680131188771972e+01,
             -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01,
             -2.400758277161838e+00, -2.549732539343734e+00,
              4.374664141464968e+00,  2.938163982698783e+00];
  const d = [ 7.784695709041462e-03,  3.224671290700398e-01,
              2.445134137142996e+00,  3.754408661907416e+00];
  const plow = 0.02425, phigh = 1 - plow;
  if (p < plow) {
    const q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= phigh) {
    const q = p - 0.5; const r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}


function renderMendResult(r) {
  const sumSlot = $('#mendSummary');
  const detSlot = $('#mendResultSlot');
  sumSlot.innerHTML = '';
  detSlot.innerHTML = '';
  if (!r) return;
  if (r.error) {
    detSlot.appendChild(el('div', {
      style: { color: 'var(--bad)', padding: '14px',
               background: 'var(--panel-2)',
               border: '1px solid var(--rule)', borderRadius: '4px' },
      text: r.error,
    }));
    return;
  }

  // Build summary cells (top-of-page)
  if (r.mode === 'dyad') {
    sumSlot.appendChild(summaryCell('mode', 'Dyad test'));
    sumSlot.appendChild(summaryCell('candidates tested', r.n_total));
    sumSlot.appendChild(summaryCell('consistent', r.n_consistent));
    sumSlot.appendChild(summaryCell('inconsistent (Mendelian errors)', r.n_inconsistent,
                                     r.n_inconsistent > 0 ? 'fail' : 'good'));
    sumSlot.appendChild(summaryCell('consistency P-value',
                                     fmt(r.consistency_p),
                                     pSeverity(r.consistency_p, $('#mendAlpha').value)));
    sumSlot.appendChild(summaryCell('informative (parent het)', r.n_informative,
                                     null,
                                     r.n_informative > 0
                                       ? `${r.n_zero}× allele 0, ${r.n_one}× allele 1`
                                       : ''));
    if (r.n_informative > 0) {
      sumSlot.appendChild(summaryCell('transmission P-value (binomial 0.5)',
                                       fmt(r.transmission_p),
                                       pSeverity(r.transmission_p, $('#mendAlpha').value)));
    }
  } else if (r.mode === 'triad') {
    sumSlot.appendChild(summaryCell('mode', 'Triad test'));
    sumSlot.appendChild(summaryCell('candidates tested', r.n_total));
    sumSlot.appendChild(summaryCell('consistent', r.n_consistent));
    sumSlot.appendChild(summaryCell('inconsistent', r.n_inconsistent,
                                     r.n_inconsistent > 0 ? 'fail' : 'good'));
    sumSlot.appendChild(summaryCell('consistency P-value',
                                     fmt(r.consistency_p),
                                     pSeverity(r.consistency_p, $('#mendAlpha').value)));
    if (r.n_het_het >= 5) {
      sumSlot.appendChild(summaryCell('het × het class size', r.n_het_het));
      sumSlot.appendChild(summaryCell('het×het χ² P (1:2:1)', fmt(r.chi2_p),
                                       pSeverity(r.chi2_p, $('#mendAlpha').value)));
    }
  } else if (r.mode === 'all_dyads') {
    sumSlot.appendChild(summaryCell('mode', 'Cohort dyad test'));
    sumSlot.appendChild(summaryCell('hub', r.hub));
    sumSlot.appendChild(summaryCell('dyads tested', r.n_dyads));
    sumSlot.appendChild(summaryCell('combined P (Stouffer)',
                                     fmt(r.combined_p),
                                     pSeverity(r.combined_p, $('#mendAlpha').value)));
  } else if (r.mode === 'all_triads') {
    sumSlot.appendChild(summaryCell('mode', 'Cohort triad test'));
    sumSlot.appendChild(summaryCell('parent₁', r.parent1));
    sumSlot.appendChild(summaryCell('parent₂', r.parent2));
    sumSlot.appendChild(summaryCell('triads tested', r.n_triads));
    sumSlot.appendChild(summaryCell('combined P (Stouffer)',
                                     fmt(r.combined_p),
                                     pSeverity(r.combined_p, $('#mendAlpha').value)));
  }

  // Detailed table
  if (r.mode === 'dyad' || r.mode === 'triad') {
    const tbl = el('table', { class: 'data-table',
                               style: { marginTop: '12px' } });
    const thead = el('thead');
    const tr = el('tr');
    const cols = r.mode === 'dyad'
      ? ['Candidate','Parent','Offspring','Status']
      : ['Candidate','Parent 1','Parent 2','Offspring','Expected support','Status'];
    cols.forEach(c => tr.appendChild(el('th', { text: c })));
    thead.appendChild(tr);
    tbl.appendChild(thead);
    const tbody = el('tbody');
    r.detail.slice(0, 50).forEach(row => {
      const t = el('tr');
      t.appendChild(el('td', { class: 'sample-id', text: row.candidate }));
      if (r.mode === 'dyad') {
        t.appendChild(el('td', { text: row.p_kt }));
        t.appendChild(el('td', { text: row.o_kt }));
      } else {
        t.appendChild(el('td', { text: row.p1_kt }));
        t.appendChild(el('td', { text: row.p2_kt }));
        t.appendChild(el('td', { text: row.o_kt }));
        t.appendChild(el('td', {
          text: row.expected.map((p, i) => p > 0
                  ? (['0/0','0/1','1/1'][i] + '(' + (p*100).toFixed(0) + '%)') : null)
                  .filter(Boolean).join(' '),
          style: { color: 'var(--ink-dim)', fontSize: '10px' }
        }));
      }
      const tdSt = el('td');
      tdSt.appendChild(el('span', {
        class: 'status-pill-cell ' + (row.status === 'pass' ? 'pass' : 'fail'),
        text: row.status.toUpperCase(),
      }));
      t.appendChild(tdSt);
      tbody.appendChild(t);
    });
    if (r.detail.length > 50) {
      const tr2 = el('tr');
      tr2.appendChild(el('td', {
        colspan: cols.length,
        text: `… ${r.detail.length - 50} more rows omitted (export to TSV for full set)`,
        style: { textAlign: 'center', color: 'var(--ink-dim)',
                 fontStyle: 'italic', padding: '8px' },
      }));
      tbody.appendChild(tr2);
    }
    tbl.appendChild(tbody);
    detSlot.appendChild(tbl);
  } else if (r.mode === 'all_dyads' || r.mode === 'all_triads') {
    const tbl = el('table', { class: 'data-table',
                               style: { marginTop: '12px' } });
    const thead = el('thead');
    const tr = el('tr');
    ['Pair','N total','N inconsistent','Consistency P','χ² P (where applicable)']
      .forEach(c => tr.appendChild(el('th', { text: c })));
    thead.appendChild(tr);
    tbl.appendChild(thead);
    const tbody = el('tbody');
    r.summary.forEach(row => {
      const t = el('tr');
      t.appendChild(el('td', { class: 'sample-id', text: row.pair }));
      t.appendChild(el('td', { class: 'num', text: String(row.n_total) }));
      t.appendChild(el('td', { class: 'num', text: String(row.n_inconsistent) }));
      const tdP = el('td', { class: 'num', text: fmt(row.consistency_p) });
      const sev = pSeverity(row.consistency_p, $('#mendAlpha').value);
      if (sev === 'fail') tdP.style.color = 'var(--bad)';
      else if (sev === 'warn') tdP.style.color = 'var(--warn)';
      else if (sev === 'good') tdP.style.color = 'var(--good)';
      t.appendChild(tdP);
      t.appendChild(el('td', { class: 'num',
        text: row.chi2_p === undefined ? '—' : fmt(row.chi2_p) }));
      tbody.appendChild(t);
    });
    tbl.appendChild(tbody);
    detSlot.appendChild(tbl);
  }
}

function summaryCell(label, value, severity = null, sub = '') {
  const cell = el('div', { class: 'mend-summary-cell' },
    el('div', { class: 'lbl', text: label }),
    el('div', { class: 'val', text: String(value) }),
    sub ? el('div', { class: 'sub', text: sub }) : null,
  );
  if (severity === 'fail')      cell.style.borderColor = 'rgba(224,85,92,0.55)';
  else if (severity === 'warn') cell.style.borderColor = 'rgba(232,196,76,0.55)';
  else if (severity === 'good') cell.style.borderColor = 'rgba(95,212,154,0.55)';
  return cell;
}

function pSeverity(p, alphaSelect) {
  if (!Number.isFinite(p)) return null;
  let alpha;
  if (alphaSelect === 'bonferroni') {
    // approximate: divide by number of inversions tested
    alpha = 0.05 / Math.max(1, getMendCandidates().length);
  } else {
    alpha = parseFloat(alphaSelect);
  }
  if (p < alpha)        return 'fail';      // significantly inconsistent
  if (p < alpha * 5)    return 'warn';
  return 'good';
}


/* Mendelian-tab event wiring */
$('#mendMode').addEventListener('change', () => updateMendModeUI());
$('#mendParent1').addEventListener('change', e => state.mend.parent1 = e.target.value);
$('#mendParent2').addEventListener('change', e => state.mend.parent2 = e.target.value);
$('#mendOffspring').addEventListener('change', e => state.mend.offspring = e.target.value);
$('#mendInvSubset').addEventListener('change', e => state.mend.inv_subset = e.target.value);
$('#mendAlpha').addEventListener('change', e => state.mend.alpha = e.target.value);
$('#mendRunBtn').addEventListener('click', runMendelianTest);
$('#mendResetBtn').addEventListener('click', () => {
  state.mend.last_result = null;
  $('#mendSummary').innerHTML = '';
  $('#mendResultSlot').innerHTML = '';
});
$('#mendExportBtn').addEventListener('click', () => {
  const r = state.mend.last_result;
  if (!r) { alert('Run a test first.'); return; }
  const tsv = mendResultToTsv(r);
  const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'mendelian_' + r.mode + '_' + Date.now() + '.tsv';
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1500);
});

function mendResultToTsv(r) {
  if (r.mode === 'dyad' || r.mode === 'triad') {
    const cols = r.mode === 'dyad'
      ? ['candidate','parent_kt','offspring_kt','status']
      : ['candidate','p1_kt','p2_kt','o_kt','expected_support','status'];
    const lines = [cols.join('\t')];
    r.detail.forEach(row => {
      const cells = r.mode === 'dyad'
        ? [row.candidate, row.p_kt, row.o_kt, row.status]
        : [row.candidate, row.p1_kt, row.p2_kt, row.o_kt,
           row.expected.map((p, i) => ['0/0','0/1','1/1'][i] + ':' + p).join(';'),
           row.status];
      lines.push(cells.join('\t'));
    });
    return lines.join('\n') + '\n';
  } else {
    const cols = ['pair','n_total','n_inconsistent','consistency_p','chi2_p'];
    const lines = [cols.join('\t')];
    r.summary.forEach(row => {
      lines.push([row.pair, row.n_total, row.n_inconsistent,
                  fmt(row.consistency_p),
                  row.chi2_p === undefined ? '' : fmt(row.chi2_p)].join('\t'));
    });
    return lines.join('\n') + '\n';
  }
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §9b Compatibility tab — breeding-partner finder.

   Given a focal individual + a target offspring karyotype + an inversion
   candidate (or a chromosome's worth of candidates), find every other
   individual in the cohort whose karyotype, when crossed with the focal,
   could produce the target karyotype in offspring.

   Mendelian compatibility table:
     focal kt   target kt   →   compatible partner kt
       0/0         0/0      →   {0/0, 0/1}      (any partner with 0 allele)
       0/0         0/1      →   {0/1, 1/1}      (any partner with 1 allele)
       0/0         1/1      →   {1/1, 0/1}      (BUT: 0/1 only gives 50% chance;
                                                     1/1 mandatory if you want 100%
                                                     yield. We separate "guaranteed"
                                                     from "possible".)
       0/1         0/0      →   {0/0, 0/1}
       0/1         0/1      →   any of {0/0, 0/1, 1/1}
       0/1         1/1      →   {0/1, 1/1}
       1/1         0/0      →   {0/0, 0/1}      (0/1 only gives 50% chance)
       1/1         0/1      →   {0/0, 0/1}
       1/1         1/1      →   {0/1, 1/1}

   For each candidate partner, we report:
     - guaranteed-yield: producing the target with probability 1.0
     - possible-yield:   probability < 1.0 but > 0
     - none:             probability = 0 (not listed)

   Per-chrom mode: compute compatibility across every inversion on the
   chosen chromosome and only list partners that satisfy the constraints
   for ALL of them simultaneously.

   Sex-aware mode: filter partners to require opposite sex.
   ╚══════════════════════════════════════════════════════════════════════╝ */

// Compatibility table: focal kt → target kt → {guaranteed: [], possible: []}
const COMPAT_RULES = {
  '0/0': {
    '0/0': { guaranteed: ['0/0'], possible: ['0/1'] },
    '0/1': { guaranteed: ['1/1'], possible: ['0/1'] },
    '1/1': { guaranteed: [],      possible: ['0/1'] },
  },
  '0/1': {
    '0/0': { guaranteed: [],      possible: ['0/0','0/1'] },
    '0/1': { guaranteed: [],      possible: ['0/0','0/1','1/1'] },
    '1/1': { guaranteed: [],      possible: ['0/1','1/1'] },
  },
  '1/1': {
    '0/0': { guaranteed: [],      possible: ['0/1'] },
    '0/1': { guaranteed: ['0/0'], possible: ['0/1'] },
    '1/1': { guaranteed: ['1/1'], possible: ['0/1'] },
  },
};

// Compute the offspring karyotype probabilities for parental cross.
// Returns {0/0:p, 0/1:p, 1/1:p}.
function offspringDist(p1, p2) {
  if (!p1 || !p2 || p1 === 'NA' || p2 === 'NA') return null;
  const probs = { '0/0': 0, '0/1': 0, '1/1': 0 };
  const a = p1.split('/').map(Number);
  const b = p2.split('/').map(Number);
  for (const aa of a) for (const bb of b) {
    const sum = aa + bb;
    const k = sum === 0 ? '0/0' : (sum === 1 ? '0/1' : '1/1');
    probs[k] += 0.25;
  }
  return probs;
}

function evaluatePartnership(focalId, partnerId, invSet, targetKt) {
  // For each inversion in invSet, compute offspring distribution and check
  // whether target is reachable. Returns aggregate verdict + per-inversion list.
  const focalK   = DEMO.karyotype_matrix[focalId]   || {};
  const partnerK = DEMO.karyotype_matrix[partnerId] || {};
  const perInv = [];
  let n_guaranteed = 0, n_possible = 0, n_impossible = 0, n_unknown = 0;
  for (const inv of invSet) {
    const fk = focalK[inv.candidate];
    const pk = partnerK[inv.candidate];
    if (!fk || !pk || fk === 'NA' || pk === 'NA') {
      perInv.push({ inv: inv.candidate, focal_kt: fk || 'NA', partner_kt: pk || 'NA',
                     prob_target: null, status: 'unknown' });
      n_unknown++; continue;
    }
    const dist = offspringDist(fk, pk);
    const p = dist[targetKt] || 0;
    let status;
    if (p >= 0.999)      { status = 'guaranteed'; n_guaranteed++; }
    else if (p > 0)      { status = 'possible';   n_possible++; }
    else                 { status = 'impossible'; n_impossible++; }
    perInv.push({ inv: inv.candidate, focal_kt: fk, partner_kt: pk,
                   prob_target: p, status });
  }
  // Aggregate verdict for the partner across the whole invSet:
  //   if any impossible → REJECT
  //   if all guaranteed → IDEAL
  //   if some guaranteed + rest possible → GOOD
  //   if all possible → POSSIBLE
  //   if mostly unknown → UNKNOWN
  let verdict;
  if (n_impossible > 0)            verdict = 'reject';
  else if (n_guaranteed === invSet.length - n_unknown && n_guaranteed > 0)
                                   verdict = 'ideal';
  else if (n_guaranteed > 0)       verdict = 'good';
  else if (n_possible > 0)         verdict = 'possible';
  else                             verdict = 'unknown';
  return { partner_id: partnerId, n_guaranteed, n_possible, n_impossible, n_unknown,
           verdict, perInv };
}

function getCompatInversionSet() {
  const scope = $('#compatInvScope').value;
  if (scope === 'single') {
    const id = $('#compatInvSingle').value;
    return DEMO.inversion_candidates_full.filter(i => i.candidate === id);
  } else if (scope === 'chrom') {
    const c = $('#compatChrom').value;
    return DEMO.inversion_candidates_full.filter(i => i.chromosome === c);
  } else {
    return DEMO.inversion_candidates_full;
  }
}

function runCompatibilitySearch() {
  const focalId  = $('#compatFocal').value;
  const targetKt = $('#compatTarget').value;
  const invSet   = getCompatInversionSet();
  const sexAware = $('#compatSexAware').checked;
  const excludeKin = $('#compatExcludeKin').checked;
  const excludeAmb = $('#compatExcludeAmbig').checked;

  state.compat = {
    focal: focalId, target: targetKt,
    scope: $('#compatInvScope').value,
    inv_single: $('#compatInvSingle').value,
    chrom: $('#compatChrom').value,
    sex_aware: sexAware,
    exclude_kin: excludeKin,
    exclude_amb: excludeAmb,
    last_results: null,
  };

  // Build candidate-partner list (all individuals except focal)
  let candidates = DEMO.individuals.filter(p => p !== focalId);

  // Sex filter
  let sex_data_available = true;
  let sex_filter_applied = false;
  if (sexAware) {
    const focalSex = (DEMO.sex || {})[focalId];
    if (!focalSex || focalSex === '?' || focalSex === 'unknown') {
      sex_data_available = false;
    } else {
      candidates = candidates.filter(p => {
        const ps = (DEMO.sex || {})[p];
        if (!ps || ps === '?') return false;
        return ps !== focalSex;
      });
      sex_filter_applied = true;
    }
  }

  // Exclude close kin (PO/FS) — uses the network_edges as a proxy
  if (excludeKin) {
    const kinList = new Set();
    DEMO.network_edges.forEach(e => {
      if (e.class === 'strong_po' || e.class === 'possible_po') {
        if (e.a === focalId) kinList.add(e.b);
        if (e.b === focalId) kinList.add(e.a);
      }
    });
    candidates = candidates.filter(p => !kinList.has(p));
  }

  // Score each candidate
  const results = candidates.map(p => evaluatePartnership(focalId, p, invSet, targetKt));

  // Optionally exclude partners with too many unknowns
  // (skipped for now — let user see all results)

  // Sort: ideal > good > possible > unknown > reject
  const order = { ideal: 0, good: 1, possible: 2, unknown: 3, reject: 4 };
  results.sort((a, b) => {
    if (order[a.verdict] !== order[b.verdict]) return order[a.verdict] - order[b.verdict];
    return b.n_guaranteed - a.n_guaranteed;
  });

  state.compat.last_results = {
    focal: focalId, target: targetKt, invSet,
    sex_data_available, sex_filter_applied,
    results,
  };
  renderCompatibilityResults();
}

function renderCompatibilityResults() {
  const sumSlot = $('#compatSummary');
  const resSlot = $('#compatResultSlot');
  sumSlot.innerHTML = ''; resSlot.innerHTML = '';
  const r = state.compat.last_results;
  if (!r) return;

  // Summary text block
  const focalSex = (DEMO.sex || {})[r.focal] || '?';
  const sexNote = state.compat.sex_aware
    ? (r.sex_data_available
        ? ` (filtering for opposite-sex partners; focal is ${focalSex})`
        : ' (sex-aware enabled but no sex data for focal — falling back to non-sex-aware)')
    : '';
  const n_total = r.results.length;
  const n_ideal = r.results.filter(x => x.verdict === 'ideal').length;
  const n_good  = r.results.filter(x => x.verdict === 'good').length;
  const n_poss  = r.results.filter(x => x.verdict === 'possible').length;
  const n_rej   = r.results.filter(x => x.verdict === 'reject').length;
  sumSlot.appendChild(el('div', { class: 'compat-summary', html:
    `<b>${r.focal}</b> ${sexBadgeHtml(r.focal)} · target offspring `
    + `<b>${r.target}</b>${sexNote}<br/>`
    + `Tested across <b>${r.invSet.length}</b> inversion candidate(s). `
    + `<b>${n_total}</b> partners considered: `
    + `<b>${n_ideal}</b> ideal · <b>${n_good}</b> good · `
    + `<b>${n_poss}</b> possible-only · <b>${n_rej}</b> rejected.`
  }));

  // Results table
  const tbl = el('table', { class: 'data-table' });
  const thead = el('thead');
  const tr = el('tr');
  ['Partner','Sex','Verdict','Guaranteed','Possible','Impossible','Per-inversion preview']
    .forEach(h => tr.appendChild(el('th', { text: h })));
  thead.appendChild(tr); tbl.appendChild(thead);
  const tbody = el('tbody');
  r.results.forEach(p => {
    const t = el('tr');
    t.appendChild(el('td', { class: 'sample-id', text: p.partner_id }));
    t.appendChild((function(){
      const td = el('td'); td.innerHTML = sexBadgeHtml(p.partner_id);
      return td;
    })());
    {
      const td = el('td');
      const cls = p.verdict === 'ideal'    ? 'pass'
                : p.verdict === 'good'     ? 'pass'
                : p.verdict === 'possible' ? 'warn'
                : p.verdict === 'unknown'  ? 'weak'
                                           : 'fail';
      td.appendChild(el('span', {
        class: p.verdict === 'unknown' ? 'tier-pill weak' : 'status-pill-cell ' + cls,
        text: p.verdict.toUpperCase(),
      }));
      t.appendChild(td);
    }
    t.appendChild(el('td', { class: 'num', text: String(p.n_guaranteed) }));
    t.appendChild(el('td', { class: 'num', text: String(p.n_possible) }));
    t.appendChild(el('td', { class: 'num', text: String(p.n_impossible),
                             style: p.n_impossible > 0 ? { color: 'var(--bad)' } : {} }));
    // Mini preview: first 5 inversions with per-inv kt × kt → outcome
    const prev = p.perInv.slice(0, 5).map(pi =>
      pi.inv + ': ' + pi.focal_kt + '×' + pi.partner_kt
        + (pi.status === 'guaranteed' ? ' ✓✓'
            : pi.status === 'possible' ? ' ~'
            : pi.status === 'impossible' ? ' ✗'
            : ' ?')
    ).join('   ');
    t.appendChild(el('td', { text: prev,
      style: { fontSize: '9px', color: 'var(--ink-dim)' } }));
    tbody.appendChild(t);
  });
  tbl.appendChild(tbody);
  resSlot.appendChild(tbl);

  // Recommendation footnote
  if (n_ideal > 0) {
    resSlot.appendChild(el('div', {
      class: 'ie-conclusion tier-strong',
      style: { marginTop: '12px' },
      html: `<div class="verdict">${n_ideal} IDEAL PARTNER${n_ideal>1?'S':''}</div>`
          + `Crossing <b>${r.focal}</b> with any of the IDEAL-tier partners `
          + `produces the <b>${r.target}</b> target with 100% certainty across all `
          + `${r.invSet.length} tested inversion(s). These are the recommended crosses.`
    }));
  } else if (n_good > 0) {
    resSlot.appendChild(el('div', {
      class: 'ie-conclusion tier-moderate',
      style: { marginTop: '12px' },
      html: `<div class="verdict">${n_good} GOOD PARTNER${n_good>1?'S':''}</div>`
          + `These produce the target with probability 1.0 for some inversions and `
          + `> 0 for others (no rejections). Pick based on which inversions you want `
          + `guaranteed vs. possible.`
    }));
  } else if (n_poss > 0) {
    resSlot.appendChild(el('div', {
      class: 'ie-conclusion tier-weak',
      style: { marginTop: '12px' },
      html: '<div class="verdict">NO IDEAL PARTNERS</div>'
          + `Only "possible" partners exist (target reachable with probability `
          + `< 1.0). Crossing in numbers may yield the target by chance, but no `
          + `partner guarantees it.`
    }));
  } else {
    resSlot.appendChild(el('div', {
      class: 'ie-conclusion tier-conflict',
      style: { marginTop: '12px' },
      html: '<div class="verdict">NO COMPATIBLE PARTNERS</div>'
          + `Every candidate produced an impossible cross for at least one tested `
          + `inversion. Try relaxing the inversion scope or the target karyotype.`
    }));
  }
}

function sexBadgeHtml(ind) {
  const s = (DEMO.sex || {})[ind] || '?';
  const cls = s === 'F' ? 'female' : s === 'M' ? 'male' : 'unk';
  return '<span class="sex-pill ' + cls + '">' + s + '</span>';
}

function exportCompatibilityTsv() {
  const r = state.compat.last_results;
  if (!r) { alert('Run a compatibility search first.'); return; }
  const cols = ['partner_id','partner_sex','verdict','n_guaranteed','n_possible','n_impossible','n_unknown'];
  const lines = [
    '# Compatibility search',
    '# Focal: ' + r.focal + '  Sex: ' + ((DEMO.sex || {})[r.focal] || '?'),
    '# Target offspring karyotype: ' + r.target,
    '# Inversion set: ' + r.invSet.length + ' candidate(s)',
    '# Sex-aware: ' + state.compat.sex_aware
      + (state.compat.sex_aware && !r.sex_data_available
          ? ' (no sex data — fallback to non-sex-aware)' : ''),
    '# Exclude close kin: ' + state.compat.exclude_kin,
    '#',
    cols.join('\t'),
  ];
  r.results.forEach(p => {
    lines.push([p.partner_id, (DEMO.sex || {})[p.partner_id] || '?',
                p.verdict, p.n_guaranteed, p.n_possible,
                p.n_impossible, p.n_unknown].join('\t'));
  });
  // Per-inversion detail at the bottom
  lines.push('');
  lines.push('# Per-inversion detail');
  lines.push(['partner_id','inv_id','focal_kt','partner_kt','prob_target','status'].join('\t'));
  r.results.forEach(p => {
    p.perInv.forEach(pi => {
      lines.push([p.partner_id, pi.inv, pi.focal_kt, pi.partner_kt,
                  pi.prob_target === null ? '' : pi.prob_target.toFixed(3),
                  pi.status].join('\t'));
    });
  });
  downloadTsv('compatibility_' + r.focal + '_' + r.target.replace('/','-')
              + '_' + Date.now() + '.tsv',
              lines.join('\n') + '\n');
}

function populateCompatSelectors() {
  // Focal picker
  const focalSel = $('#compatFocal');
  focalSel.innerHTML = '';
  DEMO.individuals.forEach(i => focalSel.appendChild(el('option', { value: i, text: i })));
  // Single-inversion picker — first 50 by candidate id (more is overwhelming)
  const invSingleSel = $('#compatInvSingle');
  invSingleSel.innerHTML = '';
  DEMO.inversion_candidates_full.slice(0, 50).forEach(i =>
    invSingleSel.appendChild(el('option', { value: i.candidate,
      text: i.candidate + ' — ' + i.chromosome + ' ' + i.start_mb.toFixed(1)
            + '-' + i.end_mb.toFixed(1) + 'Mb' })));
  // Chromosome picker
  const chromSel = $('#compatChrom');
  chromSel.innerHTML = '';
  DEMO.chromosomes.forEach(c => chromSel.appendChild(el('option', { value: c, text: c })));
  // Default to Chr28 to align with the inversion atlas focus
  chromSel.value = 'Chr28';
}

function updateCompatScopeUI() {
  const scope = $('#compatInvScope').value;
  $('#compatSingleRow').style.display = (scope === 'single') ? '' : 'none';
  $('#compatChromRow').style.display  = (scope === 'chrom')  ? '' : 'none';
}

// State container for compatibility tab
state.compat = {
  focal: 'Ind_001', target: '0/1',
  scope: 'chrom', inv_single: null, chrom: 'Chr28',
  sex_aware: false, exclude_kin: true, exclude_amb: true,
  last_results: null,
};

// Compatibility-tab event wiring
$('#compatInvScope').addEventListener('change', updateCompatScopeUI);
$('#compatRunBtn').addEventListener('click', runCompatibilitySearch);
$('#compatResetBtn').addEventListener('click', () => {
  state.compat.last_results = null;
  $('#compatSummary').innerHTML = '';
  $('#compatResultSlot').innerHTML = '';
});
$('#compatExportBtn').addEventListener('click', exportCompatibilityTsv);


/* ╔══════════════════════════════════════════════════════════════════════╗
   §10 Inspector / Stats (right column)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderInspector() {
  const pair = state.inspector_pair;
  $('#insSelection').innerHTML = '<a href="#" style="color: var(--accent); text-decoration: none;">' +
    pair.a + '</a> <span class="x">×</span> <a href="#" style="color: var(--accent); text-decoration: none;">' +
    pair.b + '</a>';

  const stats = DEMO.pairwise_stats[pair.b] || DEMO.pairwise_stats['Ind_044'];
  const grid = $('#insStatsGrid');
  grid.innerHTML = '';
  function addRow(lbl, val, type = '') {
    grid.appendChild(el('div', { class: 'ins-stat-label', text: lbl }));
    grid.appendChild(el('div', { class: 'ins-stat-value' + (type ? ' ' + type : ''),
                                  text: val }));
  }
  addRow('k0 (IBS0)', fmt(stats.k0));
  addRow('Kinship (φ)', fmt(stats.kinship));
  addRow('k1 (IBS1)', fmt(stats.k1));
  addRow('IBS0', fmt(stats.IBS0));
  addRow('k2 (IBS2)', fmt(stats.k2));
  addRow('PO distance', stats.PO_distance.toFixed(0));
  grid.appendChild(el('div', { class: 'ins-stat-label', text: '' }));
  grid.appendChild(el('div', { class: 'ins-stat-label', text: '' }));
  grid.appendChild(el('div', { class: 'ins-stat-label', text: 'Relationship class' }));
  let cls = 'tag good';
  if (stats.relationship_class.includes('possible')) cls = 'tag warn';
  else if (stats.relationship_class.includes('ambiguous')) cls = 'tag fail';
  grid.appendChild(el('div', {
    class: 'ins-stat-value ' + cls,
    text: stats.relationship_class,
    style: { textAlign: 'right' },
  }));

  // PASS/WARN/FAIL counter row
  const pwfRoot = $('#passWarnFail');
  pwfRoot.innerHTML = '';
  const cells = [
    { lbl: 'PASS',  val: stats.PASS,  cls: 'pass' },
    { lbl: 'WARN',  val: stats.WARN,  cls: 'warn' },
    { lbl: 'FAIL',  val: stats.FAIL,  cls: 'fail' },
    { lbl: 'TOTAL', val: stats.TOTAL, cls: '' },
  ];
  cells.forEach(c => {
    pwfRoot.appendChild(el('div', { class: 'pwf-cell ' + c.cls },
      el('div', { class: 'pwf-label', text: c.lbl }),
      el('div', { class: 'pwf-value', text: String(c.val) })
    ));
  });

  // Bar
  const total = stats.TOTAL || 1;
  const ppct = (stats.PASS / total * 100);
  const wpct = (stats.WARN / total * 100);
  const fpct = (stats.FAIL / total * 100);
  $('#pwfBar').innerHTML =
    `<div class="pwf-bar-pass" style="width:${ppct}%"></div>` +
    `<div class="pwf-bar-warn" style="width:${wpct}%"></div>` +
    `<div class="pwf-bar-fail" style="width:${fpct}%"></div>`;
  $('#pwfPcts').innerHTML =
    `<span>${ppct.toFixed(1)}%</span><span>${wpct.toFixed(1)}%</span><span>${fpct.toFixed(1)}%</span>`;

  // Mendelian-check sub-panel
  $('#mendCheckPair').textContent = '(' + pair.a + ' × ' + pair.b + ')';
  const mtbl = $('#mendCheckTable');
  const tbody = mtbl.querySelector('tbody');
  tbody.innerHTML = '';
  // Three rows mirroring the mockup: ChrXX with Status + Conflicts count
  // Use the per-individual stats for an at-a-glance Status assignment.
  const chrSamples = [
    { chr: 'Chr03', status: stats.WARN > 0 ? 'warn' : 'pass',
      conflicts: stats.WARN > 0 ? 2 : 0 },
    { chr: 'Chr17', status: stats.FAIL > 0 ? 'fail' : 'pass',
      conflicts: stats.FAIL > 0 ? Math.max(5, stats.FAIL * 2) : 0 },
    { chr: 'Chr28', status: 'pass', conflicts: 0 },
    { chr: '…',     status: null,   conflicts: '…' },
  ];
  chrSamples.forEach(row => {
    const tr = el('tr');
    tr.appendChild(el('td', { text: row.chr }));
    if (row.status) {
      const tdSt = el('td');
      tdSt.appendChild(el('span', {
        class: 'status-pill-cell ' + row.status,
        text: row.status.toUpperCase(),
      }));
      tr.appendChild(tdSt);
    } else {
      tr.appendChild(el('td', { text: '' }));
    }
    tr.appendChild(el('td', { class: 'num', text: String(row.conflicts) }));
    tbody.appendChild(tr);
  });
  $('#mendCheckTotal').textContent = String(
    stats.WARN * 2 + Math.max(5, stats.FAIL * 2 || 0)
  );
}

$('#insSwapBtn').addEventListener('click', () => {
  const p = state.inspector_pair;
  state.inspector_pair = { a: p.b, b: p.a };
  renderInspector();
});


/* ╔══════════════════════════════════════════════════════════════════════╗
   §11 Loaded files list (right column)
   ╚══════════════════════════════════════════════════════════════════════╝ */

function renderLoadedFiles() {
  const root = $('#loadedFilesList');
  root.innerHTML = '';
  const ICONS = { res: '🧬', beagle: '📁', prune: 'in', inv: '⌗' };
  for (const f of state.loaded_files) {
    root.appendChild(el('div', { class: 'lf-row' },
      el('div', { class: 'lf-icon', text: f.kind === 'res'    ? '◢'
                                          : f.kind === 'beagle' ? '⏷'
                                          : f.kind === 'prune'  ? '⌬'
                                          : '⌗' }),
      el('div', { class: 'lf-path', text: f.path }),
      el('div', { class: 'lf-status' + (f.loaded ? '' : ' unloaded'),
                   text: f.loaded ? 'Loaded' : 'Not loaded' }),
    ));
  }
  $('#loadedAt').textContent = 'Loaded: ' + state.loaded_at;
}


/* ╔══════════════════════════════════════════════════════════════════════╗
   §12 Initial render boot
   ╚══════════════════════════════════════════════════════════════════════╝ */

// Ancestry palette pill — fill the dots
(function fillAncestryPalette() {
  const root = $('#ancestryDots');
  root.innerHTML = '';
  DEMO.ancestry_palette.forEach(c => {
    root.appendChild(el('div', { class: 'dot', style: { background: c } }));
  });
})();

// Chromosome filter dropdown
(function fillChrFilter() {
  const sel = $('#chrFilter');
  DEMO.chromosomes.forEach(c => sel.appendChild(el('option', { value: c, text: c })));
  // Default: Chr28 (matches the screenshot)
  sel.value = 'Chr28';
  state.selected_chromosome = 'Chr28';
  sel.addEventListener('change', e => {
    state.selected_chromosome = e.target.value;
    renderInversionTables();
  });
})();

renderPopTree();
renderPopSummary();
renderNetwork();
renderKaryotypeTable('#karyoTableSlotInline');
renderKaryotypeTable('#karyoTableSlotFull',
                     { rows: DEMO.individuals,
                       columns: ['Chr01','Chr02','Chr03','Chr04','Chr05',
                                  'Chr06','Chr07','Chr08', null,
                                  'Chr17','Chr18', null, 'Chr28'] });
renderInversionTables();
populateMendSelectors();
updateMendModeUI();
populateCompatSelectors();
updateCompatScopeUI();
renderInspector();
renderLoadedFiles();

// Help button — open a quick orientation overlay (simple alert for now)
$('#helpBtn').addEventListener('click', () => {
  alert(
    'Relatedness Atlas v0.1 — Family / Individual Evidence Hub\n\n' +
    '• Left column:  Population Browser tree (search to filter; click any\n' +
    '   node to scope the view).\n' +
    '• Center column: Network → Karyotypes → Inversions → Mendelian sub-tabs.\n' +
    '   The Mendelian tab is the segregation tester — pick a dyad or triad,\n' +
    '   pick which inversions to test against, and run.\n' +
    '• Right column: Inspector / Stats — pairwise IBS coefficients for the\n' +
    '   currently selected pair, plus the Mendelian-check breakdown by\n' +
    '   chromosome.\n\n' +
    'Status pills at top: load .res, beagle, prune file, inversion candidates.\n' +
    'Theme toggle: ☀ / 📜 / 🌙 cycles through light, academic, dark.'
  );
});

console.log('[RA] Relatedness Atlas v0.1 ready.');
console.log('[RA] state =', state);

})();
