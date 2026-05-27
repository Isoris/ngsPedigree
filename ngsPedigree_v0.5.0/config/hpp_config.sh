#!/usr/bin/env bash
# hpp_config.sh — project paths and runtime knobs for HPP.
#
# Source from scripts/* or from a SLURM job script. Variables are
# intentionally minimal at MVP 1; MVP 2+ will add joint VCF and
# variant-master paths.

# Project base (set by the caller; default to LANTA layout used by
# ngsPedigree Stages 1/2).
: "${BASE:?BASE not set — export BASE=/path/to/project/root}"

# HPP output dir.
HPP_OUTDIR="${BASE}/results/ngsPedigree/04_hpp"

# Stage 3 inheritance maps (placeholder layout for now; swap when real).
STAGE3_DYAD_MAP_DIR="${BASE}/results/ngsPedigree/03_inheritance_map/dyads"
STAGE3_TRIAD_MAP_DIR="${BASE}/results/ngsPedigree/03_inheritance_map/triads"
STAGE3_PARENT_PHASE_DIR="${BASE}/results/ngsPedigree/03_inheritance_map/parent_phase"

# MODULE_CONSERVATION inputs (consumed at MVP 2+).
VARIANT_MASTER_TSV="${BASE}/results/MODULE_CONSERVATION/step16/variant_master_scored.tsv"
JOINT_VCF="${BASE}/results/MODULE_CONSERVATION/step03/joint.vcf.gz"

# Reference FASTA.
REFERENCE_FASTA="${BASE}/references/fClaHyb_Gar_LG.fa"

# Optional KBC inputs (consumed at MVP 5).
KBC_TABLE_B_TSV="${BASE}/results/catfish-variant-analysis/NN_kbc/kbc_variant_arrangement_assignments.tsv"
INVERSION_KARYOTYPE_TSV="${BASE}/results/inversion_atlas/karyotypes_pcangsd_k3.tsv"

# Damaging-variant tier (KBC §1.8): T1 | T2 | T3.
HPP_DAMAGING_TIER="${HPP_DAMAGING_TIER:-T1}"

# Bronze-segment policy: include | exclude (HANDOFF open question #2).
HPP_BRONZE_POLICY="${HPP_BRONZE_POLICY:-include}"

export HPP_OUTDIR STAGE3_DYAD_MAP_DIR STAGE3_TRIAD_MAP_DIR
export STAGE3_PARENT_PHASE_DIR VARIANT_MASTER_TSV JOINT_VCF REFERENCE_FASTA
export KBC_TABLE_B_TSV INVERSION_KARYOTYPE_TSV
export HPP_DAMAGING_TIER HPP_BRONZE_POLICY
