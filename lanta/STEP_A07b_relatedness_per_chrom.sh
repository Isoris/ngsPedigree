#!/usr/bin/env bash
###############################################################################
# STEP_A07b_relatedness_per_chrom.sh — per-chromosome ngsRelate dispatch
#
# Sibling to the existing STEP_A07_relatedness.sh (whole-genome). This script
# does the SAME ngsRelate workflow but produces one .res per region-file (RF),
# which in this pipeline corresponds to one .res per chromosome.
#
# What this dispatcher does:
#   - reads chunk_rf.list (the same list STEP_A04 used to make per-RF BEAGLEs)
#   - dispatches a SLURM array job with one task per RF
#
# What it does NOT do (intentionally — these belong to the WG run only):
#   - NAToRA culling
#   - greedy first-degree pruning
#   - 3-panel relatedness figure
#   - .results file consumed by structure_pruned
#
# Outputs land in: ${THIN_DIR}/06_relatedness/per_chrom/<rf_name>.res
#
# Consumes the SAME inputs as STEP_A07:
#   - per-RF BEAGLEs from STEP_A04 at:
#     ${THIN_DIR}/04_beagle_byRF_majmin/thin_${RELATE_THIN}/<rf>.thin_${RELATE_THIN}.beagle.gz
#   - global MAF estimates from ${GLOBAL_DIR}/02_snps/catfish.*.mafs.gz
#   - the same SAMPLE_LIST as STEP_A07
#
# Frequencies are pulled from the GLOBAL MAF estimates and aligned to each
# per-RF BEAGLE's site order separately. Allele frequencies are population-
# level and do not change per chromosome; only the alignment is per-chromosome.
# Re-estimating MAFs per chromosome would be circular (using the same data
# to estimate freqs and run relatedness).
#
# Called by: bash steps/STEP_A07b_relatedness_per_chrom.sh
# Then run:  sbatch --array=1-${N_RF}%8 slurm/SLURM_A07b_relatedness_per_chrom.sh
#
# Adapted from: STEP_A07_relatedness.sh (Quentin Andres, MS_Inversions pipeline)
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTDIR="${THIN_DIR}/06_relatedness/per_chrom"
mkdir -p "${OUTDIR}"

CHUNK_LIST="${GLOBAL_DIR}/chunk_rf.list"
[[ -s "$CHUNK_LIST" ]] || { echo "[ERROR] Missing chunk_rf.list: $CHUNK_LIST" >&2; exit 1; }
N_RF=$(wc -l < "${CHUNK_LIST}")

ARGFILE="${OUTDIR}/06_relatedness_per_chrom.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t06_relatedness_per_chrom"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "chunk_rf_list\t${CHUNK_LIST}"
  echo -e "n_rf\t${N_RF}"
  echo -e "sample_list\t${SAMPLE_LIST}"
  echo -e "n_samples\t${N_SAMPLES}"
  echo -e "relate_thin\t${RELATE_THIN}"
  echo -e "outdir\t${OUTDIR}"
} > "$ARGFILE"

echo "[$(timestamp)] Per-chromosome ngsRelate dispatch"
echo "[$(timestamp)]   ${N_RF} region files in ${CHUNK_LIST}"
echo "[$(timestamp)]   output: ${OUTDIR}/<rf_name>.res"
echo ""
echo "  sbatch --array=1-${N_RF}%8 \$(dirname \$0)/../slurm/SLURM_A07b_relatedness_per_chrom.sh"
echo ""
echo "When all array tasks complete, run ngsPedigree Stage 2:"
echo "  python scripts/STEP_PED_02_per_chromosome_qc.py \\"
echo "      --per-chrom-dir ${OUTDIR} \\"
echo "      --stage1-outdir <ngsPedigree Stage 1 output dir> \\"
echo "      --outdir <Stage 2 output dir>"
