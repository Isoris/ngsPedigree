#!/usr/bin/env bash
###############################################################################
# SLURM_A07b_relatedness_per_chrom.sh — array job, one ngsRelate per RF
#
# Dispatched by STEP_A07b_relatedness_per_chrom.sh. One SLURM array task per
# region-file (chromosome). Body of the job is the same ngsRelate logic as
# STEP_A07's "1) ngsRelate" section, sliced to one RF.
#
# Reads:  chunk_rf.list (one line per RF; full path or basename)
# Writes: ${THIN_DIR}/06_relatedness/per_chrom/<rf_name>.res
#
# The frequency-file alignment logic is copied verbatim from STEP_A07 — same
# global MAF estimates, same alignment recipe, just sliced to one RF's BEAGLE.
###############################################################################
#SBATCH --job-name=ngsrel_per_chrom
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/ngsrel_per_chrom.%A_%a.out
#SBATCH --error=logs/ngsrel_per_chrom.%A_%a.err

set -euo pipefail
source "$(dirname "$0")/../config.sh"

CHUNK_LIST="${GLOBAL_DIR}/chunk_rf.list"
[[ -s "$CHUNK_LIST" ]] || { echo "[ERROR] Missing chunk_rf.list" >&2; exit 1; }

# Resolve RF name from the array task index. Handles both forms:
#   - full path:    /scratch/.../regions/LG01.rf  → strip dir + .rf
#   - bare token:   LG01                          → use as-is
RF_LINE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "${CHUNK_LIST}")
RF_NAME=$(basename "${RF_LINE}" .rf)

OUTDIR="${THIN_DIR}/06_relatedness/per_chrom"
mkdir -p "${OUTDIR}"

WORKDIR="${OUTDIR}/_work_${RF_NAME}"
mkdir -p "${WORKDIR}"

BEAGLE="${THIN_DIR}/04_beagle_byRF_majmin/thin_${RELATE_THIN}/${RF_NAME}.thin_${RELATE_THIN}.beagle.gz"
RES="${OUTDIR}/${RF_NAME}.res"
MAF_DIR="${GLOBAL_DIR}/02_snps"

[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing BEAGLE: $BEAGLE" >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "[ERROR] Missing SAMPLE_LIST: $SAMPLE_LIST" >&2; exit 1; }

if [[ -s "$RES" && "${FORCE:-0}" -eq 0 ]]; then
  echo "[$(timestamp)] [SKIP] ${RES} exists"
  rm -rf "${WORKDIR}"
  exit 0
fi

# ===========================================================================
# Build per-chromosome freq file (logic copied from STEP_A07's ngsRelate block)
# ===========================================================================
echo "[$(timestamp)] [${RF_NAME}] Building frequency file from global MAFs..."

zcat ${MAF_DIR}/catfish.*.mafs.gz \
  | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1"_"$2, $6}' \
  | sort -k1,1 > "${WORKDIR}/all_mafs_freq.tmp"

zcat "$BEAGLE" | tail -n +2 | cut -f1 > "${WORKDIR}/beagle_sites.tmp"
NSITES=$(wc -l < "${WORKDIR}/beagle_sites.tmp")

awk 'NR==FNR {freq[$1]=$2; next} {if($1 in freq) print freq[$1]; else print "NA"}' \
  "${WORKDIR}/all_mafs_freq.tmp" "${WORKDIR}/beagle_sites.tmp" \
  > "${WORKDIR}/freq_for_ngsrelate.txt"

# ===========================================================================
# Run ngsRelate (same flags as STEP_A07's WG run)
# ===========================================================================
echo "[$(timestamp)] [${RF_NAME}] Running ngsRelate (${N_SAMPLES} samples, ${NSITES} sites)..."

ngsRelate \
  -G "$BEAGLE" \
  -f "${WORKDIR}/freq_for_ngsrelate.txt" \
  -n "${N_SAMPLES}" \
  -O "$RES" \
  -p "${SLURM_CPUS_PER_TASK}" \
  -m 1 \
  -z "$SAMPLE_LIST"

# Cleanup
rm -rf "${WORKDIR}"

echo "[$(timestamp)] [${RF_NAME}] [DONE] ${RES} ($(wc -l < "$RES") rows)"
