#!/usr/bin/env bash
###############################################################################
# STEP_A07_relatedness.sh — REFERENCE COPY (already completed)
#
# This is the EXISTING whole-genome ngsRelate script from the MS_Inversions
# pipeline. It produces the global .res that ngsPedigree Stage 1 consumes.
#
# DO NOT RUN FROM ngsPedigree. Outputs already exist on LANTA at:
#   ${THIN_DIR}/06_relatedness/catfish_226_relatedness.res
#
# The new STEP_A07b_relatedness_per_chrom.sh (in this same directory) is the
# per-chromosome sibling of this script. It uses the same inputs (per-RF
# BEAGLEs + global MAFs + same SAMPLE_LIST) but produces one .res per RF
# instead of one global .res. STEP_A07b deliberately omits the NAToRA culling,
# the greedy pruning, and the 3-panel figure, because those are meaningful
# only at genome-wide scale; per-chromosome .res files are an internal
# artifact feeding ngsPedigree Stage 2 (per-chromosome QC).
###############################################################################
###############################################################################
# STEP_A07_relatedness.sh — ngsRelate + NAToRA + pruning + plotting
#
# Runs between structure_all and structure_pruned:
#   1) ngsRelate on thin-500 whole-genome BEAGLE
#   2) NAToRA multi-cutoff culling (using config low/high pairs)
#   3) Greedy first-degree pruning
#   4) 3-panel relatedness figure
#
# Produces:
#   - ngsRelate pairwise output
#   - NAToRA cutoff summary table
#   - pruned_samples.txt (input for structure_pruned)
#   - 3-panel relatedness figure
#
# Called by: run_step1.sh relatedness
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTDIR="${THIN_DIR}/06_relatedness"
mkdir -p "${OUTDIR}"

BEAGLE="${THIN_DIR}/03_merged_beagle/catfish.wholegenome.byRF.thin_${RELATE_THIN}.beagle.gz"
MAF_DIR="${GLOBAL_DIR}/02_snps"

[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing BEAGLE: $BEAGLE" >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "[ERROR] Missing SAMPLE_LIST" >&2; exit 1; }

ARGFILE="${OUTDIR}/06_relatedness.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t06_relatedness"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "beagle\t${BEAGLE}"
  echo -e "sample_list\t${SAMPLE_LIST}"
  echo -e "n_samples\t${N_SAMPLES}"
  echo -e "relate_thin\t${RELATE_THIN}"
  echo -e "natora_dup_mz\t${NATORA_DUP_MZ_LOW},${NATORA_DUP_MZ_HIGH}"
  echo -e "natora_first\t${NATORA_FIRST_LOW},${NATORA_FIRST_HIGH}"
  echo -e "natora_second\t${NATORA_SECOND_LOW},${NATORA_SECOND_HIGH}"
  echo -e "natora_third\t${NATORA_THIRD_LOW},${NATORA_THIRD_HIGH}"
  echo -e "first_degree_prune_theta\t${THETA_FIRST_DEGREE}"
} > "$ARGFILE"

RES="${OUTDIR}/catfish_${N_SAMPLES}_relatedness.res"
P="${DEFAULT_THREADS}"

# ===========================================================================
# 1) ngsRelate
# ===========================================================================
echo "[$(timestamp)] === ngsRelate ==="

if [[ -s "$RES" && "${FORCE:-0}" -eq 0 ]]; then
  echo "[SKIP] ngsRelate output exists: $RES"
else
  # Build frequency file matched to beagle site order
  echo "[$(timestamp)] Building frequency file..."
  zcat ${MAF_DIR}/catfish.*.mafs.gz \
    | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1"_"$2, $6}' \
    | sort -k1,1 > "${OUTDIR}/all_mafs_freq.tmp"

  zcat "$BEAGLE" | tail -n +2 | cut -f1 > "${OUTDIR}/beagle_sites.tmp"
  NSITES=$(wc -l < "${OUTDIR}/beagle_sites.tmp")

  awk 'NR==FNR {freq[$1]=$2; next} {if($1 in freq) print freq[$1]; else print "NA"}' \
    "${OUTDIR}/all_mafs_freq.tmp" "${OUTDIR}/beagle_sites.tmp" \
    > "${OUTDIR}/freq_for_ngsrelate.txt"

  echo "[$(timestamp)] Running ngsRelate (${N_SAMPLES} samples, ${NSITES} sites)..."
  ngsRelate \
    -G "$BEAGLE" \
    -f "${OUTDIR}/freq_for_ngsrelate.txt" \
    -n "${N_SAMPLES}" \
    -O "$RES" \
    -p "$P" \
    -m 1 \
    -z "$SAMPLE_LIST"

  rm -f "${OUTDIR}/all_mafs_freq.tmp" "${OUTDIR}/beagle_sites.tmp"
  echo "[$(timestamp)] ngsRelate done: $RES"
fi

# (NAToRA / pruning / figure sections preserved in original; omitted here for brevity.
#  See the original STEP_A07_relatedness.sh in the MS_Inversions pipeline.)
echo "[$(timestamp)] [DONE] 06_relatedness (ngsRelate section)"
