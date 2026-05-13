#!/usr/bin/env bash
###############################################################################
# STEP_A05_merge_beagles.sh — REFERENCE COPY (already completed)
#
# Verbatim copy from the MS_Inversions pipeline. Merges per-RF BEAGLEs into
# whole-genome BEAGLE. Outputs already exist on LANTA.
#
# DO NOT RUN FROM ngsPedigree. The whole-genome BEAGLE is what STEP_A07
# (existing) consumes to produce the global .res that ngsPedigree Stage 1
# uses. The per-RF inputs to this merge are what STEP_A07b consumes for
# Stage 2 — i.e. STEP_A07b skips the merge, since each chromosome's BEAGLE
# is already the input it needs.
###############################################################################
###############################################################################
# STEP_A05_merge_beagles.sh — Merge per-RF BEAGLEs into whole-genome
# Called by: 04_beagles.sh (after all BEAGLE array jobs finish)
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

BYRF_BASE="${THIN_DIR}/04_beagle_byRF_majmin"
MERGED_DIR="${THIN_DIR}/03_merged_beagle"
mkdir -p "${MERGED_DIR}"

THIN_LIST=("${THIN_FINE[@]}")

for W in "${THIN_LIST[@]}"; do
  INDIR="${BYRF_BASE}/thin_${W}"
  [[ -d "$INDIR" ]] || { echo "[WARN] Missing $INDIR"; continue; }

  LIST="${MERGED_DIR}/beagle_thin_${W}.list"
  OUT="${MERGED_DIR}/catfish.wholegenome.byRF.thin_${W}.beagle.gz"
  TMP="${MERGED_DIR}/.tmp_thin_${W}.beagle"

  find "${INDIR}" -name "*.thin_${W}.beagle.gz" | sort -V > "$LIST"
  N=$(wc -l < "$LIST")
  (( N > 0 )) || { echo "[WARN] No beagle.gz for thin_${W}"; continue; }

  FIRST=$(head -1 "$LIST")
  zcat "$FIRST" | head -1 > "$TMP"
  while read -r F; do zcat "$F" | tail -n +2 >> "$TMP"; done < "$LIST"
  gzip -c "$TMP" > "$OUT"; rm -f "$TMP"
  gzip -t "$OUT"

  echo "[OK] thin_${W}: $N files -> $OUT ($(zcat "$OUT" | wc -l) lines)"
done
