# SLURM template for producing the ngsRelate `.res` input

If the 226-cohort `.res` does not yet exist (or only the NAToRA-pruned 81 exists),
this is the canonical command to produce it on LANTA.

ngsPedigree consumes the **standard 23-column ngsRelate output**. The classifier
needs `theta`, `IBS0`, `KING`, `R0`, `R1`, `J9`, `nSites` at minimum. Default
ngsRelate output already contains all of these.

## Prerequisites

- Genotype likelihoods in BEAGLE format for all 226 samples (not the
  NAToRA-pruned 81 — we explicitly want the full cohort here, since
  pruning removes the very pairs we're trying to find).
- `samples.txt` with 226 sample IDs in the same order as the BEAGLE columns.
- `freqs.txt` — per-site allele frequency estimates (typically from ANGSD
  `-doMaf`).

## SLURM job (LANTA)

```bash
#!/bin/bash
#SBATCH --job-name=ngsrelate_226
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/ngsrelate_226.%j.out
#SBATCH --error=logs/ngsrelate_226.%j.err

set -euo pipefail

module load ngsRelate

RUN_ID="cohort_226_full_v1"
OUT_DIR="data/cohort/relatedness/ngsrelate/${RUN_ID}"
mkdir -p "${OUT_DIR}"

BEAGLE="path/to/cohort_226.beagle.gz"
FREQS="path/to/cohort_226.mafs.gz.freqs"   # extracted from .mafs.gz column 6
SAMPLES="path/to/samples_226.txt"
N=226

# Copy samples sidecar to run folder
cp "${SAMPLES}" "${OUT_DIR}/samples.txt"

# Run ngsRelate
ngsRelate \
    -G "${BEAGLE}" \
    -f "${FREQS}" \
    -n ${N} \
    -p ${SLURM_CPUS_PER_TASK} \
    -O "${OUT_DIR}/relatedness.res"

echo "ngsRelate done. Output:"
ls -la "${OUT_DIR}/"
head -1 "${OUT_DIR}/relatedness.res"
wc -l "${OUT_DIR}/relatedness.res"
```

## Expected output

- `relatedness.res`: TSV with header, ~25,425 rows for N=226 samples
  (which is N*(N-1)/2 pairs).
- 23 columns including `a, b, nSites, J7, J8, J9, rab, Fa, Fb, theta,
  inbreed_a, inbreed_b, 2of3_IDB, FDiff, loglh, nIter, coverage, IBS0,
  IBS1, IBS2, R0, R1, KING`.
- Run time: typically 1–3 hours on 32 cores at 9× coverage.

## Then run ngsPedigree

```bash
python scripts/STEP_PED_01_annotate_relationships.py \
    --res "${OUT_DIR}/relatedness.res" \
    --samples "${OUT_DIR}/samples.txt" \
    --outdir "${OUT_DIR}/ngspedigree/" \
    --run-id "${RUN_ID}"
```

## Sanity checks before running ngsPedigree

```bash
# Confirm the column set
head -1 "${OUT_DIR}/relatedness.res"
# Should include: theta, IBS0, KING, R0, R1, J9 — at minimum

# Confirm the row count
wc -l "${OUT_DIR}/relatedness.res"
# Should be 25,426 (header + 25,425 pairs) for N=226

# Confirm samples sidecar matches indices
wc -l "${OUT_DIR}/samples.txt"
# Should be 226
```

If any of these fail, ngsPedigree will refuse to run and tell you what's wrong.
