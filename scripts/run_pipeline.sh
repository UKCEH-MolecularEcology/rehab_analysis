#!/bin/bash
# Run the REHAB Snakemake pipeline locally (no SLURM on this host).
#
# Usage: ./run_pipeline.sh [snakemake extra args...]
# Example: ./run_pipeline.sh -n              # dry-run
#          ./run_pipeline.sh                 # full run
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CORES="${SMK_CORES:-40}"

export TMPDIR="${REPO_DIR}/tmp"
mkdir -p "$TMPDIR"

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake -s "${REPO_DIR}/workflow/Snakefile" \
    --configfile "${REPO_DIR}/config/config.yaml" \
    --use-conda --conda-prefix /prj/DECODE/ea_biofilm_results/conda_envs --conda-frontend conda \
    --rerun-incomplete \
    --cores "${CORES}" -rp "$@"
