#!/bin/bash
# Parallel FASTQ download for PRJEB34634 (REHAB) via fastq-dl (Singularity container).
#
# Usage: ./download_fastq.sh <run_list_file> <outdir> [parallel_jobs]
#   run_list_file : one ENA run accession per line (e.g. cut from metadata/rehab_metagenomes_list.csv)
#   outdir        : directory to download fastq.gz files into
#   parallel_jobs : number of concurrent fastq-dl processes (default: 16)
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_LIST="$1"
mkdir -p "$2"
OUTDIR="$(cd "$2" && pwd)"
JOBS="${3:-16}"

TMPDIR="${REPO_DIR}/tmp"
SIF="${REPO_DIR}/containers/fastq-dl.sif"
LOGDIR="${REPO_DIR}/tmp/logs/download"

mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR="${TMPDIR}/singularity_cache"
export TMPDIR="$TMPDIR"

download_one() {
    acc="$1"
    singularity exec -B "${TMPDIR}:/tmp" -B "${OUTDIR}:${OUTDIR}" "${SIF}" \
        fastq-dl -a "$acc" -o "$OUTDIR" --cpus 2 \
        > "${LOGDIR}/${acc}.log" 2>&1
}
export -f download_one
export TMPDIR SIF OUTDIR LOGDIR

cat "$RUN_LIST" | xargs -P "$JOBS" -I{} bash -c 'download_one "$@"' _ {}
