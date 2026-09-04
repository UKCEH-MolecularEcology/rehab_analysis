#!/bin/bash
# Parallel FASTQ download for PRJEB34634 (REHAB) via fastq-dl (Singularity
# container), with integrity verification (gzip validity + size + md5sum
# against ENA's own read_run metadata report) and automatic retry of any
# run that fails verification.
#
# Usage: ./download_fastq.sh [--verify-only] <run_list_file> <outdir> [parallel_jobs] [metadata_tsv]
#   --verify-only : skip the initial bulk download entirely; just verify
#                   what's already on disk and re-download (with --force)
#                   only the runs that fail verification. Use this to fix
#                   up an already-mostly-complete download without
#                   re-attempting every run.
#   run_list_file : one ENA run accession per line (e.g. cut from metadata/rehab_metagenomes_list.csv)
#   outdir        : directory to download fastq.gz files into
#   parallel_jobs : number of concurrent fastq-dl / verify processes (default: 16)
#   metadata_tsv  : ENA read_run file report with fastq_bytes/fastq_md5 columns
#                   (default: data/PRJEB34634/metadata/PRJEB34634_read_run_metadata.tsv)
set -euo pipefail

VERIFY_ONLY=0
if [ "${1:-}" = "--verify-only" ]; then
    VERIFY_ONLY=1
    shift
fi

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_LIST="$1"
mkdir -p "$2"
OUTDIR="$(cd "$2" && pwd)"
JOBS="${3:-16}"
METADATA_TSV="${4:-${REPO_DIR}/data/PRJEB34634/metadata/PRJEB34634_read_run_metadata.tsv}"
MAX_RETRIES="${MAX_RETRIES:-5}"

TMPDIR="${REPO_DIR}/tmp"
SIF="${REPO_DIR}/containers/fastq-dl.sif"
LOGDIR="${REPO_DIR}/tmp/logs/download"

mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"

export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_CACHEDIR="${TMPDIR}/singularity_cache"
export TMPDIR="$TMPDIR"

download_one() {
    acc="$1"
    force_flag="${2:-}"
    singularity exec -B "${TMPDIR}:/tmp" -B "${OUTDIR}:${OUTDIR}" "${SIF}" \
        fastq-dl -a "$acc" -o "$OUTDIR" --cpus 2 ${force_flag} \
        > "${LOGDIR}/${acc}.log" 2>&1
}
export -f download_one
export TMPDIR SIF OUTDIR LOGDIR

if [ "$VERIFY_ONLY" -eq 1 ]; then
    echo "[download] --verify-only: skipping bulk download, checking what's already in ${OUTDIR}..."
else
    echo "[download] Downloading $(wc -l < "$RUN_LIST") runs into ${OUTDIR} (${JOBS} parallel)..."
    cat "$RUN_LIST" | xargs -P "$JOBS" -I{} bash -c 'download_one "$@"' _ {}
    echo "[download] Initial download complete. Verifying against ${METADATA_TSV}..."
fi

BAD_LIST="${TMPDIR}/bad_runs.txt"
attempt=0
while true; do
    if python3 "${REPO_DIR}/scripts/verify_downloads.py" "$METADATA_TSV" "$OUTDIR" "$BAD_LIST" --jobs "$JOBS"; then
        echo "[download] All runs verified OK (gzip + size + md5)."
        break
    fi

    n_bad=$(wc -l < "$BAD_LIST")
    attempt=$((attempt + 1))
    if [ "$attempt" -gt "$MAX_RETRIES" ]; then
        echo "[download] ERROR: ${n_bad} run(s) still failing verification after ${MAX_RETRIES} retries. See ${BAD_LIST}" >&2
        exit 1
    fi

    echo "[download] Attempt ${attempt}/${MAX_RETRIES}: re-downloading ${n_bad} failed run(s) with --force..."
    cat "$BAD_LIST" | xargs -P "$JOBS" -I{} bash -c 'download_one "$@" --force' _ {}
done
