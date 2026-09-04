#!/bin/bash
# Waits for the running download_fastq.sh (PID given as $1) to finish, then
# automatically runs concatenation and sample-manifest generation.
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WAIT_PID="$1"

echo "$(date) Waiting for download process PID ${WAIT_PID} to finish..."
while kill -0 "$WAIT_PID" 2>/dev/null; do
    sleep 60
done
echo "$(date) Download process finished (or was not found). Verifying completeness..."

EXPECTED=$(tail -n +2 "${REPO_DIR}/metadata/rehab_metagenomes_list.csv" | cut -d',' -f4 | sort -u | wc -l)
GOT=$(ls "${REPO_DIR}/data/PRJEB34634/fastq"/*_1.fastq.gz 2>/dev/null | wc -l)
echo "$(date) Expected ${EXPECTED} runs, found ${GOT} R1 files downloaded."

echo "$(date) Running concat_fastq.py..."
python "${REPO_DIR}/scripts/concat_fastq.py" \
    -m "${REPO_DIR}/metadata/rehab_metagenomes_list.csv" \
    -f "${REPO_DIR}/data/PRJEB34634/fastq" \
    -o "${REPO_DIR}/data/concatenated_fastq"
echo "$(date) concat_fastq.py exit code: $?"

echo "$(date) Running generate_sample_table.sh..."
"${REPO_DIR}/scripts/generate_sample_table.sh" \
    "${REPO_DIR}/data/concatenated_fastq" \
    "${REPO_DIR}/config/samples.tsv"
echo "$(date) generate_sample_table.sh exit code: $?"

echo "$(date) Chain complete."
