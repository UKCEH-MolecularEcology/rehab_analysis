#!/bin/bash
# Usage: ./generate_sample_table.sh <input_dir> <output_file>
# Example: ./generate_sample_table.sh data/concatenated_fastq config/samples.tsv

INPUT_DIR="$1"
OUTPUT_FILE="$2"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

echo -e "Sample_ID\tsR1\tsR2" > "$OUTPUT_FILE"

for r1 in "$INPUT_DIR"/*_1.fastq.gz; do
    base=$(basename "$r1" _1.fastq.gz)
    r2="$INPUT_DIR/${base}_2.fastq.gz"

    if [[ -f "$r2" ]]; then
        echo -e "${base}\t${r1}\t${r2}" >> "$OUTPUT_FILE"
    else
        echo "Warning: Missing pair for ${r1}" >&2
    fi
done

echo "Sample table written to $OUTPUT_FILE"
