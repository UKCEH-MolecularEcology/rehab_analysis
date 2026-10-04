import sys

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 3:
    print("Usage: python filter_contigs.py <input.fasta> <output.fasta>")
    sys.exit(1)

# Get input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

MIN_LEN = 1500

try:
    with open(input_file, "r") as in_f, open(output_file, "w") as out_f:
        header = None
        sequence = ""
        for line in in_f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # Flush the previous record, if any.
                if header is not None and len(sequence) >= MIN_LEN:
                    out_f.write(f"{header}\n{sequence}\n")
                header = line
                sequence = ""
            else:
                sequence += line.strip()
        # Flush the final record.
        if header is not None and len(sequence) >= MIN_LEN:
            out_f.write(f"{header}\n{sequence}\n")
    print("Processing complete")
except Exception as e:
    print(f"An error occurred: {e}")
