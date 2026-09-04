import sys

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 3:
    print("Usage: python filter_contigs.py <input.fasta> <output.fasta>")
    sys.exit(1)

# Get input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # Open input and output files
    with open(input_file, "r") as in_f, open(output_file, "w") as out_f:
        # Initialize variables to track header and sequence
        header = ""
        sequence = ""
        # Iterate over each line in the input file
        for line in in_f:
            # Strip whitespace from the line
            line = line.strip()
            # Check if the line is a header line
            if line.startswith(">s"):
                # Write the previous sequence to the output file
                if len(sequence) >= 1500:
                    out_f.write(f">{header}\n{sequence}\n")
                # Reset the header and sequence variables
                header = line
                sequence = ""
            else:
                # Append the line to the current sequence
                sequence += line
        # Write the last sequence to the output file
        if len(sequence) >= 1500:
            out_f.write(f">{header}\n{sequence}\n")
    print("Processing complete")
except Exception as e:
    print(f"An error occurred: {e}")

