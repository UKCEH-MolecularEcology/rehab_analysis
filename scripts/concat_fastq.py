import pandas as pd
import os
import subprocess
import re
import argparse


def concatenate_fastqs_by_sample(metadata_csv, fastq_dir, output_dir):
    """
    Concatenates _1.fastq.gz and _2.fastq.gz files for each sample.

    Reads metadata from a CSV file, identifies all fastq files belonging to each
    sample, and concatenates them into single sample-specific files.

    Args:
        metadata_csv (str): Path to the metadata CSV file.
        fastq_dir (str): Directory where the individual fastq.gz files are located.
        output_dir (str): Directory where the concatenated fastq.gz files
                          should be saved.
    """

    if not os.path.exists(fastq_dir):
        print(f"Error: FASTQ directory '{fastq_dir}' not found. Please ensure it exists and contains your FASTQ files.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    try:
        df = pd.read_csv(metadata_csv)
    except FileNotFoundError:
        print(f"Error: Metadata CSV file '{metadata_csv}' not found.")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: Metadata CSV file '{metadata_csv}' is empty.")
        return
    except Exception as e:
        print(f"Error reading metadata CSV: {e}")
        return

    # Ensure required columns exist
    if 'samplename' not in df.columns or 'fastq_ftp' not in df.columns:
        print("Error: Metadata CSV must contain 'samplename' and 'fastq_ftp' columns.")
        return

    # Group URLs by samplename
    sample_fastq_map = {}
    for index, row in df.iterrows():
        samplename = str(row['samplename'])  # Ensure samplename is string
        fastq_urls_str = row['fastq_ftp']

        if pd.isna(fastq_urls_str):  # Skip if fastq_ftp is NaN
            print(f"Warning: Skipping {samplename} due to missing fastq_ftp URL.")
            continue

        # Split multiple URLs by semicolon
        urls = [url.strip() for url in fastq_urls_str.split(';') if url.strip()]

        r1_files = []
        r2_files = []

        for url in urls:
            # Extract filename from the URL (handles both FTP path and simple filename)
            filename = os.path.basename(url)

            # Check for _1.fastq.gz and _2.fastq.gz pattern
            if re.search(r'_1\.fastq\.gz$', filename):
                r1_files.append(os.path.join(fastq_dir, filename))
            elif re.search(r'_2\.fastq\.gz$', filename):
                r2_files.append(os.path.join(fastq_dir, filename))
            else:
                print(f"Warning: File '{filename}' for sample '{samplename}' does not match _1.fastq.gz or _2.fastq.gz pattern. Skipping.")

        if samplename not in sample_fastq_map:
            sample_fastq_map[samplename] = {'r1': [], 'r2': []}

        sample_fastq_map[samplename]['r1'].extend(r1_files)
        sample_fastq_map[samplename]['r2'].extend(r2_files)

    # Perform concatenation for each sample
    for samplename, file_lists in sample_fastq_map.items():
        print(f"\nProcessing sample: {samplename}")

        # Process R1 files
        if file_lists['r1']:
            output_r1_path = os.path.join(output_dir, f"{samplename}_1.fastq.gz")
            existing_r1_files = [f for f in file_lists['r1'] if os.path.exists(f)]

            if not existing_r1_files:
                print(f"  No existing R1 files found on disk for {samplename}. Skipping R1 concatenation.")
            else:
                # SANITY CHECK: List the files being concatenated
                print(f"    Concatenating R1 files: {', '.join(existing_r1_files)}")
                try:
                    # Concatenate using subprocess.run with 'cat'
                    command = ['cat'] + existing_r1_files
                    with open(output_r1_path, 'wb') as outfile:  # 'wb' for binary write
                        subprocess.run(command, check=True, stdout=outfile)
                    print(f"  Successfully created {output_r1_path}")
                except subprocess.CalledProcessError as e:
                    print(f"  Error concatenating R1 files for {samplename}: {e}")
                except Exception as e:
                    print(f"  An unexpected error occurred during R1 concatenation for {samplename}: {e}")
        else:
            print(f"  No R1 files specified in metadata for {samplename}.")

        # Process R2 files
        if file_lists['r2']:
            output_r2_path = os.path.join(output_dir, f"{samplename}_2.fastq.gz")
            existing_r2_files = [f for f in file_lists['r2'] if os.path.exists(f)]

            if not existing_r2_files:
                print(f"  No existing R2 files found on disk for {samplename}. Skipping R2 concatenation.")
            else:
                # SANITY CHECK: List the files being concatenated
                print(f"    Concatenating R2 files: {', '.join(existing_r2_files)}")
                try:
                    command = ['cat'] + existing_r2_files
                    with open(output_r2_path, 'wb') as outfile:
                        subprocess.run(command, check=True, stdout=outfile)
                    print(f"  Successfully created {output_r2_path}")
                except subprocess.CalledProcessError as e:
                    print(f"  Error concatenating R2 files for {samplename}: {e}")
                except Exception as e:
                    print(f"  An unexpected error occurred during R2 concatenation for {samplename}: {e}")
        else:
            print(f"  No R2 files specified in metadata for {samplename}.")

    print("\nAll samples processed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate _1.fastq.gz and _2.fastq.gz files for each sample based on metadata.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-m', '--metadata-csv',
        type=str,
        required=True,
        help="Path to the metadata CSV file containing 'samplename' and 'fastq_ftp' columns."
    )
    parser.add_argument(
        '-f', '--fastq-dir',
        type=str,
        default="./fastq",
        help="Directory where the individual fastq.gz files are located."
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        default=".",
        help="Directory where the concatenated fastq.gz files should be saved."
    )

    args = parser.parse_args()

    concatenate_fastqs_by_sample(
        args.metadata_csv,
        args.fastq_dir,
        args.output_dir
    )
