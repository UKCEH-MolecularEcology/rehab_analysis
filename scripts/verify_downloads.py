#!/usr/bin/env python3
"""
Verify downloaded ENA FASTQ files against the ENA read_run metadata report:
checks R1 and R2 both exist, are valid gzip, and match the expected md5sum
(and byte size as a cheap pre-check before the more expensive md5).

Usage:
    verify_downloads.py <metadata_tsv> <fastq_dir> <bad_list_output> [--jobs N]

<metadata_tsv> is an ENA portal API read_run file report (must have
run_accession, fastq_bytes, fastq_md5 columns -- see
data/PRJEB34634/metadata/PRJEB34634_read_run_metadata.tsv).

Writes one bad run_accession per line to <bad_list_output>, plus a reason
per accession (mismatched/missing/corrupt) to stderr.
"""
import argparse
import csv
import gzip
import hashlib
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


def md5sum(path, chunk_size=8 * 1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(8 * 1024 * 1024):
                pass
        return True
    except Exception:
        return False


def verify_one(run_accession, fastq_dir, expected_bytes, expected_md5s):
    """Returns (run_accession, ok, reason)."""
    reasons = []
    for i, (exp_bytes, exp_md5) in enumerate(zip(expected_bytes, expected_md5s), start=1):
        path = os.path.join(fastq_dir, f"{run_accession}_{i}.fastq.gz")
        if not os.path.isfile(path):
            reasons.append(f"R{i} missing")
            continue
        actual_bytes = os.path.getsize(path)
        if exp_bytes and actual_bytes != exp_bytes:
            reasons.append(f"R{i} size mismatch (expected {exp_bytes}, got {actual_bytes})")
            continue
        if not gzip_ok(path):
            reasons.append(f"R{i} corrupt gzip")
            continue
        if exp_md5:
            actual_md5 = md5sum(path)
            if actual_md5 != exp_md5:
                reasons.append(f"R{i} md5 mismatch (expected {exp_md5}, got {actual_md5})")
    if reasons:
        return run_accession, False, "; ".join(reasons)
    return run_accession, True, ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata_tsv")
    parser.add_argument("fastq_dir")
    parser.add_argument("bad_list_output")
    parser.add_argument("--jobs", type=int, default=16)
    args = parser.parse_args()

    runs = []
    with open(args.metadata_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            acc = row["run_accession"]
            bytes_field = row.get("fastq_bytes", "")
            md5_field = row.get("fastq_md5", "")
            expected_bytes = [int(b) if b else None for b in bytes_field.split(";")] if bytes_field else [None, None]
            expected_md5s = md5_field.split(";") if md5_field else [None, None]
            # pad to 2 (some rows could be single-end, but this dataset is paired-end)
            while len(expected_bytes) < 2:
                expected_bytes.append(None)
            while len(expected_md5s) < 2:
                expected_md5s.append(None)
            runs.append((acc, expected_bytes[:2], expected_md5s[:2]))

    total = len(runs)
    bad = []
    checked = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(verify_one, acc, args.fastq_dir, exp_bytes, exp_md5s): acc
            for acc, exp_bytes, exp_md5s in runs
        }
        for future in as_completed(futures):
            acc, ok, reason = future.result()
            checked += 1
            if not ok:
                bad.append(acc)
                print(f"BAD: {acc}: {reason}", file=sys.stderr)
            if checked % 200 == 0:
                print(f"[verify] {checked}/{total} checked, {len(bad)} bad so far", file=sys.stderr)

    with open(args.bad_list_output, "w") as out:
        for acc in sorted(bad):
            out.write(acc + "\n")

    print(f"[verify] done: {total} checked, {len(bad)} bad. List written to {args.bad_list_output}", file=sys.stderr)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
