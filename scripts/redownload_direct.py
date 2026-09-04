#!/usr/bin/env python3
"""
Re-download specific ENA runs directly via their fastq_ftp URLs (from the
ENA read_run metadata report), verifying md5 immediately after each file.

Bypasses fastq-dl's own ENA/SRA API lookup entirely -- useful when that
lookup is flaky ("not found on any provider") for runs whose exact FTP
location and checksum we already have cached locally.

Usage:
    redownload_direct.py <metadata_tsv> <fastq_dir> <accessions_file> [--jobs N] [--max-attempts N]
"""
import argparse
import csv
import hashlib
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


def md5sum(path, chunk_size=8 * 1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one_file(url, dest, expected_md5, max_attempts):
    if not url.startswith(("http://", "https://", "ftp://")):
        url = "https://" + url
    for attempt in range(1, max_attempts + 1):
        tmp_dest = dest + ".part"
        result = subprocess.run(
            ["curl", "-sS", "-f", "-L", "--retry", "3", "--retry-delay", "5", "-o", tmp_dest, url],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"  attempt {attempt}: curl failed for {url}: {result.stderr.strip()}", file=sys.stderr)
            continue
        if expected_md5:
            actual_md5 = md5sum(tmp_dest)
            if actual_md5 != expected_md5:
                print(f"  attempt {attempt}: md5 mismatch for {dest} (expected {expected_md5}, got {actual_md5})", file=sys.stderr)
                os.remove(tmp_dest)
                continue
        os.replace(tmp_dest, dest)
        return True
    return False


def redownload_run(acc, fastq_dir, fastq_ftp, fastq_md5, max_attempts):
    urls = fastq_ftp.split(";") if fastq_ftp else []
    md5s = fastq_md5.split(";") if fastq_md5 else []
    while len(md5s) < len(urls):
        md5s.append(None)

    if len(urls) != 2:
        return acc, False, f"expected 2 fastq_ftp URLs, got {len(urls)}"

    for i, (url, md5) in enumerate(zip(urls, md5s), start=1):
        dest = os.path.join(fastq_dir, f"{acc}_{i}.fastq.gz")
        ok = download_one_file(url, dest, md5, max_attempts)
        if not ok:
            return acc, False, f"R{i} failed after {max_attempts} attempts"
    return acc, True, ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata_tsv")
    parser.add_argument("fastq_dir")
    parser.add_argument("accessions_file")
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--max-attempts", type=int, default=3)
    args = parser.parse_args()

    with open(args.accessions_file) as fh:
        wanted = {line.strip() for line in fh if line.strip()}

    lookup = {}
    with open(args.metadata_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["run_accession"] in wanted:
                lookup[row["run_accession"]] = (row.get("fastq_ftp", ""), row.get("fastq_md5", ""))

    missing = wanted - lookup.keys()
    if missing:
        print(f"[redownload] WARNING: {len(missing)} accession(s) not found in metadata TSV: {sorted(missing)}", file=sys.stderr)

    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(redownload_run, acc, args.fastq_dir, ftp, md5, args.max_attempts): acc
            for acc, (ftp, md5) in lookup.items()
        }
        for future in as_completed(futures):
            acc, ok, reason = future.result()
            if ok:
                print(f"[redownload] OK: {acc}")
            else:
                print(f"[redownload] FAILED: {acc}: {reason}", file=sys.stderr)
                failed.append(acc)

    if failed:
        print(f"[redownload] {len(failed)}/{len(lookup)} runs still failed: {sorted(failed)}", file=sys.stderr)
        sys.exit(1)
    print(f"[redownload] All {len(lookup)} runs downloaded and md5-verified successfully.")


if __name__ == "__main__":
    main()
