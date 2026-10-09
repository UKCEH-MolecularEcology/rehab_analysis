#!/usr/bin/env python3
"""One row per antiSMASH region of a MAG, from its <mag>.regionNNN.gbk files.
Columns: Genome, region, contig, start, end, length, products, contig_edge"""
import argparse, glob, os, re
from Bio import SeqIO

p = argparse.ArgumentParser()
p.add_argument("--mag", required=True); p.add_argument("--dir", required=True); p.add_argument("--out", required=True)
a = p.parse_args()

cols = ["Genome", "region", "contig", "start", "end", "length", "products", "contig_edge"]
rows = []
for gbk in sorted(glob.glob(os.path.join(a.dir, "*.region*.gbk"))):
    region = re.sub(r"\.gbk$", "", os.path.basename(gbk))
    for rec in SeqIO.parse(gbk, "genbank"):
        meta = rec.annotations.get("structured_comment", {}).get("antiSMASH-Data", {})
        for f in rec.features:
            if f.type != "region":
                continue
            start = int(meta.get("Orig. start", int(f.location.start)))
            end = int(meta.get("Orig. end", int(f.location.end)))
            rows.append([a.mag, region, rec.id, start, end, end - start,
                         ";".join(f.qualifiers.get("product", [])),
                         ";".join(f.qualifiers.get("contig_edge", []))])
with open(a.out, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
print(a.mag, len(rows), "regions")
