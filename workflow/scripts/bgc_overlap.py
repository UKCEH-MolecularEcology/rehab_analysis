#!/usr/bin/env python3
"""Match antiSMASH regions and GECCO clusters by MAG + contig + coordinate overlap (plain Python, no pandas).

bgc_overlap.tsv          one row per overlapping antiSMASH-GECCO pair, plus one row per region/cluster
                         without a partner; status = both | antismash_only | gecco_only
bgc_overlap_summary.tsv  counts per status and per antiSMASH product (with GECCO types of the matched clusters)
"""
import argparse, collections, csv, sys

p = argparse.ArgumentParser()
p.add_argument("--antismash", required=True)   # antismash_regions_all.tsv
p.add_argument("--gecco", required=True)       # gecco_clusters_all.tsv
p.add_argument("--out", required=True)
p.add_argument("--summary", required=True)
p.add_argument("--min_overlap_bp", type=int, default=1)
a = p.parse_args()

def read(path, need):
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in need if c not in (r.fieldnames or [])]
        if missing:
            sys.exit("{}: missing columns {} (found {})".format(path, missing, r.fieldnames))
        return list(r)

AS = read(a.antismash, ["Genome", "region", "contig", "start", "end", "products", "contig_edge"])
GC = read(a.gecco, ["Genome", "sequence_id", "cluster_id", "start", "end", "type", "average_p"])

# GECCO coordinates are 1-based inclusive; antiSMASH original coordinates 0-based half-open -> use half-open
for r in AS:
    r["_s"], r["_e"] = int(r["start"]), int(r["end"])
for r in GC:
    r["_s"], r["_e"] = int(r["start"]) - 1, int(r["end"])

gc_by = collections.defaultdict(list)
for r in GC:
    gc_by[(r["Genome"], r["sequence_id"])].append(r)

cols = ["status", "Genome", "contig",
        "antismash_region", "as_start", "as_end", "as_products", "as_contig_edge",
        "gecco_cluster", "gc_start", "gc_end", "gc_type", "gc_average_p",
        "overlap_bp", "overlap_frac_antismash", "overlap_frac_gecco"]
rows, gc_matched = [], set()
for r in AS:
    hits = []
    for g in gc_by.get((r["Genome"], r["contig"]), []):
        ov = min(r["_e"], g["_e"]) - max(r["_s"], g["_s"])
        if ov >= a.min_overlap_bp:
            hits.append((g, ov))
    base = {"Genome": r["Genome"], "contig": r["contig"], "antismash_region": r["region"],
            "as_start": r["start"], "as_end": r["end"], "as_products": r["products"],
            "as_contig_edge": r["contig_edge"]}
    if not hits:
        rows.append(dict(base, status="antismash_only"))
    for g, ov in hits:
        gc_matched.add(id(g))
        rows.append(dict(base, status="both", gecco_cluster=g["cluster_id"], gc_start=g["start"], gc_end=g["end"],
                         gc_type=g["type"], gc_average_p=g["average_p"], overlap_bp=ov,
                         overlap_frac_antismash=round(ov / max(1, r["_e"] - r["_s"]), 4),
                         overlap_frac_gecco=round(ov / max(1, g["_e"] - g["_s"]), 4)))
for g in GC:
    if id(g) not in gc_matched:
        rows.append({"status": "gecco_only", "Genome": g["Genome"], "contig": g["sequence_id"],
                     "gecco_cluster": g["cluster_id"], "gc_start": g["start"], "gc_end": g["end"],
                     "gc_type": g["type"], "gc_average_p": g["average_p"]})

with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="", extrasaction="ignore")
    w.writeheader(); w.writerows(rows)

# summary
as_ids_both = {r["antismash_region"] for r in rows if r["status"] == "both"}
gc_ids_both = {(r["Genome"], r["gecco_cluster"]) for r in rows if r["status"] == "both"}
n_as, n_gc = len(AS), len(GC)
by_product = collections.defaultdict(lambda: [0, 0, collections.Counter()])
for r in AS:
    k = r["products"] or "NA"
    by_product[k][0] += 1
    if r["region"] in as_ids_both:
        by_product[k][1] += 1
for r in rows:
    if r["status"] == "both":
        by_product[r["as_products"] or "NA"][2][r["gc_type"] or "NA"] += 1

with open(a.summary, "w") as fh:
    fh.write("section\tkey\tvalue\textra\n")
    fh.write("totals\tantismash_regions\t{}\t\n".format(n_as))
    fh.write("totals\tgecco_clusters\t{}\t\n".format(n_gc))
    fh.write("totals\tantismash_regions_with_gecco_overlap\t{}\t{:.1%}\n".format(len(as_ids_both), len(as_ids_both) / max(1, n_as)))
    fh.write("totals\tgecco_clusters_with_antismash_overlap\t{}\t{:.1%}\n".format(len(gc_ids_both), len(gc_ids_both) / max(1, n_gc)))
    for k, (n, nb, gt) in sorted(by_product.items(), key=lambda x: -x[1][0]):
        fh.write("antismash_product\t{}\t{}\t{} with GECCO overlap; GECCO types: {}\n".format(
            k, n, nb, ";".join("{}={}".format(t, c) for t, c in gt.most_common())))
print("{} antiSMASH regions, {} GECCO clusters, {} pairs/rows written".format(n_as, n_gc, len(rows)))
