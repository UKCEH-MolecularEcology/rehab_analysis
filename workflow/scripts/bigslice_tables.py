#!/usr/bin/env python3
"""Tables from a BiG-SLiCE (2.x) result database (latest clustering run):
  gcf_membership.tsv  one row per antiSMASH region: Genome, region, gcf, membership_value (distance to the GCF
                      centre; larger = more distant), on_contig_edge, length_nt, bigslice_class, + antiSMASH products
  gcf_summary.tsv     one row per GCF: n_regions, n_genomes, n_on_contig_edge, mean/min membership_value,
                      products, classes, representative (closest) region
  gcf_sample_rel_abund.tsv (optional) GCF x sample: summed relative abundance of the MAGs carrying the GCF
                      (each MAG counted once per GCF), from a CoverM MAG x sample matrix"""
import argparse, os, sqlite3
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--db", required=True)
p.add_argument("--regions", required=True, help="antismash_regions_all.tsv")
p.add_argument("--coverm", default="", help="CoverM MAG x sample relative abundance matrix (optional)")
p.add_argument("--outdir", required=True)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)

con = sqlite3.connect(a.db)
run_id = con.execute("SELECT max(run_id) FROM clustering").fetchone()[0]
if run_id is None:
    raise SystemExit("No clustering run found in " + a.db)
thr = con.execute("SELECT threshold FROM clustering WHERE run_id=?", (run_id,)).fetchone()[0]
print("BiG-SLiCE run", run_id, "threshold", thr)

mem = pd.read_sql_query("""
    SELECT b.id AS bgc_id, b.orig_folder, b.orig_filename, b.on_contig_edge, b.length_nt,
           g.id_in_run AS gcf, m.membership_value
    FROM gcf_membership m
    JOIN gcf g ON g.id = m.gcf_id
    JOIN clustering c ON c.id = g.clustering_id
    JOIN bgc b ON b.id = m.bgc_id
    WHERE m.rank = 0 AND c.run_id = ?""", con, params=(run_id,))
cls = pd.read_sql_query("""
    SELECT bc.bgc_id, cc.name || ':' || cs.name AS bigslice_class
    FROM bgc_class bc JOIN chem_subclass cs ON cs.id = bc.chem_subclass_id
    JOIN chem_class cc ON cc.id = cs.class_id""", con)
cls = cls.groupby("bgc_id")["bigslice_class"].apply(lambda s: ";".join(sorted(set(s)))).reset_index()

mem["Genome"] = mem["orig_folder"].str.rstrip("/").map(os.path.basename)
mem["region"] = mem["orig_filename"].str.replace(r"\.gbk$", "", regex=True)
mem["gcf"] = "GCF_" + mem["gcf"].astype(str)
mem = mem.merge(cls, on="bgc_id", how="left")

reg = pd.read_csv(a.regions, sep="\t")[["Genome", "region", "contig", "start", "end", "products"]]
mem = mem.merge(reg, on=["Genome", "region"], how="left")
cols = ["Genome", "region", "gcf", "membership_value", "on_contig_edge", "length_nt",
        "bigslice_class", "products", "contig", "start", "end"]
mem[cols].sort_values(["gcf", "membership_value"]).to_csv(os.path.join(a.outdir, "gcf_membership.tsv"), sep="\t", index=False)

join = lambda s: ";".join(sorted(set(";".join(s.dropna().astype(str)).split(";")) - {""}))
summ = mem.groupby("gcf").agg(
    n_regions=("region", "size"), n_genomes=("Genome", "nunique"),
    n_on_contig_edge=("on_contig_edge", lambda s: int(s.fillna(0).astype(int).sum())),
    mean_membership_value=("membership_value", "mean"), min_membership_value=("membership_value", "min"),
    products=("products", join), classes=("bigslice_class", join))
rep = mem.sort_values("membership_value").groupby("gcf").first()
summ["representative"] = rep["Genome"] + "/" + rep["region"]
summ.sort_values("n_genomes", ascending=False).to_csv(os.path.join(a.outdir, "gcf_summary.tsv"), sep="\t")
print("{} regions in {} GCFs from {} MAGs".format(len(mem), summ.shape[0], mem["Genome"].nunique()))

if a.coverm and os.path.exists(a.coverm):
    ab = pd.read_csv(a.coverm, sep="\t", index_col=0)
    pairs = mem[["gcf", "Genome"]].drop_duplicates()
    pairs = pairs[pairs["Genome"].isin(ab.index)]
    gcf_ab = ab.loc[pairs["Genome"]].set_axis(pairs["gcf"].values).groupby(level=0).sum()
    gcf_ab.index.name = "gcf"
    gcf_ab.to_csv(os.path.join(a.outdir, "gcf_sample_rel_abund.tsv"), sep="\t")
    print("GCF x sample matrix:", gcf_ab.shape)
else:
    print("CoverM matrix not found -- skipping GCF x sample abundance")
