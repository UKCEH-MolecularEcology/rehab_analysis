"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/bin_taxqual.smk --use-conda --cores 4 -rp
Latest modification: [2026-10-04] collect_bins removed -- rules.galah.output
    already represents the fully-dereplicated 7-binner ensemble (see
    rules/dereplicate.smk), so there's nothing left to merge in on top of
    it (the old version bolted raw, non-deduplicated concoct_bins back on
    top of dRep's output, which no longer makes sense now that concoct is
    one of 7 binners already going through DAS_Tool -> Rosella refine ->
    Galah, not merged in afterwards).
Purpose: Taxonomy and quality of bins
"""


############################################
rule taxqual:
    input:
        os.path.join(RESULTS_DIR, "bins/gtdbtk_final"),
        os.path.join(RESULTS_DIR, "bins/checkm2/quality_report.tsv")
    output:
        touch("status/bin_taxqual.done")


############################################
# localrules:


############################################
# GTDBTK taxonomy
rule gtdbtk:
    input:
        rules.galah.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "bins/gtdbtk_final"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB on MAGs"
    shell:
        "(date && "
        "export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output} --skip_ani_screen && "
        "date) &> >(tee {log})"

# checkm_db (downloads the CheckM2 database) now lives in
# rules/dereplicate.smk -- it's needed earlier in file-processing order
# there too (checkm2_dastool/checkm2_refine), and dereplicate.smk loads
# before this file in workflow/Snakefile's "binning" STEPS block.

# Checking bin quality
# NOTE: previously called a submodules/checkm2/bin/checkm2 binary path that
# was never actually populated (submodules/checkm2 is listed in
# .gitmodules but was never `git submodule add`-ed -- same class of bug as
# an earlier MagicLamp submodule issue fixed in this repo). This rule
# would have failed the moment it actually ran. Fixed to call `checkm2`
# directly via the conda env's own PATH, matching
# rules/dereplicate.smk's checkm2_dastool/checkm2_refine.
rule checkm_final:
    input:
        drep=rules.galah.output.final,
        db=rules.checkm_db.output[0]
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2/quality_report.tsv")
    conda:
        "/hdd0/susbus/tools/conda_envs/checkm2"
    log:
        os.path.join(RESULTS_DIR, "logs/checkm/checkm.out.log")
    threads:
        config["checkm"]["threads"]
    params:
        ext=config["checkm"]["extension"],
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    message:
        "Running Final Checkm on dereplicated output"
    shell:
        "(date && mkdir -p $(dirname $(dirname {output.tsv}))/checkm2_tmp && chmod -R 777 $(dirname $(dirname {output.tsv}))/checkm2_tmp && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --tmpdir $(dirname $(dirname {output.tsv}))/checkm2_tmp --threads {threads} -x {params.ext} --input {input.drep} --output-directory $(dirname {output.tsv}) --force && "
        "date) &> >(tee {log})"
