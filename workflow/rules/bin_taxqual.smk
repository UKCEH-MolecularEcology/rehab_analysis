"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/bin_taxqual.smk --use-conda --cores 4 -rp
Latest modification:
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
rule collect_bins:
    input:
        all=os.path.join(RESULTS_DIR, "bins/finalbins"),
        semibin=os.path.join(RESULTS_DIR, "bins/semibin/finalbins"),
        concoct=os.path.join(RESULTS_DIR,"bins/concoct_bins")
    output:
        bins=directory(os.path.join(RESULTS_DIR, "bins/all_finalbins"))
    message:
        "putting all the bins in one place"
    shell:
        "(date && mkdir {output.bins} && "
        "cp -v {input.all}/*.fa {output.bins}/. && "
        "cp -v {input.semibin}/*.fa {output.bins}/. && "
        "cp -v {input.concoct}/*.fa {output.bins}/. && "
        "date) &> >(tee {log})"

# GTDBTK taxonomy
rule gtdbtk:
    input:
        rules.collect_bins.output.bins
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

# install checkm database
rule checkm_db:
    output:
        os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_db.log")
    conda:
        "checkm2"
#        os.path.join(ENV_DIR, "checkm.yaml")
    message:
        "Downloading the checkm2 database"
    shell:
        "(date && checkm2 database --download --path $(dirname $(dirname {output})) && date) &> >(tee {log})"
   
# Checking bin quality
rule checkm_final:
    input:
        drep=rules.collect_bins.output.bins,
        db=rules.checkm_db.output[0]
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2/quality_report.tsv")
    conda:
        "/hdd0/susbus/tools/conda_envs/checkm2" # "checkm2"
#        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm/checkm.out.log")
    threads:
        config["checkm"]["threads"]
    params:
        ext=config["checkm"]["extension"],
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd"),
        checkm2=os.path.join(SUBMODULES, "checkm2/bin/checkm2")
    message:
        "Running Final Checkm on dereplicated output"
    shell:
        "(date && mkdir -p $(dirname $(dirname {output.tsv}))/checkm2_tmp && chmod -R 777 $(dirname $(dirname {output.tsv}))/checkm2_tmp && "
        "export CHECKM2DB={params.db} && "
        "{params.checkm2} predict --tmpdir $(dirname $(dirname {output.tsv}))/checkm2_tmp  --threads {threads} -x {params.ext} --input {input.drep} --output-directory $(dirname {output.tsv}) --force && "
        "date) &> >(tee {log})"
