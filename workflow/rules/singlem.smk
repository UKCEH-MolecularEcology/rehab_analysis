"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-02-12]
Run: snakemake -s workflow/rules/singlem.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run singlem on raw reads
"""


############################################
rule singlem:
    input:
        os.path.join(RESULTS_DIR, "singlem/combined_singlem_otu.csv"),
        os.path.join(RESULTS_DIR, "singlem/combined_singlem_relab.csv")
    output:
        touch("status/singlem.done")


############################################
localrules: setup_singlem_db


############################################
# Download singleM database
rule setup_singlem_db:
    output:
        db=directory(os.path.join(DB_DIR, "singlem")), 
        dummy=os.path.join(DB_DIR, "singlem/db.done")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.singlem.db.log")
    conda:
        "singlem"
#        os.path.join(ENV_DIR, "singlem.yaml")
    message:
        "Setup: download singleM database"
    shell:
        "(date && mkdir -p {output.db} && "
        "singlem data --output-directory {output.db} && "
        "touch {output.dummy} && date) &> >(tee {log})"

rule run_singlem:
    input:
        in1=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR1"], 
        in2=lambda wildcards: SAMPLES.loc[wildcards.sid, "sR2"],
        dummy=os.path.join(DB_DIR, "singlem/db.done")
    output:
        profile=os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_profile.tsv"),
        table=os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_otu.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/singlem/{sid}.log")
    conda:
        "singlem"
#        os.path.join(ENV_DIR, "singlem.yaml")
    threads:
        config["singlem"]["threads"]
    params:
        metapackage=config["singlem"]["db"]
    message:
        "Running singlem on: {wildcards.sid}"
    shell:
        # --metapackage (not the SINGLEM_METAPACKAGE_PATH env var) is required:
        # the env var routes through singlem's "acquire default backpack"
        # path, which hard-enforces an exact metapackage schema version match
        # against this installed SingleM/zenodo_backpack build (5.4.0) --
        # rejecting both the old S3.2.1 db and this newer S6.5.0 one. The
        # --metapackage flag uses the given metapackage directly, no version
        # gate. See config["singlem"]["db"]'s comment for why this specific
        # metapackage copy.
        "(date && "
        "singlem pipe -1 {input[0]} -2 {input[1]} -p {output.profile} --otu-table {output.table} --threads {threads} --metapackage {params.metapackage} && "
        "date) &> >(tee {log})"

rule summarise_singlem:
    input:
        table=expand(os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_otu.csv"), sid=SAMPLES.index), 
        profile=expand(os.path.join(RESULTS_DIR, "singlem/{sid}_singlem_profile.tsv"), sid=SAMPLES.index)
    output:
        df_otu=os.path.join(RESULTS_DIR, "singlem/combined_singlem_otu.csv"),
        df_relab=os.path.join(RESULTS_DIR, "singlem/combined_singlem_relab.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/single/combine.log")
    conda:
         "singlem"
#        os.path.join(ENV_DIR, "singlem.yaml")
    message:
        "Combined all singlem outputs"
    shell:
        "(date && singlem summarise --input-otu-tables {input.table} --output-otu-table {output.df_otu} && "
        "singlem summarise --input-otu-tables {input.table} --input-taxonomic-profiles {input.profile} --output-species-by-site-relative-abundance {output.df_relab} && "
        "date) &> >(tee {log})"
