"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-10-20]
Run: snakemake -s workflow/rules/assembly.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run SYLPH on reads to profile prokaryotes, viruses and fungi
"""


############################################
rule sylph:
    input:
        expand(os.path.join(RESULTS_DIR, "sylph/profiling_{profiles}.tsv"), profiles=["prokaryotes", "viruses", "fungi"])
    output:
        touch("status/sylph.done")


############################################
localrules: download_sylph_db


############################################
# Download SYLPH DBs
rule download_sylph_db:
    output:
        gtdb=os.path.join(DB_DIR, "sylph/gtdb-r220-c200-dbv1.syldb"),
        imgvr=os.path.join(DB_DIR, "sylph/imgvr_c200_v0.3.0.syldb"),
        fundb=os.path.join(DB_DIR, "sylph/fungi-refseq-2024-07-25-c200-v0.3.syldb")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.sylph.db.log")
    params:
        gtdb_url=config["sylph"]["gtdb_url"],
        imgvr_url=config["sylph"]["imgvr_url"],
        fundb_url=config["sylph"]["fundb_url"]
    message:
        "Setup: download SYLPH databases"
    shell:
        "(date && "
        "wget -O {output.gtdb} {params.gtdb_url} --no-check-certificate && "
        "wget -O {output.imgvr} {params.imgvr_url} --no-check-certificate && "
        "wget -O {output.fundb} {params.fundb_url} --no-check-certificate && "
        "date) &> >(tee {log})"

# Running on prokaryote db
rule prok_sylph:
    input:
        sr1=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"), sid=SAMPLES.index),
        sr2=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"), sid=SAMPLES.index),
        gtdb=os.path.join(DB_DIR, "sylph/gtdb-r220-c200-dbv1.syldb")
    output:
        os.path.join(RESULTS_DIR, "sylph/profiling_prokaryotes.tsv")
    conda:
        os.path.join(ENV_DIR, "sylph.yaml")
    threads:
        config['sylph']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/prok_sylph.log")
    message:
        "Running SYLPH on all samples for PROKARYOTES"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && "
        "sylph profile {input.gtdb} -1 {input.sr1} -2 {input.sr2} -t {threads} > {output} && "
        "date) &> >(tee {log})"

rule viral_sylph:
    input:
        sr1=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"), sid=SAMPLES.index),
        sr2=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"), sid=SAMPLES.index),
        imgvr=os.path.join(DB_DIR, "sylph/imgvr_c200_v0.3.0.syldb")
    output:
        os.path.join(RESULTS_DIR, "sylph/profiling_viruses.tsv")
    conda:
        os.path.join(ENV_DIR, "sylph.yaml")
    threads:
        config['sylph']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/viral_sylph.log")
    message:
        "Running SYLPH on all samples for VIRUSES"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && "
        "sylph profile {input.imgvr} -1 {input.sr1} -2 {input.sr2} -t {threads} > {output} && "
        "date) &> >(tee {log})"  

rule fungal_sylph:
    input:
        sr1=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"), sid=SAMPLES.index),
        sr2=expand(os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"), sid=SAMPLES.index),
        fundb=os.path.join(DB_DIR, "sylph/fungi-refseq-2024-07-25-c200-v0.3.syldb")
    output:
        os.path.join(RESULTS_DIR, "sylph/profiling_fungi.tsv")
    conda:
        os.path.join(ENV_DIR, "sylph.yaml")
    threads:
        config['sylph']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/fungal_sylph.log")
    message:
        "Running SYLPH on all samples for FUNGI"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && "
        "sylph profile {input.fundb} -1 {input.sr1} -2 {input.sr2} -t {threads} > {output} && "
        "date) &> >(tee {log})"

