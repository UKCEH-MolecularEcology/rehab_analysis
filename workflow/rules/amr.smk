"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-12-25]
Run: snakemake -s workflow/rules/amr.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run AMR on concatenated assembly
"""


############################################
rule amr:
    input:
        os.path.join(RESULTS_DIR, "amr/rgi.txt")
    output:
        touch("status/amr.done")


############################################
localrules: download_rgi_db, setup_rgi_db


############################################
# Download RGI data
rule download_rgi_db:
    output:
        archive=temp(os.path.join(DB_DIR, "rgi/card-data.tar.bz2")),
        json=os.path.join(DB_DIR, "rgi/card.json")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.rgi.db.log")
    params:
        db_url=config["rgi"]["db_url"]
    message:
        "Setup: download RGI data"
    shell:
        "(date && "
        "wget -O {output.archive} {params.db_url} --no-check-certificate && "
        "tar -C $(dirname {output.archive}) -xvf {output.archive} && "
        "date) &> >(tee {log})"

# Setup RGI: load required DB
# NOTE: to make sure that the same DB is used for all targets
rule setup_rgi_db:
    input:
        os.path.join(DB_DIR, "rgi/card.json")
    output:
        "status/rgi_setup.done"
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.setup.log")
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Setup: load RGI DB"
    shell:
        "(date && "
        "rgi clean --local && "
        "rgi load --card_json {input} --local && "
        "rgi database --version --local && "
        "touch {output} && date) &> >(tee {log})"

# Run RGI: Assembly (DNA), chunked
# NOTE: RGI on the whole concatenated assembly (~48G) runs prodigal single-threaded
# for days (killed after 23h on 2026-10-16), so split into RGI_CHUNKS pieces and run
# them in parallel, then merge the per-chunk tables into amr/rgi.txt.
RGI_CHUNKS = [f"{i:03d}" for i in range(1, config["rgi"]["chunks"] + 1)]

rule rgi_split:
    input:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta")
    output:
        temp(expand(os.path.join(RESULTS_DIR, "amr/chunks/cat_assembly_filter.part_{chunk}.fasta"), chunk=RGI_CHUNKS))
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.split.log")
    threads:
        config["filter_length"]["threads"]
    params:
        bin=config["seqkit"]["bin"],
        parts=config["rgi"]["chunks"],
        outdir=os.path.join(RESULTS_DIR, "amr/chunks")
    message:
        "Annotation: RGI: split assembly into chunks"
    shell:
        "(date && "
        "{params.bin} split2 -j {threads} -p {params.parts} -O {params.outdir} {input} && "
        "date) &> >(tee {log})"

rule rgi_chunk:
    input:
        fna=os.path.join(RESULTS_DIR, "amr/chunks/cat_assembly_filter.part_{chunk}.fasta"),
        db=os.path.join(DB_DIR, "rgi/card.json"),
        setup="status/rgi_setup.done" # NOTE: to make sure that the same DB is used for all targets
    output:
        txt=temp(os.path.join(RESULTS_DIR, "amr/chunks/rgi.part_{chunk}.txt")),
        json=temp(os.path.join(RESULTS_DIR, "amr/chunks/rgi.part_{chunk}.json"))
    log:
        os.path.join(RESULTS_DIR, "logs/rgi/rgi.part_{chunk}.log")
    threads:
        config["rgi"]["chunk_threads"]
    params:
        alignment_tool="DIAMOND",
        prefix=lambda wildcards, output: os.path.splitext(output.txt)[0]
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Annotation: RGI: chunk {wildcards.chunk}"
    shell:
        "(date && "
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run
        "(rgi main --input_sequence {input.fna} --output_file {params.prefix} --local -a {params.alignment_tool} --clean --low_quality -n {threads} || "
        "rgi main --input_sequence {input.fna} --output_file {params.prefix} --local -a {params.alignment_tool} --clean --low_quality -n {threads}) && "
        "date) &> {log}"

rule annotation_rgi:
    input:
        expand(os.path.join(RESULTS_DIR, "amr/chunks/rgi.part_{chunk}.txt"), chunk=RGI_CHUNKS)
    output:
        txt=os.path.join(RESULTS_DIR, "amr/rgi.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.annotation.log")
    message:
        "Annotation: RGI: merge chunks"
    shell:
        "(date && "
        "head -n 1 {input[0]} > {output.txt} && "
        "for f in {input}; do tail -n +2 $f >> {output.txt}; done && "
        "date) &> {log}"
