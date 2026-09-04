"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-10-02]
Run: snakemake -s workflow/rules/antismash.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run antiSMASH on MAGs
"""

FINAL_MAGS=glob_wildcards(os.path.join(RESULTS_DIR, "bins/all_finalbins/{mags}.fa")).mags


############################################
rule antismash:
    input:
        expand(os.path.join(RESULTS_DIR, "antismash/{genome}/{genome}.zip"), genome=FINAL_MAGS)
    output:
        touch("status/antismash.done")


############################################
localrules: 


############################################
rule download_antismash_db:
    output:
        done=os.path.join(RESULTS_DIR, "antismash/db_downloaded.done")
    log:
        out=os.path.join(RESULTS_DIR, "logs/antismash/setup.log")
    conda:
        os.path.join(ENV_DIR, "antismash.yaml")
    message:
        "Downloading DBs for antiSMASH"
    shell:
        "(date && download-antismash-databases && "
        "touch {output.done} && date) &> >(tee {log})"

# prodigal
rule bin_prodigal:
    input:
        os.path.join(RESULTS_DIR, "bins/all_finalbins/{genome}.fa")
    output:
#        protein = os.path.join(RESULTS_DIR, "Prodigal/{genome}.faa"),
        gbk=os.path.join(RESULTS_DIR, "Prodigal/{genome}.gbk")
    params:
        prefix = "{genome}"
    conda:
        os.path.join(ENV_DIR, "metathermo.yaml")
    threads:
        config["prodigal"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/prodigal/gbk_{genome}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prodigal/gbk_{genome}.err.log")
    message:
        "Let's first create a GBK for {wildcards.genome}"
    shell:
        """
        (date &&
        prodigal -f gbk -i {input} -o {output.gbk} -q &&
        date) 2> {log.err} > {log.out}
        """

# antiSMASH 
rule bins_antismash:
    input:
        mags=os.path.join(RESULTS_DIR, "bins/all_finalbins/{genome}.fa"),
        gbk=os.path.join(RESULTS_DIR, "Prodigal/{genome}.gbk"),
        dummy=os.path.join(RESULTS_DIR, "antismash/db_downloaded.done")
    output:
        os.path.join(RESULTS_DIR, "antismash/{genome}/{genome}.zip")
    conda:
        os.path.join(ENV_DIR, "antismash.yaml")
    threads:
        config['antismash']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/antismash.{genome}.log")
    message:
        "Running antiSMASH for all MAGS"
    shell:
        "(date && antismash --cpus {threads} --genefinding-tool none --minimal {input.gbk} --output-dir $(dirname {output}) && "
        "date) &> >(tee {log})"

# "(date && antismash --cpus {threads} --genefinding-tool prodigal --asf --cb-knownclusters --clusterhmmer {input.ma
# "(date && antismash --cpus {threads} --genefinding-tool prodigal --fullhmmer --pfam2go --asf --cb-knownclusters --clusterhmmer {input.mags} --output-dir $(dirname {output}) && "

