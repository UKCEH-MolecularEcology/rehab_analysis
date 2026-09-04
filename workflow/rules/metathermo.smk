"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-07-01]
Run: snakemake -s workflow/rules/metathermo.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run meta-thermo on MAGs
"""

MAGS=glob_wildcards(os.path.join(RESULTS_DIR, "bins/all_finalbins/{mags}.fa")).mags


############################################
rule metathermo:
    input:
        expand(os.path.join(RESULTS_DIR, "Prodigal/{genome}.faa"), genome = MAGS),
        expand(os.path.join(RESULTS_DIR, "metathermo/table_OGT.txt"))
    output:
        touch("status/metathermo.done")


############################################
localrules: 


############################################
rule mag_prodigal:
    input:
        os.path.join(RESULTS_DIR, "bins/all_finalbins/{genome}.fa")
    output:
        protein = os.path.join(RESULTS_DIR, "Prodigal/{genome}.faa"),
        gff = os.path.join(RESULTS_DIR, "Prodigal/{genome}.gff"),
        nucl = os.path.join(RESULTS_DIR, "Prodigal/{genome}.ffn")
    params:
        prefix = "{genome}"
    conda:
        os.path.join(ENV_DIR, "metathermo.yaml")
    threads:
        config["prodigal"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/prodigal/{genome}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prodigal/{genome}.err.log")
    message:
        "Let's first annotate all the genomes"
    shell:
        """
        (date && 
        prodigal -a {output.protein} -d {output.nucl} -f gff -i {input} -q -o {output.gff} && 
        date) 2> {log.err} > {log.out}
        """

rule AAproportion:
    input:
        prot=os.path.join(RESULTS_DIR, "Prodigal/{genome}.faa"),
        script = os.path.join(SRC_DIR,"getFrequency.pl")
    output:
        os.path.join(RESULTS_DIR, "Proportion/{genome}.txt")
    conda:
        os.path.join(ENV_DIR, "metathermo.yaml")
    message:
        "Calulate proportion of IVYWREL"
    shell:
        """
        perl {input.script} {input.prot} {output}
        """

rule concatResults:
    input:
        expand(os.path.join(RESULTS_DIR, "Proportion/{genome}.txt"), genome = MAGS)
    output:
        os.path.join(RESULTS_DIR, "metathermo/table_OGT.txt")
    shell:
        "cat {input} > {output}"
 
