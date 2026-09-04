"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/annotation.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run gene calling via PRODIGAL on contigs
"""


############################################
rule annotation:
    input:
        expand(os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.{filetype}"), sid=SAMPLES.index, filetype=["faa", "gff"])
    output:
        touch("status/annotation.done")


############################################
# localrules:


############################################
# Annotating the genes in assemblies
rule prodigal:
    # barrier: don't start gene calling for ANY sample until assembly has
    # finished for ALL samples -- see megahit's barrier comment in
    # rules/assembly.smk for the phase-ordering rationale.
    input:
        FASTA=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        barrier="status/assembly.done"
    output:
        FAA=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa"),
        GFF=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.gff")
    conda:
        os.path.join(ENV_DIR, "prodigal.yaml")
    threads:
        config['prodigal']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/prodigal.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running PRODIGAL on {wildcards.sid}"
    shell:
        "(date && prodigal -a {output.FAA} -p meta -i {input.FASTA} -f gff -o {output.GFF} && date) &> >(tee {log})"

