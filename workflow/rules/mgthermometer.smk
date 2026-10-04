"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Purpose: Optional module -- estimate each final dereplicated MAG's optimal
         growth temperature (OGT) proxy via the IVYWREL amino-acid
         frequency method (Zeldovich et al. 2007), using UKCEH's own
         getFrequency.pl script (workflow/scripts/getFrequency.pl, carried
         over from metag_analyses). Adapted from
         https://github.com/michoug/MAGsGeneration's mgthermometer.smk +
         its annotations.smk (only the pieces MGthermometer needs: a
         per-MAG Prodigal/pyrodigal protein prediction).

Same pooled-assembly adaptation as the other optional modules: operates
once on the final Galah-dereplicated MAG set (rules.galah.output.final)
rather than per-sample.
"""


############################################
rule mgthermometer:
    input:
        os.path.join(RESULTS_DIR, "MGthermometer/table_OGT.txt")
    output:
        touch("status/mgthermometer.done")


############################################
checkpoint mgthermometer_genomes:
    input:
        rules.galah.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "MGthermometer/intermediate/genomes"))
    log:
        os.path.join(RESULTS_DIR, "logs/mgthermometer_genomes.log")
    message:
        "Collecting final dereplicated MAGs for MGthermometer"
    shell:
        "(date && mkdir -p {output} && cp {input}/*.fa {output}/ && date) &> {log}"


rule mgthermometer_annotation:
    input:
        fasta=os.path.join(RESULTS_DIR, "MGthermometer/intermediate/genomes/{mags}.fa")
    output:
        protein=os.path.join(RESULTS_DIR, "MGthermometer/annotation/{mags}.faa")
    log:
        os.path.join(RESULTS_DIR, "logs/mgthermometer_annotation/{mags}.log")
    threads:
        4
    message:
        "Prodigal protein prediction for {wildcards.mags}"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e annotations pyrodigal -a {output.protein} -j {threads} -i {input.fasta} && \
        date) &> {log}
        """


rule mgthermometer_run:
    input:
        fasta=rules.mgthermometer_annotation.output.protein,
        script=os.path.join(SRC_DIR, "getFrequency.pl")
    output:
        proportion=os.path.join(RESULTS_DIR, "MGthermometer/intermediate/proportion_{mags}.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/MGthermometer/{mags}.log")
    message:
        "Calculating IVYWREL proportion for {wildcards.mags}"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e mgthermometer perl {input.script} {input.fasta} {output.proportion} && \
        date) &> {log}
        """


def count_mags(wildcards):
    ck_output = checkpoints.mgthermometer_genomes.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "MGthermometer/intermediate/proportion_{mags}.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule mgthermometer_concat:
    input:
        count_mags
    output:
        os.path.join(RESULTS_DIR, "MGthermometer/table_OGT.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/mgthermometer_concat.log")
    message:
        "Concatenating IVYWREL/OGT results for all MAGs"
    shell:
        "(date && cat {input} > {output} && date) &> {log}"
