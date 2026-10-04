"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Purpose: Optional module -- per-MAG synteny tracking across all 147 samples
         via SynTracker (https://github.com/leylabmpi/SynTracker), using
         each sample's own per-sample assembly as the BLAST target pool
         and the final dereplicated MAGs as query references. Adapted
         from https://github.com/michoug/MAGsGeneration's syntracker.smk.

Unlike the other optional modules, this one does NOT need the pooled-
assembly adaptation on its target-pool side -- SynTracker is specifically
about tracking a MAG's synteny across individual samples, so it needs the
147 *per-sample* assemblies (results/assembly/{sid}/{sid}.fasta), not the
single pooled assembly. Only the MAG side (final dereplicated genomes)
collapses to one pooled checkpoint instead of per-sample, same as
mgthermometer.smk's own genome checkpoint.
"""


############################################
rule syntracker:
    input:
        os.path.join(RESULTS_DIR, "syntracker/intermediate/blastDB"),
        os.path.join(RESULTS_DIR, "syntracker/syntracker_finalized.done")
    output:
        touch("status/syntracker.done")


############################################
localrules: syntracker_install, cp_assemblies


############################################
rule syntracker_install:
    output:
        dbs=directory(os.path.join(DB_DIR, "SynTracker")),
        done=touch(os.path.join(DB_DIR, "SynTracker/syntracker.done"))
    log:
        os.path.join(RESULTS_DIR, "logs/syntracker_install.log")
    message:
        "Installing SynTracker"
    shell:
        """
        (date && rm -rf "$(dirname {output.dbs})/SynTracker_clone" && \
        git clone https://github.com/leylabmpi/SynTracker.git "$(dirname {output.dbs})/SynTracker_clone" &> /dev/null && \
        rm -rf {output.dbs} && mv "$(dirname {output.dbs})/SynTracker_clone" {output.dbs} && \
        date) &> {log}
        """


rule cp_assemblies:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"), sid=SAMPLES.index)
    output:
        dir=temp(directory(os.path.join(RESULTS_DIR, "syntracker/intermediate/assemblies")))
    log:
        os.path.join(RESULTS_DIR, "logs/cp_assemblies.log")
    message:
        "Collecting per-sample assemblies as the SynTracker BLAST target pool"
    shell:
        "(date && mkdir -p {output.dir} && cp {input} {output.dir}/ && date) &> {log}"


rule syntracker_blast:
    input:
        syntracker=rules.syntracker_install.output.dbs,
        assembly=rules.cp_assemblies.output.dir
    output:
        dir=directory(os.path.join(RESULTS_DIR, "syntracker/intermediate/blastDB")),
        done=touch(os.path.join(RESULTS_DIR, "syntracker/intermediate/blastDB.done"))
    log:
        os.path.join(RESULTS_DIR, "logs/syntracker_blast.log")
    threads:
        4
    message:
        "Building SynTracker's BLAST DB from all per-sample assemblies"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.dir} && \
        pixi r -e syntracker python {input.syntracker}/syntracker_makeDB.py -target {input.assembly} -out {output.dir} && \
        date) &> {log}
        """


############################################
checkpoint syntracker_genomes:
    input:
        rules.galah.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "syntracker/intermediate/genomes"))
    log:
        os.path.join(RESULTS_DIR, "logs/syntracker_genomes.log")
    message:
        "Collecting final dereplicated MAGs for SynTracker"
    shell:
        "(date && mkdir -p {output} && cp {input}/*.fa {output}/ && date) &> {log}"


def get_syn_mags(wildcards):
    ck_output = checkpoints.syntracker_genomes.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR,
            "syntracker/intermediate/syntracker_output_{mags}/summary_output/synteny_scores_per_region.csv",
        ),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule syntracker_run:
    input:
        syntracker=rules.syntracker_install.output.dbs,
        blastdb=rules.syntracker_blast.output.dir,
        blastdb_done=rules.syntracker_blast.output.done,
        mags=os.path.join(RESULTS_DIR, "syntracker/intermediate/genomes/{mags}.fa")
    output:
        done=os.path.join(
            RESULTS_DIR,
            "syntracker/intermediate/syntracker_output_{mags}/summary_output/synteny_scores_per_region.csv",
        ),
        folder=directory(os.path.join(RESULTS_DIR, "syntracker/intermediate/syntracker_output_{mags}"))
    log:
        os.path.join(RESULTS_DIR, "logs/syntracker/syntracker_run_{mags}.log")
    threads:
        5
    message:
        "Running SynTracker for {wildcards.mags} against all per-sample assemblies"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        mkdir -p {output.folder}_mag && \
        cp {input.mags} {output.folder}_mag/ && \
        pixi r -e syntracker python {input.syntracker}/syntracker.py -ref {output.folder}_mag -mode new -blastDB {input.blastdb} -out {output.folder} -cores {threads} -length 2500 && \
        rm -r {output.folder}_mag && \
        date) &> {log}
        """


rule syntracker_finalize:
    input:
        get_syn_mags
    output:
        score_100=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_100_regions.csv"),
        score_200=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_200_regions.csv"),
        score_40=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_40_regions.csv"),
        score_60=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_60_regions.csv"),
        score_80=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_80_regions.csv"),
        score_all=os.path.join(RESULTS_DIR, "syntracker/avg_synteny_scores_all_regions.csv"),
        score_per=os.path.join(RESULTS_DIR, "syntracker/synteny_scores_per_region.csv"),
        done=os.path.join(RESULTS_DIR, "syntracker/syntracker_finalized.done")
    log:
        os.path.join(RESULTS_DIR, "logs/syntracker_finalize.log")
    message:
        "Finalizing SynTracker output across all MAGs"
    shell:
        r"""
        set -euo pipefail
        (date
        out_dir=$(dirname "{output.score_all}")
        mkdir -p "$out_dir"
        for csv in avg_synteny_scores_100_regions.csv \
            avg_synteny_scores_40_regions.csv \
            avg_synteny_scores_60_regions.csv \
            avg_synteny_scores_80_regions.csv \
            avg_synteny_scores_200_regions.csv \
            avg_synteny_scores_all_regions.csv \
            synteny_scores_per_region.csv; do
            first=1
            dest_file="${{out_dir}}/${{csv}}"
            for file in {input}; do
                src_file="$(dirname "$file")/${{csv}}"
                if [ -f "$src_file" ]; then
                    if [ $first -eq 1 ]; then
                        cat "$src_file" > "$dest_file"
                        first=0
                    else
                        tail -n +2 "$src_file" >> "$dest_file"
                    fi
                fi
            done
        done
        touch {output.done}
        date) &> {log}
        """
