"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-09-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 24 -rp
Latest modification:
Purpose: Run amrscan (H.S. Gweon's read-based AMR gene/variant detector, via
         the submodules/snake_amrscan submodule) on preprocessed reads.

Third, complementary AMR-detection approach alongside rules/amr.smk (RGI on
assembled contigs) and rules/bgc_amr.smk's curated-KEGG-marker summary
(eggNOG on assembled contigs): amrscan works directly on reads, giving a
gene/variant-level result independent of assembly quality.

amrscan's own tool (submodules/snake_amrscan/submodules/resscan) needs
bwa+diamond+samtools+pandas+numpy together (see resscan/environment.yml).
The pre-existing 'amrscan' conda env on this host has bwa/diamond/samtools
but not pandas/numpy, so PATH is stitched from that env's bin dir plus the
base miniforge3 python (which has both) -- reusing what's already installed
rather than rebuilding either.
"""


############################################
rule amrscan:
    input:
        expand(os.path.join(RESULTS_DIR, "amrscan", "{sid}", "{sid}_varscan.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "amrscan", "merged_homscan.tsv"),
        os.path.join(RESULTS_DIR, "amrscan", "merged_varscan.tsv")
    output:
        touch("status/amrscan.done")


############################################
RESSCAN_SCRIPT = os.path.join(SUBMODULES, "snake_amrscan", "submodules", "resscan", "resscan", "resscan.py")
AMRSCAN_MERGE_SCRIPT = os.path.join(SUBMODULES, "snake_amrscan", "scripts", "merge_amrscan_tables.py")


rule run_amrscan:
    """Run resscan (amrscan) on a sample's trimmed paired-end reads."""
    input:
        r1=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_1.fq.gz"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_2.fq.gz"),
    output:
        varscan=os.path.join(RESULTS_DIR, "amrscan", "{sid}", "{sid}_varscan.tsv"),
        homscan=os.path.join(RESULTS_DIR, "amrscan", "{sid}", "{sid}_homscan.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/amrscan/{sid}_amrscan.log")
    threads:
        config["amrscan"]["threads"]
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    params:
        script=RESSCAN_SCRIPT,
        db=config["amrscan"]["db"],
        python_bin=config["amrscan"]["python_bin"],
        env_bin=config["amrscan"]["env_bin"],
    message:
        "Running amrscan (resscan) on {wildcards.sid}"
    shell:
        """
        mkdir -p "$(dirname {output.varscan})" "$(dirname {log})"
        (date && cd "$(dirname {params.script})/.." && \
        PATH="{params.env_bin}:$PATH" {params.python_bin} {params.script} \
            -i {input.r1},{input.r2} --card-db-dir {params.db} \
            -o "$(dirname {output.varscan})" -t {threads} --overwrite && \
        date) &> {log}
        """

rule merge_amrscan_results:
    """Merge amrscan (homscan/varscan) results across all samples."""
    input:
        homscans=expand(os.path.join(RESULTS_DIR, "amrscan", "{sid}", "{sid}_homscan.tsv"), sid=SAMPLES.index),
        varscans=expand(os.path.join(RESULTS_DIR, "amrscan", "{sid}", "{sid}_varscan.tsv"), sid=SAMPLES.index)
    output:
        merged_homscan=os.path.join(RESULTS_DIR, "amrscan", "merged_homscan.tsv"),
        merged_varscan=os.path.join(RESULTS_DIR, "amrscan", "merged_varscan.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/amrscan/merge_amrscan.log")
    message:
        "Merging amrscan results from all samples."
    script:
        AMRSCAN_MERGE_SCRIPT
