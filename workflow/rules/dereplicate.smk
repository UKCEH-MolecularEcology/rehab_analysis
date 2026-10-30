"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/dereplicate.smk --use-conda --cores 4 -rp
Latest modification: [2026-10-30] Replaced DAS_Tool (dasTool +
    filter_dastool_bins + checkm2_dastool) with Binette (binette_run +
    checkm2_binette), ported from michoug/MAGsGeneration#37; the
    prepare_dasTool contig2bin tables are reused unchanged.
    Earlier: [2026-10-04] Upgraded from a 2-way (MetaBAT2+CONCOCT)
    to a 7-way DAS_Tool ensemble (adds Rosella, TaxVAMB, COMEBin, SemiBin2,
    MetaCAT from rules/mags_generation.smk), added a Rosella refine pass on
    the DAS_Tool consensus bins, and switched dereplication from dRep to
    Galah -- all adapted from https://github.com/michoug/MAGsGeneration,
    matching its own choices (except CheckM2, reused from this project's
    existing setup -- see the checkm2 rules below -- rather than via pixi).
Purpose: Post-processing bins
"""


############################################
rule dereplicate:
    input:
        os.path.join(RESULTS_DIR, "bins/finalbins")
    output:
        touch("status/dereplicate.done")


############################################
# localrules:


# (ass_cat_drep removed 2026-10-10: it re-concatenated all 147 samples'
# assemblies from scratch as DAS_Tool's -c contig reference, duplicating
# cat_ass_mmseqs2+cat_filter_length's work but WITHOUT the mmseqs2 dedup
# or length filter those apply -- a mismatch, not just redundant work,
# since all 7 binners actually binned against the deduplicated+filtered
# CAT_ASSEMBLY_FILTER. dasTool below now references CAT_ASSEMBLY_FILTER
# directly, matching what every binner's scaffolds2bin table was computed
# against.)


# Preparing all 7 binners' outputs as DAS_Tool scaffolds2bin tables.
# metabat/concoct: unchanged from the original 2-way ensemble.
# rosella/taxvamb/metacat: converted via Fasta_to_Contig2Bin.sh (same
#   helper DAS_Tool itself ships) directly from each tool's own per-bin
#   fasta output directory -- more robust than parsing each tool's
#   internal cluster-assignment file format, and avoids a correctness
#   issue found in the reference repo's own prepare_dasTool for MetaCAT
#   (it fed MetaCAT's .done marker file -- an empty touch-file, not real
#   data -- to its scaffolds2bin conversion; using the actual Bins/
#   directory here instead).
# comebin/semibin: perl-reformatted directly from each tool's own
#   contig/cluster assignment table, matching the reference repo's exact
#   logic (not independently re-verified against real COMEBin/SemiBin2
#   output yet, since binning hasn't run in this pipeline -- flagging
#   this assumption for whoever first sees real output from this step).
rule prepare_dasTool:
    input:
        metabat=rules.metabat2.output,
        concoct=rules.concoct.output,
        rosella=rules.rosella.output.dir,
        taxvamb=rules.taxvamb.output,
        comebin=rules.comebin.output,
        semibin=rules.semibin.output,
        metacat=rules.metacat_cluster.output.folder,
    output:
        metabatout=os.path.join(RESULTS_DIR, "bins/dastool/metabat_das.tsv"),
        concoctout=os.path.join(RESULTS_DIR, "bins/dastool/concoct_das.tsv"),
        rosellaout=os.path.join(RESULTS_DIR, "bins/dastool/rosella_das.tsv"),
        taxvambout=os.path.join(RESULTS_DIR, "bins/dastool/taxvamb_das.tsv"),
        comebinout=os.path.join(RESULTS_DIR, "bins/dastool/comebin_das.tsv"),
        semibinout=os.path.join(RESULTS_DIR, "bins/dastool/semibin_das.tsv"),
        metacatout=os.path.join(RESULTS_DIR, "bins/dastool/metacat_das.tsv"),
    params:
        value="bin",
        src=os.path.join(SRC_DIR, "Fasta_to_Contig2Bin.sh")
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool.err.log")
    message:
        "Preparing files for DasTool (7-way ensemble)"
    shell:
        """
        (date && mkdir -p $(dirname {output.metabatout}) && \
        {params.src} -e fa -i {input.metabat} > {output.metabatout} && \
        perl -pe 's/^(\S+)[^\t]*(\t.*)$/$1$2/' {output.metabatout} > t.txt && mv t.txt {output.metabatout} && \
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > t.txt && mv t.txt {output.metabatout} && \
        perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 | sed 's/>s/s/g' > {output.concoctout} && \
        {params.src} -e fna -i {input.rosella} > {output.rosellaout} && \
        perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} > t2.txt && mv t2.txt {output.rosellaout} && \
        {params.src} -e fna -i $(dirname {input.taxvamb})/bins > {output.taxvambout} && \
        perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} > t3.txt && mv t3.txt {output.taxvambout} && \
        perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} > {output.comebinout} && \
        perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 > {output.semibinout} && \
        {params.src} -e fasta -i {input.metacat} > {output.metacatout} && \
        perl -pe 's/metacat\\./{params.value}_metacat_/g' {output.metacatout} > t4.txt && mv t4.txt {output.metacatout} && \
        date) 2> {log.err} > {log.out}
        """


rule binette_run:
    # Binette (michoug/MAGsGeneration#37) replaces DAS_Tool as the
    # consensus step: it takes the same per-binner contig2bin tables
    # prepared above, scores/merges candidate bins with its own internal
    # CheckM2 run, and writes only its final non-redundant selection --
    # so the old filter_dastool_bins step (DAS_Tool's --write_bins also
    # wrote every rejected candidate) is no longer needed.
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        concoct=rules.prepare_dasTool.output.concoctout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        taxvamb=rules.prepare_dasTool.output.taxvambout,
        comebin=rules.prepare_dasTool.output.comebinout,
        semibin=rules.prepare_dasTool.output.semibinout,
        metacat=rules.prepare_dasTool.output.metacatout,
        cont=CAT_ASSEMBLY_FILTER,
        # checkm_db is defined further down this file, so reference its output by path
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    output:
        quality=os.path.join(RESULTS_DIR, "bins/binette/final_bins_quality_reports.tsv"),
        bins=directory(os.path.join(RESULTS_DIR, "bins/binette/final_bins"))
    log:
        out=os.path.join(RESULTS_DIR, "logs/binette.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/binette.err.log")
    threads:
        config["binette"]["threads"]
    params:
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    message:
        "Running Binette (7-way ensemble)"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf $(dirname {output.quality}) && \
        pixi r -e binette binette -b {input.metabat} {input.concoct} {input.rosella} {input.taxvamb} {input.comebin} {input.semibin} {input.metacat} \
            -c {input.cont} -o $(dirname {output.quality}) -t {threads} --checkm2_db {params.db} \
            --min_completeness 50 --max_length 20000000 --verbose && \
        date) 2> {log.err} > {log.out}
        """


############################################
# CheckM2 (reusing this project's existing setup -- see
# rules/bin_taxqual.smk's checkm_final for the final pass, after Galah
# dereplication; the two rules below are the pre-refine and
# post-refine/pre-Galah passes) -- deliberately NOT via pixi, since this
# is already configured and working against an existing database.
# checkm_db lives here (rather than bin_taxqual.smk, which uses it too)
# since this file loads first in workflow/Snakefile's "binning" STEPS
# block and needs it earlier.
rule checkm_db:
    output:
        os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_db.log")
    conda:
        "checkm2"
    message:
        "Downloading the checkm2 database"
    shell:
        "(date && checkm2 database --download --path $(dirname $(dirname {output})) && date) &> >(tee {log})"

rule checkm2_binette:
    # Binette already runs CheckM2 internally to score and select its final
    # bins, so re-running `checkm2 predict` on them would redo the same work.
    # Reformat Binette's own report into CheckM2's native quality_report.tsv
    # schema (Name/Completeness/Contamination) for rosella refine below.
    input:
        tsv=rules.binette_run.output.quality
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2_binette/quality_report.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_binette.log")
    params:
        script=os.path.join(config["work_dir"], "scripts/binette_quality_to_checkm2.py")
    message:
        "Reformatting Binette's quality report into CheckM2's native schema"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        mkdir -p $(dirname {output.tsv}) && \
        pixi r -e binette python {params.script} --binette-tsv {input.tsv} --out {output.tsv} && \
        date) &> {log}
        """


############################################
# Rosella refine -- a second pass over the Binette consensus bins using
# Rosella's own refine mode (density-based re-clustering, informed by
# CheckM2 quality), matching the reference repo's refine.smk.
rule rosella_refine:
    input:
        cont=CAT_ASSEMBLY_FILTER,
        bam=EXISTING_BAMS,
        check=rules.checkm2_binette.output.tsv,
        bins=rules.binette_run.output.bins
    output:
        done=touch(os.path.join(RESULTS_DIR, "bins/refine/rosella_refine.done")),
        dir=directory(os.path.join(RESULTS_DIR, "bins/refine"))
    log:
        os.path.join(RESULTS_DIR, "logs/rosella_refine.log")
    threads:
        config["rosella"]["threads"]
    params:
        bins_dir=os.path.join(RESULTS_DIR, "bins/binette/final_bins")
    message:
        "Refining Binette consensus bins with Rosella"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf {output.dir} && mkdir -p {output.dir} && \
        pixi r -e rosella rosella refine --min-bin-size 200000 --checkm-results {input.check} \
            -d {params.bins_dir} --max-contamination 15 -x fa -b {input.bam} -r {input.cont} --threads {threads} -o {output.dir} && \
        date) &> {log}
        """


############################################
# CheckM2 on the refined bins (feeds Galah's dereplication, below).
rule checkm2_refine:
    input:
        bins=rules.rosella_refine.output.done,
        db=rules.checkm_db.output[0]
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2_refine/quality_report.tsv")
    conda:
        "/hdd0/susbus/tools/conda_envs/checkm2"
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_refine.log")
    threads:
        config["checkm"]["threads"]
    params:
        ext=config["checkm"]["extension"],
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd"),
        bins_dir=rules.rosella_refine.output.dir
    message:
        "Running CheckM2 on the Rosella-refined bins"
    shell:
        "(date && mkdir -p $(dirname {output.tsv}) && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --threads {threads} -x {params.ext} --input {params.bins_dir} --output-directory $(dirname {output.tsv}) --force && "
        "date) &> {log}"


############################################
# Galah dereplication, replacing dRep -- matches the reference repo's
# galah.smk. CheckM2's quality_report.tsv is reformatted to the plain
# "genome,completeness,contamination" CSV Galah expects.
rule galah_prepare:
    input:
        check=rules.checkm2_refine.output.tsv
    output:
        final=os.path.join(RESULTS_DIR, "bins/checkm2_before_galah.tsv"),
        temp=temp(os.path.join(RESULTS_DIR, "bins/checkm2_before_galah_temp.tsv"))
    log:
        os.path.join(RESULTS_DIR, "logs/galah_prepare.log")
    message:
        "Adjusting CheckM2 output for input to Galah"
    shell:
        """
        (date && cat {input.check} | sed '/^Name/d' | awk '{{print $1","$2","$3}}' > {output.temp} && \
        (echo "genome,completeness,contamination" && cat {output.temp}) > {output.final} && \
        date) &> {log}
        """

rule galah:
    input:
        check=rules.galah_prepare.output.final,
        bins_dir=rules.rosella_refine.output.dir
    output:
        final=directory(os.path.join(RESULTS_DIR, "bins/finalbins")),
        cluster=os.path.join(RESULTS_DIR, "bins/cluster_galah.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/galah.log")
    threads:
        config["galah"]["threads"]
    params:
        comp=config["galah"]["comp"],
        cont=config["galah"]["cont"]
    message:
        "Running Galah to dereplicate the refined bins"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.final} && \
        pixi r -e galah galah cluster --genome-fasta-directory {input.bins_dir} -x fa --genome-info {input.check} \
            --min-completeness {params.comp} --max-contamination {params.cont} \
            --output-representative-fasta-directory-copy {output.final} --output-cluster-definition {output.cluster} -t {threads} && \
        date) &> {log}
        """
