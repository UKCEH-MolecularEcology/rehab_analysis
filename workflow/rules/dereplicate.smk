"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/dereplicate.smk --use-conda --cores 4 -rp
Latest modification: [2026-10-04] Upgraded from a 2-way (MetaBAT2+CONCOCT)
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


############################################
# Concatetnating the assemblies for binning
rule ass_cat_drep:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly_drep.fasta")  # "/prj/DECODE/socd/results/assembly/modified_cat_assembly_drep.fasta"
    log:
        os.path.join(RESULTS_DIR, "logs/assembly/concatenation.log")
    message:
        "Concatenating all assemblies"
    shell:
        "(date && cat {input} > {output} && date) &> >(tee {log})"


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


rule dasTool:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        concoct=rules.prepare_dasTool.output.concoctout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        taxvamb=rules.prepare_dasTool.output.taxvambout,
        comebin=rules.prepare_dasTool.output.comebinout,
        semibin=rules.prepare_dasTool.output.semibinout,
        metacat=rules.prepare_dasTool.output.metacatout,
        cont=rules.ass_cat_drep.output
    output:
        os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool.err.log")
    threads:
        config["dasTool"]["threads"]
    params:
        src=config["dasTool"]["bin"]
    message:
        "Running DasTool (7-way ensemble)"
    shell:
        "(date && "
        "{params.src} -i {input.metabat},{input.concoct},{input.rosella},{input.taxvamb},{input.comebin},{input.semibin},{input.metacat} "
        "-l metabat,concoct,rosella,taxvamb,comebin,semibin,metacat "
        "--score_threshold -42 -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads} && "
        "date) 2> {log.err} > {log.out}"


rule filter_dastool_bins:
    # DAS_Tool's --write_bins writes a FASTA for every candidate bin it
    # evaluated (including SCG-split "_sub" candidates and original bins
    # that lost every contig to a better-scoring competitor), not just its
    # final non-redundant selection. Only the bin names listed in
    # das_DASTool_summary.tsv are DAS_Tool's actual answer -- feeding the
    # raw das_DASTool_bins/ directory to CheckM2/Rosella refine (as both
    # downstream rules did before this fix) lets rejected bins that still
    # score well on their own (e.g. a near-complete bin that lost its
    # contigs to a different winning bin) pass every quality gate and
    # reach Galah as a bogus "distinct" genome that still shares raw
    # contigs with the bin that actually won them. Ported from
    # michoug/MAGsGeneration#34's filter_dastool_bins (the one piece of
    # that PR adopted here -- its larger binner-selection/ranking feature
    # was left for later, still unmerged upstream as of 2026-10-07).
    input:
        summary=rules.dasTool.output
    output:
        directory(os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_bins_filtered"))
    log:
        os.path.join(RESULTS_DIR, "logs/filter_dastool_bins.log")
    params:
        ext="fa",
        src=os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_bins")
    message:
        "Keeping only DAS_Tool's selected bins"
    shell:
        # Process substitution (not a trailing pipe) so `missing` set inside
        # the loop survives into the post-loop check -- a pipe would run the
        # loop in a subshell and silently discard the counter. Hardened to a
        # hard error (matching michoug/MAGsGeneration#34's latest commits,
        # 2026-10-07) rather than a WARNING: DAS_Tool's own summary naming a
        # bin its own --write_bins didn't produce means something is
        # actually wrong, not something to quietly skip past.
        """
        (date && mkdir -p {output} && \
        missing=0 && \
        while read -r bin; do \
            if [[ -f "{params.src}/$bin.{params.ext}" ]]; then \
                cp "{params.src}/$bin.{params.ext}" {output}/; \
            else \
                echo "ERROR: expected winning bin file missing: {params.src}/$bin.{params.ext}" >&2; \
                missing=$((missing + 1)); \
            fi; \
        done < <(tail -n +2 {input.summary} | cut -f1) && \
        if [[ $missing -gt 0 ]]; then \
            echo "ERROR: $missing winning bin(s) listed in das_DASTool_summary.tsv were missing from {params.src}" >&2; \
            exit 1; \
        fi && \
        date) &> {log}
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

rule checkm2_dastool:
    input:
        bins=rules.filter_dastool_bins.output,
        db=rules.checkm_db.output[0]
    output:
        tsv=os.path.join(RESULTS_DIR, "bins/checkm2_dastool/quality_report.tsv")
    conda:
        "/hdd0/susbus/tools/conda_envs/checkm2"
    log:
        os.path.join(RESULTS_DIR, "logs/checkm2_dastool.log")
    threads:
        config["checkm"]["threads"]
    params:
        ext=config["checkm"]["extension"],
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd"),
        bins_dir=os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_bins_filtered")
    message:
        "Running CheckM2 on the DAS_Tool consensus bins (pre-refine)"
    shell:
        "(date && mkdir -p $(dirname {output.tsv}) && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --threads {threads} -x {params.ext} --input {params.bins_dir} --output-directory $(dirname {output.tsv}) --force && "
        "date) &> {log}"


############################################
# Rosella refine -- a second pass over the DAS_Tool consensus bins using
# Rosella's own refine mode (density-based re-clustering, informed by
# CheckM2 quality), matching the reference repo's refine.smk.
rule rosella_refine:
    input:
        cont=CAT_ASSEMBLY_FILTER,
        bam=EXISTING_BAMS,
        check=rules.checkm2_dastool.output.tsv,
        bins=rules.filter_dastool_bins.output
    output:
        done=touch(os.path.join(RESULTS_DIR, "bins/refine/rosella_refine.done")),
        dir=directory(os.path.join(RESULTS_DIR, "bins/refine"))
    log:
        os.path.join(RESULTS_DIR, "logs/rosella_refine.log")
    threads:
        config["rosella"]["threads"]
    params:
        bins_dir=os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_bins_filtered")
    message:
        "Refining DAS_Tool consensus bins with Rosella"
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
