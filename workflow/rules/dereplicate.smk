"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/dereplicate.smk --use-conda --cores 4 -rp
Latest modification: [2026-11-02] Replaced the Rosella refine pass with a GUNC
    chimerism filter (gunc_run + gunc_filter, MAGsGeneration#40); Binette
    now also caps contamination (binette.max_contamination).
    Earlier: [2026-10-31] Dropped COMEBin from the ensemble (CPU
    training on 15M contigs is ~4h/epoch, i.e. weeks); now 6-way.
    Earlier: [2026-10-30] Replaced DAS_Tool (dasTool +
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


# Preparing all 6 binners' outputs as DAS_Tool scaffolds2bin tables.
# metabat/concoct: unchanged from the original 2-way ensemble.
# rosella/taxvamb/metacat: converted via Fasta_to_Contig2Bin.sh (same
#   helper DAS_Tool itself ships) directly from each tool's own per-bin
#   fasta output directory -- more robust than parsing each tool's
#   internal cluster-assignment file format, and avoids a correctness
#   issue found in the reference repo's own prepare_dasTool for MetaCAT
#   (it fed MetaCAT's .done marker file -- an empty touch-file, not real
#   data -- to its scaffolds2bin conversion; using the actual Bins/
#   directory here instead).
# semibin: perl-reformatted directly from the tool's own
#   contig/cluster assignment table, matching the reference repo's exact
#   logic (not independently re-verified against real SemiBin2
#   output yet, since binning hasn't run in this pipeline -- flagging
#   this assumption for whoever first sees real output from this step).
rule prepare_dasTool:
    input:
        metabat=rules.metabat2.output,
        concoct=rules.concoct.output,
        rosella=rules.rosella.output.dir,
        taxvamb=rules.taxvamb.output,
        semibin=rules.semibin.output,
        metacat=rules.metacat_cluster.output.folder,
    output:
        metabatout=os.path.join(RESULTS_DIR, "bins/dastool/metabat_das.tsv"),
        concoctout=os.path.join(RESULTS_DIR, "bins/dastool/concoct_das.tsv"),
        rosellaout=os.path.join(RESULTS_DIR, "bins/dastool/rosella_das.tsv"),
        taxvambout=os.path.join(RESULTS_DIR, "bins/dastool/taxvamb_das.tsv"),
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
        "Preparing contig2bin tables (6-way ensemble)"
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
        db=os.path.join(DB_DIR, "CheckM2_database/uniref100.KO.1.dmnd"),
        max_contamination=config["binette"]["max_contamination"]
    message:
        "Running Binette (6-way ensemble)"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf $(dirname {output.quality}) && \
        pixi r -e binette binette -b {input.metabat} {input.concoct} {input.rosella} {input.taxvamb} {input.semibin} {input.metacat} \
            -c {input.cont} -o $(dirname {output.quality}) -t {threads} --checkm2_db {params.db} \
            --min_completeness 50 --max_contamination {params.max_contamination} --max_length 20000000 --verbose && \
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
# GUNC database. Zenodo records and md5 sums are those of GUNC 1.1.1
# (gunc/gunc_database.py), as in MAGsGeneration#40.
GUNC_DOWNLOADS = {
    "progenomes_2.1": (
        "https://zenodo.org/records/19631736/files/gunc_db_progenomes2.1.dmnd.gz",
        "bc93a855e0760aad5c4e5f2d0e26da46",
        "447c9330056b02f29f30fe81fe4af4eb",
    ),
    "progenomes_3": (
        "https://zenodo.org/records/19631883/files/gunc_db_progenomes3.dmnd.gz",
        "a26b3be496017f291eece27740d0f37f",
        "b1a9347d219a632e9ecf4652dee9bdd1",
    ),
    "gtdb_95": (
        "https://zenodo.org/records/19631804/files/gunc_db_gtdb95.dmnd.gz",
        "14ed95a2fb1360925e28b7b55df14574",
        "3b502dc047efaafb9831a5eaec7617fd",
    ),
    "gtdb_214": (
        "https://zenodo.org/records/19632326/files/gunc_db_gtdb214.dmnd.gz",
        "3933007d83a7a85e4672295ee1f9d91f",
        "ac4d6304cab3d62a396703eeb039d3b7",
    ),
}
GUNC_CFG = config["gunc"]
if GUNC_CFG["database"] not in GUNC_DOWNLOADS:
    raise ValueError(
        f"gunc.database must be one of {sorted(GUNC_DOWNLOADS)}, got {GUNC_CFG['database']!r}"
    )
GUNC_URL, GUNC_GZ_MD5, GUNC_DMND_MD5 = GUNC_DOWNLOADS[GUNC_CFG["database"]]
GUNC_DB = os.path.join(GUNC_CFG["db_dir"], os.path.basename(GUNC_URL)[: -len(".gz")])


rule gunc_download:
    output:
        GUNC_DB + ".gz"
    log:
        os.path.join(RESULTS_DIR, "logs/gunc_download.log")
    threads: 1
    params:
        url=GUNC_URL,
        md5=GUNC_GZ_MD5
    message:
        "Downloading the GUNC database"
    shell:
        """
        (date && mkdir -p $(dirname {output}) && \
        wget --tries=3 --waitretry=30 --timeout=120 -O {output}.tmp {params.url} && \
        echo "{params.md5}  {output}.tmp" | md5sum -c - && \
        mv {output}.tmp {output} && \
        date) &> {log}
        """


rule gunc_db:
    input:
        rules.gunc_download.output
    output:
        GUNC_DB
    log:
        os.path.join(RESULTS_DIR, "logs/gunc_db.log")
    threads: 1
    params:
        md5=GUNC_DMND_MD5
    message:
        "Decompressing the GUNC database"
    shell:
        """
        (date && gzip -dc {input} > {output}.tmp && \
        echo "{params.md5}  {output}.tmp" | md5sum -c - && \
        mv {output}.tmp {output} && \
        date) &> {log}
        """


# GUNC chimerism check on the Binette bins. Replaces the Rosella refine pass
# (MAGsGeneration#40): bins flagged chimeric (pass.GUNC False) or above
# gunc.max_contamination are dropped by gunc_filter below.
rule gunc_run:
    input:
        bins=rules.binette_run.output.bins,
        db=GUNC_DB
    output:
        os.path.join(RESULTS_DIR, "bins/refine/gunc/GUNC.maxCSS_level.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/gunc_run.log")
    threads:
        GUNC_CFG["threads"]
    params:
        ext=".fa"
    message:
        "Running GUNC on the Binette bins"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        outdir=$(dirname {output}) && rm -rf $outdir && mkdir -p $outdir/tmp && \
        if compgen -G "{input.bins}/*{params.ext}" > /dev/null; then \
            pixi r -e gunc gunc run --input_dir {input.bins} --file_suffix {params.ext} --db_file {input.db} \
                --threads {threads} --out_dir $outdir --temp_dir $outdir/tmp && \
            cp $outdir/GUNC.*.maxCSS_level.tsv {output}; \
        else \
            echo "WARNING: no bins in {input.bins}, writing an empty GUNC table" >&2 && \
            printf 'genome\tpass.GUNC\n' > {output}; \
        fi && \
        rm -rf $outdir/tmp && \
        date) &> {log}
        """


rule gunc_filter:
    input:
        bins=rules.binette_run.output.bins,
        quality=rules.checkm2_binette.output.tsv,
        gunc=rules.gunc_run.output
    output:
        bins=directory(os.path.join(RESULTS_DIR, "bins/refine/gunc_bins")),
        table=os.path.join(RESULTS_DIR, "bins/refine/gunc_filter.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/gunc_filter.log")
    threads: 1
    params:
        script=os.path.join(config["work_dir"], "scripts/gunc_filter_bins.py"),
        max_contamination=GUNC_CFG["max_contamination"]
    message:
        "Filtering the Binette bins on GUNC and contamination"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf {output.bins} && \
        pixi r -e gunc python {params.script} --bins {input.bins} --quality {input.quality} --gunc {input.gunc} \
            --outdir {output.bins} --table {output.table} --max-contamination {params.max_contamination} && \
        date) &> {log}
        """


############################################
# CheckM2 on the refined bins (feeds Galah's dereplication, below).
rule checkm2_refine:
    input:
        bins=rules.gunc_filter.output.bins,
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
        bins_dir=rules.gunc_filter.output.bins
    message:
        "Running CheckM2 on the GUNC-filtered bins"
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
        bins_dir=rules.gunc_filter.output.bins
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
        "Running Galah to dereplicate the GUNC-filtered bins"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.final} && \
        pixi r -e galah galah cluster --genome-fasta-directory {input.bins_dir} -x fa --genome-info {input.check} \
            --min-completeness {params.comp} --max-contamination {params.cont} \
            --output-representative-fasta-directory-copy {output.final} --output-cluster-definition {output.cluster} -t {threads} && \
        date) &> {log}
        """
