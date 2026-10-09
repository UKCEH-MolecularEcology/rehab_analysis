"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-09]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Purpose: Group antiSMASH regions of the final dereplicated MAGs into gene
         cluster families (GCFs) with BiG-SLiCE 2, and summarise GCFs per
         MAG and per sample. Ported from UKCEH-MolecularEcology/skyline's
         annotations branch (workflow/rules/bgc_families.smk), replacing
         this project's old BiG-SCAPE pass (which clustered per-sample
         antiSMASH regions -- removed along with the per-sample antiSMASH
         run it depended on, see bgc_amr.smk).

Runs via pixi (bigslice is a pip-only package, no plain bioconda recipe),
scoped the same way as the MAG-generation module's other pixi envs --
see pixi.toml and scripts/pixi_env.sh.
"""

BS = config["bgc_families"]
BS_DIR = os.path.join(RESULTS_DIR, "bgc", "bigslice")
BS_MODELS = config["bgc_families"]["models_dir"]
BS_REGIONS = os.path.join(RESULTS_DIR, "bgc", "antismash_regions_all.tsv")
BS_ANTISMASH = os.path.join(RESULTS_DIR, "bgc", "antismash", "per_mag")
# Not a strict relative-abundance matrix like skyline's CoverM output, but
# the same shape (MAG x sample) and close enough in spirit: per-bin,
# per-sample summed trimmed-mean contig coverage, already computed by
# rules/bin_coverage.smk for the final dereplicated bins. bigslice_tables.py
# treats this input as optional and skips the GCF x sample table cleanly if
# it's missing, so this is a safe best-effort wiring rather than a new
# relative-abundance pipeline built just for this.
BS_COVERAGE = os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage.txt")


############################################
rule bgc_families:
    input:
        os.path.join(BS_DIR, "tables", "gcf_summary.tsv")
    output:
        touch("status/bgc_families.done")


############################################
localrules: download_bigslice_models, bigslice_input


############################################
rule download_bigslice_models:
    output:
        done=os.path.join(os.path.dirname(BS_MODELS), "bigslice-models.done")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/download_bigslice_models.log")
    params:
        url=BS["models_url"],
        md5=BS["models_md5"],
        models=BS_MODELS
    message:
        "Downloading BiG-SLiCE HMM models"
    shell:
        """
        (date && d=$(dirname {params.models}) && mkdir -p $d && cd $d && \
        wget -c -O bigslice-models.tar.gz {params.url} && \
        echo '{params.md5}  bigslice-models.tar.gz' | md5sum -c - && \
        rm -rf {params.models} && mkdir -p {params.models} && tar -xzf bigslice-models.tar.gz -C {params.models} && \
        rm bigslice-models.tar.gz && ls {params.models} && touch {output.done} && date) &> {log}
        """


# BiG-SLiCE input folder: one dataset, one sub-folder per MAG with its
# antiSMASH region GenBank files (symlinks)
rule bigslice_input:
    input:
        regions=BS_REGIONS
    output:
        datasets=os.path.join(BS_DIR, "input", "datasets.tsv")
    params:
        ds_dir=os.path.join(BS_DIR, "input", "rehab"),
        antismash=BS_ANTISMASH
    message:
        "Preparing BiG-SLiCE input from antiSMASH regions"
    run:
        import glob
        import shutil

        mags = set()
        with open(input.regions) as fh:
            next(fh)
            for line in fh:
                mags.add(line.split("\t", 1)[0])
        shutil.rmtree(params.ds_dir, ignore_errors=True)
        n = 0
        for m in sorted(mags):
            os.makedirs(os.path.join(params.ds_dir, m), exist_ok=True)
            for g in glob.glob(os.path.join(params.antismash, m, "*.region*.gbk")):
                os.symlink(g, os.path.join(params.ds_dir, m, os.path.basename(g)))
                n += 1
        with open(output.datasets, "w") as out:
            out.write("# dataset_name\tdataset_path\ttaxonomy_path\tdescription\n")
            out.write("rehab\trehab\t\tREHAB final dereplicated MAGs (antiSMASH regions)\n")
        print("BiG-SLiCE input: {} regions from {} MAGs".format(n, len(mags)))


rule bigslice_run:
    input:
        datasets=rules.bigslice_input.output.datasets,
        models=rules.download_bigslice_models.output.done
    output:
        db=os.path.join(BS_DIR, "output", "result", "data.db")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/bigslice_run.log")
    params:
        indir=os.path.join(BS_DIR, "input"),
        outdir=os.path.join(BS_DIR, "output"),
        models=BS_MODELS,
        threshold=BS["threshold"],
        extra=BS.get("extra", "")
    threads:
        BS["threads"]
    message:
        "BiG-SLiCE GCF clustering (threshold {})".format(BS["threshold"])
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e bigslice pip show bigslice pyhmmer | grep -E '^(Name|Version)' && \
        rm -rf {params.outdir} && \
        pixi r -e bigslice bigslice -i {params.indir} {params.outdir} -t {threads} --threshold {params.threshold} \
            --program_db_folder {params.models} {params.extra} && \
        date) &> {log}
        """


# tables: region -> GCF, GCF summary, GCF x sample coverage (if the
# per-bin coverage table exists)
rule bigslice_tables:
    input:
        db=rules.bigslice_run.output.db,
        regions=BS_REGIONS
    output:
        membership=os.path.join(BS_DIR, "tables", "gcf_membership.tsv"),
        summary=os.path.join(BS_DIR, "tables", "gcf_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/bigslice_tables.log")
    params:
        script=os.path.join(SRC_DIR, "bigslice_tables.py"),
        coverage=BS_COVERAGE
    message:
        "Summarising BiG-SLiCE GCFs"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e bigslice python3 {params.script} --db {input.db} --regions {input.regions} --coverm {params.coverage} \
            --outdir $(dirname {output.summary}) && \
        date) &> {log}
        """
