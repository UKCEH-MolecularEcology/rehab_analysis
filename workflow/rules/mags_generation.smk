"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Latest modification:
Purpose: Additional binners (Rosella, TaxVAMB, COMEBin, SemiBin2, MetaCAT),
         adapted from https://github.com/michoug/MAGsGeneration, to extend
         the existing MetaBAT2+CONCOCT ensemble in rules/binning.smk before
         DAS_Tool (rules/dereplicate.smk).

Unlike the reference repo (which assembles and bins each sample
independently, then does all-vs-all read mapping across samples for
cross-sample coverage), this project already pools ALL samples' contigs
into one assembly (results/assembly/cat_assembly_filter.fasta, built by
rules/cluster.smk) and already maps all 147 samples' reads against that
one pooled reference (results/bam/{sid}/cat_assembly_{sid}.bam, built by
rules/coverage.smk's filter_mapping). That already gives full
147-sample coverage signal with no new mapping at all -- so every new
binner here is fed the SAME pooled contigs + SAME existing BAMs that
rules/binning.smk's metabat2/concoct already use, rather than
replicating the reference's per-sample/all-vs-all design (which would
not scale: ~21,600 BAM files at this project's 147-sample size).

All five new binners run via pixi (see pixi.toml / scripts/pixi_env.sh)
rather than conda: several of them (VAMB, MetaCAT) are pip/wheel
installs with no plain bioconda package, so this reuses the reference
repo's already-validated pixi manifest instead of re-deriving install
recipes from scratch. CheckM2 and GTDB-Tk are deliberately NOT pixi
envs -- see rules/dereplicate.smk and rules/bin_taxqual.smk, which
reuse this project's existing conda-based setups and databases.
"""

import os

PIXI_ENV_SCRIPT = os.path.join(config["work_dir"], "scripts/pixi_env.sh")
CAT_ASSEMBLY_FILTER = os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta")
EXISTING_BAMS = expand(os.path.join(RESULTS_DIR, "bam/{sid}/cat_assembly_{sid}.bam"), sid=SAMPLES.index)


############################################
rule mags_generation:
    input:
        os.path.join(RESULTS_DIR, "bins/rosella/rosella_bins.done"),
        os.path.join(RESULTS_DIR, "bins/taxvamb/vaevae_clusters_unsplit.tsv"),
        os.path.join(RESULTS_DIR, "bins/comebin/comebin_res/comebin_res.tsv"),
        os.path.join(RESULTS_DIR, "bins/semibin/contig_bins.tsv"),
        os.path.join(RESULTS_DIR, "bins/metacat/metacat.done"),
    output:
        touch("status/mags_generation.done")


############################################
localrules: mmseqs_gtdb_db, metacat_install


############################################
# Rosella -- density-based binner using UMAP + HDBSCAN on composition+coverage
rule rosella:
    input:
        cont=CAT_ASSEMBLY_FILTER,
        bam=EXISTING_BAMS
    output:
        done=touch(os.path.join(RESULTS_DIR, "bins/rosella/rosella_bins.done")),
        dir=directory(os.path.join(RESULTS_DIR, "bins/rosella"))
    log:
        os.path.join(RESULTS_DIR, "logs/rosella.log")
    threads:
        config["rosella"]["threads"]
    message:
        "Running Rosella on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf {output.dir} && \
        pixi r -e rosella rosella recover -b {input.bam} -r {input.cont} --threads {threads} -o {output.dir} && \
        date) &> {log}
        """


############################################
# TaxVAMB -- taxonomy-informed VAMB binning. Needs an mmseqs2 GTDB
# taxonomy classification of the pooled contigs as an extra input.
rule mmseqs_gtdb_db:
    output:
        done=touch(os.path.join(DB_DIR, "mmseqs_gtdb/mmseqs_gtdb.done"))
    log:
        os.path.join(RESULTS_DIR, "logs/mmseqs_gtdb_db.log")
    threads:
        config["mmseqs_gtdb"]["threads"]
    message:
        "Downloading the mmseqs2-formatted GTDB taxonomy database"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        mkdir -p "$(dirname {output.done})/tmp" && \
        pixi r -e mmseqs mmseqs databases GTDB "$(dirname {output.done})/GTDB_mmseqs" "$(dirname {output.done})/tmp" --threads {threads} && \
        rm -rf "$(dirname {output.done})/tmp" && \
        date) &> {log}
        """

rule contig_taxonomy_mmseqs:
    input:
        fasta=CAT_ASSEMBLY_FILTER,
        db=rules.mmseqs_gtdb_db.output.done
    output:
        taxonomy=os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_tax.tsv"),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "bins/taxonomy/tmp"))),
        dbfasta=temp(os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_fasta"))
    log:
        os.path.join(RESULTS_DIR, "logs/contig_taxonomy_mmseqs.log")
    threads:
        config["mmseqs_gtdb"]["threads"]
    message:
        "Running mmseqs2 GTDB taxonomy classification on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e mmseqs mmseqs createdb {input.fasta} {output.dbfasta} && \
        pixi r -e mmseqs mmseqs taxonomy {output.dbfasta} "$(dirname {input.db})/GTDB_mmseqs" "$(dirname {output.taxonomy})/cat_assembly" {output.tmp} --threads {threads} --tax-lineage 1 && \
        pixi r -e mmseqs mmseqs createtsv {output.dbfasta} "$(dirname {output.taxonomy})/cat_assembly" {output.taxonomy} && \
        date) &> {log}
        """

rule taxconverter_convert:
    input:
        tax=rules.contig_taxonomy_mmseqs.output.taxonomy
    output:
        temp(os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_tax_temp.tsv"))
    log:
        os.path.join(RESULTS_DIR, "logs/taxconverter.log")
    message:
        "Converting mmseqs2 taxonomy to TaxVAMB's unified format"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e taxconverter taxconverter mmseqs2 -i {input.tax} -o {output} && \
        date) &> {log}
        """

rule tax_add_missing:
    input:
        fasta=CAT_ASSEMBLY_FILTER,
        tax=rules.taxconverter_convert.output
    output:
        os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_tax_good.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/tax_add_missing.log")
    message:
        "Filling in 'unknown' taxonomy for any contig mmseqs2 didn't classify"
    run:
        import pandas as pd

        contigs = []
        with open(input.fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    contigs.append(line[1:].strip().split()[0])
        tax_df = pd.read_csv(input.tax, sep="\t", header=None, names=["contigs", "predictions"])
        all_df = pd.DataFrame({"contigs": contigs})
        merged = all_df.merge(tax_df, on="contigs", how="left")
        merged["predictions"] = merged["predictions"].fillna("unknown")
        merged.to_csv(output[0], sep="\t", index=False, header=True)

rule vamb_length:
    input:
        CAT_ASSEMBLY_FILTER
    output:
        temp(os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_filter_vamb.fasta"))
    log:
        os.path.join(RESULTS_DIR, "logs/vamb_length.log")
    threads:
        4
    message:
        "Applying VAMB's own >=1501bp contig length requirement"
    shell:
        """
        (date && {config[seqkit][bin]} seq -j {threads} --remove-gaps -o {output} -m 1501 {input} && date) &> {log}
        """

rule vamb_filter_taxonomy:
    input:
        tax=rules.tax_add_missing.output,
        fasta=rules.vamb_length.output
    output:
        temp(os.path.join(RESULTS_DIR, "bins/taxonomy/cat_assembly_tax_vamb.tsv"))
    log:
        os.path.join(RESULTS_DIR, "logs/vamb_filter_taxonomy.log")
    message:
        "Subsetting taxonomy to only the >=1501bp contigs VAMB will see"
    run:
        import pandas as pd

        ids = []
        with open(input.fasta, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].strip())
        df = pd.read_csv(input.tax, sep="\t")
        filtered_df = df[df.iloc[:, 0].isin(ids)]
        filtered_df.to_csv(output[0], sep="\t", index=False, header=True)

rule vamb_depth:
    # Separate from metabat2_mapping's own depth file (rules/binning.smk):
    # --noIntraDepthVariance drops the per-sample variance columns VAMB
    # doesn't want, matching the reference repo's clean_coverage rule.
    input:
        bam=EXISTING_BAMS
    output:
        temp=temp(os.path.join(RESULTS_DIR, "bins/taxonomy/vamb_depth_temp.txt")),
        final=os.path.join(RESULTS_DIR, "bins/taxonomy/vamb_depth_cov.txt")
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    threads:
        config["vamb"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/vamb_depth.log")
    message:
        "Computing VAMB-specific coverage depth (no intra-depth variance)"
    shell:
        """
        (date && jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output.temp} {input.bam} && \
        cat {output.temp} | awk '{{if ($2>1500) print $0 }}' | cut -f -1,4- > {output.final} && \
        date) &> {log}
        """

rule vamb_coverage:
    input:
        rules.vamb_depth.output.final
    output:
        os.path.join(RESULTS_DIR, "bins/taxonomy/vamb_aemb_cov.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/vamb_coverage.log")
    message:
        "Reformatting depth file to VAMB's expected coverage table"
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")
        df.rename(columns={"contigName": "contigname"}, inplace=True)
        df = df.loc[:, (df != 0).any(axis=0)]
        df.to_csv(output[0], sep="\t", index=False, header=True)

rule taxvamb:
    input:
        contig=rules.vamb_length.output,
        tax=rules.vamb_filter_taxonomy.output,
        cov=rules.vamb_coverage.output
    output:
        os.path.join(RESULTS_DIR, "bins/taxvamb/vaevae_clusters_unsplit.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/taxvamb.log")
    threads:
        config["vamb"]["threads"]
    message:
        "Running TaxVAMB on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf "$(dirname {output})" && \
        pixi r -e vamb vamb bin taxvamb --outdir "$(dirname {output})" -m 1501 --fasta {input.contig} --taxonomy {input.tax} --minfasta 500000 --abundance_tsv {input.cov} -p {threads} -o && \
        date) &> {log}
        """


############################################
# COMEBin -- contrastive-learning-based binner
rule comebin:
    input:
        cont=CAT_ASSEMBLY_FILTER,
        bam=EXISTING_BAMS
    output:
        os.path.join(RESULTS_DIR, "bins/comebin/comebin_res/comebin_res.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/comebin.log")
    threads:
        config["comebin"]["threads"]
    params:
        bam_dir=os.path.join(RESULTS_DIR, "bam_flat")
    message:
        "Running COMEBin on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        mkdir -p {params.bam_dir} && \
        for b in {input.bam}; do ln -sf "$b" {params.bam_dir}/; ln -sf "$b.bai" {params.bam_dir}/ 2>/dev/null || true; done && \
        rm -rf "$(dirname $(dirname {output}))" && \
        set +e && \
        pixi r -e comebin run_comebin.sh -a {input.cont} -o "$(dirname $(dirname {output}))" -t {threads} -p {params.bam_dir} && \
        status=$? && set -e && \
        if [[ ! -s {output} ]]; then \
            mkdir -p "$(dirname {output})" && : > {output} && \
            echo "COMEBin produced no bins; created empty output file {output}"; \
        elif [[ $status -ne 0 ]]; then exit $status; fi && \
        date) &> {log}
        """


############################################
# SemiBin2 -- self-supervised binner. Reuses the pooled assembly + existing
# BAMs; single_easy_bin (not multi_easy_bin) since coverage already spans
# all 147 samples via one pooled reference, not per-sample assemblies.
rule semibin:
    input:
        cont=CAT_ASSEMBLY_FILTER,
        bam=EXISTING_BAMS
    output:
        os.path.join(RESULTS_DIR, "bins/semibin/contig_bins.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/semibin.log")
    threads:
        config["semibin"]["threads"]
    message:
        "Running SemiBin2 on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && cd {config[work_dir]} && \
        rm -rf "$(dirname {output})" && \
        pixi r -e semibin SemiBin2 single_easy_bin -i {input.cont} -b {input.bam} -o "$(dirname {output})" -t {threads} --no-recluster --minfasta 200000 && \
        touch {output} && \
        date) &> {log}
        """


############################################
# MetaCAT -- seed-gene + coverage/composition clustering binner
rule metacat_install:
    output:
        done=touch(os.path.join(DB_DIR, "metacat/metacat.installed"))
    log:
        os.path.join(RESULTS_DIR, "logs/metacat_install.log")
    message:
        "Installing MetaCAT (pip wheel, via pixi)"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e metacat MetaCAT --version && \
        date) &> {log}
        """

rule metacat_cov:
    input:
        db=rules.metacat_install.output.done,
        bam=EXISTING_BAMS
    output:
        os.path.join(RESULTS_DIR, "bins/metacat/metacat.cov")
    log:
        os.path.join(RESULTS_DIR, "logs/metacat_cov.log")
    threads:
        config["metacat"]["threads"]
    message:
        "Running MetaCAT coverage"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        mkdir -p "$(dirname {output})" && \
        pixi r -e metacat MetaCAT coverage -b {input.bam} -ti {threads} -tc {threads} -o {output} && \
        date) &> {log}
        """

rule metacat_seed:
    input:
        db=rules.metacat_install.output.done,
        cont=CAT_ASSEMBLY_FILTER
    output:
        os.path.join(RESULTS_DIR, "bins/metacat/metacat.seed")
    log:
        os.path.join(RESULTS_DIR, "logs/metacat_seed.log")
    threads:
        config["metacat"]["threads"]
    message:
        "Running MetaCAT seed"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e metacat MetaCAT seed -f {input.cont} -o {output} -t {threads} && \
        date) &> {log}
        """

rule metacat_cluster:
    input:
        db=rules.metacat_install.output.done,
        cont=CAT_ASSEMBLY_FILTER,
        seed=rules.metacat_seed.output,
        cov=rules.metacat_cov.output
    output:
        folder=directory(os.path.join(RESULTS_DIR, "bins/metacat/cat_assembly_metacat/Bins")),
        done=touch(os.path.join(RESULTS_DIR, "bins/metacat/metacat.done"))
    log:
        os.path.join(RESULTS_DIR, "logs/metacat_cluster.log")
    threads:
        config["metacat"]["threads"]
    message:
        "Running MetaCAT clustering"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.folder} && mkdir -p {output.folder} && \
        pixi r -e metacat MetaCAT cluster -t {threads} -f {input.cont} -s {input.seed} -c {input.cov} -o {output.folder}/cat_assembly_metacat && \
        date) &> {log}
        """
