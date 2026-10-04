"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Purpose: Optional module -- recover and quality-assess EUKARYOTIC bins that
         the main (prokaryote-focused) DAS_Tool/Rosella-refine/Galah chain
         in rules/dereplicate.smk would otherwise discard. Adapted from
         https://github.com/michoug/MAGsGeneration's eukaryotes_filter.smk.

Same pooled-assembly adaptation as rules/mags_generation.smk: the
reference repo loops this whole module per-sample; here there is only
one pooled bin set per binner, so every rule below runs once rather
than per-sample. Matches the reference's own choice of which 5 raw
binner outputs to scan for eukaryotic signal (MetaBAT2, Rosella,
CONCOCT, SemiBin2, TaxVAMB -- COMEBin/MetaCAT were not part of its
eukaryote-candidate pool either).

Two tracks feed the final eukaryotic MAG set, exactly as in the
reference:
  1. DeepMicroClass: classify every contig in every RAW candidate bin
     (from the 5 binners above, before DAS_Tool/refine) as eukaryotic
     or prokaryotic; keep bins where >=80% of sequence is eukaryotic.
  2. REMAG: a binner specifically designed to recover eukaryotic genomes
     missed by prokaryote-focused binning, run directly on the pooled
     assembly + its coverage; keep bins >=2.5Mb.
Both tracks' candidates get BUSCO (eukaryota lineage) quality assessment
and a final Galah dereplication pass.
"""


############################################
rule eukaryotes_filter:
    input:
        os.path.join(RESULTS_DIR, "bins/euk/remag_bins/done"),
        os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags_deepmicroclass.txt"),
        os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags_busco_done.txt"),
        os.path.join(RESULTS_DIR, "bins/euk/Eukaryotic_mags_dereplicated/done.txt"),
    output:
        touch("status/eukaryotes_filter.done")


############################################
localrules: busco_db_euk


############################################
rule busco_db_euk:
    output:
        dir=directory(config["busco_euk"]["db"]),
        done=touch(os.path.join(config["busco_euk"]["db"], "busco_euk.done"))
    log:
        os.path.join(RESULTS_DIR, "logs/busco_db_euk.log")
    message:
        "Reusing/extending the existing eukaryote BUSCO lineage cache"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        mkdir -p {output.dir} && \
        pixi r -e busco busco --download all --download_path {output.dir} && \
        date) &> {log}
        """


############################################
# Collect every raw candidate bin (pre-DAS_Tool/refine) from the 5 binners
# the reference repo scans for eukaryotic signal into one flat folder.
rule get_all_mags_euk:
    input:
        contigs=CAT_ASSEMBLY_FILTER,
        metabat=rules.metabat2.output,
        rosella=rules.rosella.output.dir,
        concoct=rules.concoct_bins.output,
        semibin=rules.semibin.output,
        taxvamb=rules.taxvamb.output,
    output:
        dir=directory(os.path.join(RESULTS_DIR, "bins/euk/candidates")),
        done=touch(os.path.join(RESULTS_DIR, "bins/euk/candidates/done.txt"))
    log:
        os.path.join(RESULTS_DIR, "logs/get_all_mags_euk.log")
    message:
        "Collecting raw candidate bins from all binners for eukaryote screening"
    shell:
        """
        set -euo pipefail
        (date
        mkdir -p {output.dir}
        cp {input.metabat}/*.fa {output.dir}/ 2>/dev/null || true
        for f in {input.rosella}/*.fna; do [ -e "$f" ] || continue; cp "$f" {output.dir}/rosella_$(basename "$f"); done
        for f in {input.concoct}/*.fa; do [ -e "$f" ] || continue; cp "$f" {output.dir}/concoct_$(basename "$f"); done
        for f in "$(dirname {input.taxvamb})"/bins/*.fna; do [ -e "$f" ] || continue; cp "$f" {output.dir}/taxvamb_$(basename "$f" .fna).fa; done
        for f in "$(dirname {input.semibin})"/output_bins/*.fa.gz; do [ -e "$f" ] || continue; cp "$f" {output.dir}/; done
        for f in {output.dir}/*.gz; do [ -e "$f" ] || continue; gunzip -f "$f"; done
        for f in {output.dir}/*.fna; do [ -e "$f" ] || continue; mv "$f" "${{f%.fna}}.fa"; done
        date) &> {log}
        """


rule get_length_euk:
    input:
        dir=rules.get_all_mags_euk.output.dir,
        done=rules.get_all_mags_euk.output.done
    output:
        txt=temp(os.path.join(RESULTS_DIR, "bins/euk/candidates_stats.txt")),
        final=os.path.join(RESULTS_DIR, "bins/euk/candidates_length_filtered_list.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/get_length_euk.log")
    threads:
        4
    message:
        "Length-filtering candidate bins (>=2.5Mb, >=2500 contigs) before DeepMicroClass"
    shell:
        """
        (date && {config[seqkit][bin]} stats -j {threads} {input.dir}/*.fa -T > {output.txt} && \
        awk '$5>=2500000 && $7>=2500' {output.txt} | cut -f1 | tail -n +2 > {output.final} && \
        date) &> {log}
        """


checkpoint get_filtered_mags_euk:
    input:
        list=rules.get_length_euk.output.final,
        dir=rules.get_all_mags_euk.output.dir
    output:
        files=directory(os.path.join(RESULTS_DIR, "bins/euk/length_filtered"))
    log:
        os.path.join(RESULTS_DIR, "logs/get_filter_euk.log")
    message:
        "Copying length-filtered candidates for DeepMicroClass"
    shell:
        """
        set -euo pipefail
        (date
        mkdir -p {output.files}
        cat {input.list} | while read -r line; do
            base=$(basename "$line")
            cp "{input.dir}/$base" {output.files}/$base
        done
        date) &> {log}
        """


def get_file_names_dmc(wildcards):
    ck_output = checkpoints.get_filtered_mags_euk.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "bins/euk/deepmicroclass/collated/{mags}_percentage.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule deepmicroclass:
    input:
        os.path.join(RESULTS_DIR, "bins/euk/length_filtered/{mags}.fa")
    output:
        os.path.join(RESULTS_DIR, "bins/euk/deepmicroclass/{mags}.fa_pred_one-hot_hybrid.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/deepmicroclass/{mags}.log")
    threads:
        3
    message:
        "Running DeepMicroClass on {wildcards.mags}"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e deepmicroclass DeepMicroClass predict --input {input} --output_dir $(dirname {output}) --device cpu --cpu_thread {threads} && \
        date) &> {log}
        """


rule deepmicroclass_percentage:
    input:
        pred=rules.deepmicroclass.output
    output:
        percent=os.path.join(RESULTS_DIR, "bins/euk/deepmicroclass/collated/{mags}_percentage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/deepmicroclass_collate_{mags}.log")
    message:
        "Getting the eukaryotic fraction for {wildcards.mags}"
    run:
        import pandas as pd

        df = pd.read_csv(input.pred, sep="\t", header=0)
        choice_labels = {1: "Eukaryote", 2: "Eukaryote", 3: "Prokaryote", 4: "Prokaryote", 5: "Prokaryote"}
        df["pred"] = df["best_choice"].map(choice_labels)
        count = pd.DataFrame(df.groupby("pred")["pred"].count())
        count["percent"] = (count["pred"] / count["pred"].sum()) * 100
        count_trans = count.transpose()
        count_trans["bin"] = os.path.basename(input.pred).replace("_pred_one-hot_hybrid.tsv", "")
        count_trans.drop(labels="pred").to_csv(output.percent, sep="\t", header=True)


rule concatenate_deepmicroclass:
    input:
        get_file_names_dmc
    output:
        collated=os.path.join(RESULTS_DIR, "bins/euk/deepmicroclass/collated_percentages.txt"),
        final_list=os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags_deepmicroclass.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/concatenate_deepmicroclass.log")
    message:
        "Collecting all DeepMicroClass percentages; flagging bins >=80% eukaryotic"
    run:
        import pandas as pd

        dfs = [pd.read_csv(f, sep="\t", header=0) for f in input]
        concatenated = pd.concat(dfs, ignore_index=True)
        concatenated.to_csv(output.collated, sep="\t", index=False)
        eukaryotic_bins = concatenated[concatenated["Eukaryote"] >= 80]
        eukaryotic_bins[["bin"]].to_csv(output.final_list, sep="\t", index=False, header=False)


############################################
# REMAG -- a second, independent track for recovering eukaryotic genomes,
# run directly on the pooled assembly + its existing MetaBAT2-style depth.
rule remag:
    input:
        contig=CAT_ASSEMBLY_FILTER,
        tsv=rules.metabat2_mapping.output
    output:
        folder=directory(os.path.join(RESULTS_DIR, "bins/euk/remag")),
        file=os.path.join(RESULTS_DIR, "bins/euk/remag/bins.csv"),
        bins=directory(os.path.join(RESULTS_DIR, "bins/euk/remag/bins"))
    log:
        os.path.join(RESULTS_DIR, "logs/remag.log")
    threads:
        config["remag"]["threads"]
    message:
        "Running REMAG on the pooled assembly"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        pixi r -e remag remag --fasta {input.contig} --coverage {input.tsv} -o {output.folder} -t {threads} && \
        date) &> {log}
        """


rule remag_bins_flat:
    input:
        rules.remag.output.file
    output:
        folder=directory(os.path.join(RESULTS_DIR, "bins/euk/remag_bins")),
        done=os.path.join(RESULTS_DIR, "bins/euk/remag_bins/done")
    log:
        os.path.join(RESULTS_DIR, "logs/remag_concat.log")
    message:
        "Flattening REMAG bins into one folder"
    shell:
        """
        set -euo pipefail
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        for f in "$(dirname {input})/bins/"*.fa; do
            [ -e "$f" ] || continue
            cp "$f" {output.folder}/$(basename "$f")
        done
        touch {output.done}
        date) &> {log}
        """


rule get_length_remag:
    input:
        rules.remag_bins_flat.output.folder
    output:
        txt=temp(os.path.join(RESULTS_DIR, "bins/euk/remag_stats.txt")),
        final=os.path.join(RESULTS_DIR, "bins/euk/remag_length_filtered_list.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/get_length_remag.log")
    threads:
        4
    shell:
        """
        (date && {config[seqkit][bin]} stats -j {threads} {input}/*.fa -T > {output.txt} && \
        awk '$5>=2500000' {output.txt} | cut -f1 | tail -n +2 > {output.final} && \
        date) &> {log}
        """


############################################
# Combine both tracks (DeepMicroClass->80%-eukaryotic hits + REMAG
# length-filtered candidates) into one final eukaryotic MAG pool.
checkpoint get_eukaryotic_mags:
    input:
        folder=rules.get_filtered_mags_euk.output.files,
        deep=rules.concatenate_deepmicroclass.output.final_list,
        remag_list=rules.get_length_remag.output.final,
        remag_dir=rules.remag_bins_flat.output.folder
    output:
        folder=directory(os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags")),
        done=os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags/done.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/get_eukaryotic_mags.log")
    message:
        "Combining DeepMicroClass and REMAG eukaryotic MAG candidates"
    shell:
        """
        set -euo pipefail
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        cat {input.deep} | while read -r line; do
            base=$(basename "$line")
            cp {input.folder}/$base {output.folder}/$base 2>/dev/null || true
        done
        while read -r line; do
            base=$(basename "$line")
            cp "{input.remag_dir}/$base" {output.folder}/remag_$base 2>/dev/null || true
        done < {input.remag_list}
        touch {output.done}
        date) &> {log}
        """


def get_file_eukaryotes(wildcards):
    ck_output = checkpoints.get_eukaryotic_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "bins/euk/busco/{mags}/{mags}.done"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule busco_euk:
    input:
        folder=os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags/{mags}.fa"),
        db=rules.busco_db_euk.output.dir
    output:
        folder=directory(os.path.join(RESULTS_DIR, "bins/euk/busco/{mags}")),
        done=os.path.join(RESULTS_DIR, "bins/euk/busco/{mags}/{mags}.done")
    log:
        os.path.join(RESULTS_DIR, "logs/busco/{mags}.log")
    threads:
        10
    message:
        "Running BUSCO on eukaryotic MAG {wildcards.mags}"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.folder} && \
        pixi r -e busco busco -i {input.folder} --out_path {output.folder} -m genome --auto-lineage-euk -c {threads} -f --metaeuk --download_path {input.db} --offline \
            || echo "BUSCO failed (possibly no genes found) -- continuing anyway" && \
        touch {output.done} && \
        date) &> {log}
        """


rule busco_done:
    input:
        get_file_eukaryotes
    output:
        os.path.join(RESULTS_DIR, "bins/euk/eukaryotic_mags_busco_done.txt")
    message:
        "All BUSCO runs on eukaryotic MAGs are done"
    shell:
        "touch {output}"


rule busco_multiqc_euk:
    input:
        rules.busco_done.output
    output:
        html=os.path.join(RESULTS_DIR, "qc/busco_euk/multiqc_report.html"),
        stats=os.path.join(RESULTS_DIR, "qc/busco_euk/multiqc_data/multiqc_general_stats.txt"),
        source=os.path.join(RESULTS_DIR, "qc/busco_euk/multiqc_data/multiqc_sources.txt"),
        folder=directory(os.path.join(RESULTS_DIR, "qc/busco_euk"))
    log:
        os.path.join(RESULTS_DIR, "logs/multiqc_busco_euk.log")
    message:
        "Running MultiQC on eukaryotic BUSCO results"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        rm -rf {output.folder} && mkdir -p {output.folder}/multiqc_data && \
        pixi r -e multiqc multiqc $(dirname $(dirname {input[0]}))/busco --verbose --ignore auto_lineage/ -o {output.folder} --force || true && \
        if [ ! -s {output.html} ]; then printf '%s\\n' '<!doctype html>' '<html><body><p>No BUSCO results for eukaryotic MAGs.</p></body></html>' > {output.html}; fi && \
        if [ ! -s {output.stats} ]; then printf "Sample\\tbusco-complete\\tbusco-complete_single_copy\\tbusco-fragmented\\tbusco-missing\\n" > {output.stats}; fi && \
        if [ ! -s {output.source} ]; then printf "Sample Name\\tSource\\n" > {output.source}; fi && \
        date) &> {log}
        """


rule quality_eukaryotes:
    input:
        stats=rules.busco_multiqc_euk.output.stats,
        source=rules.busco_multiqc_euk.output.source,
        mags_done=rules.get_eukaryotic_mags.output.done
    output:
        os.path.join(RESULTS_DIR, "qc/busco_euk.csv")
    log:
        os.path.join(RESULTS_DIR, "logs/quality_eukaryotes.log")
    message:
        "Extracting quality metrics for eukaryotic MAGs from the MultiQC report"
    run:
        from pathlib import Path
        import pandas as pd

        df = pd.read_csv(input.stats, sep="\t")
        src = pd.read_csv(input.source, sep="\t")
        busco_total = df.get("busco-complete_single_copy", 0) + df.get("busco-fragmented", 0) + df.get("busco-missing", 0)
        if "busco-complete_duplicated" in df.columns:
            duplicated = df["busco-complete_duplicated"]
        else:
            duplicated = (df.get("busco-complete", 0) - df.get("busco-complete_single_copy", 0)).clip(lower=0)
        src["sample_key"] = src["Sample Name"].str.strip()
        src["genome"] = src["Source"].str.extract(r"/busco/([^/]+)/", expand=False)
        sample_to_genome = src.dropna(subset=["genome"]).drop_duplicates(subset="sample_key").set_index("sample_key")["genome"]
        df["sample_key"] = df["Sample"].str.strip()
        df["genome"] = df["sample_key"].map(sample_to_genome)
        fallback = df["sample_key"].str.replace(r"^short_summary\.specific\.\w+_odb\d+\.", "", regex=True).str.replace(r"\.fa(?:\.txt)?$", "", regex=True)
        df["genome"] = df["genome"].fillna(fallback)
        df["completeness"] = (df.get("busco-complete", 0) + df.get("busco-fragmented", 0)) * 100 / busco_total.replace(0, pd.NA)
        df["contamination"] = duplicated * 100 / busco_total.replace(0, pd.NA)
        result = df[["genome", "completeness", "contamination"]].sort_values("completeness", ascending=False).drop_duplicates(subset="genome", keep="first")
        mags_folder = Path(input.mags_done).parent
        all_genomes = {p.stem for p in mags_folder.glob("*.fa")}
        missing_genomes = sorted(all_genomes - set(result["genome"].astype(str)))
        if missing_genomes:
            result = pd.concat([result, pd.DataFrame({"genome": missing_genomes, "completeness": 0.0, "contamination": 0.0})], ignore_index=True)
        result.sort_values("genome").to_csv(output[0], index=False)


rule galah_euk:
    input:
        bins=rules.get_eukaryotic_mags.output.folder,
        busco=rules.quality_eukaryotes.output
    output:
        final=directory(os.path.join(RESULTS_DIR, "bins/euk/Eukaryotic_mags_dereplicated")),
        done=os.path.join(RESULTS_DIR, "bins/euk/Eukaryotic_mags_dereplicated/done.txt"),
        cluster=os.path.join(RESULTS_DIR, "bins/euk/cluster_galah_euk.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/galah_euk.log")
    threads:
        10
    message:
        "Dereplicating eukaryotic MAGs with Galah (skips if BUSCO table is empty)"
    shell:
        """
        (date && source {PIXI_ENV_SCRIPT} && \
        if [[ $(wc -l < "{input.busco}") -le 1 ]]; then \
            echo "No BUSCO quality rows available; skipping galah clustering." && \
            mkdir -p {output.final} && : > {output.cluster}; \
        else \
            pixi r -e galah galah cluster --genome-fasta-directory {input.bins} --genome-info {input.busco} -x fa \
                --output-representative-fasta-directory-copy {output.final} -t {threads} --output-cluster-definition {output.cluster}; \
        fi && \
        touch {output.done} && \
        date) &> {log}
        """
