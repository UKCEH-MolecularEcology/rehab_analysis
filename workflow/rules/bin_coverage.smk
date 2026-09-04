"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/bin_coverage.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Getting coverage for all bins
"""


############################################
rule bin_coverage:
    input:
        os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_cov.txt"),
        os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage.txt")
    output:
        touch("status/bin_coverage.done")


############################################
# localrules: 


############################################
# CONCATENATING MAGS
rule concatenate_mags:
    input:
        mags=os.path.join(RESULTS_DIR, "bins/all_finalbins")    # os.path.join(RESULTS_DIR, "bins/finalbins")
    output:
        os.path.join(RESULTS_DIR, "concat/cat_mags.fa")
    message:
        "Concatenating all dereplicated MAGS"
    shell:
        "(date && cat {input}/*.fa > {output} && "
        "sed -i 's/>>s/>s/g' {output} && date)"

rule mags_mapping_all:
    input:
        read1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        read2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"),
        cont=rules.concatenate_mags.output
        #idx=rules.mapping_index.output
    output:
        temp(os.path.join(RESULTS_DIR, "bam/cat_mags.fa.{sid}_filtered.R1.fq.bam"))
    threads:
        config["mapping"]["threads"]
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    params:
        tmp=os.path.join(RESULTS_DIR, "tmpdir")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/cat_mags_{sid}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/cat_mags_{sid}.err.log")
    message:
        "Running coverm to produce sorted bams for: {wildcards.sid}"
    shell:
        "(date && "
        "TMPDIR={params.tmp} coverm make -r {input.cont} -t {threads} -o $(dirname {output}) --discard-unmapped -c {input.read1} {input.read2} && "
        "date) 2> {log.err} > {log.out}"

# COVERAGE
rule cat_mag_coverage:
    input:
        os.path.join(RESULTS_DIR, "bam/cat_mags.fa.{sid}_filtered.R1.fq.bam")
    output:
        final=os.path.join(RESULTS_DIR, "bins/coverage/cat_mags_cov_{sid}.txt")
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")        
    threads:
        config["mapping"]["threads"]
    priority: 
        50
    log:
        out=os.path.join(RESULTS_DIR, "logs/coverage/cat_mag_coverage_{sid}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/coverage/cat_mag_coverage_{sid}.err.log")
    message:
        "Estimating coverage for concatenated MAGs: {wildcards.sid}"
    shell:
        "(date && coverm contig -b {input} -m trimmed_mean -t {threads} -o {output.final} && date) 2> {log.err} > {log.out}"

rule concat_cat_mag_coverage:
    input:
        files = expand(os.path.join(RESULTS_DIR, "bins/coverage/cat_mags_cov_{sid}.txt"), sid= SAMPLES.index),
        script=os.path.join(SRC_DIR, "mergeTables.pl")
    output:
        cov=os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_cov.txt") 
    log:
        out=os.path.join(RESULTS_DIR, "logs/concat_cat_mag_coverage.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concat_cat_mag_coverage.err.log")
    message:
        "Concatenating the coverage files"
    shell:
        "(date && perl {input.script} {input.files} {output.cov} && date) 2> {log.err} > {log.out}"

rule get_contigs:
    input:
        bins=os.path.join(RESULTS_DIR, "bins/all_finalbins"),
        src=os.path.join(SRC_DIR, "get_all_contigs.pl")
    output:
        contigs=os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/coverage/get_contigs.log")
    conda:
        os.path.join(ENV_DIR, "perl.yaml")
    message:
        "Getting list of all contigs from all final bins"
    shell:
        "(date && perl {input.src} {input.bins} > {output.contigs} && date) &> {log}"

rule merge_contigs_coverage:
    input:
        contigs=rules.get_contigs.output.contigs,
        cov=rules.concat_cat_mag_coverage.output.cov
    output:
        merged=os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_all_contigs_coverage.txt"),
        grouped=os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/coverage/merge_contigs_coverage.log")
    message:
        "Merging the contigs and coverage files for finalbins"
    run:
        import pandas as pd
        
        # reading in the contigs list file
        contigs=pd.read_csv(input.contigs, sep="\t")

        # reading in the coverage file
        coverage=pd.read_csv(input.cov, sep="\t")

        # Merge the two dataframes on the '0_Contig' column
        merged_data=pd.merge(contigs, coverage, on='0_Contig', how='inner')

        # Print the merged data
        print(merged_data)

        # Save the merged data to a new file
        merged_data.to_csv(output.merged, sep="\t", index=False) 

        # Drop '0_Contig' and 'length' columns
        merged_data=merged_data.drop(['0_Contig', 'length'], axis=1)

        # Summarize based on 'bin' column
        summary_df=merged_data.groupby('bin').sum().reset_index()

        # Print the summary DataFrame
        print(summary_df)

        # Save the summary DataFrame to a new CSV file
        summary_df.to_csv(output.grouped, sep="\t", index=False)
