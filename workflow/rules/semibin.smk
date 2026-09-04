"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/semibin.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run semibin on concatenated assemblies
"""


############################################
rule semibin_all:
    input:
        os.path.join(RESULTS_DIR, "bins/semibin/bins"),
        os.path.join(RESULTS_DIR, "dummy/unzip.semibin.done"),
        os.path.join(RESULTS_DIR, "bins/semibin/finalbins")
    output:
        touch("status/semibin.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Adding filename to assemblies for easy tracking
rule modify_semibin_fasta:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "semibin/concatenated.fa.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/semibin/modify_fasta.log")
    conda:
        os.path.join(ENV_DIR, "semibin.yaml")
    message:
        "Adding filename to fasta headers"
    shell:
        "(date && SemiBin2 concatenate_fasta --input-fasta {input} --output $(dirname {output}) && date) &> >(tee {log})"

rule semibin_filter_length:
    input:
        rules.modify_semibin_fasta.output[0]
    output:
        os.path.join(RESULTS_DIR, "semibin/concatenated_filter.fasta")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    threads:
        config['filter_length']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/semibin/semibin_length.log")
    message:
        "Removing contigs less than 1.5K"
    shell:
        "(date && seqkit seq -j {threads} -o {output} -m 1499 {input} && date) &> >(tee {log})"

rule semibin_mapping_index:
    input:
        rules.semibin_filter_length.output
    output:
        os.path.join(RESULTS_DIR,"semibin/concatenated_filter.fasta.sa")
    log:
        os.path.join(RESULTS_DIR, "logs/semibin/semibin_mapping.bwa.index.log")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    message:
        "Mapping: BWA index for semibin assembly mapping"
    shell:
        "(date && bwa index {input} && date) &> >(tee {log})"    

rule semibin_mapping:
    input:
        read1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        read2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"),
        cont=rules.semibin_filter_length.output,
        idx=rules.semibin_mapping_index.output
    output:
        os.path.join(RESULTS_DIR,"semibin/bam/{sid}/concatenated_filter_{sid}.bam")
    threads:
        config["mapping"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/semibin/concatenated_filter_{sid}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/semibin/concatenated_filter_{sid}.err.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running bwa to produce sorted bams: {wildcards.sid}"
    shell:
        "(date && bwa mem -t {threads} {input.cont} {input.read1} {input.read2} | samtools sort -@{threads} -o {output} - && " 
        "samtools index {output} && date) 2> {log.err} > {log.out}"

rule semibin:
    input:
        contig=rules.semibin_filter_length.output,
        bam=expand(os.path.join(RESULTS_DIR,"semibin/bam/{sid}/concatenated_filter_{sid}.bam"), sid=SAMPLES.index)
    output:
        bins=directory(os.path.join(RESULTS_DIR, "bins/semibin/bins"))
    conda:
        os.path.join(ENV_DIR, "semibin.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/semibin/semibin.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/semibin/semibin.err.log")
    params:
        env=config["semibin"]["env"]
    threads:
        config["semibin"]["threads"]
    message:
        "Running SemiBin2 on concatenated assembly"
    shell:
        "(date && SemiBin2 multi_easy_bin -i {input.contig} -b {input.bam} -o {output.bins} -t {threads} && date) 2> {log.err} > {log.out}"

rule unzip_semibin:
    input:
        rules.semibin.output.bins
    output:
        dummy=os.path.join(RESULTS_DIR, "dummy/unzip.semibin.done")
    log:
        out=os.path.join(RESULTS_DIR, "logs/semibin/semibin.unzip.out"),
        err=os.path.join(RESULTS_DIR, "logs/semibin/semibin.unzip.err")
    message:
        "Unzipping the fasta files"
    shell:
        "(date && cd {input} && "
        """find . -type f -name "*.fa.gz" -exec sh -c 'gunzip -c "$1" > "${1%.gz}"' _ {} \\; && """
        "touch {output.dummy} && date) 2> {log.err} > {log.out}"

rule semibin_checkm:
    input:
        os.path.join(RESULTS_DIR, "bins/semibin/bins")
    output:
        os.path.join(RESULTS_DIR, "bins/semibin_dastool_checkm.tsv")
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool_checkm_semibin.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool_checkm_semibin.err.log")
    threads:
        config["checkm"]["threads"]
    message:
        "Running CheckM on semibin output"
    shell:
        """
        (date && checkm lineage_wf -t {threads} -f {output} --tab_table -x fa {input}/bins $(dirname {output})/semibin_dastool_checkm && date) 2> {log.err} > {log.out}
        """

rule semibin_drep_prepare:
    input:
        check=rules.semibin_checkm.output
    output:
        temp=temp(os.path.join(RESULTS_DIR, "bins/semibin_checkm_temp.tsv")),
        final=os.path.join(RESULTS_DIR, "bins/semibin_checkmbeforedrep.tsv")
    shell:
        """
        awk '{{print $1".fa,"$13","$14}}' {input.check} | tail -n+2 > {output.temp} && 
        (echo "genome,completeness,contamination" && cat {output.temp}) > {output.final}
        """
       
rule semibin_drep:
    input:
        check=rules.semibin_drep_prepare.output.final,
        bins=os.path.join(RESULTS_DIR, "bins/semibin/bins")
    output:
        temp=directory(os.path.join(RESULTS_DIR, "bins/drep/semibin/dereplicated_genomes")),
        final=directory(os.path.join(RESULTS_DIR, "bins/semibin/finalbins"))
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/drep/drep.semibin.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/drep/drep.semibin.err.log")
    params:
        comp=config["drep"]["comp"],
        cont=config["drep"]["cont"]
    threads:
        config["drep"]["threads"]
    message:
        "Running dRep on all semibin bins"
    shell:
        "(date && "
        "dRep dereplicate $(dirname {output.temp}) -p {threads} -comp {params.comp} -con {params.cont} --genomeInfo {input.check} -g {input.bins}/bins/*fa && "
        "cp -r {output.temp} {output.final} && "
        """cd {output.final} && for filename in *.fa; do awk -v fname="$(basename -s '.fa' "$filename")" '/^>/ {{print $0 "_" fname; next}} {{print}}' "$filename" > "temp.fa" && mv "temp.fa" "$filename"; done && """
        "date) 2> {log.err} > {log.out}"
