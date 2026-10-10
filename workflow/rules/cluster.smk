"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/cluster.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To prepare a concatenated assembly for binning
"""


############################################
rule cluster:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta.sa"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR,"bam/{sid}/cat_assembly_{sid}.bam"), sid=SAMPLES.index)
    output:
        touch("status/cluster.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Adding filename to assemblies for easy tracking
rule modify_fasta:
    # barrier: don't start for ANY sample until gene calling (Prodigal) has
    # finished for ALL samples -- see megahit's barrier comment in
    # rules/assembly.smk for the phase-ordering rationale. Also requires
    # status/eggnog.done, status/functions.done, and status/amr.done: this
    # rule is binning's entry point, and user-requested phase ordering is
    # singlem -> eggnog -> functions (incl. AMR) -> binning -> antismash
    # (see emapper's barrier in rules/eggnog.smk and pre_antismash_barrier
    # in rules/bgc_amr.smk for the rest of the chain). Safe to add this
    # late (2026-10-10): binning hadn't started yet (0/7 binners) when
    # this barrier was tightened, so there's no completed work to
    # invalidate -- see the singlem->eggnog comment on rule emapper for
    # why a hard barrier added to an already-in-progress pipeline is
    # risky in general.
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        barrier="status/annotation.done",
        eggnog_barrier="status/eggnog.done",
        functions_barrier="status/functions.done",
        amr_barrier="status/amr.done"
    output:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/modify_fasta/{sid}.log")
    message:
        "Adding filename to fasta headers: {wildcards.sid}"
    shell:
        """(date && filename={wildcards.sid} && awk -v filename="$filename" '/^>/ {{ sub(">", ">" filename ":"); }} 1' {input.fasta} > {output} && date) &> >(tee {log})"""

## Concatetnating the assemblies for binning 
#rule ass_cat:
#    input:
#        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta"), sid=SAMPLES.index)
#    output:
#        os.path.join(RESULTS_DIR, "assembly/cat_assembly.fasta")
#    log:
#        os.path.join(RESULTS_DIR, "logs/assembly/concatenation.log")
#    message:
#        "Concatenating all assemblies"
#    shell:
#        "(date && cat {input} > {output} && date) &> >(tee {log})"
#
## Clustering the assemblies for binning
#rule cat_mmseqs2:
#    input:
#        rules.ass_cat.output
#    output:
#        os.path.join(RESULTS_DIR, "mmseqs/cat_assembly_rep_seq.fasta")
#    log:
#        os.path.join(RESULTS_DIR, "logs/mmseqs_cluster.log")
#    conda:
#        os.path.join(ENV_DIR, "mmseqs2.yaml")
#    threads:
#        config["mmseqs2"]["threads"]
#    params:
#        tmp=os.path.join(RESULTS_DIR, "tmpdir"),
#        cov=config["mmseqs2"]["cov"],
#        cov_mode=config["mmseqs2"]["cov_mode"],
#        min_id=config["mmseqs2"]["min_id"]
#    shell:
#        "(date && "
#        "mmseqs easy-linclust --threads {threads} --force-reuse 0 --cov-mode {params.cov_mode} --min-seq-id {params.min_id} -c {params.cov} {input} $(echo {output} | sed 's,_rep_seq.fasta,,g') {params.tmp} && "
#        "date) &> >(tee {log})"

# Clustering the assemblies for binning
rule ass_mmseqs2:
    input:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta")
    output:
        os.path.join(RESULTS_DIR, "mmseqs/{sid}/{sid}_rep_seq.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/mmseqs/{sid}_cluster.log")
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    threads:
        config["mmseqs2"]["threads"]
    params:
        tmp=os.path.join(RESULTS_DIR, "tmpdir"),
        cov=config["mmseqs2"]["cov"],
        cov_mode=config["mmseqs2"]["cov_mode"],
        min_id=config["mmseqs2"]["min_id"]
    shell:
        "(date && "
        "mmseqs easy-linclust --threads {threads} --force-reuse 0 --cov-mode {params.cov_mode} --min-seq-id {params.min_id} -c {params.cov} {input} $(echo {output} | sed 's,_rep_seq.fasta,,g') {params.tmp} && "
        "date) &> >(tee {log})"

# Concatetnating the non-redundant assemblies for binning
rule cat_ass_mmseqs2:
    input:
        expand(os.path.join(RESULTS_DIR, "mmseqs/{sid}/{sid}_rep_seq.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/assembly/concatenation.log")
    message:
        "Concatenating all assemblies"
    shell:
        "(date && cat {input} > {output} && date) &> >(tee {log})"

rule cat_filter_length:
    # seqkit instead of the old hand-rolled filter_contigs.py (2026-10-10):
    # multithreaded via -j, and -m filters purely by sequence length with
    # no header parsing at all -- the whole class of header-format bugs
    # filter_contigs.py existed to work around (">{SAMPLE_ID}:{contig}"
    # vs. bare headers) just doesn't apply to seqkit.
    input:
        rules.cat_ass_mmseqs2.output[0]
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta")
    threads:
        config['filter_length']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/filter/length.log")
    params:
        bin=config["seqkit"]["bin"]
    message:
        "Removing contigs less than 1.5K"
    shell:
        "(date && {params.bin} seq -j {threads} -m 1500 -o {output} {input} && date) &> >(tee {log})"

rule cat_bin_mapping_index:
    input:
        rules.cat_filter_length.output
    output:
        os.path.join(RESULTS_DIR,"assembly/cat_assembly_filter.fasta.sa")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping/mapping.bwa.index.log")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    message:
        "Mapping: BWA index for assembly mapping"
    shell:
        "(date && bwa index {input} && date) &> >(tee {log})"    

rule filter_mapping:
    input:
        read1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        read2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"),
        cont=rules.cat_filter_length.output,
        idx=rules.cat_bin_mapping_index.output
    output:
        os.path.join(RESULTS_DIR,"bam/{sid}/cat_assembly_{sid}.bam")
    threads:
        config["mapping"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/cat_assembly_{sid}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/cat_assembly_{sid}.err.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running bwa to produce sorted bams: {wildcards.sid}"
    shell:
        "(date && bwa mem -t {threads} {input.cont} {input.read1} {input.read2} | samtools sort -@{threads} -o {output} - && " 
        "samtools index {output} && date) 2> {log.err} > {log.out}"
