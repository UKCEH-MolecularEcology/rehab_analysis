"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/bin_coverage.smk --use-conda --cores 4 -rp
Latest modification: [2026-10-10] Replaced concatenate_mags -> mags_mapping_all
    (147 fresh bwa+coverm mapping jobs) -> cat_mag_coverage (147 more) with
    a single coverm genome --genome-definition call reusing EXISTING_BAMS
    (already computed for binning against the pooled cat_assembly_filter.fasta
    in rules/mags_generation.smk) -- CoverM can compute per-genome coverage
    directly from a BAM aligned to a superset reference plus a contig->genome
    map, no second alignment pass needed ("When specified, no read mapping
    algorithm is undertaken" -- confirmed via `coverm genome --full-help`).
Purpose: Getting coverage for all bins
"""


############################################
rule bin_coverage:
    input:
        os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage.txt")
    output:
        touch("status/bin_coverage.done")


############################################
localrules: genome_definition


############################################
rule get_contigs:
    input:
        bins=rules.galah.output.final,
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


rule genome_definition:
    # CoverM's --genome-definition wants "genome_name<TAB>contig", no
    # header -- get_contigs' own output is already "bin\t0_Contig\tlength"
    # (see get_all_contigs.pl), so this is just dropping the header and
    # the length column.
    input:
        contigs=rules.get_contigs.output.contigs
    output:
        os.path.join(RESULTS_DIR, "bins/coverage/genome_definition.tsv")
    message:
        "Building the bin->contig genome-definition file for CoverM"
    shell:
        "tail -n +2 {input.contigs} | cut -f1,2 > {output}"


rule coverm_genome_coverage:
    input:
        bam=EXISTING_BAMS,
        ref=CAT_ASSEMBLY_FILTER,
        definition=rules.genome_definition.output
    output:
        raw=temp(os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage_raw.txt"))
    conda:
        os.path.join(ENV_DIR, "coverm.yaml")
    threads:
        config["coverm_genome"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/coverage/coverm_genome_coverage.log")
    message:
        "Computing per-MAG coverage across all samples (reusing the existing binning BAMs)"
    shell:
        "(date && coverm genome -b {input.bam} --reference {input.ref} --genome-definition {input.definition} "
        "-m trimmed_mean -t {threads} -o {output.raw} && date) &> {log}"


rule clean_genome_coverage:
    # CoverM names each BAM's column after its filename
    # ("cat_assembly_{sid}.bam Trimmed Mean") -- strip that down to the
    # bare sample id, matching the old mergeTables.pl-based table's column
    # naming, and rename the genome column to "bin" to match every other
    # finalbins table in this project.
    input:
        rules.coverm_genome_coverage.output.raw
    output:
        os.path.join(RESULTS_DIR, "bins/coverage/all_finalbins_contigs_coverage.txt")
    message:
        "Cleaning up CoverM's per-MAG coverage column names"
    run:
        import re
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")
        df = df.rename(columns={df.columns[0]: "bin"})
        df.columns = [
            re.sub(r"^cat_assembly_(.+)\.bam Trimmed Mean$", r"\1", c) if c != "bin" else c
            for c in df.columns
        ]
        df.to_csv(output[0], sep="\t", index=False)
