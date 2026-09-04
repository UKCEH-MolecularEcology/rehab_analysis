"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/eggnog.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run eggnog on proteins
"""


############################################
rule eggnog:
    input:
        expand(os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.annotations"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}_antibiotic_biosynthesis.txt"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "eggnog/merged_antibiotic_biosynthesis.txt")
    output:
        touch("status/eggnog.done")


############################################
localrules: download_eggnogDB


############################################
# EggNOG database installation
rule download_eggnogDB:
    output:
        done=os.path.join(RESULTS_DIR, "eggnog/db_download.done")
    log:
        out=os.path.join(RESULTS_DIR, "logs/setup.eggnog_DB.log")
    params:
        path=config["eggnog"]["db"]
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    message:
        "Setup: EggNOG database"
    shell:
        "(date && mkdir -p {params.path} && "
        "download_eggnog_data.py -y --data_dir {params.path} && "
        "touch {output} && date) &> >(tee {log})"

# EGGNOG mapping to annotations
rule emapper:
    input:
        dummy=os.path.join(RESULTS_DIR, "eggnog/db_download.done"),
        fasta=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa")
    output:
        os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.seed_orthologs")
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    threads:
        config["eggnog"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/eggnog/emapper_{sid}.log")
    params:
        itype=config["eggnog"]["itype"],
        genepred=config["eggnog"]["genepred"],
        db=config["eggnog"]["db"]
    message:
        "Running EggNog-mapper on {wildcards.sid}"
    shell:
        "(date && mkdir -p $(dirname {output}) && "
        "emapper.py -m diamond --data_dir {params.db} --itype {params.itype} --no_file_comments --cpu {threads} -i {input.fasta} -o {wildcards.sid} --output_dir $(dirname {output}) && "
        "date) &> >(tee {log})"

# Final annotations
rule emapper_final:
    input:
        rules.emapper.output[0]
    output:
        os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.annotations")
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    threads:
        config["eggnog"]["final_threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/eggnog/{sid}.final_emapper.out.log")
    params:
        itype=config["eggnog"]["itype"],
        genepred=config["eggnog"]["genepred"],
        db=config["eggnog"]["db"]
    message:
        "Running EggNog-mapper annotations on {wildcards.sid}"
    shell:
        "(date && "
        "emapper.py --data_dir {params.db} --annotate_hits_table {input} --no_file_comments -o $(echo {output} | sed 's/.emapper.annotations//g' ) --cpu {threads} --dbmem && "
        "date) &> >(tee {log})"


# ANTIBIOTIC BIOSYNTHESIS GENE SEARCH
rule identify_abx_biosynthesis:
    input:
        annotation=os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.annotations"),
        gene_cov=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_coverage.txt")
    output:
        kegg=os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}_antibiotic_biosynthesis.txt")
    params:
        keywords=config["kegg_ids"]
    message:
        "Searching for Antibiotic biosynthesis genes in {wildcards.sid}"
    run:
        import pandas as pd
        import re
        import os
        
        # Define column names based on eggnog output
        column_names = [
            'query', 'seed_ortholog', 'evalue', 'score', 
            'eggNOG_OGs', 'max_annot_lvl', 'COG_category', 
            'Description', 'Preferred_name', 'GOs', 
            'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 
            'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 
            'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        ]

        # Read the eggNOG annotations file
        annotations_df = pd.read_csv(input.annotation, sep='\t', comment='#', names=column_names, on_bad_lines = 'skip')
        print(annotations_df.head())

        # Create an empty list to store results
        results = []

        # Reading in the keywords
        keywords=params.keywords
        print(keywords)

        # Get sample name from wildcards
        sample_name = wildcards.sid  # Accessing wildcards here is correct

        # Iterate through the annotations DataFrame to find matches
        for index, row in annotations_df.iterrows():
            if pd.notnull(row['KEGG_Pathway']):
                if re.search(keywords, row['KEGG_Pathway']):
                    contig_id = row['query']
                    matching_keywords = re.findall(keywords, row['KEGG_Pathway'])
                    for keyword in matching_keywords:
                        results.append((contig_id, sample_name, keyword))

        # Create a DataFrame from the results
        results_df = pd.DataFrame(results, columns=['contig_id', 'sample_name', 'keyword'])
        print(results_df.head())

        # Read the gene coverage data
        coverage_df = pd.read_csv(input.gene_cov, sep='\s+', header=None, usecols=[1, 2], names=['contig_id', 'coverage'])
        print(coverage_df.head())

        # Merge the results with the coverage data
        merged_df = pd.merge(results_df, coverage_df, on='contig_id', how='left')

        # Save to a file
        merged_df.to_csv(output.kegg, sep='\t', index=False)

rule merge_abx_biosynthesis:
    input:
        expand(os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}_antibiotic_biosynthesis.txt"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "eggnog/merged_antibiotic_biosynthesis.txt")
    run:
        import pandas as pd

        # Create an empty list to store DataFrames
        dfs = []

        # Read each input file and append the DataFrame to the list
        for file in input[0]:
            df = pd.read_csv(file, sep='\t')
            dfs.append(df)

        # Concatenate all DataFrames in the list
        merged_df = pd.concat(dfs, ignore_index=True)

        # Save the merged DataFrame to a file
        merged_df.to_csv(output[0], sep='\t', index=False)

