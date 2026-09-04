"""
Author: Susheel Bhanu BUSI & Amy Thorpe
Affiliation: Molecular Ecology group, UKCEH
Date: [2024-07-01]
Run: snakemake -s workflow/rules/microtraint.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run microtrait on MAGs
"""

MAGS=glob_wildcards(os.path.join(DATA_DIR, "mags/{mags}.fa")).mags


############################################
rule microtrait:
    input:
        os.path.join(RESULTS_DIR, "microtrait/genomeset_results.rds")  
    output:
        touch("status/microtrait.done")


############################################
localrules: 


############################################
# rule install_microtrait:
#     output:
#         done = os.path.join(RESULTS_DIR, "microtrait/microtrait.installed")
#     log:
#         out = os.path.join(RESULTS_DIR, "logs/setup.microtrait.log")
#     conda:
#         "/hdd0/susbus/tools/conda_envs/r_env"
# #        os.path.join(ENV_DIR, "microtrait.yaml")
#     message:
#         "installing microtrait"
#     script:
#         os.path.join(SRC_DIR, "install_microtrait.R")

rule run_microtrait:
    input:
        mags = os.path.join(RESULTS_DIR, "bins/all_finalbins")
#        done = os.path.join(RESULTS_DIR, "microtrait/microtrait.installed")
    output:
        rds = os.path.join(RESULTS_DIR, "microtrait/genomeset_results.rds")
    conda:
        "/hdd0/susbus/tools/conda_envs/r_env"
#        os.path.join(ENV_DIR, "microtrait.yaml")
    log:
        out = os.path.join(RESULTS_DIR, "logs/microtrait.log")
    message:
        "running microtrait"
    script:
        os.path.join(SRC_DIR, "run_microtrait.R")
