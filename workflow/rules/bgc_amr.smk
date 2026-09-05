"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-09-04]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 64 -rp
Latest modification:
Purpose: KEGG/eggNOG antibiotic marker summary, BGC prediction (GECCO +
         antiSMASH on per-sample assemblies) and clustering (BiG-SCAPE 2),
         and average maximal growth rate (gRodon2). Adapted from the
         SOCD project's singlem_bgc_Snakefile_no_metadata
         (/prj/DECODE/socd/results/singlem_bgc_Snakefile_no_metadata).

Distinct from rules/antismash.smk: that rule runs antiSMASH on dereplicated
MAGs (post-binning); this one runs GECCO/antiSMASH/BiG-SCAPE directly on
each sample's assembly (pre-binning), matching the SOCD comparison of
AMR (KEGG markers) against BGC calls per sample.

GECCO, BiG-SCAPE and antiSMASH here reuse the pre-existing named conda
envs on this host (/home/susbus/miniforge3/envs/{gecco,bigscape,antismash})
rather than rebuilding them via `--use-conda`, per project convention of
reusing what's already installed. gRodon2 has no existing install, so it
runs via a Singularity image pulled into containers/ (project-local, not
$HOME or /tmp -- see config["work_dir"]).
"""

import json
import re

############################################
# NOTE: samples are NOT pre-filtered by whether eggNOG/coverage output
# already exists on disk (a prior version did this via a parse-time
# os.path.isfile check, copied from the SOCD original -- but that project
# ran this step in a separate snakemake invocation *after* eggNOG/coverage
# were already fully computed by an earlier pipeline run, so the check
# meant something there). Since this repo runs the whole pipeline
# (assembly -> annotation -> eggnog/coverage -> bgc_amr) in a single
# invocation, a parse-time filesystem check always sees an empty result at
# the start -- which made every expand(..., sid=SAMPLES.index) resolve to an
# empty input list, so Snakemake treated those rules as having no
# prerequisites and ran them immediately, before antiSMASH/eggNOG/coverage
# had produced anything. Using SAMPLES.index directly lets Snakemake's own
# DAG resolution require the real per-sample eggnog/coverage/antismash
# outputs as normal rule inputs, enforcing the correct order.
BGC_OUTDIR = os.path.join(RESULTS_DIR, "amr_bgc")

GECCO_BIN = config["gecco"]["bin"]
BIGSCAPE_BIN = config["bigscape"]["bin"]
ANTISMASH_BIN = config["antismash_assembly"]["bin"]
SEQKIT_BIN = config["seqkit"]["bin"]

THREADS_ANTISMASH_ASSEMBLY = config["antismash_assembly"]["threads"]
ANTISMASH_ASSEMBLY_MIN_CONTIG_LEN = config["antismash_assembly"]["min_contig_len"]
THREADS_BIGSCAPE = config["bigscape"]["threads"]
PFAM_DB = config["bigscape"]["pfam_db"]

GRODON_IMAGE = config["grodon"]["image"]
GRODON_DIR = os.path.join(BGC_OUTDIR, "growth_rate")
GRODON_SIF = os.path.join(config["work_dir"], "containers", "grodon2.sif")
GRODON_SCRIPT = os.path.join(GRODON_DIR, "scripts", "run_gRodon.R")


############################################
rule bgc_amr:
    input:
        os.path.join(BGC_OUTDIR, "abx_bgc", "abx_bgc_master_results.tsv"),
        expand(os.path.join(BGC_OUTDIR, "gecco", "{sid}", "{sid}.genes.tsv"), sid=SAMPLES.index),
        os.path.join(BGC_OUTDIR, "bigscape_results"),
        os.path.join(BGC_OUTDIR, "antismash_summaries", "antismash_master_results.tsv"),
        os.path.join(GRODON_DIR, "growth_rate_master_results.tsv")
    output:
        touch("status/bgc_amr.done")


############################################
# 1. Antibiotic Marker Analysis (KEGG / EggNOG), curated reference table
############################################
rule summarise_sample_abx:
    """Cross-reference eggNOG annotations with the curated KEGG ABX/BGC
    marker reference (resources/kegg_abx_bgc.txt) + gene coverage.
    Complements rules/eggnog.smk's identify_abx_biosynthesis, which matches
    KEGG_Pathway by keyword regex; this instead joins on exact KO id against
    a small curated marker set (biosynthesis + resistance genes)."""
    input:
        annot=os.path.join(RESULTS_DIR, "eggnog", "{sid}", "{sid}.emapper.annotations"),
        cov=os.path.join(RESULTS_DIR, "coverage", "{sid}", "{sid}_gene_coverage.txt"),
        ref=config["abx_bgc"]["kegg_ref"]
    output:
        tsv=temp(os.path.join(BGC_OUTDIR, "abx_bgc", "{sid}_summary.tsv"))
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        """
        mkdir -p "$(dirname {output.tsv})"
        awk -v ref="{input.ref}" -v cov_file="{input.cov}" -v sname="{wildcards.sid}" '
            BEGIN {{ FS="\\t"; OFS="\\t" }}
            FILENAME == ref {{
                if (FNR > 1) meta[$1] = $2 "\\t" $3 "\\t" $4 "\\t" $5
                next
            }}
            FILENAME == cov_file {{
                split($0, a, /[[:space:]]+/)
                coverage[a[2]] = a[3]
                next
            }}
            /^#/ {{ next }}
            {{
                query = $1; ko_field = $12
                if (ko_field == "-") next
                split(ko_field, kos, ",")
                for (i in kos) {{
                    clean_ko = kos[i]
                    sub(/^ko:/, "", clean_ko)
                    if (clean_ko in meta) {{
                        val = (query in coverage) ? coverage[query] : "0"
                        print sname, query, clean_ko, val, meta[clean_ko]
                    }}
                }}
            }}
        ' "{input.ref}" "{input.cov}" "{input.annot}" > {output.tsv}
        """

rule master_abx_summary:
    """Merge all per-sample ABX summaries and add a header."""
    input:
        expand(os.path.join(BGC_OUTDIR, "abx_bgc", "{sid}_summary.tsv"), sid=SAMPLES.index)
    output:
        os.path.join(BGC_OUTDIR, "abx_bgc", "abx_bgc_master_results.tsv")
    shell:
        """
        echo -e "Sample\\tQuery\\tKO\\tCoverage\\tType\\tFunction\\tClassification\\tGroup" > {output}
        cat {input} >> {output}
        """


############################################
# 2. BGC Prediction (GECCO + antiSMASH, per-sample assembly)
############################################
rule run_gecco:
    """Predict BGCs with GECCO.

    barrier: don't start for ANY sample until gene calling (Prodigal) has
    finished for ALL samples -- see megahit's barrier comment in
    rules/assembly.smk for the phase-ordering rationale."""
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly", "{sid}", "{sid}.fasta"),
        barrier="status/annotation.done"
    output:
        os.path.join(BGC_OUTDIR, "gecco", "{sid}", "{sid}.genes.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs", "gecco", "{sid}_gecco.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    params:
        bin=GECCO_BIN
    shell:
        """
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        (date && {params.bin} run --genome {input.fasta} -o $(dirname {output}) && date) &> {log}
        """

rule pre_antismash_barrier:
    """antiSMASH is the heaviest, most contig-count-sensitive step in the
    pipeline; run it strictly last, after every other step has fully
    finished for every sample -- not just after gene calling (unlike
    GECCO/eggNOG/coverage/binning/etc., which only wait on annotation.done
    and may still be running when antiSMASH would otherwise start)."""
    input:
        "status/preprocessing.done",
        "status/taxonomy.done",
        "status/assembly.done",
        "status/annotation.done",
        "status/coverage.done",
        "status/functions.done",
        "status/amr.done",
        "status/amrscan.done",
        "status/binning.done",
        "status/seed.done",
    output:
        touch("status/pre_antismash_barrier.done")


rule filter_assembly_min_len:
    """Drop contigs shorter than ANTISMASH_ASSEMBLY_MIN_CONTIG_LEN before
    antiSMASH ever sees them -- a real multi-gene BGC can't fit intact on a
    much shorter contig, but antiSMASH's hmmsearch cost scales with the
    genes annotated on them regardless.

    barrier: don't start for ANY sample until every other step (not just
    gene calling) has finished for ALL samples -- see pre_antismash_barrier."""
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly", "{sid}", "{sid}.fasta"),
        barrier="status/pre_antismash_barrier.done"
    output:
        fasta=os.path.join(BGC_OUTDIR, "antismash_filtered_fasta", "{sid}.fna")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    params:
        min_len=ANTISMASH_ASSEMBLY_MIN_CONTIG_LEN,
        bin=SEQKIT_BIN
    shell:
        """
        mkdir -p "$(dirname {output.fasta})"
        {params.bin} seq -m {params.min_len} {input.fasta} > {output.fasta}
        """

rule prep_antismash_gff:
    """Re-prefix Prodigal's GFF seqids to match the assembly fasta's headers
    ("{sid}:contig_..." vs. Prodigal's bare "contig_..."), and drop any GFF
    entries for contigs filter_assembly_min_len removed -- antiSMASH's
    cluster-detection hmmsearch runs against GFF-annotated genes, not a raw
    scan of contig nucleotides."""
    input:
        gff=os.path.join(RESULTS_DIR, "prodigal", "{sid}", "{sid}.gff"),
        fasta_filtered=rules.filter_assembly_min_len.output.fasta,
    output:
        gff3=os.path.join(BGC_OUTDIR, "antismash_gff", "{sid}.gff3")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    run:
        prefix = f"{wildcards.sid}:"
        keep_ids = set()
        with open(input.fasta_filtered) as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].split()[0]
                    bare_id = header[len(prefix):] if header.startswith(prefix) else header
                    keep_ids.add(bare_id)

        os.makedirs(os.path.dirname(output.gff3), exist_ok=True)
        with open(input.gff) as fh, open(output.gff3, "w") as out:
            for line in fh:
                if line.startswith("#"):
                    out.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0] in keep_ids:
                    fields[0] = prefix + fields[0]
                    out.write("\t".join(fields) + "\n")

rule prep_antismash_gff_all:
    """Barrier: every prep_antismash_gff job finishes before any
    run_antismash job starts."""
    input:
        expand(rules.prep_antismash_gff.output.gff3, sid=SAMPLES.index)
    output:
        touch(os.path.join(BGC_OUTDIR, "antismash_gff", ".all_prepped"))

rule run_antismash:
    """Predict BGCs with antiSMASH on the per-sample assembly.

    Uses --genefinding-gff3 to reuse this project's existing Prodigal calls
    (no re-running gene-finding), --minimal --cb-knownclusters (only
    knownclusterblast hits are consumed downstream), and runs against
    filter_assembly_min_len's output (contigs >= the configured minimum
    only). --reuse-results is applied automatically when a prior run's JSON
    is already present, so an interrupted run can resume."""
    input:
        fasta=rules.filter_assembly_min_len.output.fasta,
        gff3=rules.prep_antismash_gff.output.gff3,
        gff_done=rules.prep_antismash_gff_all.output[0],
    output:
        directory(os.path.join(BGC_OUTDIR, "antismash", "{sid}"))
    log:
        os.path.join(RESULTS_DIR, "logs", "antismash_assembly", "{sid}_antismash.log")
    params:
        db=config["antismash_assembly"]["db"],
        bin=ANTISMASH_BIN
    threads:
        THREADS_ANTISMASH_ASSEMBLY
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        """
        mkdir -p "$(dirname {log})"
        REUSE_ARG=""
        if [ -f "{output}/{wildcards.sid}.json" ]; then
            REUSE_ARG="--reuse-results {output}/{wildcards.sid}.json"
        fi
        (date && {params.bin} -c {threads} --genefinding-gff3 {input.gff3} --databases {params.db} --minimal --cb-knownclusters $REUSE_ARG --output-dir {output} {input.fasta} && date) &> {log}
        """


############################################
# 2b. BGC summary table (antiSMASH regions + coverage + knownclusterblast)
############################################
rule summarise_sample_bgc:
    """Extract one row per BGC region: type, coordinates, contig coverage,
    and the top knownclusterblast hit (if any), from antiSMASH's own JSON."""
    input:
        antismash_dir=os.path.join(BGC_OUTDIR, "antismash", "{sid}"),
        depth=os.path.join(RESULTS_DIR, "coverage", "{sid}", "{sid}_depth.txt"),
    output:
        tsv=temp(os.path.join(BGC_OUTDIR, "antismash_summaries", "{sid}_bgc_summary.tsv"))
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    run:
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        contig_depth = {}
        with open(input.depth) as fh:
            next(fh)
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                contig = fields[0].split(":", 1)[-1]
                contig_depth[contig] = fields[2]

        json_path = os.path.join(input.antismash_dir, f"{wildcards.sid}.json")
        with open(json_path) as fh:
            data = json.load(fh)

        rows = []
        for record in data.get("records", []):
            areas = record.get("areas", [])
            if not areas:
                continue
            # antiSMASH strips the ":" from "{sid}:contig_..." headers and
            # concatenates the remainder onto the sample id with no
            # separator -- strip the known sample prefix back off to match
            # contig_depth's bare keys.
            contig = record["id"]
            if contig.startswith(wildcards.sid):
                contig = contig[len(wildcards.sid):]
            coverage = contig_depth.get(contig, "NA")

            kc_by_region = {
                res["region_number"]: res
                for res in record.get("modules", {})
                                  .get("antismash.modules.clusterblast", {})
                                  .get("knowncluster", {})
                                  .get("results", [])
            }

            for i, area in enumerate(areas):
                region_number = i + 1
                bgc_id = f"Region_{region_number}"
                bgc_type = ";".join(area.get("products", []))

                best_hit, similarity = "None", 0
                kc_result = kc_by_region.get(region_number)
                if kc_result and kc_result.get("total_hits", 0) > 0 and kc_result.get("ranking"):
                    hit_info, score_info = kc_result["ranking"][0]
                    best_hit = f"{hit_info.get('accession', 'NA')}: {hit_info.get('description', 'NA')}"
                    similarity = score_info.get("similarity", 0)

                rows.append([
                    wildcards.sid, bgc_id, contig, bgc_type,
                    str(area.get("start", "NA")), str(area.get("end", "NA")),
                    str(coverage), best_hit, str(similarity),
                ])

        with open(output.tsv, "w") as out:
            for row in rows:
                out.write("\t".join(row) + "\n")

rule master_bgc_summary:
    """Merge all per-sample BGC summaries and add a header."""
    input:
        expand(os.path.join(BGC_OUTDIR, "antismash_summaries", "{sid}_bgc_summary.tsv"), sid=SAMPLES.index)
    output:
        os.path.join(BGC_OUTDIR, "antismash_summaries", "antismash_master_results.tsv")
    shell:
        """
        echo -e "Sample\\tBGC_ID\\tContig\\tType\\tStart\\tEnd\\tMean_Coverage\\tBest_Hit\\tSimilarity" > {output}
        cat {input} >> {output}
        """


############################################
# 3. BiG-SCAPE 2 -- BGC clustering across samples
############################################
rule gather_bigscape_input:
    """Symlink all antiSMASH region GBKs into a flat input directory.
    Rebuilds the directory from scratch each time (rm -rf first) so stale
    entries from an earlier partial run don't accumulate."""
    input:
        expand(os.path.join(BGC_OUTDIR, "antismash", "{sid}"), sid=SAMPLES.index)
    output:
        outdir=directory(os.path.join(BGC_OUTDIR, "bigscape_input"))
    log:
        os.path.join(RESULTS_DIR, "logs", "gather_bigscape_input.log")
    shell:
        """
        mkdir -p "$(dirname {log})"
        (rm -rf {output.outdir}
        mkdir -p {output.outdir}
        find {BGC_OUTDIR}/antismash/ -name "*.region*.gbk" \
            -exec ln -sf {{}} {output.outdir}/ \\;) &> {log}
        """

rule run_bigscape:
    """Cluster BGCs into Gene Cluster Families with BiG-SCAPE 2."""
    input:
        gbk_dir=os.path.join(BGC_OUTDIR, "bigscape_input")
    output:
        results=directory(os.path.join(BGC_OUTDIR, "bigscape_results"))
    log:
        os.path.join(RESULTS_DIR, "logs", "run_bigscape.log")
    threads:
        THREADS_BIGSCAPE
    params:
        bin=BIGSCAPE_BIN,
        pfam=PFAM_DB
    shell:
        """
        mkdir -p "$(dirname {log})"
        (rm -rf {output.results} && {params.bin} cluster --input-dir {input.gbk_dir} --output-dir {output.results} --mibig-version 3.1 --cores {threads} --gcf-cutoffs 0.3,0.4,0.5 --pfam-path {params.pfam} --mix --force-gbk --label REHAB_Project) &> {log}
        """


############################################
# 4. Average Maximal Growth Rate (gRodon2, codon usage bias)
############################################
# gRodon2 predicts community-average maximal growth rate from codon usage
# bias between highly-expressed genes (ribosomal proteins) and bulk gene
# content. Reuses Prodigal's existing .faa headers (exact nucleotide
# coordinates) to recover CDS sequences from the existing assembly FASTA --
# no new gene predictions. "Highly expressed" genes are identified by
# string-matching "ribosomal protein" in eggNOG's Description column.
# Runs via the published gRodon2 Docker image under Singularity (no conda
# recipe covers Bioconductor + CRAN together) -- image and all temp files
# live under this repo's containers/ and tmp/, never $HOME or /tmp.
rule pull_grodon_sif:
    output:
        sif=GRODON_SIF
    log:
        os.path.join(RESULTS_DIR, "logs", "pull_grodon_sif.log")
    shell:
        """
        mkdir -p "$(dirname {output.sif})" "$(dirname {log})"
        SINGULARITY_TMPDIR="{config[work_dir]}/tmp" SINGULARITY_CACHEDIR="{config[work_dir]}/tmp/singularity_cache" \
            singularity pull --force {output.sif} docker://{GRODON_IMAGE} &> {log}
        """

rule write_grodon_script:
    """Write the gRodon2 R driver script as a real Snakemake output, so it
    self-heals if ever missing rather than depending on a hand-placed file."""
    output:
        script=GRODON_SCRIPT
    run:
        os.makedirs(os.path.dirname(output.script), exist_ok=True)
        with open(output.script, "w") as f:
            f.write(r'''#!/usr/bin/env Rscript
# Run gRodon2 (metagenome mode) on one sample's extracted CDS + highly-
# expressed (ribosomal protein) gene set, and write the result as a
# single-row TSV.
# Usage: run_gRodon.R <cds.fna> <highly_expressed.txt> <sample_id> <output.tsv> <depth.tsv>
suppressMessages(library(Biostrings))
suppressMessages(library(gRodon))

args <- commandArgs(trailingOnly = TRUE)
cds_path   <- args[1]
he_path    <- args[2]
sample_id  <- args[3]
out_path   <- args[4]
depth_path <- args[5]

genes <- readDNAStringSet(cds_path)
he_ids <- readLines(he_path)
he_ids <- he_ids[he_ids != ""]
highly_expressed <- names(genes) %in% he_ids

depth_table <- read.delim(depth_path, header = FALSE, col.names = c("gene_id", "depth"))
depths <- depth_table$depth
names(depths) <- depth_table$gene_id
depth_of_coverage <- depths[names(genes)]

result <- predictGrowth(genes, highly_expressed, mode = "metagenome_v2",
                         depth_of_coverage = depth_of_coverage)

df <- data.frame(
  sample            = sample_id,
  d                 = result$d,
  LowerCI           = result$LowerCI,
  UpperCI           = result$UpperCI,
  CUBHE             = result$CUBHE,
  CUB               = result$CUB,
  dCUB              = result$dCUB,
  ConsistencyHE     = result$ConsistencyHE,
  CPB               = result$CPB,
  GC                = result$GC,
  GCdiv             = result$GCdiv,
  nHE               = result$nHE,
  FilteredSequences = result$FilteredSequences
)
write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
''')

rule extract_grodon_input:
    """Build gRodon2's inputs for {sid} without any new gene-calling or
    read-mapping: nucleotide CDS for every Prodigal-called gene (recovered
    from the assembly FASTA via Prodigal's .faa header coordinates), a
    highly-expressed (ribosomal protein) flag list from eggNOG, and a
    per-gene mean-depth table reusing the existing coverage output."""
    input:
        faa=os.path.join(RESULTS_DIR, "prodigal", "{sid}", "{sid}.faa"),
        assembly=os.path.join(RESULTS_DIR, "assembly", "{sid}", "{sid}.fasta"),
        annot=os.path.join(RESULTS_DIR, "eggnog", "{sid}", "{sid}.emapper.annotations"),
        coverage=os.path.join(RESULTS_DIR, "coverage", "{sid}", "{sid}_gene_coverage.txt"),
    output:
        cds=os.path.join(GRODON_DIR, "{sid}", "{sid}_cds.fna"),
        he=os.path.join(GRODON_DIR, "{sid}", "{sid}_highly_expressed.txt"),
        depth=os.path.join(GRODON_DIR, "{sid}", "{sid}_depth.tsv"),
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    run:
        os.makedirs(os.path.dirname(output.cds), exist_ok=True)

        he_ids = set()
        with open(input.annot) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                gene_id, description = fields[0], fields[7]
                if "ribosomal protein" in description.lower():
                    he_ids.add(gene_id)
        with open(output.he, "w") as out:
            out.write("\n".join(sorted(he_ids)) + "\n")

        with open(input.coverage) as fh, open(output.depth, "w") as out:
            for line in fh:
                fields = line.split()
                if len(fields) < 3:
                    continue
                gene_id, depth = fields[1], fields[2]
                out.write(f"{gene_id}\t{depth}\n")

        seqs = {}
        name, parts = None, []
        with open(input.assembly) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(parts)
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            if name is not None:
                seqs[name] = "".join(parts)

        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

        with open(input.faa) as fh, open(output.cds, "w") as out:
            for line in fh:
                if not line.startswith(">"):
                    continue
                fields = line[1:].strip().split(" # ")
                gene_id = fields[0]
                start, end, strand = int(fields[1]), int(fields[2]), fields[3]
                contig = re.sub(r"_\d+$", "", gene_id)
                chrom = f"{wildcards.sid}:{contig}"
                seq = seqs.get(chrom)
                if seq is None:
                    continue
                cds = seq[start - 1:end]
                if strand == "-1":
                    cds = cds.translate(comp)[::-1]
                out.write(f">{gene_id}\n{cds}\n")

rule run_gRodon:
    """Run gRodon2 (metagenome_v2 mode) on {sid}'s recovered CDS set,
    depth-weighted using the per-gene coverage table."""
    input:
        cds=rules.extract_grodon_input.output.cds,
        he=rules.extract_grodon_input.output.he,
        depth=rules.extract_grodon_input.output.depth,
        sif=GRODON_SIF,
        script=rules.write_grodon_script.output.script,
    output:
        tsv=os.path.join(GRODON_DIR, "{sid}", "{sid}_growth_rate.tsv"),
    log:
        os.path.join(RESULTS_DIR, "logs", "{sid}_run_gRodon.log"),
    params:
        data_dir=lambda wc, input: os.path.dirname(input.cds),
        script_dir=os.path.dirname(GRODON_SCRIPT),
        tmpdir=os.path.join(config["work_dir"], "tmp"),
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        """
        mkdir -p "$(dirname {log})"
        (SINGULARITY_TMPDIR={params.tmpdir} singularity exec \
            --bind {params.data_dir}:/data \
            --bind {params.script_dir}:/scripts \
            {input.sif} \
            Rscript /scripts/run_gRodon.R \
                /data/{wildcards.sid}_cds.fna \
                /data/{wildcards.sid}_highly_expressed.txt \
                {wildcards.sid} \
                /data/{wildcards.sid}_growth_rate.tsv \
                /data/{wildcards.sid}_depth.tsv) &> {log}
        """

rule master_growth_rate_summary:
    """Concatenate every sample's gRodon2 growth-rate estimate into one table."""
    input:
        expand(rules.run_gRodon.output.tsv, sid=SAMPLES.index)
    output:
        tsv=os.path.join(GRODON_DIR, "growth_rate_master_results.tsv")
    run:
        header = None
        rows = []
        for f in input:
            with open(f) as fh:
                lines = [l.rstrip("\n") for l in fh if l.strip()]
            if header is None:
                header = lines[0]
            rows.extend(lines[1:])
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w") as out:
            out.write(header + "\n")
            out.write("\n".join(rows) + "\n")
