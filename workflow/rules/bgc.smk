"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2026-10-09]
Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 64 -rp
Purpose: Predict biosynthetic gene clusters (BGCs) in the final dereplicated
         MAGs with antiSMASH and GECCO, batched, collated into one table
         per tool. Ported from UKCEH-MolecularEcology/skyline's annotations
         branch (workflow/rules/bgc.smk), which runs this on a large
         pre-existing dereplicated-MAG collection as a separate downstream
         pipeline. Adapted here for the pooled-assembly architecture, where
         the MAG set doesn't exist until rules.galah runs in this SAME
         Snakemake invocation: skyline's parse-time `glob_wildcards` over
         an already-populated MAG directory becomes a checkpoint here.

Replaces the old per-sample-assembly antiSMASH pass in bgc_amr.smk (too
slow at 147 samples) -- see that file's docstring. antiSMASH uses the
official container image (databases bundled, no separate DB setup step),
GECCO reuses this project's existing named conda env. Both run in batches
of MAGs with internal parallelism (xargs -P), not one Snakemake job per
MAG, to keep job count reasonable; a batch only fails outright if more
than `bgc.max_failed_frac` of its MAGs fail, so a handful of genuinely
bad MAGs don't abort the whole run (tracked per-MAG via a STATUS file,
also letting reruns skip MAGs that already finished).
"""

import glob

BGC = config["bgc"]
BGC_DIR = os.path.join(RESULTS_DIR, "bgc")
BGC_SIF = os.path.join(config["work_dir"], "containers", "antismash_standalone.sif")


############################################
rule bgc:
    input:
        os.path.join(BGC_DIR, "antismash_regions_all.tsv"),
        os.path.join(BGC_DIR, "gecco_clusters_all.tsv"),
        os.path.join(BGC_DIR, "bgc_overlap.tsv")
    output:
        touch("status/bgc.done")


############################################
localrules: pull_antismash_image, collate_antismash, collate_gecco, bgc_overlap, bgc_mags


############################################
# Dynamic MAG discovery + batching -- the final dereplicated MAG set only
# exists once rules.galah has run in this same invocation, so (unlike
# skyline's parse-time glob over a pre-existing directory) this needs a
# checkpoint. `barrier` preserves "antiSMASH runs strictly last" (see
# bgc_amr.smk's pre_antismash_barrier) now that GECCO/antiSMASH moved here.
checkpoint bgc_mags:
    input:
        mags=rules.galah.output.final,
        barrier=rules.pre_antismash_barrier.output
    output:
        directory(os.path.join(BGC_DIR, "mags"))
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/bgc_mags.log")
    message:
        "Snapshotting the final dereplicated MAG set for BGC prediction"
    shell:
        """
        (date && mkdir -p {output} && \
        for f in {input.mags}/*.fa; do ln -sf "$(readlink -f "$f")" {output}/; done && \
        date) &> {log}
        """


def all_mags(wildcards):
    ck = checkpoints.bgc_mags.get(**wildcards).output[0]
    return sorted(glob_wildcards(os.path.join(ck, "{mag}.fa")).mag)


def all_batches(wildcards):
    mags = all_mags(wildcards)
    size = BGC["batch_size"]
    return {
        "batch{:04d}".format(i // size): mags[i:i + size]
        for i in range(0, len(mags), size)
    }


def batch_mags(wildcards):
    return all_batches(wildcards)[wildcards.batch]


def batch_mag_paths(wildcards):
    ck = checkpoints.bgc_mags.get(**wildcards).output[0]
    return [os.path.join(ck, m + ".fa") for m in batch_mags(wildcards)]


############################################ antiSMASH
# official image incl. databases (needs internet to pull); build into
# <sif>.tmp, rename when complete -- an interrupted pull never leaves a
# half-written file at the real path.
rule pull_antismash_image:
    output:
        sif=BGC_SIF
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/pull_antismash_image.log")
    message:
        "Pulling antiSMASH image ({})".format(BGC["antismash"]["image"])
    shell:
        """
        (date && mkdir -p "$(dirname {output.sif})" && \
        SINGULARITY_TMPDIR="{config[work_dir]}/tmp" SINGULARITY_CACHEDIR="{config[work_dir]}/tmp/singularity_cache" \
            singularity pull {output.sif}.tmp docker://{BGC[antismash][image]} && \
        mv {output.sif}.tmp {output.sif} && \
        date) &> {log}
        """


# antiSMASH on a batch of MAGs, several in parallel; finished MAGs (own
# STATUS=ok already present) are skipped on reruns.
rule antismash_batch:
    input:
        fa=batch_mag_paths,
        sif=BGC_SIF
    output:
        done=os.path.join(BGC_DIR, "antismash/batches/{batch}.done")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/antismash_{batch}.log")
    params:
        outdir=os.path.join(BGC_DIR, "antismash/per_mag"),
        pairs=lambda wc: " ".join(
            "{} {}".format(m, p) for m, p in zip(batch_mags(wc), batch_mag_paths(wc))
        ),
        mags=lambda wc: " ".join(batch_mags(wc)),
        parallel=BGC["antismash"]["parallel"],
        cpus=BGC["antismash"]["cpus_per_mag"],
        extra=BGC["antismash"].get("extra", ""),
        runner=os.path.join(SRC_DIR, "run_antismash_one.sh"),
        regions=os.path.join(SRC_DIR, "antismash_regions.py"),
        max_fail=BGC.get("max_failed_frac", 0.05)
    threads:
        BGC["antismash"]["parallel"] * BGC["antismash"]["cpus_per_mag"]
    message:
        "antiSMASH on {wildcards.batch}"
    shell:
        """
        (date && mkdir -p {params.outdir} && \
        echo {params.pairs} | xargs -n 2 -P {params.parallel} sh -c \
            'bash {params.runner} "$0" "$1" {params.outdir} {input.sif} {params.cpus} {params.regions} {config[work_dir]} {params.extra}' && \
        n=0; f=0; for m in {params.mags}; do n=$((n+1)); grep -q '^ok' {params.outdir}/$m/STATUS || f=$((f+1)); done; \
        echo "antiSMASH failed for $f of $n MAGs"; \
        if [ $(awk -v f=$f -v n=$n 'BEGIN{{print (f > {params.max_fail} * n) ? 1 : 0}}') -eq 1 ]; then exit 1; fi && \
        touch {output.done} && date) &> {log}
        """


# all MAGs -> one table of regions
rule collate_antismash:
    input:
        lambda wc: expand(
            os.path.join(BGC_DIR, "antismash/batches/{batch}.done"),
            batch=sorted(all_batches(wc))
        )
    output:
        os.path.join(BGC_DIR, "antismash_regions_all.tsv")
    params:
        outdir=os.path.join(BGC_DIR, "antismash/per_mag")
    message:
        "Collating antiSMASH regions for all MAGs"
    run:
        mags = all_mags(wildcards)
        header, failed = None, []
        with open(output[0], "w") as out:
            for m in mags:
                f = os.path.join(params.outdir, m, m + ".regions.tsv")
                if not os.path.exists(f):
                    failed.append(m)
                    continue
                with open(f) as fh:
                    h = fh.readline()
                    if header is None:
                        header = h
                        out.write(h)
                    for line in fh:
                        out.write(line)
        with open(output[0].replace(".tsv", ".failed_mags.txt"), "w") as fh:
            fh.write("\n".join(failed) + ("\n" if failed else ""))
        print("antiSMASH: {} MAGs without results (see *.failed_mags.txt)".format(len(failed)))


############################################ GECCO
rule gecco_batch:
    input:
        fa=batch_mag_paths
    output:
        done=os.path.join(BGC_DIR, "gecco/batches/{batch}.done")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/gecco_{batch}.log")
    params:
        outdir=os.path.join(BGC_DIR, "gecco/per_mag"),
        pairs=lambda wc: " ".join(
            "{} {}".format(m, p) for m, p in zip(batch_mags(wc), batch_mag_paths(wc))
        ),
        parallel=BGC["gecco"]["parallel"],
        cpus=BGC["gecco"]["cpus_per_mag"],
        extra=BGC["gecco"].get("extra", ""),
        bin=BGC["gecco"]["bin"],
        mags=lambda wc: " ".join(batch_mags(wc)),
        max_fail=BGC.get("max_failed_frac", 0.05)
    threads:
        BGC["gecco"]["parallel"] * BGC["gecco"]["cpus_per_mag"]
    message:
        "GECCO on {wildcards.batch}"
    shell:
        # one output folder per MAG; a MAG counts as done once its
        # .clusters.tsv exists (GECCO writes it last)
        """
        (date && {params.bin} --version && mkdir -p {params.outdir} && \
        echo {params.pairs} | xargs -n 2 -P {params.parallel} sh -c '\
        d={params.outdir}/$0; if ls $d/*.clusters.tsv >/dev/null 2>&1; then echo "skip $0"; exit 0; fi; \
        rm -rf $d $d.FAILED && {params.bin} run --genome $1 --output-dir $d --jobs {params.cpus} {params.extra} > $d.log 2>&1 \
        && rm -f $d.log && echo "done $0" || {{ mv $d.log $d.FAILED; echo "FAILED $0"; }}' && \
        n=0; f=0; for m in {params.mags}; do n=$((n+1)); [ -e {params.outdir}/$m.FAILED ] && f=$((f+1)); done; \
        echo "GECCO failed for $f of $n MAGs"; \
        if [ $(awk -v f=$f -v n=$n 'BEGIN{{print (f > {params.max_fail} * n) ? 1 : 0}}') -eq 1 ]; then exit 1; fi && \
        touch {output.done} && date) &> {log}
        """


# all MAGs -> one table of clusters (Genome column added)
rule collate_gecco:
    input:
        lambda wc: expand(
            os.path.join(BGC_DIR, "gecco/batches/{batch}.done"),
            batch=sorted(all_batches(wc))
        )
    output:
        os.path.join(BGC_DIR, "gecco_clusters_all.tsv")
    params:
        outdir=os.path.join(BGC_DIR, "gecco/per_mag")
    message:
        "Collating GECCO clusters for all MAGs"
    run:
        mags = all_mags(wildcards)
        header = None
        with open(output[0], "w") as out:
            for m in mags:
                for f in sorted(glob.glob(os.path.join(params.outdir, m, "*.clusters.tsv"))):
                    with open(f) as fh:
                        h = fh.readline()
                        if header is None:
                            header = h
                            out.write("Genome\t" + h)
                        for line in fh:
                            out.write(m + "\t" + line)


############################################ antiSMASH vs GECCO
# regions/clusters matched by MAG + contig + coordinate overlap
rule bgc_overlap:
    input:
        antismash=rules.collate_antismash.output[0],
        gecco=rules.collate_gecco.output[0]
    output:
        pairs=os.path.join(BGC_DIR, "bgc_overlap.tsv"),
        summary=os.path.join(BGC_DIR, "bgc_overlap_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/bgc/bgc_overlap.log")
    params:
        script=os.path.join(SRC_DIR, "bgc_overlap.py"),
        min_bp=BGC.get("overlap_min_bp", 1)
    message:
        "Matching antiSMASH regions and GECCO clusters"
    shell:
        """
        (date && python3 {params.script} --antismash {input.antismash} --gecco {input.gecco} \
            --out {output.pairs} --summary {output.summary} --min_overlap_bp {params.min_bp} && \
        cat {output.summary} && date) &> {log}
        """
