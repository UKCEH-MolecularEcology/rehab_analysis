#!/bin/bash
# Run antiSMASH on one MAG inside the antiSMASH container, keep only the compact outputs,
# summarise regions. Always exits 0; the outcome is written to <outdir>/<mag>/STATUS (ok | failed).
# Usage: run_antismash_one.sh <mag> <fasta> <outdir> <sif> <cpus> <regions_script> <bind_root> [extra antiSMASH args...]
#
# <bind_root> is explicitly --bind-mounted into the container at the same
# path on both sides: unlike some HPC singularity/apptainer setups, this
# host doesn't auto-bind arbitrary paths (only /tmp), so <fasta>, <outdir>,
# and <regions_script> all need to live under it for the container to see
# them -- in practice this is the pipeline's whole work_dir.
mag=$1; fa=$2; outdir=$3; sif=$4; cpus=$5; regions=$6; bind_root=$7; shift 7; extra="$*"
d=$outdir/$mag
if [ -s $d/STATUS ] && grep -q '^ok' $d/STATUS; then echo "skip $mag"; exit 0; fi
rm -rf $d && mkdir -p $(dirname $d)
if singularity exec --bind $bind_root:$bind_root $sif antismash --cpus $cpus --genefinding-tool prodigal-m --allow-long-headers \
     --output-dir $d --output-basename $mag $extra $fa > $outdir/$mag.antismash.log 2>&1; then
  # keep results JSON + per-region GenBank; drop the HTML report (large, per MAG)
  find $d -mindepth 1 -maxdepth 1 ! -name "$mag.json" ! -name '*.region*.gbk' -exec rm -rf {} +
  singularity exec --bind $bind_root:$bind_root $sif python3 $regions --mag $mag --dir $d --out $d/$mag.regions.tsv \
    && echo ok > $d/STATUS && rm -f $outdir/$mag.antismash.log && echo "done $mag" && exit 0
fi
mkdir -p $d && { echo failed; tail -5 $outdir/$mag.antismash.log; } > $d/STATUS
echo "FAILED $mag"; exit 0
