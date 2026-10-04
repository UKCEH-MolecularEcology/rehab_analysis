# REHAB Analysis

Metagenomic taxonomy, assembly, and AMR analysis pipeline for the **REHAB** study
([PRJEB34634](https://www.ebi.ac.uk/ena/browser/view/PRJEB34634)) — a UKCEH project
investigating the transmission of antibiotic-resistant bacteria and antibiotic
resistance genes across farm animals, human/animal sewage, sewage treatment works,
and rivers.

Snakemake workflow structure is adapted from
[UKCEH-MolecularEcology/metag_analyses](https://github.com/UKCEH-MolecularEcology/metag_analyses).

## Data

PRJEB34634 comprises **2,076 sequencing runs** (paired-end Illumina, ~2.1 TB) which
map, via `metadata/rehab_metagenomes_list.csv`, to **147 real biological samples**
(site/matrix/distance/season combinations, e.g. `DID-DW-100-AUT` = Didcot,
Downstream, 100 m, Autumn), each split across 7–18 sequencing runs/lanes that need
concatenating.

> Note: ENA's own per-run metadata (`sample_accession`, `run_alias`, `library_name`)
> is **not** informative for this project — every run shares the same BioSample
> (`SAMEA5989477`). The real sample identity only exists in
> `metadata/rehab_metagenomes_list.csv`.

Files in `metadata/`:
- `rehab_metagenomes_list.csv` — run → samplename mapping + fastq FTP URLs (used by `scripts/concat_fastq.py`)
- `EBI_data.csv` / `EBI_REHAB.csv` (identical) — pre-existing ENA assemblies for a 106-sample subset (project PRJEB45055, `SEQUENCE_ASSEMBLY`)
- `All_chemistry_data.txt` — water chemistry covariates per site/season

`data/PRJEB34634/metadata/PRJEB34634_read_run_metadata.tsv` holds the full ENA
`read_run` file report (run/experiment/sample accessions, instrument, byte sizes,
md5 checksums, etc.) pulled directly from the ENA Portal API.

## Tools

All new tool installs for this repo (not already provided by the shared conda
environments in `/hdd0/susbus/tools/conda_envs`) are Singularity containers stored
under `containers/` (gitignored — rebuild with the pull commands below). All
tool/Singularity temp files go to `tmp/` (gitignored), never `/tmp` or `$HOME`.

```bash
export SINGULARITY_TMPDIR=containers/../tmp
export SINGULARITY_CACHEDIR=tmp/singularity_cache
singularity pull containers/fastq-dl.sif docker://quay.io/biocontainers/fastq-dl:4.0.1--pyhdfd78af_0
```

> **Known issue**: the shared `eggnog-mapper=2.1.9` conda env at
> `/prj/DECODE/ea_biofilm_results/conda_envs/da281a3695a41196430106915101b45c_`
> was solved against Python 3.14, which removed the stdlib `distutils` module
> that eggnog-mapper 2.1.9 still imports directly — every `emapper.py`
> invocation (including the `download_eggnogDB` setup rule) failed with
> `ModuleNotFoundError: No module named 'distutils'`. Fixed by installing
> `setuptools` into that env (`pip install setuptools`), which ships a
> `distutils` compatibility shim — confirmed `emapper.py` now imports and
> runs correctly. This is an env-level fix (not tracked in this repo); if the
> shared env cache is ever rebuilt from its `.yaml` alone, this will need
> reapplying.

## 1. Download raw FASTQ

```bash
cd /prj/DECODE/rehab
cut -d',' -f4 metadata/rehab_metagenomes_list.csv | tail -n +2 > tmp/run_accessions.txt
./scripts/download_fastq.sh tmp/run_accessions.txt data/PRJEB34634/fastq 16
```

Downloads run in parallel (16 concurrent `fastq-dl` processes by default) via the
Singularity container; per-run logs land in `tmp/logs/download/`. Safe to re-run —
`fastq-dl` skips files it has already fetched successfully.

## 2. Concatenate runs into per-sample FASTQ

```bash
python scripts/concat_fastq.py \
    -m metadata/rehab_metagenomes_list.csv \
    -f data/PRJEB34634/fastq \
    -o data/concatenated_fastq
```

## 3. Generate the sample manifest

```bash
./scripts/generate_sample_table.sh data/concatenated_fastq config/samples.tsv
```

## 4. Run the Snakemake pipeline

No SLURM on this host — run locally:

```bash
./scripts/run_pipeline.sh -n   # dry-run
./scripts/run_pipeline.sh      # full run
```

`config/config.yaml` controls which pipeline `steps` run. Currently enabled:
`preprocessing`, `taxonomy`, `assembly`, `annotation` (Prodigal), `coverage`,
`functions` (eggNOG + KEGG antibiotic-biosynthesis marker search + MagicLamp),
`amr` (RGI), `amrscan`, `bgc_amr`, `binning` (mmseqs2 dedup → 7-binner
ensemble → DAS_Tool → Rosella refine → Galah dereplication → GTDB-Tk/CheckM2
quality — see **MAG generation** below), `mgthermometer`, `syntracker`,
`seed`. See `workflow/rules/*.smk` for the full set mirrored from
`metag_analyses`; each addition is verified with `snakemake -n` before being
enabled.

Two pre-existing bugs in the mirrored `metag_analyses` rules were fixed here
(both blocked `binning` outright, not something introduced by this repo):
- `rules/cluster.smk`'s `cat_ass_mmseqs2` referenced a non-existent path
  (`mmseqs/{sid}/{sid}_modified.fasta`); fixed to use `ass_mmseqs2`'s real
  output (`mmseqs/{sid}/{sid}_rep_seq.fasta` — the correct path was already
  present as a dead commented-out alternative on the same line).
- `rules/semibin.smk`'s `unzip_semibin` shell command had unescaped
  `${1%.gz}` / `{}` (find's placeholder), which Snakemake's shell
  `.format()` misinterprets; escaped as `${{1%.gz}}` / `{{}}`.

## MAG generation (`binning`)

Adopts [michoug/MAGsGeneration](https://github.com/michoug/MAGsGeneration)'s
multi-binner ensemble approach, adapted to this pipeline's existing
**pooled/merged co-assembly** architecture: all 147 samples' deduplicated
contigs are already concatenated into one `results/assembly/cat_assembly_filter.fasta`,
and all 147 samples' reads are already mapped against that single reference
(`results/bam/{sid}/cat_assembly_{sid}.bam`). Every new binner below
consumes that same pooled assembly + existing BAMs directly — there is no
per-sample or scoped-grouping variant, since the coverage-signal problem
the reference repo's own per-sample design solves (needing cross-sample
coverage for good binning) is already solved by the pooled design.

**Binners** (`workflow/rules/mags_generation.smk`, `workflow/rules/binning.smk`):
MetaBAT2, CONCOCT, Rosella, TaxVAMB (taxonomy-informed VAMB, using an
mmseqs2 GTDB classification of the pooled contigs), COMEBin, SemiBin2,
MetaCAT — 7 binners total. MetaBinner and a plain (non-taxonomy) SemiBin
run are disabled for now (see the commented-out lines in `binning.smk`).

**Dereplication** (`workflow/rules/dereplicate.smk`, `bin_taxqual.smk`):
DAS_Tool ensembles all 7 binners' outputs (`--score_threshold -42`) →
Rosella refine (density/UMAP re-clustering informed by CheckM2 quality) →
Galah dereplication (replacing dRep) → the pipeline's pre-existing
GTDB-Tk/CheckM2 setup runs unchanged on the final set
(`results/bins/finalbins`).

**Package management**: several of the new binners (VAMB/TaxVAMB, MetaCAT)
are pip/wheel installs with no plain bioconda package, so this module uses
[pixi](https://pixi.sh) (`pixi.toml`/`pixi.lock` at the repo root) instead
of conda — scoped to just these new tools, not a replacement for the
conda/Singularity environments used everywhere else in the pipeline. CheckM2
and GTDB-Tk are deliberately *not* in `pixi.toml` — they reuse the
already-configured conda env + databases (`config.yaml`'s `checkm2`/`gtdbtk`
sections). `scripts/pixi_env.sh` must be sourced before any `pixi` command
(including from rule shell blocks) — it keeps pixi's binary, cache, and
resolved environments entirely project-local
(`tools/pixi_home/`, `.pixi/`, both gitignored — ~27 GB resolved), never
`$HOME` or `/tmp`.

**Databases**: new DBs needed by this module (mmseqs2 GTDB taxonomy,
MetaCAT) go under `/hdd0/susbus/databases/` alongside the project's
existing shared DBs (CheckM2, GTDB-Tk, BUSCO, SingleM, etc.), not inside
this repo.

### Optional sub-modules

Adapted from the same reference repo's optional steps:

- **MGthermometer** (`workflow/rules/mgthermometer.smk`, standard step) —
  per-MAG optimal growth temperature (OGT) proxy via the IVYWREL
  amino-acid frequency method, using pyrodigal protein prediction +
  `workflow/scripts/getFrequency.pl` (carried over from `metag_analyses`).
- **SynTracker** (`workflow/rules/syntracker.smk`, standard step) — per-MAG
  synteny tracking across all 147 samples. The one piece of this module
  that does *not* use the pooled assembly: it needs each sample's own
  per-sample assembly (`results/assembly/{sid}/{sid}.fasta`) as the BLAST
  target pool, since tracking synteny *across* samples requires keeping
  them distinct.
- **Eukaryotic MAG filtering** (`workflow/rules/eukaryotes_filter.smk`,
  **opt-in** — add `"eukaryotic_mags"` to `config.yaml`'s `steps` list to
  enable) — recovers eukaryotic bins that the prokaryote-focused
  DAS_Tool/Rosella-refine/Galah chain above would otherwise discard, via
  two independent tracks: DeepMicroClass classification of every raw
  candidate bin (≥80% eukaryotic sequence kept) and REMAG (a binner
  designed specifically for eukaryotic genome recovery) run directly on
  the pooled assembly. Both tracks' candidates get BUSCO (`--auto-lineage-euk`,
  reusing the existing `eukaryota_odb10` cache at
  `/hdd0/susbus/databases/busco_downloads`) and a final Galah
  dereplication pass, matching the main MAG chain.

### A note on Snakemake barriers vs. priority hints

When ordering steps against an already-partially-complete pipeline run, a
*hard* input-file barrier (e.g. requiring a new marker file as input) can
force Snakemake to invalidate and re-run expensive already-completed jobs —
the new marker file is always newer by mtime, and this holds even under
`--rerun-triggers mtime`. Where this pipeline needs step B to merely run
*after* step A completes rather than strictly *depend on* A's output (e.g.
SingleM before eggNOG), a soft `priority:` hint is used instead of a real
barrier, specifically to avoid invalidating completed work. See the
relevant rule files' comments for details.

## Read-based AMR detection (`amrscan`)

Adds [UKCEH-MolecularEcology/snake_amrscan](https://github.com/UKCEH-MolecularEcology/snake_amrscan)
as a git submodule (`submodules/snake_amrscan`), wired in via
`workflow/rules/amrscan.smk`. Runs H.S. Gweon's
[resscan](https://github.com/hsgweon/resscan) tool directly on trimmed reads
(gene/variant-level AMR calls, independent of assembly quality) — a third,
complementary approach alongside `rules/amr.smk` (RGI on assembled contigs)
and `bgc_amr.smk`'s curated-KEGG-marker summary.

> **Known issue**: `snake_amrscan`'s own nested submodule (`hsgweon/resscan`)
> is pinned at a commit that was never pushed to its public GitHub remote, so
> `git submodule update --init --recursive` fails for that inner submodule on
> a fresh clone. Worked around here by copying the working tree directly from
> the known-good checkout at `/prj/DECODE/salisbury_plain_results/snake_amrscan/submodules/resscan`
> (untracked, not a real git submodule) — re-run that copy if
> `submodules/snake_amrscan/submodules/resscan/resscan/resscan.py` is ever
> missing after a fresh clone.

Reuses the pre-existing `amrscan` conda env
(`/prj/DECODE/conda_envs/amrscan`, has bwa/diamond/samtools) stitched with
the base miniforge3 python (has pandas/numpy, which that env lacks) — see
`config/config.yaml`'s `amrscan` section.

## BGC / AMR-marker / growth-rate step (`bgc_amr`)

Adapted from the SOCD project's `singlem_bgc_Snakefile_no_metadata`
(`/prj/DECODE/socd/results/`). Runs per-sample (not per-MAG, unlike
`rules/antismash.smk`), scoped to samples with eggNOG + coverage already
computed:

- **ABX/KEGG markers** — cross-references eggNOG KO annotations against a
  curated antibiotic biosynthesis/resistance reference
  (`resources/kegg_abx_bgc.txt`), joined with gene coverage. Complements
  (does not replace) `rules/eggnog.smk`'s pathway-keyword-based
  `identify_abx_biosynthesis`.
- **BGC prediction** — GECCO + antiSMASH (contig-length-filtered,
  `--genefinding-gff3` reusing existing Prodigal calls, `--minimal
  --cb-knownclusters` only) on each sample's assembly.
- **BGC clustering** — BiG-SCAPE 2 across all samples' predicted regions.
- **Growth rate** — gRodon2 (codon usage bias → average maximal growth
  rate), via a Singularity image (no conda recipe covers Bioconductor +
  CRAN together).

GECCO, antiSMASH, and BiG-SCAPE reuse the pre-existing named conda
environments on this host (`/home/susbus/miniforge3/envs/{gecco,bigscape,
antismash}`) rather than being rebuilt — see `config/config.yaml`'s
`gecco`/`antismash_assembly`/`bigscape` sections to repoint them elsewhere.

## Repo layout

```
config/      config.yaml, samples.tsv (generated), schema-validated
data/        ENA metadata + raw/concatenated FASTQ (gitignored, huge)
metadata/    REHAB sample tracking + chemistry data (tracked in git)
resources/   pipeline reference data (e.g. kegg_abx_bgc.txt), tracked in git
scripts/     download / concat / manifest / pipeline-runner scripts, plus
             pixi_env.sh (sourced by MAG-generation rules before any `pixi` call)
workflow/    Snakefile, rules/, envs/, scripts/ (mirrored from metag_analyses,
             plus rules/bgc_amr.smk adapted from SOCD and
             rules/{mags_generation,dereplicate,bin_taxqual,eukaryotes_filter,
             mgthermometer,syntracker}.smk adapted from michoug/MAGsGeneration)
schemas/     config & sample-sheet JSON schemas
containers/  Singularity images for newly-installed tools (gitignored)
pixi.toml    pixi manifest for the MAG-generation module's new binners
             (VAMB, MetaCAT, etc.) -- scoped to that module only, everything
             else still uses conda/Singularity
pixi.lock    pinned pixi resolution, tracked for reproducibility
tools/       pixi's own binary/cache/resolved envs (gitignored, ~27 GB)
tmp/         scratch space for all tool/Singularity/conda/pixi temp files (gitignored)
```
