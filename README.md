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

`config/config.yaml` controls which pipeline `steps` run (preprocessing, taxonomy,
assembly, AMR, binning, ...). See `workflow/rules/*.smk` for the full set mirrored
from `metag_analyses`.

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
scripts/     download / concat / manifest / pipeline-runner scripts
workflow/    Snakefile, rules/, envs/, scripts/ (mirrored from metag_analyses,
             plus rules/bgc_amr.smk adapted from the SOCD project)
schemas/     config & sample-sheet JSON schemas
containers/  Singularity images for newly-installed tools (gitignored)
tmp/         scratch space for all tool/Singularity/conda temp files (gitignored)
```
