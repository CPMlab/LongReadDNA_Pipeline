# Operations

Running, monitoring and troubleshooting a real cohort run.

## Before starting

```bash
# 1. Sample sheet points at readable BAMs
cut -f3 samples.tsv | tail -n+2 | xargs -I{} ls -l {} | head

# 2. Reference paths in config.yaml exist
grep -E "reference_fasta|vep_cache|_container|_jar|hmf_resources" config.yaml

# 3. Disk space — roughly 2-3 TB per patient
df -h .

# 4. Dry run
snakemake -n --cores 4
```

## Running

```bash
# SLURM (recommended)
sbatch slurm/run_snakemake.sh

# Single patient
bash slurm/run_patient.sh PT001 96

# Interactive, if the session will stay alive
bash slurm/run_snakemake.sh
```

`slurm/run_snakemake.sh` reads these environment variables, so you rarely need to edit it:

| Variable | Default | Purpose |
|---|---|---|
| `CONDA_ENV` | `long_read_pipeline` | conda environment to activate |
| `SINGULARITY_CACHEDIR` | `$HOME/singularity_cache` | container cache |
| `TMP_BASE` | `$HOME/tmp` | scratch directory for the run |
| `SNAKE_CORES` | 384 | `--cores` |
| `SNAKE_JOBS` | 10 | `--jobs` |

## Monitoring

```bash
# Scheduler log
tail -f snakemake_fullpower.o<JOBID>
tail -f snakemake_fullpower.e<JOBID>

# What is left to do
snakemake -n --cores 4 | tail -30

# Per-rule logs
ls -lt results/PT001/logs/ | head
tail -f results/PT001/logs/severus_PT001_multimode.log

# Files produced so far
find results/ -name "*.vcf.gz" -o -name "*.tsv" | wc -l
```

If the run died and left a lock behind:

```bash
snakemake --unlock
```

## Where the results are

| What | Path |
|---|---|
| Germline variants | `results/<p>/annotation/<p>.NORMAL.germline.vep.vcf.gz` |
| Somatic variants | `results/<p>/annotation/<p>.<TUMOR>.somatic.vep.vcf.gz` |
| Annotated SVs in cancer genes | `results/<p>/sv/<p>.sv.annotsv_intogenCCG.tsv` |
| SV circos plot | `results/<p>/sv/circos_<p>/<p>.sv.circos.svg` |
| CNV (SAVANA) | `results/<p>/cnv/savana_<p>_<TUMOR>/<p>.<TUMOR>.classified.somatic.vcf` |
| CNV (Wakhan) | `results/<p>/cnv/wakhan_<p>_<TUMOR>/<p>.<TUMOR>.copynumbers_segments.bed` |
| Purity / ploidy (PURPLE) | `results/<p>/cnv/purple_<p>_<TUMOR>/purity_ploidy.tsv` |
| Allele-specific CN (PURPLE) | `results/<p>/cnv/purple_<p>_<TUMOR>/purple/<p>.<TUMOR>.purple.cnv.somatic.tsv` |
| HRD prediction | `results/<p>/biomarkers/<p>.<TUMOR>_chord_prediction.txt` |
| Mutational signatures | `results/<p>/biomarkers/<p>.<TUMOR>.mut_sigs.tsv` |
| MSI score | `results/<p>/biomarkers/<p>.<TUMOR>.owl-scores.txt` |
| TMB | `results/<p>/biomarkers/<p>.<TUMOR>.tmb_estimate.json` (and `.gencode_coding.json`) |
| Methylation | `results/<p>/methylation/<p>.<TYPE>.cpg.combined.bed.gz` |
| Differential methylation | `results/<p>/dmr/<p>.<TUMOR>_vs_NORMAL.annotated_DMR.tsv.gz` |

Three purity/ploidy estimates are produced — SAVANA, Wakhan and PURPLE. They disagree on low-purity
samples; PURPLE is usually the one to quote, with the other two as a sanity check.

## Troubleshooting

### Out of memory

```bash
snakemake --cores 48          # instead of 96
```

or lower `threads` / `threads_low` in `config.yaml`. The heavy steps are alignment, DeepSomatic,
Severus and the HMFtools jars (`purple_java_mem`, default 60G).

### Out of disk

```bash
du -sh results/*/ | sort -h | tail
# Intermediates that can go once the run is finished:
rm -rf results/*/snv/clair3_*/tmp*
```

### Container problems

```bash
module load Program/singularity-4.2.2     # if your cluster uses modules
rm -rf "$SINGULARITY_CACHEDIR"/*          # stale cache
singularity exec resource/container/chord.sif ls /opt/chord
```

### conda environment scripts fail with `bad interpreter`

Environments hardcode absolute paths in shebangs, so moving or renaming a home directory breaks
every console script in them. Call the module directly:

```bash
/path/to/env/bin/python -m snakemake --version
```

### Rules rerun unexpectedly

Usually a timestamp change on a config or environment file:

```bash
snakemake --rerun-triggers mtime --cores 48
```

### A single rule keeps failing

Read its log first — every rule writes one:

```bash
cat results/PT001/logs/<rule>_PT001_<sample>.log
```

Then rerun just that target with `--printshellcmds` to see the exact command.
