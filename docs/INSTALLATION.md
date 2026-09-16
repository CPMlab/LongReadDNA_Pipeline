# Installation

## Contents

1. [Requirements](#1-requirements)
2. [Conda environment](#2-conda-environment)
3. [Container images](#3-container-images)
4. [Reference data](#4-reference-data)
5. [HMFtools (PURPLE) setup](#5-hmftools-purple-setup)
6. [Verifying the installation](#6-verifying-the-installation)
7. [Notes and troubleshooting](#7-notes-and-troubleshooting)

## 1. Requirements

- **OS**: Linux (RHEL 8/9, Ubuntu 20.04+)
- **Memory**: 64 GB minimum, 128 GB+ recommended (500 GB+ for the full cohort runs we do)
- **Disk**: 2–3 TB of working space per patient, plus ~60 GB for reference data and containers
- **CPU**: 48+ cores recommended
- **Java**: 11 or later, required by the HMFtools jars (Amber/Cobalt/PURPLE)
- **Singularity / Apptainer**: version 3 or later

On a module-based cluster, load Singularity first:

```bash
module load Program/singularity-4.2.2    # adjust to your cluster
```

### Conda / Mamba

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

## 2. Conda environment

```bash
cd LONG_READ_DNA_WGS
mamba env create -f env.yaml
conda activate long_read_pipeline

snakemake --version
samtools --version
```

## 3. Container images

```bash
mkdir -p resource/container

# --- variant calling / SV / CNV ---
singularity pull resource/container/clair3_latest.sif      docker://hkubal/clair3:latest
singularity pull resource/container/deepsomatic_1.8.0.sif  docker://google/deepsomatic:1.8.0
singularity pull resource/container/savana_1.3.4.sif       docker://quay.io/biocontainers/savana
singularity pull resource/container/wakhan_latest.sif      docker://mkolmogo/wakhan:dev_c717baa

# --- biomarkers ---
# Digests pinned to the versions used by HiFi-somatic-WDL.
singularity pull resource/container/chord.sif \
  docker://scwatts/chord@sha256:9f6aa44ffefe3f736e66a0e2d7941d4f3e1cc6d848a9a11a17e85a6525e63a77
singularity pull resource/container/somatic_r_tools.sif \
  docker://quay.io/pacbio/somatic_r_tools@sha256:68dc04908a37e26b30dc9795fa6cc0e85a238c8695afe805ad164a071193fb48
singularity pull resource/container/owl.sif \
  docker://quay.io/pacbio/owl@sha256:753b83abe1fb5d8c4f1e2e4ef200bfbecf1a342827ab4974a0d271911675461d
singularity pull resource/container/tmb_calculator.sif \
  docker://quay.io/pacbio/tmb_calculator@sha256:93f89b7f2777bb27fc7e8ba5fb0b54a56c132c7b5a5b3f01b95696b5b0b3b63b
```

Expect roughly 17 GB in total.

Build these on **node-local disk, not on a network filesystem**. Converting a Docker image to SIF
unpacks its whole root filesystem first, which is hundreds of thousands of small files; on Lustre or
NFS that crawls, and a large image such as `somatic_r_tools` can sit for over an hour with no
progress and then produce a corrupt SIF. Set both variables to a local path and copy the finished
`.sif` to shared storage afterwards:

```bash
LOCAL=/tmp/$USER/singularity
export SINGULARITY_CACHEDIR=$LOCAL/cache
export SINGULARITY_TMPDIR=$LOCAL/tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

singularity pull "$LOCAL/somatic_r_tools.sif" docker://quay.io/pacbio/somatic_r_tools@sha256:...
mv "$LOCAL/somatic_r_tools.sif" resource/container/
```

Check each image afterwards. A truncated build still looks like a normal file, and only fails when
a rule tries to run it:

```bash
singularity exec resource/container/chord.sif ls /opt/chord/extractSigPredictHRD.R
singularity exec resource/container/somatic_r_tools.sif ls /app/mutational_pattern.R
singularity exec resource/container/owl.sif owl --version
singularity exec resource/container/tmb_calculator.sif ls /opt/venv/bin/calculate_tmb.py
```

`bad superblock for squashfs image` means the pull did not finish; delete the file and pull again.

## 4. Reference data

### 4.1 Reference genome (required)

```bash
mkdir -p resource/ref
wget -O resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"
gunzip resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz
samtools faidx resource/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta
```

### 4.2 VEP cache (required)

```bash
mkdir -p resource/vep_cache
wget -O resource/vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz \
  "https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz"
```

### 4.3 AnnotSV annotations (required)

```bash
mkdir -p resource/annotsv
wget -O resource/annotsv/annotsv_cache.tar.gz \
  "https://www.lbgi.fr/~geoffroy/Annotations/AnnotSV_annotations_3.4.tar.gz"
```

`scripts/install_annotsv_annotations.sh` unpacks it into the layout AnnotSV expects.

### 4.4 Other resources

```bash
# Severus VNTR BED
mkdir -p resource/severus
wget -O resource/severus/human_GRCh38_no_alt_analysis_set.trf.bed \
  "https://github.com/KolmogorovLab/Severus/raw/main/resources/human_GRCh38_no_alt_analysis_set.trf.bed"

# IntOGen Compendium of Cancer Genes — download from https://www.intogen.org
mkdir -p resource/intogen_genelist

# Mitelman database MCGENE dump (used to flag known fusions in the circos plot)
# https://mitelmandatabase.isb-cgc.org  -> place the MCGENE table at resource/circos/mitel
mkdir -p resource/circos
```

`resources/chr.bed`, `resources/hg38.bed` and `resources/hg38_cytoband.tsv` are already in this
repository; no download needed.

## 5. HMFtools (PURPLE) setup

PURPLE runs from local jars rather than a container.

```bash
mkdir -p resource/purple
# Amber, Cobalt and PURPLE jars (HMFtools releases)
#   https://github.com/hartwigmedical/hmftools
# The pipeline was validated with:
#   amber-4.0-jar-with-dependencies.gamma1000.jar
#   cobalt-1.16.0-jar-with-dependencies.jar
#   purple-4.0-jar-with-dependencies.jar

# GRCh38 resource bundle (~9 GB), from the HMF resource release:
#   hmf_pipeline_resources.38_v2.0.0--3.tar.gz
```

Place all four files in `resource/purple/` and check the paths in `config.yaml`
(`amber_jar`, `cobalt_jar`, `purple_jar`, `hmf_resources_tarball`). The workflow unpacks the bundle
once into `results/shared/hmf_resources/` and reuses it for every sample.

## 6. Verifying the installation

```bash
conda activate long_read_pipeline

pbmm2 --version
samtools --version
bcftools --version
hiphase --version
severus --version
aligned_bam_to_cpg_scores --version
java -version
Rscript -e "library(DSS); library(annotatr); sessionInfo()"

singularity run resource/container/clair3_latest.sif /opt/bin/run_clair3.sh --version
singularity run resource/container/deepsomatic_1.8.0.sif run_deepsomatic --version

# Finally, a dry run against your own samples.tsv
snakemake -n --cores 4
```

## 7. Notes and troubleshooting

- **Memory**: alignment and variant calling are the heavy steps; lower `threads` in `config.yaml`
  if jobs are killed by the scheduler.
- **Runtime**: a full cohort run takes more than a day.
- **Disk**: intermediates are large — see the input data scale section in the README.
- **Conda conflicts**: recreate with `mamba env create --force -f env.yaml`.
- **Broken console scripts after a home directory move**: conda environments hardcode absolute paths
  in shebangs. Call the module directly instead — `<env>/bin/python -m snakemake`.
- **Container errors**: check the Singularity version and that `--bind` covers your working directory.
