# LONG_READ_DNA_WGS

A Snakemake workflow for somatic cancer genome analysis from **PacBio HiFi whole-genome sequencing**
(tumor/normal). It covers small variants, structural variants, copy number, purity/ploidy,
CpG methylation and differential methylation, and clinical biomarkers (HRD, MSI, TMB,
mutational signatures) — organised per patient.

The analysis design follows PacBio's
**[HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL)**; this repository
re-implements it in Snakemake for a SLURM + conda + Singularity cluster, and extends it to handle
several tumor samples per patient (for example primary and metastasis) in a single run.

---

## How it works

You describe your cohort in one `samples.tsv` — patient, sample type, BAM path — and the workflow
does the rest. Every sample type other than `NORMAL` is treated as a tumor sample, so the same code
runs a plain tumor/normal pair and a primary + metastasis patient without any change.

![workflow rulegraph](docs/figures/rulegraph.png)

### Pipeline steps

| Stage | Tool | Rule |
|---|---|---|
| Alignment | pbmm2 | `mapping` |
| Coverage QC | mosdepth | `mosdepth` |
| Germline small variants | Clair3 | `clair3` |
| Somatic small variants | DeepSomatic | `deepsomatic` |
| Phasing | HiPhase | `hiphase_normal`, `hiphase_tumor` |
| VCF normalisation | bcftools | `normalize_normal_vcf`, `normalize_tumor_vcf` |
| Variant annotation | VEP | `vep_annotate_normal`, `vep_annotate_somatic` |
| Structural variants | Severus (multimode) | `severus_multimode` |
| SV filtering | tabix, SVpack, mate-BND recovery | `tabix_filter`, `svpack`, `recover_mate_bnd` |
| SV annotation | AnnotSV + IntOGen cancer genes | `annotsv`, `sv_intogen` |
| Copy number | SAVANA, Wakhan | `savana_sv`, `wakhan_cnv` |
| Purity / ploidy / allele-specific CN | Amber → Cobalt → PURPLE (HMFtools) | `purple_amber`, `purple_cobalt`, `purple` |
| HRD prediction | CHORD | `chord_hrd` |
| Mutational signatures | MutationalPatterns | `mutational_signature` |
| MSI | owl | `owl_msi_profile`, `owl_msi_score` |
| TMB | tmb-calculator | `tmb_estimate` |
| Methylation | pb-CpG-tools | `cpg_methylation` |
| Differential methylation | DSS + annotatr | `dss_dmr`, `annotate_dmr` |
| Visualisation | circos (bundled script) | `circosplot`, `collect_sample_circos` |

Severus runs in multimode: the normal sample is the control and every tumor sample of that patient
is analysed together, so SVs shared between primary and metastasis stay comparable.

## Quick start

```bash
# 1. Environment (details in docs/INSTALLATION.md)
mamba env create -f env.yaml
conda activate long_read_pipeline

# 2. Describe your samples
cp samples.example.tumor_normal.tsv samples.tsv    # or samples.example.multi_tumor.tsv
vi samples.tsv

# 3. Point config.yaml at your reference data and containers
vi config.yaml

# 4. Check the plan before running anything
snakemake -n --cores 4

# 5. Run (SLURM)
sbatch slurm/run_snakemake.sh

# ...or a single patient
bash slurm/run_patient.sh PT001 96
```

## Sample sheet

```tsv
patient_id	sample_type	bam_path
PT001	NORMAL	/path/to/PT001_N.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_1.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_2.hifi_reads.bam
```

- `NORMAL` is reserved for the control sample and is required for every patient.
- Any other label (`TUMOR`, `PRIMARY`, `META`, `RELAPSE`, …) is treated as a tumor sample.
- Repeating a `patient_id` + `sample_type` (one line per SMRT cell) merges those BAMs during alignment.
- To add a patient, add rows — existing results are untouched and only the new patient is processed.
- `samples.tsv` is gitignored because it holds patient identifiers; only the anonymised examples are tracked.

## Output layout

```
results/
└── PT001/
    ├── mapping/         aligned BAMs
    ├── qc/              mosdepth coverage
    ├── snv/             Clair3 germline, DeepSomatic somatic
    ├── phasing/         HiPhase BAM/VCF, normalised VCFs
    ├── annotation/      VEP-annotated germline and somatic VCFs
    ├── sv/              Severus → SVpack → AnnotSV, IntOGen hits, circos
    ├── cnv/             SAVANA, Wakhan, PURPLE (purity_ploidy.tsv)
    ├── biomarkers/      CHORD HRD, mutational signatures, MSI, TMB
    ├── methylation/     pb-CpG-tools 5mC
    ├── dmr/             DSS differential methylation (tumor vs normal)
    └── logs/
```

## Input data scale

Input is **unaligned BAM** (`*.hifi_reads.bam`) straight from the instrument. From our breast cancer
WGS cohort (~30x): roughly 200–370 GB per SMRT cell, so about **1 TB of input for one
normal + tumor patient**, and 2–3 TB of working space once intermediates and results are included.
Worth checking before a run — the alignment step is where a full disk usually stops the pipeline.

## Reference data to prepare

Only code and small reference files are tracked here. The following go under `resource/`
(gitignored) — see [docs/INSTALLATION.md](docs/INSTALLATION.md) for how to fetch them.

| config key | Content | Approx. size |
|---|---|---|
| `reference_fasta` | GRCh38 no-alt analysis set | ~3 GB |
| `vep_cache` | VEP 112 RefSeq cache | ~26 GB |
| `annotsv_cache` | AnnotSV annotations | ~5 GB |
| `clair3/deepsomatic/savana/wakhan_container` | Variant calling images | ~11 GB |
| `chord/somatic_r_tools/owl/tmb_container` | Biomarker images | ~6 GB |
| `amber/cobalt/purple_jar`, `hmf_resources_tarball` | HMFtools jars + GRCh38 resources | ~9 GB |
| `vntr_bed` | Severus VNTR BED | ~7 MB |
| `svpack_match_vcf`, `reference_gff` | SVpack control VCF, GFF3 | ~700 MB |
| `compendium_file` | IntOGen Compendium Cancer Genes | ~1 MB |
| `mitelman_mcgene` | Mitelman database MCGENE dump | ~4 MB |

## Validation

Run end to end on a breast cancer long-read cohort: three tumor/normal pairs and two patients with
primary + metastasis samples. These were previously two separate copies of the pipeline
(`TUMOR/NORMAL` and `PRIMARY/META`); they are merged here into one, and both configurations were
confirmed to build a complete DAG after the merge.

| Configuration | Jobs |
|---|---|
| 1 normal + 1 tumor | 28 |
| 1 normal + primary + metastasis | 40 |
| the above, with PURPLE and biomarkers | 57 |

## Not included

Compared with HiFi-somatic-WDL, this repository does not implement CNVkit segmentation, the summary
HTML report, seqkit/csvtk alignment statistics, tumor-only mode, or ClairS as an alternative somatic
caller (DeepSomatic only).

## Credits

- Workflow design: [HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL) (PacBio) — please
  follow its licence and citations when using this workflow.
- `scripts/svpack.py` from PacBio [svpack](https://github.com/PacificBiosciences/svpack).
- `scripts/circosplot.py`, `scripts/DSS_tumor_normal.R`, `scripts/annotatr_dmr.R` are part of this repository.
- IntOGen Compendium of Cancer Genes and the Mitelman Database are redistributed under their own terms and
  are not bundled here.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — sample sheet, running specific stages, output structure
- [docs/OPERATIONS.md](docs/OPERATIONS.md) — running, monitoring, troubleshooting
- [docs/INSTALLATION.md](docs/INSTALLATION.md) — environment, containers, reference data
