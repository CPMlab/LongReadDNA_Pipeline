# LONG_READ_DNA_WGS

Snakemake workflow for somatic analysis of PacBio HiFi whole-genome sequencing, tumor against
matched normal.

It covers alignment, germline and somatic small variants, structural variants, copy number with
purity and ploidy, CpG methylation and differential methylation, and the biomarkers derived from
those (HRD, MSI, TMB, mutational signatures). Output is organised per patient.

A patient may have more than one tumor sample — primary and metastasis, a relapse, a second biopsy.
They are analysed against the same normal in one run, and structural variants are called across all
of them at once. The cohort is described in a single TSV, so adding a patient or a sample means
adding rows rather than editing the workflow.

Tool choices and filtering strategy were developed against PacBio's
[HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL).

---

## Sample types

`NORMAL` is the control. Every other sample type is treated as a tumor sample, and the label is free
text: `TUMOR`, `PRIMARY`, `META`, `RELAPSE` all work, and the label becomes part of the output file
names. Each patient needs exactly one `NORMAL`.

Stages run at one of three levels:

| Level | Stages |
|---|---|
| Per sample, normal included | alignment, coverage QC, Clair3 germline calling, CpG methylation |
| Per tumor sample | DeepSomatic, tumor phasing, VCF normalisation, VEP somatic annotation, SAVANA, Wakhan, PURPLE, CHORD, mutational signatures, MSI, TMB, DSS differential methylation |
| Per patient, covering all tumors together | normal phasing and annotation, Severus multimode, SV filtering, AnnotSV, IntOGen intersection, circos |

For a patient with `NORMAL`, `PRIMARY` and `META`, that means alignment and germline calling run
three times; somatic calling, copy number, purity, the biomarkers and differential methylation run
twice, once per tumor; and Severus runs once with the normal as control and both tumors as targets.
That single somatic SV VCF carries a genotype column per tumor, so breakpoints shared between
primary and metastasis can be separated from private ones without a second comparison. The circos
step splits it again per tumor and also writes a combined fusion table.

A plain tumor/normal pair is the same thing with one tumor sample; nothing needs to be configured
differently.

If two tumors should not be compared jointly, give them separate `patient_id` values instead of
separate sample types. Each then needs its own `NORMAL` row, and SVs are called separately.

## Workflow

![workflow rulegraph](docs/figures/rulegraph.png)

| Stage | Tool | Rule |
|---|---|---|
| Alignment | pbmm2 | `mapping` |
| Coverage | mosdepth | `mosdepth` |
| Germline small variants | Clair3 | `clair3` |
| Somatic small variants | DeepSomatic | `deepsomatic` |
| Phasing | HiPhase | `hiphase_normal`, `hiphase_tumor` |
| VCF normalisation | bcftools | `normalize_normal_vcf`, `normalize_tumor_vcf` |
| Variant annotation | VEP | `vep_annotate_normal`, `vep_annotate_somatic` |
| Structural variants | Severus, multimode | `severus_multimode` |
| SV filtering | tabix, SVpack, mate-BND recovery | `tabix_filter`, `svpack`, `recover_mate_bnd` |
| SV annotation | AnnotSV, IntOGen cancer genes | `annotsv`, `sv_intogen` |
| Copy number | SAVANA, Wakhan | `savana_sv`, `wakhan_cnv` |
| Purity, ploidy, allele-specific CN | Amber, Cobalt, PURPLE | `purple_amber`, `purple_cobalt`, `purple` |
| HRD | CHORD | `chord_hrd` |
| Mutational signatures | MutationalPatterns | `mutational_signature` |
| MSI | owl | `owl_msi_profile`, `owl_msi_score` |
| TMB | tmb-calculator | `tmb_estimate` |
| Methylation | pb-CpG-tools | `cpg_methylation` |
| Differential methylation | DSS, annotatr | `dss_dmr`, `annotate_dmr` |
| Circos | bundled script | `circosplot`, `collect_sample_circos` |

Three purity and ploidy estimates are produced, from SAVANA, Wakhan and PURPLE. They disagree on
low-purity samples; PURPLE is the one we quote, with the other two as a check.

## Running

```bash
mamba env create -f env.yaml
conda activate long_read_pipeline

cp samples.example.tumor_normal.tsv samples.tsv    # or samples.example.multi_tumor.tsv
vi samples.tsv                                     # patient, sample type, BAM path
vi config.yaml                                     # reference data and container paths

snakemake -n --cores 4                             # dry run
sbatch slurm/run_snakemake.sh                      # full cohort
bash slurm/run_patient.sh PT001 96                 # one patient
```

Sample sheet format:

```tsv
patient_id	sample_type	bam_path
PT001	NORMAL	/path/to/PT001_N.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_1.hifi_reads.bam
PT001	TUMOR	/path/to/PT001_T_2.hifi_reads.bam
```

Repeating a patient and sample type, one row per SMRT cell, merges those BAMs during alignment.
`samples.tsv` is gitignored because it carries patient identifiers; only the anonymised examples are
tracked.

## Output

```
results/
├── shared/hmf_resources/   HMFtools reference bundle, unpacked once
└── PT001/
    ├── mapping/            aligned BAMs
    ├── qc/                 mosdepth coverage
    ├── snv/                Clair3 germline, DeepSomatic somatic
    ├── phasing/            HiPhase BAM and VCF, normalised VCFs
    ├── annotation/         VEP-annotated germline and somatic VCFs
    ├── sv/                 Severus, SVpack, AnnotSV, IntOGen hits, circos
    ├── cnv/                SAVANA, Wakhan, PURPLE
    ├── biomarkers/         HRD, mutational signatures, MSI, TMB
    ├── methylation/        5mC from pb-CpG-tools
    ├── dmr/                DSS differentially methylated regions
    └── logs/               one log per rule
```

## Input and disk

Input is unaligned BAM (`*.hifi_reads.bam`) as delivered by the instrument. In our breast cancer
cohort at roughly 30x, one SMRT cell is 200–370 GB, so a normal plus a two-cell tumor is about 1 TB
of input. Allow 2–3 TB of working space per patient for intermediates and results; alignment is
where a run usually stops if the disk fills.

## Reference data

Code and small reference files are tracked here. The rest goes under `resource/`, which is
gitignored — see [docs/INSTALLATION.md](docs/INSTALLATION.md).

| config key | Content | Size |
|---|---|---|
| `reference_fasta` | GRCh38 no-alt analysis set | ~3 GB |
| `vep_cache` | VEP 112 RefSeq cache | ~26 GB |
| `annotsv_cache` | AnnotSV annotations | ~5 GB |
| `clair3`, `deepsomatic`, `savana`, `wakhan` containers | variant calling images | ~11 GB |
| `chord`, `somatic_r_tools`, `owl`, `tmb` containers | biomarker images | ~6 GB |
| `amber_jar`, `cobalt_jar`, `purple_jar`, `hmf_resources_tarball` | HMFtools | ~9 GB |
| `vntr_bed` | Severus VNTR BED | ~7 MB |
| `svpack_match_vcf`, `reference_gff` | SVpack control VCF, GFF3 | ~700 MB |
| `compendium_file` | IntOGen Compendium of Cancer Genes | ~1 MB |
| `mitelman_mcgene` | Mitelman database MCGENE dump | ~4 MB |

## Validation

Run end to end on a breast cancer cohort of three tumor/normal pairs and two patients with primary
and metastasis samples. The two configurations were maintained as separate copies of the pipeline
until they were merged here; after merging, both build a complete DAG:

| Configuration | Jobs |
|---|---|
| normal + tumor | 28 |
| normal + primary + metastasis | 40 |
| the same, with PURPLE and the biomarkers | 57 |

## Not included

CNVkit segmentation, the summary HTML report, seqkit and csvtk alignment statistics, tumor-only
mode, and ClairS as an alternative somatic caller. Somatic calling is DeepSomatic only.

## Credits

Reference workflow: [HiFi-somatic-WDL](https://github.com/PacificBiosciences/HiFi-somatic-WDL)
(PacBio). Follow its licence, and cite the tools it lists, when publishing results from this
workflow.

`scripts/svpack.py` comes from PacBio [svpack](https://github.com/PacificBiosciences/svpack).
`scripts/circosplot.py`, `scripts/DSS_tumor_normal.R` and `scripts/annotatr_dmr.R` are part of this
repository. The IntOGen Compendium of Cancer Genes and the Mitelman Database are used under their
own terms and are not redistributed here.

## Documentation

- [docs/USAGE.md](docs/USAGE.md) — sample sheet, running individual stages, output structure
- [docs/OPERATIONS.md](docs/OPERATIONS.md) — running, monitoring, troubleshooting
- [docs/INSTALLATION.md](docs/INSTALLATION.md) — environment, containers, reference data
