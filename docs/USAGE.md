# Usage

How to describe your cohort, run part or all of the workflow, and find the results.

## Sample sheet

`samples.tsv` is the only file describing your data. Three columns, tab separated:

```tsv
patient_id	sample_type	bam_path
PT101	NORMAL	/data/PT101_N.hifi_reads.bam
PT101	PRIMARY	/data/PT101_primary_1.hifi_reads.bam
PT101	PRIMARY	/data/PT101_primary_2.hifi_reads.bam
PT101	META	/data/PT101_meta.hifi_reads.bam
PT102	NORMAL	/data/PT102_N.hifi_reads.bam
PT102	TUMOR	/data/PT102_T.hifi_reads.bam
```

| Column | Meaning |
|---|---|
| `patient_id` | Patient identifier; becomes the output directory name |
| `sample_type` | `NORMAL` for the control, anything else is treated as a tumor sample |
| `bam_path` | Absolute path to the unaligned HiFi BAM (`*.hifi_reads.bam`) |

Rules:

- Every patient needs exactly one `NORMAL` sample type.
- Tumor labels are free text — `TUMOR`, `PRIMARY`, `META`, `RELAPSE` all work. The workflow reads
  them from this file and builds its wildcard constraints accordingly, so no code changes are needed.
- Repeating the same `patient_id` + `sample_type` on several rows (one per SMRT cell) merges those
  BAMs during alignment.
- Two starting points are provided: `samples.example.tumor_normal.tsv` and
  `samples.example.multi_tumor.tsv`.

## Running

```bash
# Everything in samples.tsv
snakemake --configfile config.yaml --cores 48

# Dry run first — always worth it
snakemake --configfile config.yaml --cores 48 -n

# On SLURM
sbatch slurm/run_snakemake.sh

# One patient only (reads its sample types from samples.tsv)
bash slurm/run_patient.sh PT101 96
```

### Running part of the workflow

Ask for the files you want and Snakemake works out what to run:

```bash
# Stop after alignment
snakemake --configfile config.yaml --cores 48 \
    results/PT101/mapping/PT101.NORMAL.aligned.bam \
    results/PT101/mapping/PT101.PRIMARY.aligned.bam \
    results/PT101/mapping/PT101.META.aligned.bam

# Only the annotated structural variants for one patient
snakemake --configfile config.yaml --cores 48 \
    results/PT101/sv/PT101.sv.annotsv_intogenCCG.tsv

# Only the biomarkers for one tumor sample
snakemake --configfile config.yaml --cores 16 \
    results/PT101/biomarkers/PT101.PRIMARY_chord_prediction.txt \
    results/PT101/biomarkers/PT101.PRIMARY.owl-scores.txt \
    results/PT101/biomarkers/PT101.PRIMARY.tmb_estimate.json
```

## Output structure

```
results/
├── shared/
│   └── hmf_resources/          unpacked HMFtools reference bundle (shared by all samples)
└── PT101/
    ├── mapping/                PT101.<TYPE>.aligned.bam
    ├── qc/                     PT101.<TYPE>.mosdepth.summary.txt, regions.bed.gz
    ├── snv/                    Clair3 germline, DeepSomatic somatic VCFs
    ├── phasing/                HiPhase BAM/VCF, normalized VCFs
    ├── annotation/             PT101.NORMAL.germline.vep.vcf.gz
    │                           PT101.<TUMOR>.somatic.vep.vcf.gz
    ├── sv/                     severus_PT101/, svpack_PT101/, recovered_PT101/
    │                           PT101.sv.annotsv_intogenCCG.tsv
    │                           circos_PT101/PT101.sv.circos.svg
    ├── cnv/                    savana_PT101_<TUMOR>/, wakhan_PT101_<TUMOR>/
    │                           purple_PT101_<TUMOR>/purity_ploidy.tsv
    ├── biomarkers/             PT101.<TUMOR>_chord_prediction.txt      (HRD)
    │                           PT101.<TUMOR>.mut_sigs.tsv              (signatures)
    │                           PT101.<TUMOR>.owl-scores.txt            (MSI)
    │                           PT101.<TUMOR>.tmb_estimate.json         (TMB)
    ├── methylation/            PT101.<TYPE>.cpg.combined.bed.gz / .bw
    ├── dmr/                    PT101.<TUMOR>_vs_NORMAL.DMR.tsv
    │                           PT101.<TUMOR>_vs_NORMAL.annotated_DMR.tsv.gz
    └── logs/                   one log per rule invocation
```

## Patients with more than one tumor sample

`NORMAL` is the control; every other sample type is a tumor. Adding a second tumor (metastasis,
relapse, a second biopsy) needs no configuration change — just extra rows.

What changes when a patient has `NORMAL` + `PRIMARY` + `META` instead of `NORMAL` + `TUMOR`:

```
results/PT101/
├── snv/         PT101.PRIMARY.somatic.vcf.gz      PT101.META.somatic.vcf.gz
├── annotation/  PT101.PRIMARY.somatic.vep.vcf.gz  PT101.META.somatic.vep.vcf.gz
├── cnv/         savana_PT101_PRIMARY/  savana_PT101_META/
│                wakhan_PT101_PRIMARY/  wakhan_PT101_META/
│                purple_PT101_PRIMARY/  purple_PT101_META/
├── biomarkers/  PT101.PRIMARY_chord_prediction.txt   PT101.META_chord_prediction.txt
│                PT101.PRIMARY.mut_sigs.tsv           PT101.META.mut_sigs.tsv
│                PT101.PRIMARY.owl-scores.txt         PT101.META.owl-scores.txt
│                PT101.PRIMARY.tmb_estimate.json      PT101.META.tmb_estimate.json
├── dmr/         PT101.PRIMARY_vs_NORMAL.DMR.tsv      PT101.META_vs_NORMAL.DMR.tsv
└── sv/          severus_PT101/            <- ONE call covering both tumors
                 PT101.sv.annotsv_intogenCCG.tsv
                 circos_PT101/             <- one plot per tumor + combined table
```

Everything that is tumor-specific is produced once per tumor sample. Structural variants are the
exception: Severus runs in multimode with the normal as control and both tumors as targets, so a
single somatic SV VCF covers the patient, with one genotype column per tumor. Shared and private
breakpoints therefore stay directly comparable, and the circos step splits that VCF per tumor while
also writing a combined fusion table (`<patient>_fusion_calls_combined.tsv`).

If you would rather analyze the tumors completely independently, give them different `patient_id`
values instead of different sample types — but then each one needs its own `NORMAL` row and the SVs
are no longer called jointly.

## Adding a patient

Add the rows to `samples.tsv` and rerun. Existing results are left alone and only the new patient is
processed.

## Common errors

| Message | Cause |
|---|---|
| workflow stops saying a patient has no NORMAL sample | That patient has no `NORMAL` row in `samples.tsv` |
| workflow stops saying a required column is missing | Check the header line of `samples.tsv` |
| `MissingInputException` on a BAM | `bam_path` is wrong or not readable |
| Job killed by the scheduler | Lower `threads` in `config.yaml`, or request more memory |

Logs live in `results/<patient>/logs/`, one file per rule, for example
`results/PT101/logs/severus_PT101_multimode.log`.
