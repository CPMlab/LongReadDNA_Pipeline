# Somatic cancer genome workflow for PacBio HiFi long reads, driven by a TSV sample sheet
# Usage: snakemake --cores N [targets] --configfile config.yaml

import os
import re
import pandas as pd
from os.path import join
from collections import defaultdict

# Configuration
configfile: "config.yaml"

# Load sample information from the TSV sample sheet
def load_samples(samples_file):
    """Read the sample sheet and return {patient: {sample_type: [bam, ...]}}."""
    df = pd.read_csv(samples_file, sep='\t')
    
    # Required columns
    required_cols = ['patient_id', 'sample_type', 'bam_path']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' is missing from {samples_file}.")
    
    # Group samples by patient
    samples_dict = defaultdict(lambda: defaultdict(list))
    
    for _, row in df.iterrows():
        patient = row['patient_id']
        sample_type = row['sample_type']
        bam_path = row['bam_path']
        
        samples_dict[patient][sample_type].append(bam_path)
    
    return dict(samples_dict)

# Load samples
SAMPLES_DATA = load_samples(config["samples_file"])

# Patients
PATIENTS = list(SAMPLES_DATA.keys())

# Sample types per patient
def get_tumor_samples(patient):
    """Return the tumor sample types of a patient (everything except NORMAL)."""
    sample_types = list(SAMPLES_DATA[patient].keys())
    return [s for s in sample_types if s != 'NORMAL']

def get_normal_sample(patient):
    """Return the control sample type of a patient."""
    if 'NORMAL' in SAMPLES_DATA[patient]:
        return 'NORMAL'
    else:
        raise ValueError(f"Patient {patient} has no NORMAL sample.")

# Globals
REF_FASTA = config["reference_fasta"]
VEP_CACHE = config["vep_cache"]
VNTR_BED = config["vntr_bed"]
THREADS = config["threads"]
OUTPUT_DIR = config["output_dir"]

# Patient / sample-type combinations used to build the target list
def get_all_sample_combinations():
    """All (patient, sample_type) pairs."""
    combinations = []
    for patient in PATIENTS:
        for sample_type in SAMPLES_DATA[patient].keys():
            combinations.append((patient, sample_type))
    return combinations

def get_all_tumor_combinations():
    """All (patient, tumor_sample_type) pairs."""
    combinations = []
    for patient in PATIENTS:
        for tumor_type in get_tumor_samples(patient):
            combinations.append((patient, tumor_type))
    return combinations

ALL_SAMPLE_COMBINATIONS = get_all_sample_combinations()
ALL_TUMOR_COMBINATIONS = get_all_tumor_combinations()

# Tumor sample type constraint, derived from samples.tsv
#   - tumor/normal pairs        -> "TUMOR"
#   - primary + metastasis      -> "PRIMARY|META"
# This used to be hardcoded in rules/*.smk, which meant copying the pipeline per dataset.
TUMOR_SAMPLE_TYPES = sorted({t for p in PATIENTS for t in get_tumor_samples(p)})
if not TUMOR_SAMPLE_TYPES:
    raise ValueError("samples.tsv contains no tumor sample type other than NORMAL.")
TUMOR_TYPE_CONSTRAINT = "|".join(re.escape(t) for t in TUMOR_SAMPLE_TYPES)

# Rule priority: the normal rule wins over the tumor rule
ruleorder: hiphase_normal > hiphase_tumor

# Rule files
include: "rules/mapping.smk"
include: "rules/variant_calling.smk"
include: "rules/annotation.smk"
include: "rules/structural_variants.smk"
include: "rules/visualization.smk"
include: "rules/methylation.smk"
include: "rules/copy_number.smk"
include: "rules/purple.smk"
include: "rules/biomarkers.smk"

rule all:
    input:
        # Annotated structural variants
        expand(join(OUTPUT_DIR, "{patient}", "sv", "{patient}.sv.annotsv_intogenCCG.tsv"), patient=PATIENTS),
        # Alignments
        [join(OUTPUT_DIR, patient, "mapping", f"{patient}.{sample_type}.aligned.bam") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # Coverage QC
        [join(OUTPUT_DIR, patient, "qc", f"{patient}.{sample_type}.mosdepth.summary.txt") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # Germline variants (normal sample)
        expand(join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.NORMAL.germline.vep.vcf.gz"), patient=PATIENTS),
        # Somatic variants (tumor samples)
        [join(OUTPUT_DIR, patient, "annotation", f"{patient}.{tumor_type}.somatic.vep.vcf.gz") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # SV circos (multi-sample)
        expand(join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "{patient}_fusion_calls_combined.tsv"), patient=PATIENTS),
        expand(join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "sample_files.txt"), patient=PATIENTS),
        # Methylation (all samples)
        [join(OUTPUT_DIR, patient, "methylation", f"{patient}.{sample_type}.cpg.combined.bed.gz") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "methylation", f"{patient}.{sample_type}.cpg.combined.bw") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # Copy number (SAVANA)
        [join(OUTPUT_DIR, patient, "cnv", f"savana_{patient}_{tumor_type}", f"{patient}.{tumor_type}.classified.somatic.vcf") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"savana_{patient}_{tumor_type}", f"{patient}.{tumor_type}_fitted_purity_ploidy.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Copy number (Wakhan)
        [join(OUTPUT_DIR, patient, "cnv", f"wakhan_{patient}_{tumor_type}", f"{patient}.{tumor_type}.copynumbers_segments.bed") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"wakhan_{patient}_{tumor_type}", "purity_ploidy.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Purity/ploidy and allele-specific copy number (PURPLE)
        [join(OUTPUT_DIR, patient, "cnv", f"purple_{patient}_{tumor_type}", "purity_ploidy.tsv")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"purple_{patient}_{tumor_type}", "purple", f"{patient}.{tumor_type}.purple.cnv.somatic.tsv")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # HRD prediction (CHORD)
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}_chord_prediction.txt")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Mutational signatures (MutationalPatterns)
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}.mut_sigs.tsv")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # MSI (owl)
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}.owl-scores.txt")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # TMB
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}.tmb_estimate.json")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}.tmb_estimate.gencode_coding.json")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Differential methylation (DSS)
        [join(OUTPUT_DIR, patient, "dmr", f"{patient}.{tumor_type}_vs_NORMAL.DMR.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "dmr", f"{patient}.{tumor_type}_vs_NORMAL.annotated_DMR.tsv.gz") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS] 