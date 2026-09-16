# Long-read (PacBio HiFi) 암 게놈 변이 검출 워크플로우 (TSV 입력 지원)
# 사용법: snakemake --cores N [타겟명] --configfile config.yaml

import os
import re
import pandas as pd
from os.path import join
from collections import defaultdict

# 설정 파일
configfile: "config.yaml"

# TSV 파일에서 샘플 정보 읽기
def load_samples(samples_file):
    """TSV 파일에서 샘플 정보를 읽어서 구조화된 딕셔너리로 반환"""
    df = pd.read_csv(samples_file, sep='\t')
    
    # 필수 컬럼 확인
    required_cols = ['patient_id', 'sample_type', 'bam_path']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"필수 컬럼 '{col}'이 {samples_file}에 없습니다.")
    
    # 환자별로 샘플 정보 그룹화
    samples_dict = defaultdict(lambda: defaultdict(list))
    
    for _, row in df.iterrows():
        patient = row['patient_id']
        sample_type = row['sample_type']
        bam_path = row['bam_path']
        
        samples_dict[patient][sample_type].append(bam_path)
    
    return dict(samples_dict)

# 샘플 정보 로드
SAMPLES_DATA = load_samples(config["samples_file"])

# 환자 리스트
PATIENTS = list(SAMPLES_DATA.keys())

# 각 환자별 샘플 타입 확인
def get_tumor_samples(patient):
    """환자별 종양 샘플 타입들 반환 (NORMAL 제외)"""
    sample_types = list(SAMPLES_DATA[patient].keys())
    return [s for s in sample_types if s != 'NORMAL']

def get_normal_sample(patient):
    """환자별 정상 샘플 타입 반환"""
    if 'NORMAL' in SAMPLES_DATA[patient]:
        return 'NORMAL'
    else:
        raise ValueError(f"환자 {patient}에게 NORMAL 샘플이 없습니다.")

# 전역 변수
REF_FASTA = config["reference_fasta"]
VEP_CACHE = config["vep_cache"]
VNTR_BED = config["vntr_bed"]
THREADS = config["threads"]
OUTPUT_DIR = config["output_dir"]

# 생성될 최종 결과물들 - 환자별 샘플 타입 조합 생성
def get_all_sample_combinations():
    """모든 환자의 모든 샘플 타입 조합 반환"""
    combinations = []
    for patient in PATIENTS:
        for sample_type in SAMPLES_DATA[patient].keys():
            combinations.append((patient, sample_type))
    return combinations

def get_all_tumor_combinations():
    """모든 환자의 모든 종양 샘플 타입 조합 반환"""
    combinations = []
    for patient in PATIENTS:
        for tumor_type in get_tumor_samples(patient):
            combinations.append((patient, tumor_type))
    return combinations

ALL_SAMPLE_COMBINATIONS = get_all_sample_combinations()
ALL_TUMOR_COMBINATIONS = get_all_tumor_combinations()

# 종양 샘플 타입 제약 (samples.tsv 에서 자동 유도)
#   - TUMOR/NORMAL 쌍 데이터  -> "TUMOR"
#   - PRIMARY/META 다중 종양   -> "PRIMARY|META"
# 이전에는 rules/*.smk 에 하드코딩되어 데이터셋마다 파이프라인을 복제해야 했음.
TUMOR_SAMPLE_TYPES = sorted({t for p in PATIENTS for t in get_tumor_samples(p)})
if not TUMOR_SAMPLE_TYPES:
    raise ValueError("samples.tsv 에 NORMAL 이외의 종양 샘플 타입이 하나도 없습니다.")
TUMOR_TYPE_CONSTRAINT = "|".join(re.escape(t) for t in TUMOR_SAMPLE_TYPES)

# 규칙 우선순위 설정 - normal 규칙이 tumor 규칙보다 우선
ruleorder: hiphase_normal > hiphase_tumor

# 분리된 규칙 파일들 포함
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
        # 모든 환자의 최종 SV 분석 결과
        expand(join(OUTPUT_DIR, "{patient}", "sv", "{patient}.sv.annotsv_intogenCCG.tsv"), patient=PATIENTS),
        # 매핑 결과
        [join(OUTPUT_DIR, patient, "mapping", f"{patient}.{sample_type}.aligned.bam") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # QC 결과
        [join(OUTPUT_DIR, patient, "qc", f"{patient}.{sample_type}.mosdepth.summary.txt") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # SNV 결과 (정상 샘플)
        expand(join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.NORMAL.germline.vep.vcf.gz"), patient=PATIENTS),
        # 체세포 변이 결과 (종양 샘플들)
        [join(OUTPUT_DIR, patient, "annotation", f"{patient}.{tumor_type}.somatic.vep.vcf.gz") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # SV 분석 결과 (멀티샘플 circos)
        expand(join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "{patient}_fusion_calls_combined.tsv"), patient=PATIENTS),
        expand(join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "sample_files.txt"), patient=PATIENTS),
        # 메틸화 분석 결과 (모든 샘플)
        [join(OUTPUT_DIR, patient, "methylation", f"{patient}.{sample_type}.cpg.combined.bed.gz") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "methylation", f"{patient}.{sample_type}.cpg.combined.bw") 
         for patient, sample_type in ALL_SAMPLE_COMBINATIONS],
        # CNV 분석 결과 (SAVANA - 종양 샘플들)
        [join(OUTPUT_DIR, patient, "cnv", f"savana_{patient}_{tumor_type}", f"{patient}.{tumor_type}.classified.somatic.vcf") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"savana_{patient}_{tumor_type}", f"{patient}.{tumor_type}_fitted_purity_ploidy.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Copy Number 분석 결과 (Wakhan - 종양 샘플들)
        [join(OUTPUT_DIR, patient, "cnv", f"wakhan_{patient}_{tumor_type}", f"{patient}.{tumor_type}.copynumbers_segments.bed") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"wakhan_{patient}_{tumor_type}", "purity_ploidy.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # PURPLE purity/ploidy 및 allele-specific CNV (종양 샘플들)
        [join(OUTPUT_DIR, patient, "cnv", f"purple_{patient}_{tumor_type}", "purity_ploidy.tsv")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "cnv", f"purple_{patient}_{tumor_type}", "purple", f"{patient}.{tumor_type}.purple.cnv.somatic.tsv")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # HRD 예측 (CHORD)
        [join(OUTPUT_DIR, patient, "biomarkers", f"{patient}.{tumor_type}_chord_prediction.txt")
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        # Mutational signature (MutationalPatterns)
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
        # 차등 메틸화 분석 결과 (DSS - 종양 샘플들)
        [join(OUTPUT_DIR, patient, "dmr", f"{patient}.{tumor_type}_vs_NORMAL.DMR.tsv") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS],
        [join(OUTPUT_DIR, patient, "dmr", f"{patient}.{tumor_type}_vs_NORMAL.annotated_DMR.tsv.gz") 
         for patient, tumor_type in ALL_TUMOR_COMBINATIONS] 