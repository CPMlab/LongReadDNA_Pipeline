#!/bin/bash
# 특정 환자 한 명만 분석 (samples.tsv 에서 sample_type 자동 인식)
# 사용법: bash slurm/run_patient.sh PATIENT_ID [CORES]
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "사용법: $0 PATIENT_ID [CORES]"
    echo "예시: $0 PT001 96"
    exit 1
fi

PATIENT_ID=$1
CORES=${2:-96}
SAMPLES_TSV=${SAMPLES_TSV:-samples.tsv}
OUTPUT_DIR=${OUTPUT_DIR:-results}

source ~/.bashrc
conda activate "${CONDA_ENV:-long_read_pipeline}"

export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$HOME/singularity_cache}"
mkdir -p "${SINGULARITY_CACHEDIR}"

WORK_DIR=$(pwd)
TMP_DIR="${TMP_BASE:-$HOME/tmp}/${PATIENT_ID}_$$"
mkdir -p "${TMP_DIR}"

if [ ! -f "${SAMPLES_TSV}" ]; then
    echo "❌ ${SAMPLES_TSV} 파일을 찾을 수 없습니다!"
    exit 1
fi

# 이 환자의 샘플 타입 목록 (NORMAL + 종양 타입들)
mapfile -t SAMPLE_TYPES < <(awk -F'\t' -v p="${PATIENT_ID}" 'NR>1 && $1==p {print $2}' "${SAMPLES_TSV}" | sort -u)
if [ ${#SAMPLE_TYPES[@]} -eq 0 ]; then
    echo "❌ 환자 ${PATIENT_ID} 를 ${SAMPLES_TSV} 에서 찾을 수 없습니다!"
    exit 1
fi

TUMOR_TYPES=()
for st in "${SAMPLE_TYPES[@]}"; do
    [ "$st" != "NORMAL" ] && TUMOR_TYPES+=("$st")
done

echo "=== 환자별 분석: ${PATIENT_ID} ==="
echo "샘플 타입: ${SAMPLE_TYPES[*]}"
echo "종양 타입: ${TUMOR_TYPES[*]}"
echo "시작 시간: $(date)"

# 분석 타겟 구성
TARGETS=()
for st in "${SAMPLE_TYPES[@]}"; do
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/mapping/${PATIENT_ID}.${st}.aligned.bam")
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/methylation/${PATIENT_ID}.${st}.cpg.combined.bed.gz")
done
TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/annotation/${PATIENT_ID}.NORMAL.germline.vep.vcf.gz")
TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/sv/${PATIENT_ID}.sv.annotsv_intogenCCG.tsv")
for tt in "${TUMOR_TYPES[@]}"; do
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/annotation/${PATIENT_ID}.${tt}.somatic.vep.vcf.gz")
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/cnv/savana_${PATIENT_ID}_${tt}/${PATIENT_ID}.${tt}.classified.somatic.vcf")
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/cnv/wakhan_${PATIENT_ID}_${tt}/${PATIENT_ID}.${tt}.copynumbers_segments.bed")
    TARGETS+=("${OUTPUT_DIR}/${PATIENT_ID}/dmr/${PATIENT_ID}.${tt}_vs_NORMAL.annotated_DMR.tsv.gz")
done

echo "분석 타겟: ${#TARGETS[@]}개 파일"

snakemake --unlock --directory "${WORK_DIR}"

set +e
snakemake \
    --directory "${WORK_DIR}" \
    --cores "${CORES}" \
    --jobs "${SNAKE_JOBS:-8}" \
    --latency-wait 120 \
    --rerun-incomplete \
    --keep-going \
    --use-singularity \
    --singularity-args "-B ${WORK_DIR}:${WORK_DIR} -B ${TMP_DIR}:${TMP_DIR}" \
    --default-resources "tmpdir='${TMP_DIR}'" "mem_mb=100000" \
    --resources "mem_mb=500000" "deepsomatic_slots=1" \
    --printshellcmds \
    "${TARGETS[@]}"
EXIT_CODE=$?
set -e

echo "=== 환자 ${PATIENT_ID} 분석 종료 (코드: ${EXIT_CODE}) ==="
echo "종료 시간: $(date)"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ 완료. 결과 위치: ${OUTPUT_DIR}/${PATIENT_ID}/"
    for target in "${TARGETS[@]}"; do
        [ -f "${target}" ] && echo "  ✅ ${target}" || echo "  ❌ ${target} (생성되지 않음)"
    done
else
    echo "❌ 분석 중 오류가 발생했습니다."
fi

rm -rf "${TMP_DIR}"
exit ${EXIT_CODE}
