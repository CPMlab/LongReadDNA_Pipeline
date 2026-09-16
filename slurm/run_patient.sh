#!/bin/bash
# Run a single patient; sample types are read from samples.tsv
# Usage: bash slurm/run_patient.sh PATIENT_ID [CORES]
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 PATIENT_ID [CORES]"
    echo "Example: $0 PT001 96"
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
    echo "ERROR: ${SAMPLES_TSV} not found"
    exit 1
fi

# Sample types of this patient (NORMAL plus tumor types)
mapfile -t SAMPLE_TYPES < <(awk -F'\t' -v p="${PATIENT_ID}" 'NR>1 && $1==p {print $2}' "${SAMPLES_TSV}" | sort -u)
if [ ${#SAMPLE_TYPES[@]} -eq 0 ]; then
    echo "ERROR: patient ${PATIENT_ID} not found in ${SAMPLES_TSV}"
    exit 1
fi

TUMOR_TYPES=()
for st in "${SAMPLE_TYPES[@]}"; do
    [ "$st" != "NORMAL" ] && TUMOR_TYPES+=("$st")
done

echo "=== Single-patient run: ${PATIENT_ID} ==="
echo "Sample types: ${SAMPLE_TYPES[*]}"
echo "Tumor types: ${TUMOR_TYPES[*]}"
echo "Start time: $(date)"

# Build the target list
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

echo "Targets: ${#TARGETS[@]} files"

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

echo "=== Finished ${PATIENT_ID} (exit code: ${EXIT_CODE}) ==="
echo "End time: $(date)"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Done. Results: ${OUTPUT_DIR}/${PATIENT_ID}/"
    for target in "${TARGETS[@]}"; do
        [ -f "${target}" ] && echo "  OK ${target}" || echo "  MISSING ${target}"
    done
else
    echo "The run failed."
fi

rm -rf "${TMP_DIR}"
exit ${EXIT_CODE}
