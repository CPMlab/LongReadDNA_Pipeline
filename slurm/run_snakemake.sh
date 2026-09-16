#!/bin/bash
#SBATCH -J snakemake_fullpower
#SBATCH -p debug          # 클러스터에 맞게 수정
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=300G
#SBATCH -o snakemake_fullpower.o%j
#SBATCH -e snakemake_fullpower.e%j

# 환경 설정
source ~/.bashrc

# Conda 환경 활성화
# env.yaml 로 만든 conda 환경 (기본: long_read_pipeline)
conda activate "${CONDA_ENV:-long_read_pipeline}"

# Singularity 캐시 설정 (컨테이너 동시성 문제 해결)
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$HOME/singularity_cache}"
mkdir -p ${SINGULARITY_CACHEDIR}
echo "Singularity 캐시 디렉토리: ${SINGULARITY_CACHEDIR}"

# 작업 디렉토리 설정
WORK_DIR=$(pwd)
TMP_DIR="${TMP_BASE:-$HOME/tmp}/snakemake_${SLURM_JOB_ID}"

# 임시 디렉토리 생성
mkdir -p ${TMP_DIR}

echo "=== 🚀 FULL POWER Snakemake 워크플로우 시작 🚀 ==="
echo "작업 디렉토리: ${WORK_DIR}"
echo "임시 디렉토리: ${TMP_DIR}"
echo "작업 ID: ${SLURM_JOB_ID}"
echo "할당된 노드: ${SLURM_JOB_NODELIST}"
echo "🔥 FULL POWER 리소스: 384 스레드 (4×96) + 2TB 메모리 (4×500GB) 🔥"
echo "실행 방식: mapping(384스레드) → clair3(384스레드) → severus(384스레드) 순차 풀파워"
echo "활성화된 Conda 환경: $(conda info --envs | grep '*')"
echo "사용 가능한 도구들 확인:"
echo "  - pbmm2: $(which pbmm2)"
echo "  - samtools: $(which samtools)"
echo "  - severus: $(which severus)"
echo "시작 시간: $(date)"

# Snakemake 실행 (FULL POWER 로컬 모드)
# 먼저 락 해제 (이전 작업이 비정상 종료된 경우 대비)
echo "기존 락 파일 확인 및 해제 중..."
snakemake --unlock --directory ${WORK_DIR}

echo "Snakemake 워크플로우 실행 시작..."
snakemake \
    --directory ${WORK_DIR} \
    --cores "${SNAKE_CORES:-384}" \
    --jobs "${SNAKE_JOBS:-10}" \
    --latency-wait 180 \
    --rerun-incomplete \
    --keep-going \
    --use-singularity \
    --singularity-args "-B ${WORK_DIR}:${WORK_DIR} -B ${TMP_DIR}:${TMP_DIR}" \
    --default-resources "tmpdir='${TMP_DIR}'" "mem_mb=400000" \
    --resources "mem_mb=2000000" "deepsomatic_slots=1" \
    --printshellcmds \
    --verbose \
    all

EXIT_CODE=$?

echo "=== 🎉 FULL POWER Snakemake 워크플로우 완료 🎉 ==="
echo "종료 시간: $(date)"
echo "종료 코드: ${EXIT_CODE}"

# 작업 완료 후 임시 디렉토리 정리
echo "임시 디렉토리 정리 중: ${TMP_DIR}"
rm -rf ${TMP_DIR}

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ 워크플로우가 성공적으로 완료되었습니다!"
else
    echo "❌ 워크플로우 실행 중 오류가 발생했습니다. (종료 코드: ${EXIT_CODE})"
fi

exit ${EXIT_CODE}