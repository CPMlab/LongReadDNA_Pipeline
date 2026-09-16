#!/bin/bash
#SBATCH -J snakemake_fullpower
#SBATCH -p debug          # adjust to your cluster
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=300G
#SBATCH -o snakemake_fullpower.o%j
#SBATCH -e snakemake_fullpower.e%j

# Environment
source ~/.bashrc

# Activate the conda environment
# conda environment created from env.yaml (default: long_read_pipeline)
conda activate "${CONDA_ENV:-long_read_pipeline}"

# Singularity cache (avoids concurrent-pull problems)
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$HOME/singularity_cache}"
mkdir -p ${SINGULARITY_CACHEDIR}
echo "Singularity cache: ${SINGULARITY_CACHEDIR}"

# Working directory
WORK_DIR=$(pwd)
TMP_DIR="${TMP_BASE:-$HOME/tmp}/snakemake_${SLURM_JOB_ID}"

# Create the temporary directory
mkdir -p ${TMP_DIR}

echo "=== Snakemake workflow start ==="
echo "Working directory: ${WORK_DIR}"
echo "Temporary directory: ${TMP_DIR}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Nodes: ${SLURM_JOB_NODELIST}"
echo "Resources: up to 384 threads and 2 TB of memory"
echo "Heavy rules (mapping, clair3, severus) run one at a time with all threads"
echo "Conda environment: $(conda info --envs | grep '*')"
echo "Tool check:"
echo "  - pbmm2: $(which pbmm2)"
echo "  - samtools: $(which samtools)"
echo "  - severus: $(which severus)"
echo "Start time: $(date)"

# Run Snakemake locally on the allocated node(s)
# Release a stale lock from a previous aborted run
echo "Releasing any stale lock"
snakemake --unlock --directory ${WORK_DIR}

echo "Starting the workflow"
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

echo "=== Snakemake workflow finished ==="
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"

# Clean up the temporary directory
echo "Removing temporary directory: ${TMP_DIR}"
rm -rf ${TMP_DIR}

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Workflow completed successfully."
else
    echo "Workflow failed (exit code: ${EXIT_CODE})."
fi

exit ${EXIT_CODE}