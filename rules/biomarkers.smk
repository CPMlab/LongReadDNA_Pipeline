# 바이오마커 분석: HRD (CHORD), mutational signature (MutationalPatterns),
# MSI (owl), TMB (tmb-calculator)
# 모두 HiFi-somatic-WDL 의 annotation.wdl / common.wdl / biomarker.wdl 과 동일한 호출 방식이며
# 컨테이너 이미지는 config 의 *_container 경로(Singularity .sif)를 사용한다.

BIOMARKER_DIR = join(OUTPUT_DIR, "{patient}", "biomarkers")

# 1. CHORD - 상동재조합결핍(HRD) 예측
#    체세포 SNV/indel VCF 와 Severus SV VCF 를 함께 사용한다.
#    컨테이너의 기본 SV caller 설정이 GRIDSS 이므로 manta 로 바꿔서 실행한다 (원본 WDL 과 동일).
rule chord_hrd:
    input:
        small_variant_vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.somatic.vcf.gz"),
        sv_vcf = join(OUTPUT_DIR, "{patient}", "sv", "severus_{patient}", "somatic_SVs", "severus_somatic.vcf")
    output:
        prediction = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}_chord_prediction.txt"),
        signature = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}_chord_signatures.txt")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = BIOMARKER_DIR,
        pname = "{patient}.{tumor_sample_type}"
    threads: 4
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "chord_hrd_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== CHORD HRD 예측 시작: {params.pname} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        SNV_ABS=$(realpath {input.small_variant_vcf})
        SV_ABS=$(realpath {input.sv_vcf})
        OUT_ABS=$(realpath {params.out_dir})
        echo "SNV VCF: $SNV_ABS" >> {log}
        echo "SV VCF: $SV_ABS" >> {log}

        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[chord_container]} \
            /bin/bash -c "
                set -euxo pipefail
                cd $OUT_ABS

                # 컨테이너 기본값(GRIDSS)을 manta 로 교체
                sed 's/gridss/manta/g' /opt/chord/extractSigPredictHRD.R > ./extractSigPredictHRD.R
                chmod +x ./extractSigPredictHRD.R

                ./extractSigPredictHRD.R . {params.pname} $SNV_ABS $SV_ABS 38

                rm -f ./extractSigPredictHRD.R
            " >> {log} 2>&1

        echo "=== CHORD HRD 완료: $(date) ===" >> {log}
        """

# 2. MutationalPatterns - 체세포 SNV 기반 mutational signature
rule mutational_signature:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.somatic.vcf.gz")
    output:
        mutsig = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.mut_sigs.tsv"),
        recon = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.reconstructed_sigs.tsv"),
        occurrences = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.type_occurences.tsv"),
        profile_pdf = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.mutation_profile.pdf")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = BIOMARKER_DIR,
        pname = "{patient}.{tumor_sample_type}",
        max_delta = config["mutsig_max_delta"]
    threads: 4
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "mutsig_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Mutational signature 시작: {params.pname} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        VCF_ABS=$(realpath {input.vcf})
        OUT_ABS=$(realpath {params.out_dir})

        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[somatic_r_tools_container]} \
            /bin/bash -c "
                set -euxo pipefail
                cd $OUT_ABS
                Rscript --vanilla /app/mutational_pattern.R \
                    $VCF_ABS \
                    {params.pname} \
                    {params.max_delta}
            " >> {log} 2>&1

        echo "=== Mutational signature 완료: $(date) ===" >> {log}
        """

# 3. owl - micro-satellite instability (MSI) 프로파일링
rule owl_msi_profile:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam"),
        bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam.bai")
    output:
        profile = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.owl.txt")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = BIOMARKER_DIR,
        pname = "{patient}.{tumor_sample_type}"
    threads: 2
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "owl_profile_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== owl MSI profile 시작: {params.pname} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        BAM_ABS=$(realpath {input.bam})
        OUT_ABS=$(realpath {output.profile})

        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[owl_container]} \
            /bin/bash -c "
                set -euxo pipefail

                # 컨테이너 내장 마커 BED (GRCh38) 사용
                gunzip -c /opt/owl/data/GRCh38_owl_markers.bed.gz > \$TMPDIR/owl_markers.bed 2>/dev/null || \
                    gunzip -c /opt/owl/data/GRCh38_owl_markers.bed.gz > /tmp/owl_markers.bed
                MARKERS=\$TMPDIR/owl_markers.bed
                [ -f \"\$MARKERS\" ] || MARKERS=/tmp/owl_markers.bed

                owl profile --bam $BAM_ABS \
                    --regions \$MARKERS \
                    --sample {params.pname} \
                    > $OUT_ABS
            " >> {log} 2>&1

        echo "=== owl MSI profile 완료: $(date) ===" >> {log}
        """

# 4. owl - MSI 스코어 계산
rule owl_msi_score:
    input:
        profile = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.owl.txt")
    output:
        scores = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.owl-scores.txt"),
        motif_counts = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.owl-motif-counts.txt")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = BIOMARKER_DIR,
        pname = "{patient}.{tumor_sample_type}",
        min_depth = config["msi_min_depth"]
    threads: 2
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "owl_score_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})

        echo "=== owl MSI score 시작: {params.pname} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        PROFILE_ABS=$(realpath {input.profile})
        OUT_ABS=$(realpath {params.out_dir})

        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[owl_container]} \
            /bin/bash -c "
                set -euxo pipefail
                cd $OUT_ABS
                owl score --file $PROFILE_ABS \
                    --prefix {params.pname} \
                    --min-depth {params.min_depth}
            " >> {log} 2>&1

        echo "=== owl MSI score 완료: $(date) ===" >> {log}
        """

# 5. TMB (tumor mutational burden) 추정
#    전체 게놈 기준과 Gencode CDS 영역 기준 두 가지를 산출한다.
rule tmb_estimate:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.{tumor_sample_type}.somatic.vep.vcf.gz"),
        coverage = join(OUTPUT_DIR, "{patient}", "qc", "{patient}.{tumor_sample_type}.regions.bed.gz")
    output:
        tmb_json = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.tmb_estimate.json"),
        tmb_tsv = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.tmb_estimate.tsv"),
        tmb_gencode_json = join(BIOMARKER_DIR, "{patient}.{tumor_sample_type}.tmb_estimate.gencode_coding.json")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = BIOMARKER_DIR,
        pname = "{patient}.{tumor_sample_type}",
        min_depth = config["tmb_min_depth"],
        min_coverage = config["tmb_min_coverage"],
        min_vaf = config["tmb_min_vaf"],
        max_af = config["tmb_gnomad_max_af"]
    threads: 2
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "tmb_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== TMB 추정 시작: {params.pname} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        VCF_ABS=$(realpath {input.vcf})
        COV_ABS=$(realpath {input.coverage})
        OUT_ABS=$(realpath {params.out_dir})

        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[tmb_container]} \
            /bin/bash -c "
                set -euxo pipefail
                cd $OUT_ABS

                # 전체 게놈 기준
                python /opt/venv/bin/calculate_tmb.py \
                    --vcf $VCF_ABS \
                    --coverage $COV_ABS \
                    --min-depth {params.min_depth} \
                    --min-coverage {params.min_coverage} \
                    --min-vaf {params.min_vaf} \
                    --max-af-cutoff {params.max_af} \
                    --output {params.pname}.tmb_estimate.json \
                    --debug-tsv {params.pname}.tmb_estimate.tsv

                # Gencode CDS 영역 기준 (컨테이너 내장 BED)
                python /opt/venv/bin/calculate_tmb.py \
                    --vcf $VCF_ABS \
                    --coverage $COV_ABS \
                    --min-depth {params.min_depth} \
                    --min-coverage {params.min_coverage} \
                    --min-vaf {params.min_vaf} \
                    --max-af-cutoff {params.max_af} \
                    --output {params.pname}.tmb_estimate.gencode_coding.json \
                    --debug-tsv {params.pname}.tmb_estimate.gencode_coding.tsv \
                    --region-bed /opt/gencode_46_coding.bed.gz
            " >> {log} 2>&1

        echo "=== TMB 추정 완료: $(date) ===" >> {log}
        """
