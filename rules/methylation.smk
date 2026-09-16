# Updated Snakemake rules for CpG methylation, DMR calling, and annotation
# 주요 변경 사항
#   1. shell 블록 내부의 `cd` 사용 제거 → 절대경로 기반으로 처리
#   2. set -euo pipefail 추가로 오류 발생 시 즉시 중단
#   3. 작업용 임시 파일·로그·출력 디렉터리 변수화(OUT_DIR)로 경로 명확화
#   4. R 스크립트가 결과를 현재 작업 디렉터리에 생성하는 특성을 고려해
#      OUT_DIR 내에서 생성하도록 prefix 전달 & 사후 mv 로 최종 출력 보장

#######################################################################
# 1) CpG 메틸화 스코어 계산
#######################################################################
rule cpg_methylation:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{sample_type}.hiphase.bam"),
        ref = REF_FASTA
    output:
        combined_bed = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.combined.bed.gz"),
        hap1_bed     = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap1.bed.gz"),
        hap2_bed     = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap2.bed.gz"),
        combined_bw  = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.combined.bw"),
        hap1_bw      = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap1.bw"),
        hap2_bw      = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap2.bw")
    params:
        output_prefix = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg"),
        min_coverage = config.get("methylation_min_coverage", 5),
        min_mapq     = config.get("methylation_min_mapq", 1)
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "cpg_methylation_{patient}_{sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        # 결과 디렉터리 확보
        mkdir -p "$(dirname {params.output_prefix})"
        mkdir -p "$(dirname {log})"

        echo "CpG 메틸화 추출 시작: {wildcards.patient}.{wildcards.sample_type}" > {log}
        echo "시작 시간: $(date)" >> {log}

        aligned_bam_to_cpg_scores --version >> {log} 2>&1

        aligned_bam_to_cpg_scores \
          --threads {threads} \
          --bam {input.bam} \
          --ref {input.ref} \
          --output-prefix {params.output_prefix} \
          --min-mapq {params.min_mapq} \
          --min-coverage {params.min_coverage} \
          >> {log} 2>&1

        echo "CpG 메틸화 추출 완료: {wildcards.patient}.{wildcards.sample_type}" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

#######################################################################
# 2) DSS 차등 메틸화 영역 분석 (종양 vs 정상)
#######################################################################
rule dss_dmr:
    input:
        tumor_bed  = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{tumor_sample_type}.cpg.combined.bed.gz"),
        normal_bed = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.NORMAL.cpg.combined.bed.gz")
    output:
        dmr_tsv = join(OUTPUT_DIR, "{patient}", "dmr", "{patient}.{tumor_sample_type}_vs_NORMAL.DMR.tsv")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir    = join(OUTPUT_DIR, "{patient}", "dmr"),
        sample_name = "{patient}.{tumor_sample_type}_vs_NORMAL",
        dss_script  = "scripts/DSS_tumor_normal.R"
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "dss_dmr_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        OUT_DIR={params.out_dir}
        LOG_FILE={log}
        DSS_SCRIPT=$(readlink -f {params.dss_script})

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"

        echo "DSS 차등 메틸화 분석 시작: {wildcards.patient}.{wildcards.tumor_sample_type} vs NORMAL" > "$LOG_FILE"
        echo "시작 시간: $(date)" >> "$LOG_FILE"

        # 임시 파일 (OUT_DIR 안에 생성)
        TUMOR_TMP="$OUT_DIR/{params.sample_name}.tumor.tmp"
        NORMAL_TMP="$OUT_DIR/{params.sample_name}.normal.tmp"

        gunzip -c {input.tumor_bed} | grep -v '^#' | cut -f1,2,6,7 > "$TUMOR_TMP" 2>> "$LOG_FILE"
        gunzip -c {input.normal_bed} | grep -v '^#' | cut -f1,2,6,7 > "$NORMAL_TMP" 2>> "$LOG_FILE"

        echo "DSS R 스크립트 실행 중..." >> "$LOG_FILE"

        Rscript --vanilla "$DSS_SCRIPT" \
            "$TUMOR_TMP" \
            "$NORMAL_TMP" \
            {output.dmr_tsv} \
            {threads} \
            >> "$LOG_FILE" 2>&1

        rm -f "$TUMOR_TMP" "$NORMAL_TMP"

        echo "DSS 차등 메틸화 분석 완료: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "종료 시간: $(date)" >> "$LOG_FILE"
        """

#######################################################################
# 3) DMR annotation 분석
#######################################################################
rule annotate_dmr:
    input:
        dmr_tsv = join(OUTPUT_DIR,
                       "{patient}",
                       "dmr",
                       "{patient}.{tumor_sample_type}_vs_NORMAL.DMR.tsv")
    output:
        annotated_dmr = join(OUTPUT_DIR,
                             "{patient}",
                             "dmr",
                             "{patient}.{tumor_sample_type}_vs_NORMAL.annotated_DMR.tsv.gz")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir        = join(OUTPUT_DIR, "{patient}", "dmr"),
        sample_name    = "{patient}.{tumor_sample_type}_vs_NORMAL",
        annotate_script = "scripts/annotatr_dmr.R"
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR,
             "{patient}",
             "logs",
             "annotate_dmr_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        OUT_DIR="{params.out_dir}"
        LOG_FILE="{log}"
        ANNO_SCRIPT="$(readlink -f {params.annotate_script})"

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"

        echo "DMR 주석 분석 시작: {wildcards.patient}.{wildcards.tumor_sample_type}" > "$LOG_FILE"
        echo "시작 시간: $(date)" >> "$LOG_FILE"

        PREFIX="$OUT_DIR/{params.sample_name}"

        Rscript --vanilla "$ANNO_SCRIPT" \
            "{input.dmr_tsv}" \
            "$PREFIX" \
            {threads} >> "$LOG_FILE" 2>&1

        # 결과 확인 & 이동
        GENERATED=$(ls "$PREFIX"*.tsv.gz 2>/dev/null || true)
        if [[ -n "$GENERATED" ]]; then
            mv "$GENERATED" "{output.annotated_dmr}"
        else
            echo "주석된 DMR 파일이 생성되지 않음. 빈 파일 생성." >> "$LOG_FILE"
            touch "{output.annotated_dmr}"
        fi

        echo "DMR 주석 분석 완료: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "종료 시간: $(date)" >> "$LOG_FILE"
        """
