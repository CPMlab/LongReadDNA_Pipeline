# 시각화 및 주석 관련 규칙들

# 1. circosplot 생성 (멀티샘플용)
rule circosplot:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "sv", "recovered_{patient}", "{patient}.recovered.somatic.sv.vcf.gz"),
        circos_bed_file = config["circos_config"],
        cytoband_file = config["cytoband_file"],
        mitelman_file = config["mitelman_mcgene"]
    output:
        circos_dir = directory(join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}")),
        fusion_tsv = join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "{patient}_fusion_calls_combined.tsv")
    params:
        script = "scripts/circosplot.py",
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "circosplot_{patient}.log")
    shell:
        """
        mkdir -p {output.circos_dir}
        
        python {params.script} \
            --input-vcf {input.vcf} \
            --sample-name {wildcards.patient} \
            --mitelman-mcgene {input.mitelman_file} \
            --circos-bed {input.circos_bed_file} \
            --cytoband-file {input.cytoband_file} \
            --output-svg {output.circos_dir}/{wildcards.patient}.sv.circos.svg \
            --output-png {output.circos_dir}/{wildcards.patient}.sv.circos.png \
            --output-fusion-tsv {output.fusion_tsv} \
            > {log} 2>&1 
        """

# 2. 개별 샘플 circos 파일들을 수집하는 보조 rule
rule collect_sample_circos:
    input:
        circos_dir = join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}")
    output:
        sample_files = join(OUTPUT_DIR, "{patient}", "sv", "circos_{patient}", "sample_files.txt")
    shell:
        """
        # 생성된 샘플별 파일들의 목록을 생성
        find {input.circos_dir} -name "*_*.svg" -o -name "*_*.png" -o -name "*_*.tsv" | \
        grep -v "combined" > {output.sample_files} || touch {output.sample_files}
        
        echo "Generated sample-specific circos files:" >> {output.sample_files}
        ls -la {input.circos_dir}/*_*.* 2>/dev/null >> {output.sample_files} || \
        echo "No sample-specific files found" >> {output.sample_files}
        """

# 2. AnnotSV 주석
rule annotsv:
    input:
        vcf_gz = join(OUTPUT_DIR, "{patient}", "sv", "recovered_{patient}", "{patient}.recovered.somatic.sv.vcf.gz"),
        annotsv_annotations_tar = config["annotsv_cache"]
    output:
        tsv = join(OUTPUT_DIR, "{patient}", "sv", "{patient}.sv.annotsv.tsv")
    params:
        tmp_work_dir_prefix = lambda wildcards: f"annotsv_work_{wildcards.patient}",
        tmp_vcf_uncompressed = "tmp.vcf", 
        tmp_vcf_processed = "tmp_processed.vcf",
        annotations_untar_parent_dir = "annotations_from_tar" 
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "annotsv_{patient}.log")
    shell:
        """
        LOG_FILE_PATH_ABS=$(readlink -f {log})
        INPUT_VCF_GZ_ABS=$(readlink -f {input.vcf_gz})
        INPUT_ANNOTATIONS_TAR_ABS=$(readlink -f {input.annotsv_annotations_tar})
        OUTPUT_TSV_ABS=$(readlink -f {output.tsv})

        WORK_DIR=$(mktemp -d -p "." {params.tmp_work_dir_prefix}.XXXXXX)
        trap 'echo "--- 임시 작업 디렉토리 ($WORK_DIR) 삭제 ---" >> "$LOG_FILE_PATH_ABS"; rm -rf "$WORK_DIR"' EXIT
        
        echo "=== AnnotSV 작업 시작 @ $(date) ===" > "$LOG_FILE_PATH_ABS"
        echo "WORK_DIR: $(readlink -f "$WORK_DIR")" >> "$LOG_FILE_PATH_ABS"

        ProcessedVCF="{params.tmp_vcf_processed}"
        UncompressedVCF="{params.tmp_vcf_uncompressed}"

        cd "$WORK_DIR"

        echo "--- 입력 VCF 준비 ---" >> "$LOG_FILE_PATH_ABS"
        gunzip -c "$INPUT_VCF_GZ_ABS" > "$UncompressedVCF" 2>> "$LOG_FILE_PATH_ABS"
        
        awk -F'\\t' -v OFS='\\t' '
        {{ 
            if (NR==2) {{ print "##INFO=<ID=SV_ALT,Number=1,Type=String,Description=\\"Square bracketed notation for BND event\\">" }}
        }}
        {{
            if ($0 ~ /^#/) {{ print $0; }} 
            else {{ if ($8 ~ /SVTYPE=BND/) {{ $8 = $8 ";SV_ALT=" $5; $5 = "<BND>"; }} print $0; }}
        }}' "$UncompressedVCF" > "$ProcessedVCF" 2>> "$LOG_FILE_PATH_ABS"

        AnnotationsUntarParentDirRelative="{params.annotations_untar_parent_dir}"
        mkdir -p "$AnnotationsUntarParentDirRelative"
        echo "--- AnnotSV 주석 데이터 압축 해제 ---" >> "$LOG_FILE_PATH_ABS"
        tar -xzf "$INPUT_ANNOTATIONS_TAR_ABS" -C "$AnnotationsUntarParentDirRelative" 2>> "$LOG_FILE_PATH_ABS"
        
        AnnotationsDirForCmd="$AnnotationsUntarParentDirRelative/AnnotSV"
        TmpOutputTsvInWorkDir="$(basename "$OUTPUT_TSV_ABS")"

        echo "--- AnnotSV 명령어 실행 ---" >> "$LOG_FILE_PATH_ABS"
        
        AnnotSV \
          -SVinputFile "$ProcessedVCF" \
          -annotationsDir "$AnnotationsDirForCmd" \
          -outputFile "$TmpOutputTsvInWorkDir" \
          -outputDir "." \
          -SVinputInfo 1 \
            -genomeBuild GRCh38 \
          2>> "$LOG_FILE_PATH_ABS"
        
        AnnotSV_EXIT_CODE=$?

        if [ $AnnotSV_EXIT_CODE -ne 0 ]; then
            echo "오류: AnnotSV 명령어가 실패했습니다 (종료 코드: $AnnotSV_EXIT_CODE)." >> "$LOG_FILE_PATH_ABS"
            exit 1
        fi

        if [ ! -s "$TmpOutputTsvInWorkDir" ]; then
            echo "오류: AnnotSV가 출력 파일을 생성하지 못했거나 비어있습니다." >> "$LOG_FILE_PATH_ABS"
            exit 1
        fi
        
        mv "$TmpOutputTsvInWorkDir" "$OUTPUT_TSV_ABS" >> "$LOG_FILE_PATH_ABS" 2>&1
        
        echo "--- AnnotSV 완료 ---" >> "$LOG_FILE_PATH_ABS"
        echo "최종 출력 파일: $OUTPUT_TSV_ABS" >> "$LOG_FILE_PATH_ABS"
        """

# 3. SV IntOGen 주석
rule sv_intogen:
    input:
        annotsv_tsv = join(OUTPUT_DIR, "{patient}", "sv", "{patient}.sv.annotsv.tsv"),
        compendium = config["compendium_file"]
    output:
        tsv = join(OUTPUT_DIR, "{patient}", "sv", "{patient}.sv.annotsv_intogenCCG.tsv")
    wildcard_constraints:
        patient="[^.]+"
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "sv", "prioritized_output")
    threads: 1
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "sv_intogen_{patient}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "Processing SV IntOGen for: {input.annotsv_tsv}" > {log}

        TEMP_NOQUOTE_TSV="{params.out_dir}/$(basename {input.annotsv_tsv})_noquote.tsv"
        
        # 1. Remove any quote from the file
        sed 's/"//g' {input.annotsv_tsv} > ${{TEMP_NOQUOTE_TSV}} 2>> {log}

        # 2. csvtk join and summary
        HEADERS_FOR_SUMMARY_SV=$(head -n 1 ${{TEMP_NOQUOTE_TSV}} | tr '\\t' ',' | sed 's/,$//g')

        csvtk join -t \
            ${{TEMP_NOQUOTE_TSV}} \
            {input.compendium} \
            -f "Gene_name;SYMBOL" 2>> {log} | \
        csvtk filter2 -t -f '$Annotation_mode == "split"' 2>> {log} | \
        csvtk summary -t -g "${{HEADERS_FOR_SUMMARY_SV}}" \
            -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" 2>> {log} | \
                sed 's/:collapse//g' > {output.tsv}
                
        rm -f ${{TEMP_NOQUOTE_TSV}}
        
        echo "SV IntOGen annotation complete: {output.tsv}" >> {log}
        """ 