# 매핑 및 QC 관련 규칙들

# 1. 매핑 규칙 (모든 샘플에 대해)
rule mapping:
    input:
        bam_files = lambda wildcards: SAMPLES_DATA[wildcards.patient][wildcards.sample_type],
        ref = REF_FASTA
    output:
        bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{sample_type}.aligned.bam"),
        bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{sample_type}.aligned.bam.bai")
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "mapping_{patient}_{sample_type}.log")
    run:
        import os
        import tempfile
        from datetime import datetime
        
        # 절대 경로로 변환
        log_file = os.path.abspath(log[0])
        output_bam = os.path.abspath(output.bam)
        ref_fasta = os.path.abspath(input.ref)
        
        # 출력 디렉토리 생성
        os.makedirs(os.path.dirname(output_bam), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        
        # 고유한 작업 디렉토리 생성 (절대 경로)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        work_dir = os.path.abspath(os.path.join(os.path.dirname(output_bam), f"mapping_work_{wildcards.patient}_{wildcards.sample_type}_{timestamp}"))
        os.makedirs(work_dir, exist_ok=True)
        
        # 로그 초기화
        shell(f"echo '=== 매핑 시작: {wildcards.patient}.{wildcards.sample_type} ===' > {log_file}")
        shell(f"echo '시작 시간: $(date)' >> {log_file}")
        shell(f"echo '입력 파일들: {input.bam_files}' >> {log_file}")
        shell(f"echo '작업 디렉토리: {work_dir}' >> {log_file}")
        
        try:
            # BAM 파일이 하나만 있는 경우
            if len(input.bam_files) == 1:
                shell(f"""
                    echo "단일 BAM 파일 매핑 시작" >> {log_file}
                    cd {work_dir}
                    pbmm2 align \
                        {ref_fasta} \
                        {input.bam_files[0]} \
                        {output_bam} \
                        --sample {wildcards.patient}.{wildcards.sample_type} \
                        --sort -j {threads} \
                        --preset HIFI --unmapped -A 2 \
                        --log-level INFO --log-file pbmm2.log \
                        >> {log_file} 2>&1
                
                    echo "BAM 인덱싱 시작" >> {log_file}
                    samtools index -@{threads} {output_bam} >> {log_file} 2>&1
                    echo "단일 BAM 파일 매핑 완료" >> {log_file}
                """)
        
            # BAM 파일이 여러 개인 경우
            else:
                shell(f"echo '다중 BAM 파일 매핑 시작 (총 {len(input.bam_files)}개)' >> {log_file}")
                
                temp_bams = []
                for i, bam in enumerate(input.bam_files):
                    # 고유한 임시 파일명 생성 (절대 경로)
                    temp_bam = os.path.join(work_dir, f"{wildcards.patient}_{wildcards.sample_type}_temp_bam{i+1}_{timestamp}.bam")
                    temp_bams.append(temp_bam)
                    
                    shell(f"""
                        echo "BAM 파일 {i+1}/{len(input.bam_files)} 매핑 시작: $(basename {bam})" >> {log_file}
                        cd {work_dir}
                        pbmm2 align \
                            {ref_fasta} \
                            {bam} \
                            {temp_bam} \
                            --sample {wildcards.patient}.{wildcards.sample_type} \
                            --sort -j {threads} \
                            --preset HIFI --unmapped -A 2 \
                            --log-level INFO --log-file pbmm2_{wildcards.patient}_{wildcards.sample_type}_bam{i+1}_{timestamp}.log \
                            >> {log_file} 2>&1
                        echo "BAM 파일 {i+1}/{len(input.bam_files)} 매핑 완료" >> {log_file}
                    """)
            
                # 결과 병합 및 인덱싱
                temp_bams_str = " ".join(temp_bams)
                shell(f"""
                    echo "BAM 파일 병합 시작" >> {log_file}
                    cd {work_dir}
                    samtools merge -@{threads} {output_bam} {temp_bams_str} >> {log_file} 2>&1
                    echo "BAM 인덱싱 시작" >> {log_file}
                    samtools index -@{threads} {output_bam} >> {log_file} 2>&1
                    echo "다중 BAM 파일 매핑 완료" >> {log_file}
                """)
        
        finally:
            # 작업 디렉토리 정리
            shell(f"echo '작업 디렉토리 정리: {work_dir}' >> {log_file}")
            shell(f"rm -rf {work_dir}")
            shell(f"echo '=== 매핑 완료: {wildcards.patient}.{wildcards.sample_type} ===' >> {log_file}")
            shell(f"echo '종료 시간: $(date)' >> {log_file}")

# 2. 커버리지 계산 (모든 샘플에 대해)
rule mosdepth:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{sample_type}.aligned.bam")
    output:
        summary = join(OUTPUT_DIR, "{patient}", "qc", "{patient}.{sample_type}.mosdepth.summary.txt"),
        regions_bed = join(OUTPUT_DIR, "{patient}", "qc", "{patient}.{sample_type}.regions.bed.gz")
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "mosdepth_{patient}_{sample_type}.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.summary})
        
        echo "=== MosDepth 커버리지 분석 시작: {wildcards.patient}.{wildcards.sample_type} ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 BAM: {input.bam}" >> {log}
        echo "출력 prefix: {OUTPUT_DIR}/{wildcards.patient}/qc/{wildcards.patient}.{wildcards.sample_type}" >> {log}
        
        # mosdepth 커버리지 분석
        mosdepth -t {threads} \
            -n --fast-mode \
            --by 500 \
            {OUTPUT_DIR}/{wildcards.patient}/qc/{wildcards.patient}.{wildcards.sample_type} \
            {input.bam} \
            >> {log} 2>&1
        
        echo "=== MosDepth 커버리지 분석 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """ 