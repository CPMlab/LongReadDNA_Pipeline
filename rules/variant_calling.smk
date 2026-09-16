# 변이 검출 및 위상결정 관련 규칙들

# 1. Clair3 변이 검출 (모든 샘플에 대해)
rule clair3:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{sample_type}.aligned.bam"),
        ref = REF_FASTA
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{sample_type}.clair3.small_variants.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{sample_type}.clair3.small_variants.vcf.gz.tbi")
    params:
        sample = "{patient}.{sample_type}",
        out_dir = join(OUTPUT_DIR, "{patient}", "snv", "clair3_{patient}_{sample_type}")
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "clair3_{patient}_{sample_type}.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Clair3 변이 검출 시작: {wildcards.patient}.{wildcards.sample_type} ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 BAM: {input.bam}" >> {log}
        echo "출력 디렉토리: {params.out_dir}" >> {log}

        # 임시 파일 이름 정의
        TEMP_VCF=$(mktemp {params.out_dir}/temp_clair3_{wildcards.sample_type}_XXXXXX.vcf)
        echo "임시 VCF: $TEMP_VCF" >> {log}

        # Clair3 변이 검출 (Singularity 내부에서 실행)
        echo "Clair3 Singularity 컨테이너 실행 중..." >> {log}
        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[clair3_container]} \
            /bin/bash -c "
                set -euxo pipefail

                # Clair3 실행
                /opt/bin/run_clair3.sh \
                    --bam_fn={input.bam} \
                    --ref_fn={input.ref} \
                    --threads={threads} \
                    --platform='hifi' \
                    --model_path='/opt/models/hifi_revio' \
                    --output={params.out_dir} \
                    --sample_name={params.sample}

                # LowQual 필터 제외하여 임시 VCF 파일 생성
                gunzip -c {params.out_dir}/merge_output.vcf.gz | \
                    awk 'BEGIN{{OFS=\"\\t\"}} {{if(/^#/ || \\$7 != \"LowQual\") {{print \\$0}}}}' > ${{TEMP_VCF}}
            " >> {log} 2>&1

        echo "Clair3 실행 완료, 후처리 중..." >> {log}

        # 생성된 임시 VCF 파일을 bgzip으로 압축
        bgzip -f ${{TEMP_VCF}} 2>> {log}
        mv ${{TEMP_VCF}}.gz {output.vcf} 2>> {log}

        # 최종 VCF 파인덱싱
        echo "VCF 인덱싱 중..." >> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        
        echo "=== Clair3 변이 검출 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

# 2. 체세포 변이 검출 (각 tumor sample에 대해 normal과 비교)
rule deepsomatic:
    input:
        tumor_bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam"),
        normal_bam = lambda wildcards: join(OUTPUT_DIR, wildcards.patient, "mapping", f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}.aligned.bam"),
        ref = REF_FASTA
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.somatic.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.somatic.vcf.gz.tbi")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "snv", "deepsomatic_{patient}_{tumor_sample_type}"),
        tumor_sample_name = "{patient}.{tumor_sample_type}",
        normal_sample_name = lambda wildcards: f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}"
    threads: THREADS
    resources:
        deepsomatic_slots=1  # DeepSomatic 작업 동시 실행 방지
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "deepsomatic_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "DeepSomatic 체세포 변이 검출 시작: {wildcards.patient}.{wildcards.tumor_sample_type}" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 종양 BAM: {input.tumor_bam}" >> {log}
        echo "입력 정상 BAM: {input.normal_bam}" >> {log}
        echo "참조 게놈: {input.ref}" >> {log}
        
        # DeepSomatic 실행 (기본 설정)
        singularity run \
            -B $(pwd):$(pwd) \
            {config[deepsomatic_container]} \
            run_deepsomatic \
            --model_type=PACBIO \
            --ref={input.ref} \
            --reads_tumor={input.tumor_bam} \
            --reads_normal={input.normal_bam} \
            --output_vcf={params.out_dir}/somatic.vcf.gz \
            --sample_name_tumor={params.tumor_sample_name} \
            --sample_name_normal={params.normal_sample_name} \
            --num_shards={threads} \
            --vcf_stats_report=true \
            --logging_dir={params.out_dir}/logs \
            2>&1 | tee -a {log}
        
        echo "DeepSomatic 컨테이너 실행 완료" >> {log}
        echo "후처리 시작 시간: $(date)" >> {log}
            
        # PASS 변이만 필터링
        bcftools view \
            -f PASS -Oz \
            -o {params.out_dir}/somatic_PASS.vcf.gz \
            {params.out_dir}/somatic.vcf.gz \
            2>&1 | tee -a {log}
            
        # 호모 변이(1/1)를 헤테로 변이(0/1)로 변경
        bcftools +setGT {params.out_dir}/somatic_PASS.vcf.gz -- -t q -i 'GT="1/1"' -n c:"0/1" | \
            bcftools sort -Oz -o {output.vcf} \
            2>&1 | tee -a {log}
            
        # 인덱스 생성
        tabix -p vcf {output.vcf} 2>&1 | tee -a {log}
        
        echo "DeepSomatic 완료: {wildcards.patient}.{wildcards.tumor_sample_type}" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """ 

# 3. 정상 샘플 위상 결정
rule hiphase_normal:
    input:
        bam = lambda wildcards: join(OUTPUT_DIR, wildcards.patient, "mapping", f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}.aligned.bam"),
        vcf = lambda wildcards: join(OUTPUT_DIR, wildcards.patient, "snv", f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}.clair3.small_variants.vcf.gz"),
        ref = REF_FASTA
    output:
        bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.bam"),
        bam_index = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.bam.bai"),
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.vcf.gz.tbi"),
        stats = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.stats"),
        summary = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.summary.tsv")
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "hiphase_{patient}_NORMAL_germline.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        
        echo "=== HiPhase 위상결정 시작: {wildcards.patient}.NORMAL ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 BAM: {input.bam}" >> {log}
        echo "입력 VCF: {input.vcf}" >> {log}
        echo "출력 BAM: {output.bam}" >> {log}
        
        # HiPhase 위상결정
        hiphase --bam {input.bam} \
            -t {threads} \
            --output-bam {output.bam} \
            --vcf {input.vcf} \
            --output-vcf {output.vcf} \
            -r {input.ref} \
            --stats-file {output.stats} \
            --summary-file {output.summary} \
            --ignore-read-groups \
            >> {log} 2>&1
            
        echo "BAM 인덱싱 중..." >> {log}
        # BAM 인덱스 생성
        samtools index -@{threads} {output.bam} >> {log} 2>&1
        
        echo "=== HiPhase 위상결정 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

# 4. 종양 샘플 위상 결정 (체세포 변이 포함)
rule hiphase_tumor:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam"),
        germline_vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.clair3.small_variants.vcf.gz"),
        somatic_vcf = join(OUTPUT_DIR, "{patient}", "snv", "{patient}.{tumor_sample_type}.somatic.vcf.gz"),
        ref = REF_FASTA
    output:
        bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam"),
        bam_index = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam.bai"),
        germline_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.germline_like.hiphase.vcf.gz"),
        somatic_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.hiphase.vcf.gz"),
        stats = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.hiphase.stats"),
        summary = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.hiphase.summary.tsv")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "hiphase_{patient}_{tumor_sample_type}_somatic.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        
        echo "=== HiPhase 종양 위상결정 시작: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 BAM: {input.bam}" >> {log}
        echo "입력 생식세포 VCF: {input.germline_vcf}" >> {log}
        echo "입력 체세포 VCF: {input.somatic_vcf}" >> {log}
        echo "출력 BAM: {output.bam}" >> {log}
        
        # 종양 샘플 위상 결정 (체세포 변이 포함)
        hiphase --bam {input.bam} \
            -t {threads} \
            --output-bam {output.bam} \
            --vcf {input.germline_vcf} \
            --output-vcf {output.germline_vcf} \
            --vcf {input.somatic_vcf} \
            --output-vcf {output.somatic_vcf} \
            -r {input.ref} \
            --stats-file {output.stats} \
            --summary-file {output.summary} \
            --ignore-read-groups \
            >> {log} 2>&1
            
        echo "BAM 인덱싱 중..." >> {log}
        # BAM 인덱스 생성
        samtools index -@{threads} {output.bam} >> {log} 2>&1
        
        echo "=== HiPhase 종양 위상결정 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

# 5. 정상 샘플 위상 결정된 VCF 정규화
rule normalize_normal_vcf:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.vcf.gz"),
        ref = REF_FASTA
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz.tbi")
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "normalize_vcf_{patient}_NORMAL.log")
    shell:
        """
        # 로그 디렉토리 생성
        mkdir -p $(dirname {log})
        
        echo "=== VCF 정규화 시작: {wildcards.patient}.NORMAL ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 VCF: {input.vcf}" >> {log}
        echo "출력 VCF: {output.vcf}" >> {log}
        
        # VCF 인덱싱
        echo "VCF 인덱싱 중..." >> {log}
        bcftools index --threads {threads} {input.vcf} >> {log} 2>&1
        
        # 정규화 실행
        echo "VCF 정규화 중..." >> {log}
        bcftools view {input.vcf} | \
          sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \
          bcftools norm --threads {threads} --multiallelics - \
            --output-type b --fasta-ref {input.ref} | \
          bcftools sort -Oz -o {output.vcf} \
          2>> {log}
        
        # 정규화된 VCF 인덱스 생성
        echo "정규화된 VCF 인덱싱 중..." >> {log}
        bcftools index --threads {threads} -t {output.vcf} >> {log} 2>&1
        
        echo "=== VCF 정규화 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

# 6. 종양 샘플 위상 결정된 VCF 정규화
rule normalize_tumor_vcf:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.hiphase.vcf.gz"),
        ref = REF_FASTA
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.normalized.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.normalized.vcf.gz.tbi")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "normalize_vcf_{patient}_{tumor_sample_type}_somatic.log")
    shell:
        """
        # 로그 디렉토리 생성
        mkdir -p $(dirname {log})
        
        echo "=== 종양 VCF 정규화 시작: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 VCF: {input.vcf}" >> {log}
        echo "출력 VCF: {output.vcf}" >> {log}
        
        # VCF 인덱싱
        echo "VCF 인덱싱 중..." >> {log}
        bcftools index --threads {threads} {input.vcf} >> {log} 2>&1
        
        # 정규화 실행
        echo "VCF 정규화 중..." >> {log}
        bcftools view {input.vcf} | \
          sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \
          bcftools norm --threads {threads} --multiallelics - \
            --output-type b --fasta-ref {input.ref} | \
          bcftools sort -Oz -o {output.vcf} \
          2>> {log}
        
        # 정규화된 VCF 인덱스 생성
        echo "정규화된 VCF 인덱싱 중..." >> {log}
        bcftools index --threads {threads} -t {output.vcf} >> {log} 2>&1
        
        echo "=== 종양 VCF 정규화 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """ 