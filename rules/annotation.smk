# VEP 주석 관련 규칙들

# 1. VEP 주석 (정상 샘플)
rule vep_annotate_normal:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz"),
        ref = REF_FASTA,
        vep_cache = VEP_CACHE
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.NORMAL.germline.vep.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.NORMAL.germline.vep.vcf.gz.tbi")
    params:
        tmp_dir = "vep_data_{patient}_NORMAL"
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "vep_{patient}_NORMAL.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.vcf})
        
        echo "=== VEP 주석 시작: {wildcards.patient}.NORMAL ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 VCF: {input.vcf}" >> {log}
        echo "출력 VCF: {output.vcf}" >> {log}
        echo "임시 디렉토리: {params.tmp_dir}" >> {log}
        
        # VEP 캐시 압축 해제
        echo "VEP 캐시 압축 해제 중..." >> {log}
        tar -xzf {input.vep_cache} -C {params.tmp_dir} >> {log} 2>&1
        
        # VEP 주석 추가
        echo "VEP 주석 처리 중..." >> {log}
        vep \
          --cache \
          --cache_version 112 \
          --offline \
          --dir {params.tmp_dir} \
          --fasta {input.ref} \
          --format vcf \
          --fork {threads} \
          --species homo_sapiens \
          --assembly GRCh38 \
          --symbol \
          --hgvs \
          --refseq \
          --check_existing \
          --vcf \
          --pick \
          --flag_pick_allele_gene \
          --everything \
          --compress_output bgzip \
          -i {input.vcf} \
          -o {output.vcf} \
          >> {log} 2>&1
          
        # 인덱스 생성
        echo "VCF 인덱싱 중..." >> {log}
        tabix -p vcf {output.vcf} >> {log} 2>&1
        
        # 임시 파일 정리
        echo "임시 파일 정리 중..." >> {log}
        rm -rf {params.tmp_dir} >> {log} 2>&1
        
        echo "=== VEP 주석 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """

# 2. VEP 주석 (체세포 변이)
rule vep_annotate_somatic:
    input:
        vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.somatic.normalized.vcf.gz"),
        ref = REF_FASTA,
        vep_cache = VEP_CACHE
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.{tumor_sample_type}.somatic.vep.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "annotation", "{patient}.{tumor_sample_type}.somatic.vep.vcf.gz.tbi")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        tmp_dir = "vep_data_{patient}_{tumor_sample_type}_somatic"
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "vep_{patient}_{tumor_sample_type}_somatic.log")
    shell:
        """
        # 로그 및 출력 디렉토리 생성
        mkdir -p $(dirname {log})
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.vcf})
        
        echo "=== VEP 체세포 주석 시작: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "시작 시간: $(date)" >> {log}
        echo "입력 VCF: {input.vcf}" >> {log}
        echo "출력 VCF: {output.vcf}" >> {log}
        echo "임시 디렉토리: {params.tmp_dir}" >> {log}
        
        # VEP 캐시 압축 해제
        echo "VEP 캐시 압축 해제 중..." >> {log}
        tar -xzf {input.vep_cache} -C {params.tmp_dir} >> {log} 2>&1
        
        # VEP 주석 추가
        echo "VEP 주석 처리 중..." >> {log}
        vep \
          --cache \
          --cache_version 112 \
          --offline \
          --dir {params.tmp_dir} \
          --fasta {input.ref} \
          --format vcf \
          --fork {threads} \
          --species homo_sapiens \
          --assembly GRCh38 \
          --symbol \
          --hgvs \
          --refseq \
          --check_existing \
          --vcf \
          --pick \
          --flag_pick_allele_gene \
          --everything \
          --compress_output bgzip \
          -i {input.vcf} \
          -o {output.vcf} \
          >> {log} 2>&1
          
        # 인덱스 생성
        echo "VCF 인덱싱 중..." >> {log}
        tabix -p vcf {output.vcf} >> {log} 2>&1
        
        # 임시 파일 정리
        echo "임시 파일 정리 중..." >> {log}
        rm -rf {params.tmp_dir} >> {log} 2>&1
        
        echo "=== VEP 체세포 주석 완료 ===" >> {log}
        echo "종료 시간: $(date)" >> {log}
        """ 