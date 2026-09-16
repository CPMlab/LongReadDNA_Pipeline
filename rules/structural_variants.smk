# 구조적 변이 검출 및 처리 관련 규칙들

# 1. Severus를 이용한 구조적 변이 검출 (Multimode)
rule severus_multimode:
    input:
        normal_bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.hiphase.bam"),
        tumor_bams = lambda wildcards: [join(OUTPUT_DIR, wildcards.patient, "phasing", f"{wildcards.patient}.{t}.hiphase.bam") for t in get_tumor_samples(wildcards.patient)],
        phasing_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz"),
        vntr_bed = VNTR_BED
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "sv", "severus_{patient}", "somatic_SVs", "severus_somatic.vcf")
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "sv", "severus_{patient}"),
        min_supp = 3
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "severus_{patient}_multimode.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "Severus multimode 실행 시작" > {log}
        echo "Control BAM: {input.normal_bam}" >> {log}
        echo "Target BAMs: {input.tumor_bams}" >> {log}
        
        # Severus multimode 실행 - 여러 target-bam 지원
        severus \
            --control-bam {input.normal_bam} \
            --target-bam {input.tumor_bams} \
            --phasing-vcf {input.phasing_vcf} \
            --vntr-bed {input.vntr_bed} \
            --out-dir {params.out_dir} \
            -t {threads} \
            --min-support {params.min_supp} \
            --resolve-overlaps \
            --between-junction-ins \
            --single-bp \
            >> {log} 2>&1
            
        echo "Severus multimode 완료" >> {log}
        """

# 2. tabix 필터링
rule tabix_filter:
    input:
        uncompressed_vcf = join(OUTPUT_DIR, "{patient}", "sv", "severus_{patient}", "somatic_SVs", "severus_somatic.vcf"),
        contig_bed = config["contig_bed"]
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "sv", "tabix_{patient}", "{patient}.tabix.somatic.sv.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "sv", "tabix_{patient}", "{patient}.tabix.somatic.sv.vcf.gz.tbi")
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "tabix_filter_{patient}_somatic.log")
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        
        echo "VCF 파일 인덱싱 및 필터링 시작: {input.uncompressed_vcf}" > {log}
        
        bcftools sort -Oz -o tmp.{wildcards.patient}.sorted.vcf.gz "{input.uncompressed_vcf}" >> {log} 2>&1
        
        tabix -p vcf tmp.{wildcards.patient}.sorted.vcf.gz >> {log} 2>&1
        
        bcftools view -R "{input.contig_bed}" -Oz -o "{output.vcf}" tmp.{wildcards.patient}.sorted.vcf.gz >> {log} 2>&1
        
        tabix -p vcf "{output.vcf}" >> {log} 2>&1
        
        rm -f tmp.{wildcards.patient}.sorted.vcf.gz tmp.{wildcards.patient}.sorted.vcf.gz.tbi >> {log} 2>&1
        
        echo "완료: VCF 파일 처리가 완료되었습니다. 결과 파일: {output.vcf}, {output.tbi}" >> {log}
        """

# 3. SVpack 필터링 및 주석 
rule svpack:
    input:
        filtered_sv_vcf = join(OUTPUT_DIR, "{patient}", "sv", "tabix_{patient}", "{patient}.tabix.somatic.sv.vcf.gz"),
        match_vcf = config["svpack_match_vcf"],
        reference_gff = config["reference_gff"]
    output:
        svpack_vcf_gz = join(OUTPUT_DIR, "{patient}", "sv", "svpack_{patient}", "{patient}.tabix.somatic.svpack.sv.vcf.gz"),
        svpack_tbi = join(OUTPUT_DIR, "{patient}", "sv", "svpack_{patient}", "{patient}.tabix.somatic.svpack.sv.vcf.gz.tbi")
    params:
        svpack_executable = config["svpack_script"],
        out_dir = join(OUTPUT_DIR, "{patient}", "sv", "svpack_{patient}"),
        temp_uncompressed_svpack_vcf = join(OUTPUT_DIR, "{patient}", "sv", "svpack_{patient}", "{patient}.tabix.somatic.svpack.intermediate.sv.vcf")
    threads: 1
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "svpack_{patient}_somatic.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "SVPACK 필터링 및 주석 시작: {input.filtered_sv_vcf}" > {log}
        
        # svpack 파이프라인 실행 결과를 임시 압축되지 않은 파일에 저장
        python {params.svpack_executable} filter --pass-only "{input.filtered_sv_vcf}" 2>> {log} | \
        python {params.svpack_executable} filter --min-svlen 50 - 2>> {log} | \
        python {params.svpack_executable} match -v - "{input.match_vcf}" 2>> {log} | \
        python {params.svpack_executable} consequence - "{input.reference_gff}" 2>> {log} | \
        python {params.svpack_executable} tagzygosity - 2>> {log} \
        > {params.temp_uncompressed_svpack_vcf}
        
        # 생성된 임시 파일이 비어있지 않은지 확인
        if [ ! -s "{params.temp_uncompressed_svpack_vcf}" ]; then
            echo "Error: svpack 파이프라인에서 {params.temp_uncompressed_svpack_vcf} 파일 생성 실패 또는 비어있음" >> {log}
            exit 1
        fi

        echo "임시 VCF 압축 및 인덱싱 중..." >> {log}
        # 임시 파일을 bgzip으로 압축하여 최종 출력 파일(.gz) 생성
        bgzip -f -c {params.temp_uncompressed_svpack_vcf} > {output.svpack_vcf_gz} 2>> {log}
        # 최종 압축 파일 인덱싱
        tabix -f -p vcf {output.svpack_vcf_gz} 2>> {log}
        
        # 임시 압축되지 않은 파일 삭제
        rm -f {params.temp_uncompressed_svpack_vcf}
        
        echo "SVPACK 완료: {output.svpack_vcf_gz}" >> {log}
        """ 

# 4. recover_mate_bnd 실행
rule recover_mate_bnd:
    input:
        original_vcf = join(OUTPUT_DIR, "{patient}", "sv", "tabix_{patient}", "{patient}.tabix.somatic.sv.vcf.gz"),
        filtered_vcf_gz = join(OUTPUT_DIR, "{patient}", "sv", "svpack_{patient}", "{patient}.tabix.somatic.svpack.sv.vcf.gz")
    output:
        vcf = join(OUTPUT_DIR, "{patient}", "sv", "recovered_{patient}", "{patient}.recovered.somatic.sv.vcf.gz"),
        tbi = join(OUTPUT_DIR, "{patient}", "sv", "recovered_{patient}", "{patient}.recovered.somatic.sv.vcf.gz.tbi")
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "sv", "recovered_{patient}"),
        missing_mate_list = lambda wildcards: f"missing_mate_{wildcards.patient}.txt",
        missing_mate_vcf = lambda wildcards: f"missing_mate_{wildcards.patient}.tmp.vcf.gz",
        tmp_ids = lambda wildcards: f"tmp_ids_{wildcards.patient}.txt",
        tmp_mate_ids = lambda wildcards: f"tmp_mate_ids_{wildcards.patient}.txt"
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "recover_mate_bnd_{patient}_somatic.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "=== recover_mate_bnd 시작 ===" > {log}
        echo "실행 디렉토리: $(pwd)" >> {log}
        echo "입력 original_vcf: {input.original_vcf}" >> {log}
        echo "입력 filtered_vcf_gz: {input.filtered_vcf_gz}" >> {log}

        echo "누락된 짝 ID 찾는 중..." >> {log}
        
        bcftools query -f '%ID\\t%MATE_ID\\n' {input.filtered_vcf_gz} 2>> {log} | grep -v '\\.' 2>> {log} | cut -f1 | sort > {params.tmp_ids}
        bcftools query -f '%ID\\t%MATE_ID\\n' {input.filtered_vcf_gz} 2>> {log} | grep -v '\\.' 2>> {log} | cut -f2 | sort > {params.tmp_mate_ids}

        comm -13 {params.tmp_ids} {params.tmp_mate_ids} > {params.missing_mate_list} 2>> {log}

        if [ -s "{params.missing_mate_list}" ]; then
            echo "누락된 짝 ID가 {params.missing_mate_list} 에 발견됨. 추출 및 병합 진행." >> {log}
            
            bcftools view -i ID=@$(basename {params.missing_mate_list}) "{input.original_vcf}" 2>> {log} | \
                bcftools sort -Oz -o {params.missing_mate_vcf} >> {log} 2>&1
            
            if [ ! -s "{params.missing_mate_vcf}" ]; then
                echo "오류: bcftools view 또는 sort 이후 MISSING_MATE_VCF 파일이 비어있거나 생성되지 않았습니다." >> {log}
                exit 1
            fi
            tabix -f -p vcf {params.missing_mate_vcf} >> {log} 2>&1
            
            echo "파일 병합 및 정렬 중..." >> {log}
            bcftools concat -a "{input.filtered_vcf_gz}" {params.missing_mate_vcf} 2>> {log} | \
                bcftools sort -Oz -o "{output.vcf}" >> {log} 2>&1
        else
            echo "누락된 짝 ID가 없습니다. 원본 필터링된 VCF를 최종 결과로 사용합니다." >> {log}
            cp "{input.filtered_vcf_gz}" "{output.vcf}" >> {log} 2>&1
        fi
        
        tabix -f -p vcf "{output.vcf}" >> {log} 2>&1
        
        # 임시 파일 최종 정리
        rm -f {params.tmp_ids} {params.tmp_mate_ids} {params.missing_mate_list} {params.missing_mate_vcf} {params.missing_mate_vcf}.tbi >> {log} 2>&1
        
        echo "recover_mate_bnd 완료: {output.vcf}" >> {log}
        """ 