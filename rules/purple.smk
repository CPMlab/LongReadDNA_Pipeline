# PURPLE (HMFtools) 기반 purity/ploidy 및 allele-specific CNV
# Amber (BAF) -> Cobalt (read ratio) -> Purple (purity/ploidy/CNV) 순서로 실행한다.
# HiFi-somatic-WDL 의 clonality.wdl 과 동일한 호출 방식이며, 컨테이너 대신
# 로컬 jar (config 의 amber_jar / cobalt_jar / purple_jar) 를 직접 사용한다.

PURPLE_DIR = join(OUTPUT_DIR, "{patient}", "cnv", "purple_{patient}_{tumor_sample_type}")

# 0. HMF 참조 리소스 tarball 압축 해제 (여러 샘플이 공유하므로 한 번만)
rule hmf_resources:
    input:
        tarball = config["hmf_resources_tarball"]
    output:
        res_dir = directory(join(OUTPUT_DIR, "shared", "hmf_resources"))
    log:
        join(OUTPUT_DIR, "shared", "logs", "hmf_resources.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {output.res_dir}

        echo "=== HMF 참조 리소스 압축 해제 ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        tar -xzf {input.tarball} -C {output.res_dir} 2>> {log}

        # 필수 파일 확인
        AMBER_LOCI=$(find {output.res_dir} -name "AmberGermlineSites.*.tsv.gz" | head -1)
        GC_PROFILE=$(find {output.res_dir} -name "GC_profile.*.cnp" | head -1)
        ENSEMBL_DATA=$(find {output.res_dir} -type d -name "ensembl_data" | head -1)

        if [ -z "$AMBER_LOCI" ] || [ -z "$GC_PROFILE" ] || [ -z "$ENSEMBL_DATA" ]; then
            echo "❌ HMF 리소스에서 필수 파일을 찾지 못했습니다." >> {log}
            echo "  AmberGermlineSites: $AMBER_LOCI" >> {log}
            echo "  GC_profile: $GC_PROFILE" >> {log}
            echo "  ensembl_data: $ENSEMBL_DATA" >> {log}
            exit 1
        fi

        echo "AmberGermlineSites: $AMBER_LOCI" >> {log}
        echo "GC_profile: $GC_PROFILE" >> {log}
        echo "ensembl_data: $ENSEMBL_DATA" >> {log}
        echo "=== 완료: $(date) ===" >> {log}
        """

# 1. Amber - B-allele frequency 및 contamination 추정
rule purple_amber:
    input:
        normal_bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.NORMAL.aligned.bam"),
        normal_bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.NORMAL.aligned.bam.bai"),
        tumor_bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam"),
        tumor_bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam.bai"),
        ref = REF_FASTA,
        res_dir = join(OUTPUT_DIR, "shared", "hmf_resources")
    output:
        baf_pcf = join(PURPLE_DIR, "amber", "{patient}.{tumor_sample_type}.amber.baf.pcf")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = join(PURPLE_DIR, "amber"),
        normal_name = "{patient}.NORMAL",
        tumor_name = "{patient}.{tumor_sample_type}",
        jar = config["amber_jar"],
        java_mem = config["purple_java_mem"]
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "amber_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Amber 시작: {params.tumor_name} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        AMBER_LOCI=$(find {input.res_dir} -name "AmberGermlineSites.*.tsv.gz" | head -1)
        echo "loci: $AMBER_LOCI" >> {log}

        java -Xmx{params.java_mem} -jar {params.jar} \
            -reference {params.normal_name} \
            -reference_bam {input.normal_bam} \
            -tumor {params.tumor_name} \
            -tumor_bam {input.tumor_bam} \
            -output_dir {params.out_dir} \
            -threads {threads} \
            -ref_genome {input.ref} \
            -ref_genome_version V38 \
            -loci "$AMBER_LOCI" >> {log} 2>&1

        echo "=== Amber 완료: $(date) ===" >> {log}
        """

# 2. Cobalt - read depth ratio 및 GC 보정
rule purple_cobalt:
    input:
        normal_bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.NORMAL.aligned.bam"),
        normal_bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.NORMAL.aligned.bam.bai"),
        tumor_bam = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam"),
        tumor_bai = join(OUTPUT_DIR, "{patient}", "mapping", "{patient}.{tumor_sample_type}.aligned.bam.bai"),
        ref = REF_FASTA,
        res_dir = join(OUTPUT_DIR, "shared", "hmf_resources")
    output:
        ratio_pcf = join(PURPLE_DIR, "cobalt", "{patient}.{tumor_sample_type}.cobalt.ratio.pcf")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = join(PURPLE_DIR, "cobalt"),
        normal_name = "{patient}.NORMAL",
        tumor_name = "{patient}.{tumor_sample_type}",
        jar = config["cobalt_jar"],
        java_mem = config["purple_java_mem"],
        pcf_gamma = config["cobalt_pcf_gamma"]
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "cobalt_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Cobalt 시작: {params.tumor_name} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        GC_PROFILE=$(find {input.res_dir} -name "GC_profile.*.cnp" | head -1)
        echo "gc_profile: $GC_PROFILE" >> {log}

        java -Xmx{params.java_mem} -jar {params.jar} \
            -reference {params.normal_name} \
            -reference_bam {input.normal_bam} \
            -tumor {params.tumor_name} \
            -tumor_bam {input.tumor_bam} \
            -ref_genome {input.ref} \
            -output_dir {params.out_dir} \
            -threads {threads} \
            -pcf_gamma {params.pcf_gamma} \
            -validation_stringency SILENT \
            -gc_profile "$GC_PROFILE" >> {log} 2>&1

        echo "=== Cobalt 완료: $(date) ===" >> {log}
        """

# 3. Purple - purity / ploidy / allele-specific CNV
#    DeepSomatic VCF 는 PURPLE 이 기대하는 tumor+normal 2-sample 형식이 아니므로
#    somatic VCF 없이 Amber/Cobalt 결과만으로 적합한다 (기존 수동 실행과 동일한 방식).
rule purple:
    input:
        baf_pcf = join(PURPLE_DIR, "amber", "{patient}.{tumor_sample_type}.amber.baf.pcf"),
        ratio_pcf = join(PURPLE_DIR, "cobalt", "{patient}.{tumor_sample_type}.cobalt.ratio.pcf"),
        germline_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz"),
        ref = REF_FASTA,
        res_dir = join(OUTPUT_DIR, "shared", "hmf_resources")
    output:
        purity_tsv = join(PURPLE_DIR, "purple", "{patient}.{tumor_sample_type}.purple.purity.tsv"),
        cnv_somatic = join(PURPLE_DIR, "purple", "{patient}.{tumor_sample_type}.purple.cnv.somatic.tsv"),
        cnv_gene = join(PURPLE_DIR, "purple", "{patient}.{tumor_sample_type}.purple.cnv.gene.tsv"),
        purity_ploidy = join(PURPLE_DIR, "purity_ploidy.tsv")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        amber_dir = join(PURPLE_DIR, "amber"),
        cobalt_dir = join(PURPLE_DIR, "cobalt"),
        out_dir = join(PURPLE_DIR, "purple"),
        normal_name = "{patient}.NORMAL",
        tumor_name = "{patient}.{tumor_sample_type}",
        jar = config["purple_jar"],
        java_mem = config["purple_java_mem"],
        min_purity = config["purple_min_purity"],
        max_purity = config["purple_max_purity"],
        min_ploidy = config["purple_min_ploidy"],
        max_ploidy = config["purple_max_ploidy"]
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "purple_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Purple 시작: {params.tumor_name} ===" > {log}
        echo "시작 시간: $(date)" >> {log}

        GC_PROFILE=$(find {input.res_dir} -name "GC_profile.*.cnp" | head -1)
        ENSEMBL_DATA=$(find {input.res_dir} -type d -name "ensembl_data" | head -1)
        echo "gc_profile: $GC_PROFILE" >> {log}
        echo "ensembl_data: $ENSEMBL_DATA" >> {log}

        java -Xmx{params.java_mem} -jar {params.jar} \
            -reference {params.normal_name} \
            -germline_vcf {input.germline_vcf} \
            -tumor {params.tumor_name} \
            -output_dir {params.out_dir} \
            -amber {params.amber_dir} \
            -cobalt {params.cobalt_dir} \
            -gc_profile "$GC_PROFILE" \
            -ref_genome {input.ref} \
            -ref_genome_version 38 \
            -ensembl_data_dir "$ENSEMBL_DATA" \
            -threads {threads} \
            -min_purity {params.min_purity} -max_purity {params.max_purity} \
            -min_ploidy {params.min_ploidy} -max_ploidy {params.max_ploidy} >> {log} 2>&1

        # purity / ploidy 요약 추출 (1열 purity, 5열 ploidy)
        cut -f1,5 {output.purity_tsv} | tail -n+2 > {output.purity_ploidy}

        echo "purity/ploidy: $(cat {output.purity_ploidy})" >> {log}
        echo "=== Purple 완료: $(date) ===" >> {log}
        """
