# Copy Number 및 Purity/Ploidy 분석 관련 규칙들

# 1. SAVANA SV 호출 및 CNV 분석 (각 종양 샘플에 대해)
rule savana_sv:
    input:
        tumor_bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam"),
        tumor_bam_index = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam.bai"),
        normal_bam = lambda wildcards: join(OUTPUT_DIR, wildcards.patient, "phasing", f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}.hiphase.bam"),
        normal_bam_index = lambda wildcards: join(OUTPUT_DIR, wildcards.patient, "phasing", f"{wildcards.patient}.{get_normal_sample(wildcards.patient)}.hiphase.bam.bai"),
        ref = REF_FASTA,
        ref_index = REF_FASTA + ".fai",
        phased_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz")
    output:
        savana_vcf = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}.classified.somatic.vcf"),
        cnv_segments = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}_read_counts_mnorm_log2r_segmented.tsv"),
        cnv_absolute_cn = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}_segmented_absolute_copy_number.tsv"),
        purity_ploidy = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}_fitted_purity_ploidy.tsv"),
        purity_ploidy_solutions = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}_ranked_solutions.tsv")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "cnv", "savana_{patient}_{tumor_sample_type}"),
        sample_name = "{patient}.{tumor_sample_type}",
        min_supp_reads = config.get("savana_min_supp", 3),
        min_af = config.get("savana_min_af", 0.05),
        svlength = config.get("savana_svlength", 50)
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "savana_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail
        LOG_FILE=$(readlink -f {log})
        OUT_DIR=$(readlink -f {params.out_dir})
        CONTAINER=$(readlink -f {config[savana_container]})
        TUMOR_BAM=$(readlink -f {input.tumor_bam})
        NORMAL_BAM=$(readlink -f {input.normal_bam})
        REF_FASTA=$(readlink -f {input.ref})
        PHASED_VCF=$(readlink -f {input.phased_vcf})
        CONTIGS_FILE="$OUT_DIR/contigs.txt"

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"
        
        echo "SAVANA SV 및 CNV 분석 시작: {wildcards.patient}.{wildcards.tumor_sample_type}" > "$LOG_FILE"
        echo "시작 시간: $(date)" >> "$LOG_FILE"

        # Contig 목록 생성
        echo -e "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY" > "$CONTIGS_FILE"

        singularity run -B "$PWD":"$PWD" \
            "${{CONTAINER}}" \
            savana \
              --tumour "{input.tumor_bam}" \
              --normal "{input.normal_bam}" \
              --ref "{input.ref}" \
              --threads {threads} \
              --outdir "{params.out_dir}" \
              --contigs "$CONTIGS_FILE" \
              --pb \
              --single_bnd \
              --no_blacklist \
              --sample "{params.sample_name}" \
              --length {params.svlength} \
              --snp_vcf "{input.phased_vcf}" \
              --min_support {params.min_supp_reads} \
              --min_af {params.min_af} \
              >> "$LOG_FILE" 2>&1

        echo "SAVANA 완료: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "종료 시간: $(date)" >> "$LOG_FILE"
        """

# 2. Wakhan Copy Number 분석 (Severus 결과 활용)
rule wakhan_cnv:
    input:
        tumor_bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam"),
        tumor_bam_index = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{tumor_sample_type}.hiphase.bam.bai"),
        ref = REF_FASTA,
        ref_index = REF_FASTA + ".fai",
        severus_sv_vcf = join(OUTPUT_DIR, "{patient}", "sv", "severus_{patient}", "somatic_SVs", "severus_somatic.vcf"),
        normal_germline_vcf = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.NORMAL.normalized.vcf.gz")
    output:
        wakhan_tar = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}_wakhan.tar.gz"),
        purity_ploidy = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}", "purity_ploidy.tsv"),
        copynumbers_segments = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}.copynumbers_segments.bed"),
        loh_regions = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}.loh_regions.bed"),
        cancer_genes_copynumber = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}", "{patient}.{tumor_sample_type}.cancer_genes_copynumber.bed")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir = join(OUTPUT_DIR, "{patient}", "cnv", "wakhan_{patient}_{tumor_sample_type}"),
        sample_name = "{patient}.{tumor_sample_type}",
        purity_range = config.get("wakhan_purity_range", "0.2-1.0"),
        ploidy_range = config.get("wakhan_ploidy_range", "1-6")
    threads: config["threads"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "wakhan_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail
        LOG_FILE=$(readlink -f {log})
        OUT_DIR=$(readlink -f {params.out_dir})
        CONTAINER=$(readlink -f {config[wakhan_container]})
        TUMOR_BAM=$(readlink -f {input.tumor_bam})
        REF_FASTA=$(readlink -f {input.ref})
        SEVERUS_VCF=$(readlink -f {input.severus_sv_vcf})
        NORMAL_VCF=$(readlink -f {input.normal_germline_vcf})

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"

        echo "Wakhan Copy Number 분석 시작: {wildcards.patient}.{wildcards.tumor_sample_type}" > "$LOG_FILE"
        echo "시작 시간: $(date)" >> "$LOG_FILE"

        singularity run -B "$PWD":"$PWD" \
            "${{CONTAINER}}" \
            wakhan \
              --threads {threads} \
              --target-bam "$TUMOR_BAM" \
              --reference "$REF_FASTA" \
              --genome-name "{params.sample_name}" \
              --out-dir "$OUT_DIR/{params.sample_name}_wakhan" \
              --breakpoints "$SEVERUS_VCF" \
              --loh-enable \
              --ploidy-range "{params.ploidy_range}" \
              --purity-range "{params.purity_range}" \
              --centromere /opt/wakhan/Wakhan/src/annotations/grch38.cen_coord.curated.bed \
              --pdf-enable \
              --normal-phased-vcf "$NORMAL_VCF" \
              >> "$LOG_FILE" 2>&1
        
        WAKHAN_OUT_DIR="$OUT_DIR/{params.sample_name}_wakhan"
        
        # Purity/Ploidy 결과 정리
        echo -e "folder_name\tploidy\tpurity\tconfidence" > "$OUT_DIR/folder_numbers.tsv"
        
        find "$WAKHAN_OUT_DIR" -type d -regex ".*/[0-9.]+_[0-9.]+_[0-9.]+$" | while read dir; do
            folder_name=$(basename "$dir")
            ploidy=$(echo "$folder_name" | cut -d'_' -f1)
            purity=$(echo "$folder_name" | cut -d'_' -f2)
            confidence=$(echo "$folder_name" | cut -d'_' -f3)
            echo -e "$folder_name\t$ploidy\t$purity\t$confidence"
        done | sort -t$'\t' -k4,4rn > "$OUT_DIR/temp.tsv"
        
        (head -n 1 "$OUT_DIR/folder_numbers.tsv"; cat "$OUT_DIR/temp.tsv") > {output.purity_ploidy}
        rm "$OUT_DIR/temp.tsv" "$OUT_DIR/folder_numbers.tsv"
        
        # 전체 결과 압축
        tar -czvf {output.wakhan_tar} -C "$OUT_DIR" "{params.sample_name}_wakhan" >> "$LOG_FILE" 2>&1
        
        # 최고 해법 파일 복사
        best_folder=$(head -n 2 {output.purity_ploidy} | tail -n 1 | cut -f1)
        cp -r "$WAKHAN_OUT_DIR/${{best_folder}}" "$OUT_DIR/{params.sample_name}_wakhan_best"
        rm -rf "$OUT_DIR/{params.sample_name}_wakhan_best/variation_plots"
        
        # 결과 파일 이름 변경
        mv "$OUT_DIR/{params.sample_name}_wakhan_best/bed_output/"*copynumbers_segments.bed {output.copynumbers_segments}
        mv "$OUT_DIR/{params.sample_name}_wakhan_best/bed_output/loh_regions.bed" {output.loh_regions}
        mv "$OUT_DIR/{params.sample_name}_wakhan_best/bed_output/cancer_genes_copynumber_states.bed" {output.cancer_genes_copynumber}
        
        rm -rf "$WAKHAN_OUT_DIR" "$OUT_DIR/{params.sample_name}_wakhan_best"
        
        echo "Wakhan 완료: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "종료 시간: $(date)" >> "$LOG_FILE"
        """ 