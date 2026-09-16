# Purity/ploidy and allele-specific copy number with PURPLE (HMFtools)
# Runs Amber (BAF) -> Cobalt (read ratio) -> Purple (purity/ploidy/CN) in order.
# Same invocation as clonality.wdl in HiFi-somatic-WDL, but using local jars
# (config: amber_jar / cobalt_jar / purple_jar) instead of a container.

PURPLE_DIR = join(OUTPUT_DIR, "{patient}", "cnv", "purple_{patient}_{tumor_sample_type}")

# 0. Unpack the HMF reference bundle once; every sample reuses it
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

        echo "=== Unpacking HMF reference bundle ===" > {log}
        echo "Start time: $(date)" >> {log}

        tar -xzf {input.tarball} -C {output.res_dir} 2>> {log}

        # Check that the required files are present
        AMBER_LOCI=$(find {output.res_dir} -name "AmberGermlineSites.*.tsv.gz" | head -1)
        GC_PROFILE=$(find {output.res_dir} -name "GC_profile.*.cnp" | head -1)
        ENSEMBL_DATA=$(find {output.res_dir} -type d -name "ensembl_data" | head -1)

        if [ -z "$AMBER_LOCI" ] || [ -z "$GC_PROFILE" ] || [ -z "$ENSEMBL_DATA" ]; then
            echo "ERROR: required files not found in the HMF bundle." >> {log}
            echo "  AmberGermlineSites: $AMBER_LOCI" >> {log}
            echo "  GC_profile: $GC_PROFILE" >> {log}
            echo "  ensembl_data: $ENSEMBL_DATA" >> {log}
            exit 1
        fi

        echo "AmberGermlineSites: $AMBER_LOCI" >> {log}
        echo "GC_profile: $GC_PROFILE" >> {log}
        echo "ensembl_data: $ENSEMBL_DATA" >> {log}
        echo "=== done: $(date) ===" >> {log}
        """

# 1. Amber - B-allele frequency and contamination
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

        echo "=== Amber start: {params.tumor_name} ===" > {log}
        echo "Start time: $(date)" >> {log}

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

        echo "=== Amber done: $(date) ===" >> {log}
        """

# 2. Cobalt - read depth ratio and GC correction
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

        echo "=== Cobalt start: {params.tumor_name} ===" > {log}
        echo "Start time: $(date)" >> {log}

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

        echo "=== Cobalt done: $(date) ===" >> {log}
        """

# 3. Purple - purity / ploidy / allele-specific CNV
#    The DeepSomatic VCF is not the two-sample tumor+normal format PURPLE expects, so
#    the fit uses Amber/Cobalt output only (same as the earlier manual runs).
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

        echo "=== Purple start: {params.tumor_name} ===" > {log}
        echo "Start time: $(date)" >> {log}

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

        # Extract purity (col 1) and ploidy (col 5)
        cut -f1,5 {output.purity_tsv} | tail -n+2 > {output.purity_ploidy}

        echo "purity/ploidy: $(cat {output.purity_ploidy})" >> {log}
        echo "=== Purple done: $(date) ===" >> {log}
        """
