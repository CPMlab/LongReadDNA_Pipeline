# Structural variant calling and processing rules

# 1. Structural variant calling with Severus (multimode)
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
        
        echo "Severus multimode start" > {log}
        echo "Control BAM: {input.normal_bam}" >> {log}
        echo "Target BAMs: {input.tumor_bams}" >> {log}
        
        # Severus multimode - several target BAMs at once
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
            
        echo "Severus multimode done" >> {log}
        """

# 2. tabix filtering
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
        
        echo "Indexing and filtering the VCF: {input.uncompressed_vcf}" > {log}
        
        bcftools sort -Oz -o tmp.{wildcards.patient}.sorted.vcf.gz "{input.uncompressed_vcf}" >> {log} 2>&1
        
        tabix -p vcf tmp.{wildcards.patient}.sorted.vcf.gz >> {log} 2>&1
        
        bcftools view -R "{input.contig_bed}" -Oz -o "{output.vcf}" tmp.{wildcards.patient}.sorted.vcf.gz >> {log} 2>&1
        
        tabix -p vcf "{output.vcf}" >> {log} 2>&1
        
        rm -f tmp.{wildcards.patient}.sorted.vcf.gz tmp.{wildcards.patient}.sorted.vcf.gz.tbi >> {log} 2>&1
        
        echo "Done. Result: {output.vcf}, {output.tbi}" >> {log}
        """

# 3. SVpack filtering and annotation
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
        
        echo "SVpack start: {input.filtered_sv_vcf}" > {log}
        
        # Write the svpack pipeline result to an uncompressed temporary file
        python {params.svpack_executable} filter --pass-only "{input.filtered_sv_vcf}" 2>> {log} | \
        python {params.svpack_executable} filter --min-svlen 50 - 2>> {log} | \
        python {params.svpack_executable} match -v - "{input.match_vcf}" 2>> {log} | \
        python {params.svpack_executable} consequence - "{input.reference_gff}" 2>> {log} | \
        python {params.svpack_executable} tagzygosity - 2>> {log} \
        > {params.temp_uncompressed_svpack_vcf}
        
        # Make sure the temporary file is not empty
        if [ ! -s "{params.temp_uncompressed_svpack_vcf}" ]; then
            echo "Error: the svpack pipeline did not produce {params.temp_uncompressed_svpack_vcf} (missing or empty)" >> {log}
            exit 1
        fi

        echo "Compressing and indexing the temporary VCF" >> {log}
        # bgzip the temporary file into the final output
        bgzip -f -c {params.temp_uncompressed_svpack_vcf} > {output.svpack_vcf_gz} 2>> {log}
        # Index the final compressed file
        tabix -f -p vcf {output.svpack_vcf_gz} 2>> {log}
        
        # Remove the uncompressed temporary file
        rm -f {params.temp_uncompressed_svpack_vcf}
        
        echo "SVpack done: {output.svpack_vcf_gz}" >> {log}
        """ 

# 4. Recovering missing BND mates
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
        
        echo "=== recover_mate_bnd start ===" > {log}
        echo "Working directory: $(pwd)" >> {log}
        echo "Input original_vcf: {input.original_vcf}" >> {log}
        echo "Input filtered_vcf_gz: {input.filtered_vcf_gz}" >> {log}

        echo "Looking for missing mate IDs" >> {log}
        
        bcftools query -f '%ID\\t%MATE_ID\\n' {input.filtered_vcf_gz} 2>> {log} | grep -v '\\.' 2>> {log} | cut -f1 | sort > {params.tmp_ids}
        bcftools query -f '%ID\\t%MATE_ID\\n' {input.filtered_vcf_gz} 2>> {log} | grep -v '\\.' 2>> {log} | cut -f2 | sort > {params.tmp_mate_ids}

        comm -13 {params.tmp_ids} {params.tmp_mate_ids} > {params.missing_mate_list} 2>> {log}

        if [ -s "{params.missing_mate_list}" ]; then
            echo "Missing mate IDs found in {params.missing_mate_list} ; extracting and merging." >> {log}
            
            bcftools view -i ID=@$(basename {params.missing_mate_list}) "{input.original_vcf}" 2>> {log} | \
                bcftools sort -Oz -o {params.missing_mate_vcf} >> {log} 2>&1
            
            if [ ! -s "{params.missing_mate_vcf}" ]; then
                echo "ERROR: MISSING_MATE_VCF is empty or missing after bcftools view/sort." >> {log}
                exit 1
            fi
            tabix -f -p vcf {params.missing_mate_vcf} >> {log} 2>&1
            
            echo "Merging and sorting" >> {log}
            bcftools concat -a "{input.filtered_vcf_gz}" {params.missing_mate_vcf} 2>> {log} | \
                bcftools sort -Oz -o "{output.vcf}" >> {log} 2>&1
        else
            echo "No missing mate IDs; using the filtered VCF as the final result." >> {log}
            cp "{input.filtered_vcf_gz}" "{output.vcf}" >> {log} 2>&1
        fi
        
        tabix -f -p vcf "{output.vcf}" >> {log} 2>&1
        
        # Final cleanup of temporary files
        rm -f {params.tmp_ids} {params.tmp_mate_ids} {params.missing_mate_list} {params.missing_mate_vcf} {params.missing_mate_vcf}.tbi >> {log} 2>&1
        
        echo "recover_mate_bnd done: {output.vcf}" >> {log}
        """ 