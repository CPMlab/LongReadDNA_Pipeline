# Variant calling and phasing rules

# 1. Clair3 germline variant calling (all samples)
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p {params.out_dir}

        echo "=== Clair3 start: {wildcards.patient}.{wildcards.sample_type} ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output directory: {params.out_dir}" >> {log}

        # Temporary file name
        TEMP_VCF=$(mktemp {params.out_dir}/temp_clair3_{wildcards.sample_type}_XXXXXX.vcf)
        echo "Temporary VCF: $TEMP_VCF" >> {log}

        # Clair3 inside the Singularity container
        echo "Running the Clair3 container" >> {log}
        singularity exec \
            --bind $(pwd):$(pwd) \
            {config[clair3_container]} \
            /bin/bash -c "
                set -euxo pipefail

                # Run Clair3
                /opt/bin/run_clair3.sh \
                    --bam_fn={input.bam} \
                    --ref_fn={input.ref} \
                    --threads={threads} \
                    --platform='hifi' \
                    --model_path='/opt/models/hifi_revio' \
                    --output={params.out_dir} \
                    --sample_name={params.sample}

                # Drop LowQual calls into a temporary VCF
                gunzip -c {params.out_dir}/merge_output.vcf.gz | \
                    awk 'BEGIN{{OFS=\"\\t\"}} {{if(/^#/ || \\$7 != \"LowQual\") {{print \\$0}}}}' > ${{TEMP_VCF}}
            " >> {log} 2>&1

        echo "Clair3 finished, post-processing" >> {log}

        # Compress the temporary VCF with bgzip
        bgzip -f ${{TEMP_VCF}} 2>> {log}
        mv ${{TEMP_VCF}}.gz {output.vcf} 2>> {log}

        # Index the final VCF
        echo "Indexing the VCF" >> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        
        echo "=== Clair3 done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# 2. Somatic variant calling (each tumor sample against the normal)
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
        deepsomatic_slots=1  # only one DeepSomatic job at a time
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "deepsomatic_{patient}_{tumor_sample_type}.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        echo "DeepSomatic start: {wildcards.patient}.{wildcards.tumor_sample_type}" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input tumor BAM: {input.tumor_bam}" >> {log}
        echo "Input normal BAM: {input.normal_bam}" >> {log}
        echo "Reference genome: {input.ref}" >> {log}
        
        # Run DeepSomatic with default settings
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
        
        echo "DeepSomatic container finished" >> {log}
        echo "Post-processing start time: $(date)" >> {log}
            
        # Keep PASS variants only
        bcftools view \
            -f PASS -Oz \
            -o {params.out_dir}/somatic_PASS.vcf.gz \
            {params.out_dir}/somatic.vcf.gz \
            2>&1 | tee -a {log}
            
        # Convert homozygous (1/1) calls to heterozygous (0/1)
        bcftools +setGT {params.out_dir}/somatic_PASS.vcf.gz -- -t q -i 'GT="1/1"' -n c:"0/1" | \
            bcftools sort -Oz -o {output.vcf} \
            2>&1 | tee -a {log}
            
        # Build the index
        tabix -p vcf {output.vcf} 2>&1 | tee -a {log}
        
        echo "DeepSomatic done: {wildcards.patient}.{wildcards.tumor_sample_type}" >> {log}
        echo "End time: $(date)" >> {log}
        """ 

# 3. Phasing the normal sample
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        
        echo "=== HiPhase start: {wildcards.patient}.NORMAL ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output BAM: {output.bam}" >> {log}
        
        # HiPhase phasing
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
            
        echo "Indexing the BAM" >> {log}
        # Index the BAM
        samtools index -@{threads} {output.bam} >> {log} 2>&1
        
        echo "=== HiPhase done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# 4. Phasing the tumor sample (including somatic variants)
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.bam})
        
        echo "=== HiPhase (tumor) start: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Input germline VCF: {input.germline_vcf}" >> {log}
        echo "Input somatic VCF: {input.somatic_vcf}" >> {log}
        echo "Output BAM: {output.bam}" >> {log}
        
        # Phasing the tumor sample (including somatic variants)
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
            
        echo "Indexing the BAM" >> {log}
        # Index the BAM
        samtools index -@{threads} {output.bam} >> {log} 2>&1
        
        echo "=== HiPhase (tumor) done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# 5. Normalizing the phased normal VCF
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
        # Create the log directory
        mkdir -p $(dirname {log})
        
        echo "=== VCF normalization start: {wildcards.patient}.NORMAL ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output VCF: {output.vcf}" >> {log}
        
        # Index the VCF
        echo "Indexing the VCF" >> {log}
        bcftools index --threads {threads} {input.vcf} >> {log} 2>&1
        
        # Run normalization
        echo "Normalizing the VCF" >> {log}
        bcftools view {input.vcf} | \
          sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \
          bcftools norm --threads {threads} --multiallelics - \
            --output-type b --fasta-ref {input.ref} | \
          bcftools sort -Oz -o {output.vcf} \
          2>> {log}
        
        # Index the normalized VCF
        echo "Indexing the normalized VCF" >> {log}
        bcftools index --threads {threads} -t {output.vcf} >> {log} 2>&1
        
        echo "=== VCF normalization done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# 6. Normalizing the phased tumor VCF
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
        # Create the log directory
        mkdir -p $(dirname {log})
        
        echo "=== Tumor VCF normalization start: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output VCF: {output.vcf}" >> {log}
        
        # Index the VCF
        echo "Indexing the VCF" >> {log}
        bcftools index --threads {threads} {input.vcf} >> {log} 2>&1
        
        # Run normalization
        echo "Normalizing the VCF" >> {log}
        bcftools view {input.vcf} | \
          sed -e 's/ID=AD,Number=\\./ID=AD,Number=R/' | \
          bcftools norm --threads {threads} --multiallelics - \
            --output-type b --fasta-ref {input.ref} | \
          bcftools sort -Oz -o {output.vcf} \
          2>> {log}
        
        # Index the normalized VCF
        echo "Indexing the normalized VCF" >> {log}
        bcftools index --threads {threads} -t {output.vcf} >> {log} 2>&1
        
        echo "=== Tumor VCF normalization done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """ 