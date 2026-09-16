# VEP annotation rules

# 1. VEP annotation (normal sample)
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.vcf})
        
        echo "=== VEP start: {wildcards.patient}.NORMAL ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output VCF: {output.vcf}" >> {log}
        echo "Temporary directory: {params.tmp_dir}" >> {log}
        
        # Extract the VEP cache
        echo "Extracting the VEP cache" >> {log}
        tar -xzf {input.vep_cache} -C {params.tmp_dir} >> {log} 2>&1
        
        # Run VEP
        echo "Running VEP" >> {log}
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
          
        # Build the index
        echo "Indexing the VCF" >> {log}
        tabix -p vcf {output.vcf} >> {log} 2>&1
        
        # Clean up temporary files
        echo "Cleaning up temporary files" >> {log}
        rm -rf {params.tmp_dir} >> {log} 2>&1
        
        echo "=== VEP done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """

# 2. VEP annotation (somatic variants)
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {output.vcf})
        
        echo "=== VEP (somatic) start: {wildcards.patient}.{wildcards.tumor_sample_type} ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input VCF: {input.vcf}" >> {log}
        echo "Output VCF: {output.vcf}" >> {log}
        echo "Temporary directory: {params.tmp_dir}" >> {log}
        
        # Extract the VEP cache
        echo "Extracting the VEP cache" >> {log}
        tar -xzf {input.vep_cache} -C {params.tmp_dir} >> {log} 2>&1
        
        # Run VEP
        echo "Running VEP" >> {log}
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
          
        # Build the index
        echo "Indexing the VCF" >> {log}
        tabix -p vcf {output.vcf} >> {log} 2>&1
        
        # Clean up temporary files
        echo "Cleaning up temporary files" >> {log}
        rm -rf {params.tmp_dir} >> {log} 2>&1
        
        echo "=== VEP (somatic) done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """ 