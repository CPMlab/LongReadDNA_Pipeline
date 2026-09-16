# Alignment and QC rules

# 1. Alignment (all samples)
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
        
        # Convert to an absolute path
        log_file = os.path.abspath(log[0])
        output_bam = os.path.abspath(output.bam)
        ref_fasta = os.path.abspath(input.ref)
        
        # Create the output directory
        os.makedirs(os.path.dirname(output_bam), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        
        # Unique working directory (absolute path)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        work_dir = os.path.abspath(os.path.join(os.path.dirname(output_bam), f"mapping_work_{wildcards.patient}_{wildcards.sample_type}_{timestamp}"))
        os.makedirs(work_dir, exist_ok=True)
        
        # Reset the log
        shell(f"echo '=== Alignment start: {wildcards.patient}.{wildcards.sample_type} ===' > {log_file}")
        shell(f"echo 'Start time: $(date)' >> {log_file}")
        shell(f"echo 'Input files: {input.bam_files}' >> {log_file}")
        shell(f"echo 'Working directory: {work_dir}' >> {log_file}")
        
        try:
            # Single BAM
            if len(input.bam_files) == 1:
                shell(f"""
                    echo "Aligning a single BAM" >> {log_file}
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
                
                    echo "Indexing the BAM" >> {log_file}
                    samtools index -@{threads} {output_bam} >> {log_file} 2>&1
                    echo "Single BAM alignment done" >> {log_file}
                """)
        
            # Multiple BAMs
            else:
                shell(f"echo 'Aligning multiple BAMs (total {len(input.bam_files)})' >> {log_file}")
                
                temp_bams = []
                for i, bam in enumerate(input.bam_files):
                    # Unique temporary file name (absolute path)
                    temp_bam = os.path.join(work_dir, f"{wildcards.patient}_{wildcards.sample_type}_temp_bam{i+1}_{timestamp}.bam")
                    temp_bams.append(temp_bam)
                    
                    shell(f"""
                        echo "BAM {i+1}/{len(input.bam_files)} alignment start: $(basename {bam})" >> {log_file}
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
                        echo "BAM {i+1}/{len(input.bam_files)} alignment done" >> {log_file}
                    """)
            
                # Merge results and index
                temp_bams_str = " ".join(temp_bams)
                shell(f"""
                    echo "Merging BAMs" >> {log_file}
                    cd {work_dir}
                    samtools merge -@{threads} {output_bam} {temp_bams_str} >> {log_file} 2>&1
                    echo "Indexing the BAM" >> {log_file}
                    samtools index -@{threads} {output_bam} >> {log_file} 2>&1
                    echo "Multi-BAM alignment done" >> {log_file}
                """)
        
        finally:
            # Clean up the working directory
            shell(f"echo 'Cleaning the working directory: {work_dir}' >> {log_file}")
            shell(f"rm -rf {work_dir}")
            shell(f"echo '=== Alignment done: {wildcards.patient}.{wildcards.sample_type} ===' >> {log_file}")
            shell(f"echo 'End time: $(date)' >> {log_file}")

# 2. Coverage (all samples)
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
        # Create log and output directories
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.summary})
        
        echo "=== mosdepth start: {wildcards.patient}.{wildcards.sample_type} ===" > {log}
        echo "Start time: $(date)" >> {log}
        echo "Input BAM: {input.bam}" >> {log}
        echo "Output prefix: {OUTPUT_DIR}/{wildcards.patient}/qc/{wildcards.patient}.{wildcards.sample_type}" >> {log}
        
        # mosdepth coverage
        mosdepth -t {threads} \
            -n --fast-mode \
            --by 500 \
            {OUTPUT_DIR}/{wildcards.patient}/qc/{wildcards.patient}.{wildcards.sample_type} \
            {input.bam} \
            >> {log} 2>&1
        
        echo "=== mosdepth done ===" >> {log}
        echo "End time: $(date)" >> {log}
        """ 