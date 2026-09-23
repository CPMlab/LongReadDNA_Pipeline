# Updated Snakemake rules for CpG methylation, DMR calling, and annotation
# Notable changes
#   1. No `cd` inside the shell block; absolute paths are used instead
#   2. set -euo pipefail so failures stop immediately
#   3. Temp/log/output paths go through OUT_DIR variables
#   4. The R script writes into the current directory, so a prefix is passed
#      and the result is moved into place afterwards

#######################################################################
# 1. CpG methylation scores
#######################################################################
rule cpg_methylation:
    input:
        bam = join(OUTPUT_DIR, "{patient}", "phasing", "{patient}.{sample_type}.hiphase.bam"),
        ref = REF_FASTA
    output:
        combined_bed = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.combined.bed.gz"),
        hap1_bed     = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap1.bed.gz"),
        hap2_bed     = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap2.bed.gz"),
        combined_bw  = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.combined.bw"),
        hap1_bw      = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap1.bw"),
        hap2_bw      = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg.hap2.bw")
    params:
        output_prefix = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{sample_type}.cpg"),
        min_coverage = config.get("methylation_min_coverage", 5),
        min_mapq     = config.get("methylation_min_mapq", 1)
    threads: THREADS
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "cpg_methylation_{patient}_{sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        # Make sure the output directory exists
        mkdir -p "$(dirname {params.output_prefix})"
        mkdir -p "$(dirname {log})"

        echo "CpG methylation start: {wildcards.patient}.{wildcards.sample_type}" > {log}
        echo "Start time: $(date)" >> {log}

        aligned_bam_to_cpg_scores --version >> {log} 2>&1

        aligned_bam_to_cpg_scores \
          --threads {threads} \
          --bam {input.bam} \
          --ref {input.ref} \
          --output-prefix {params.output_prefix} \
          --min-mapq {params.min_mapq} \
          --min-coverage {params.min_coverage} \
          >> {log} 2>&1

        echo "CpG methylation done: {wildcards.patient}.{wildcards.sample_type}" >> {log}
        echo "End time: $(date)" >> {log}
        """

#######################################################################
# 2. DSS differentially methylated regions (tumor vs normal)
#######################################################################
rule dss_dmr:
    input:
        tumor_bed  = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.{tumor_sample_type}.cpg.combined.bed.gz"),
        normal_bed = join(OUTPUT_DIR, "{patient}", "methylation", "{patient}.NORMAL.cpg.combined.bed.gz")
    output:
        dmr_tsv = join(OUTPUT_DIR, "{patient}", "dmr", "{patient}.{tumor_sample_type}_vs_NORMAL.DMR.tsv")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir    = join(OUTPUT_DIR, "{patient}", "dmr"),
        sample_name = "{patient}.{tumor_sample_type}_vs_NORMAL",
        dss_script  = "scripts/DSS_tumor_normal.R"
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR, "{patient}", "logs", "dss_dmr_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        OUT_DIR={params.out_dir}
        LOG_FILE={log}
        DSS_SCRIPT=$(readlink -f {params.dss_script})

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"

        echo "DSS start: {wildcards.patient}.{wildcards.tumor_sample_type} vs NORMAL" > "$LOG_FILE"
        echo "Start time: $(date)" >> "$LOG_FILE"

        # Temporary files (created inside OUT_DIR)
        TUMOR_TMP="$OUT_DIR/{params.sample_name}.tumor.tmp"
        NORMAL_TMP="$OUT_DIR/{params.sample_name}.normal.tmp"

        gunzip -c {input.tumor_bed} | grep -v '^#' | cut -f1,2,6,7 > "$TUMOR_TMP" 2>> "$LOG_FILE"
        gunzip -c {input.normal_bed} | grep -v '^#' | cut -f1,2,6,7 > "$NORMAL_TMP" 2>> "$LOG_FILE"

        echo "Running the DSS R script" >> "$LOG_FILE"

        Rscript --vanilla "$DSS_SCRIPT" \
            "$TUMOR_TMP" \
            "$NORMAL_TMP" \
            {output.dmr_tsv} \
            {threads} \
            >> "$LOG_FILE" 2>&1

        rm -f "$TUMOR_TMP" "$NORMAL_TMP"

        echo "DSS done: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "End time: $(date)" >> "$LOG_FILE"
        """

#######################################################################
# 3. DMR annotation
#######################################################################
rule annotate_dmr:
    input:
        dmr_tsv = join(OUTPUT_DIR,
                       "{patient}",
                       "dmr",
                       "{patient}.{tumor_sample_type}_vs_NORMAL.DMR.tsv")
    output:
        annotated_dmr = join(OUTPUT_DIR,
                             "{patient}",
                             "dmr",
                             "{patient}.{tumor_sample_type}_vs_NORMAL.annotated_DMR.tsv.gz")
    wildcard_constraints:
        tumor_sample_type = TUMOR_TYPE_CONSTRAINT
    params:
        out_dir        = join(OUTPUT_DIR, "{patient}", "dmr"),
        sample_name    = "{patient}.{tumor_sample_type}_vs_NORMAL",
        annotate_script = "scripts/annotatr_dmr.R"
    threads: config["threads_low"]
    log:
        join(OUTPUT_DIR,
             "{patient}",
             "logs",
             "annotate_dmr_{patient}_{tumor_sample_type}.log")
    shell:
        r"""
        set -euo pipefail

        OUT_DIR="{params.out_dir}"
        LOG_FILE="{log}"
        ANNO_SCRIPT="$(readlink -f {params.annotate_script})"

        mkdir -p "$OUT_DIR"
        mkdir -p "$(dirname "$LOG_FILE")"

        echo "DMR annotation start: {wildcards.patient}.{wildcards.tumor_sample_type}" > "$LOG_FILE"
        echo "Start time: $(date)" >> "$LOG_FILE"

        PREFIX="$OUT_DIR/{params.sample_name}"

        Rscript --vanilla "$ANNO_SCRIPT" \
            "{input.dmr_tsv}" \
            "$PREFIX" \
            {threads} >> "$LOG_FILE" 2>&1

        # Check the result and move it into place
        # annotatr_dmr.R writes a summary plus one file per genomic region; the summary is the
        # rule's declared output (the per-region files stay alongside it). Copy just the summary,
        # not a glob, or mv gets several sources and one destination and fails.
        SUMMARY="${PREFIX}_dmr_annotation_summary.tsv.gz"
        if [[ -s "$SUMMARY" ]]; then
            cp "$SUMMARY" "{output.annotated_dmr}"
        else
            echo "No annotated DMR file was produced; writing an empty one." >> "$LOG_FILE"
            touch "{output.annotated_dmr}"
        fi

        echo "DMR annotation done: {wildcards.patient}.{wildcards.tumor_sample_type}" >> "$LOG_FILE"
        echo "End time: $(date)" >> "$LOG_FILE"
        """
