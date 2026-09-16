# %%
import pycirclize
from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset
import pandas as pd
import re
import argparse
import json
import gzip
import os

def expand_info_field(df):
    # Split INFO field into individual components
    info_fields = df['INFO'].str.split(';')
    
    # Dictionary to store all key-value pairs
    info_dict = {}
    
    # Process each row
    for idx, fields in enumerate(info_fields):
        row_dict = {}
        if fields is not None:  # Handle potential NaN values
            for field in fields:
                if '=' in field:  # Only process fields with key=value format
                    key, value = field.split('=', 1)  # Split on first '=' only
                    row_dict[key] = value
        info_dict[idx] = row_dict
    
    # Create DataFrame from dictionary
    info_df = pd.DataFrame.from_dict(info_dict, orient='index')
    
    # Merge expanded columns with original DataFrame
    result_df = pd.concat([df, info_df], axis=1)
    
    return result_df

def extract_chr2_pos2(df):
    # Filter for BND entries and create a copy
    bnd_df = df[df['SVTYPE'] == 'BND'].copy()
    # If empty, return empty DataFrame
    if bnd_df.empty:
        return bnd_df
    
    # Extract chr2 and pos2 using regex
    # ALT field format example: ]chr2:123456]N
    pattern = r'[[\]]?(chr\w+):(\d+)[[\]]'
    
    # Create new columns
    bnd_df['CHR2'] = bnd_df['ALT'].str.extract(pattern)[0]
    bnd_df['POS2'] = bnd_df['ALT'].str.extract(pattern)[1].astype(int)

    # If BCSQ column exists, extract gene name (second field split by |)
    if 'BCSQ' in bnd_df.columns:
        bnd_df['GENE'] = bnd_df['BCSQ'].str.split('|').str[1]
        # If NA, replace with empty string
        bnd_df['GENE'] = bnd_df['GENE'].fillna('')
    
    return bnd_df

def create_fusion_annotation(df):
    # Create empty FUSION column
    df['FUSION'] = ''
    
    # Filter for BND records only
    bnd_records = df[df['SVTYPE'] == 'BND'].copy()
    
    for idx, row in bnd_records.iterrows():
        # Find mate record using MATEID
        mate = df[df['ID'] == row['MATE_ID']].iloc[0] if not df[df['ID'] == row['MATE_ID']].empty else None
        
        # Check if mate exists and both have GENE information
        if mate is not None and row['GENE'] != '' and mate['GENE'] != '':
            # If both genes are the same, skip
            if row['GENE'] == mate['GENE']:
                continue
            # Sort genes alphabetically and join with ::
            genes = sorted([row['GENE'], mate['GENE']])
            fusion = '::'.join(genes)
            
            # Update FUSION column for both records
            df.loc[idx, 'FUSION'] = fusion
            df.loc[df['ID'] == row['MATE_ID'], 'FUSION'] = fusion

    # Fill empty FUSION with empty string
    df['FUSION'] = df['FUSION'].fillna('')
    
    return df

# %%

def save_empty_circos(output_svg_path, output_png_path, hg38_bed_path):
    circos_empty = Circos.initialize_from_bed(hg38_bed_path, space=3)
    circos_empty.text("No SVs found", size=15)
    fig = circos_empty.plotfig(dpi=300)
    if output_svg_path:
        fig.savefig(output_svg_path, dpi=300, bbox_inches='tight')
    if output_png_path:
        fig.savefig(output_png_path, dpi=300, bbox_inches='tight')

def get_sample_names_from_vcf(vcf_file):
    """Extract the sample names from a VCF."""
    if vcf_file.endswith('.gz'):
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # Sample names come from the last VCF header line
                    columns = line.strip().split('\t')
                    if len(columns) > 9:
                        # Sample names start at column 10
                        return columns[9:]
                    else:
                        return []
    else:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # Sample names come from the last VCF header line
                    columns = line.strip().split('\t')
                    if len(columns) > 9:
                        # Sample names start at column 10
                        return columns[9:]
                    else:
                        return []
    return []

def filter_sv_for_sample(df, sample_idx):
    """Keep only the SVs of one sample (genotype other than 0/0)."""
    if f'SAMPLE_{sample_idx}' not in df.columns:
        return pd.DataFrame()
    
    sample_col = f'SAMPLE_{sample_idx}'
    
    # Take the GT field (before the first ':')
    gt_only = df[sample_col].str.split(':').str[0]
    
    # Keep real genotypes (drop 0/0, ./. and .)
    mask = (gt_only.notna()) & (~gt_only.isin(['0/0', './.', '.', '0|0']))
    
    return df[mask].copy()

def main():
    parser = argparse.ArgumentParser(description='Generate circos plot from VCF file')
    parser.add_argument('--input-vcf', required=True, help='Input VCF file (can be .vcf or .vcf.gz)')
    parser.add_argument('--sample-name', required=True, help='Sample name to be used in the plot')
    parser.add_argument('--mitelman-mcgene', required=True, help='Mitelman fusion database MCGENE file')
    parser.add_argument('--circos-bed', required=True, help='Path to hg38.bed file for Circos initialization')
    parser.add_argument('--cytoband-file', required=True, help='Path to hg38_cytoband.tsv file')
    parser.add_argument('--output-svg', required=True, help='Output path for SVG plot')
    parser.add_argument('--output-png', required=True, help='Output path for PNG plot')
    parser.add_argument('--output-fusion-tsv', required=True, help='Output path for fusion calls TSV')

    args = parser.parse_args()

    # Sample names from the VCF
    actual_sample_names = get_sample_names_from_vcf(args.input_vcf)
    
    if not actual_sample_names:
        print("No samples found in VCF file. Creating empty plot.")
        save_empty_circos(args.output_svg, args.output_png, args.circos_bed)
        pd.DataFrame(columns=['SAMPLE_NAME', 'CHROM', 'POS', 'CHR2', 'POS2', 'GENE', 'FUSION']).to_csv(
            args.output_fusion_tsv, sep='\t', index=False
        )
        return

    try:
        df = pd.read_csv(args.input_vcf, sep="\t", header=0, comment="#")
        if df.empty:
             raise pd.errors.EmptyDataError
    except pd.errors.EmptyDataError:
        save_empty_circos(args.output_svg, args.output_png, args.circos_bed)
        pd.DataFrame(columns=['SAMPLE_NAME', 'CHROM', 'POS', 'CHR2', 'POS2', 'GENE', 'FUSION']).to_csv(
            args.output_fusion_tsv, sep='\t', index=False
        )
        print(f"Input VCF {args.input_vcf} is empty. Empty plots and fusion table created.")
        return

    # Define the standard VCF columns
    vcf_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    
    # Read the number of columns from the dataframe
    num_cols = len(df.columns)
    
    # Generate column names dynamically
    # The first 9 are standard, the rest are sample columns
    if num_cols > 9:
        sample_cols = [f"SAMPLE_{i}" for i in range(num_cols - 9)]
        df.columns = vcf_cols + sample_cols
    else:
        df.columns = vcf_cols[:num_cols]

    # Expand INFO field once for all samples
    df = expand_info_field(df)
    
    # Collects the fusion calls of every sample
    all_fusion_results = []
    
    # Process each sample separately
    for sample_idx, sample_name in enumerate(actual_sample_names):
        print(f"Processing sample: {sample_name} (index: {sample_idx})")
        
        # Keep this sample's SVs
        sample_df = filter_sv_for_sample(df, sample_idx)
        
        if sample_df.empty:
            print(f"No variants found for sample {sample_name}")
            continue
            
        # Extract BNDs
        bnd_only = extract_chr2_pos2(sample_df)
        
        if bnd_only.empty:
            print(f"No BND entries found for sample {sample_name}")
            continue
            
        # Annotate fusions
        bnd_only = create_fusion_annotation(bnd_only)
        bnd_only = bnd_only[abs(bnd_only['POS2'] - bnd_only['POS']) > 100000]
        
        if bnd_only.empty:
            print(f"No BND entries left after filtering for sample {sample_name}")
            continue
            
        # Include the sample name in the output file name
        base_svg = os.path.splitext(args.output_svg)[0]
        base_png = os.path.splitext(args.output_png)[0]
        base_tsv = os.path.splitext(args.output_fusion_tsv)[0]
        
        sample_svg = f"{base_svg}_{sample_name}.svg"
        sample_png = f"{base_png}_{sample_name}.png"
        sample_tsv = f"{base_tsv}_{sample_name}.tsv"
        
        # Draw the circos plot
        circos = Circos.initialize_from_bed(args.circos_bed, space=3)
        circos.text(f"{sample_name} Translocations (GRCh38)", size=15)
        circos.add_cytoband_tracks((95, 100), args.cytoband_file)
        cytoband_df = pd.read_csv(args.cytoband_file, sep="\t")
        
        # Mitelman database
        mitelman_db = pd.read_csv(args.mitelman_mcgene, sep="\t")
        mitelman_db = mitelman_db[mitelman_db['Gene'].str.contains("::")]
        mitelman_db['Fusion'] = mitelman_db['Gene'].str.split("::").apply(lambda x: "::".join(sorted(x)))
        mitelman_db = mitelman_db.groupby('Fusion').filter(lambda x: len(x) >= 3)
        unique_mitelman_fusion = mitelman_db['Fusion'].unique()

        # Filter chromosomes
        bnd_only = bnd_only[bnd_only['CHROM'].isin(cytoband_df['#chrom'])]
        bnd_only = bnd_only[bnd_only['CHR2'].isin(cytoband_df['#chrom'])]

        if bnd_only.empty:
            print(f"No BND entries left after chromosome filtering for sample {sample_name}")
            continue

        # Gene labels
        for sector in circos.sectors:
            sector.text(sector.name, size=10)
            bnd_chr = bnd_only[bnd_only['CHROM'] == sector.name]
            bnd_chr = bnd_chr[bnd_chr['FUSION'] != '']
            if not bnd_chr.empty:
                bnd_track = sector.get_track('cytoband')
                label_pos = bnd_chr['POS']
                labels = bnd_chr['GENE']
                bnd_track.xticks(
                    label_pos,
                    labels,
                    label_orientation="vertical",
                    outer = False,
                    show_bottom_line=True,
                    label_size=10,
                    line_kws=dict(ec="grey"),
                )

        # Link colour
        def get_link_color(fusion):
            if fusion == '':
                return 'lightgrey'
            elif fusion in unique_mitelman_fusion:
                return 'red'
            else:
                return 'black'

        # Draw links
        bnd_only.apply(lambda row: circos.link(
            (row['CHROM'], row['POS'], row['POS']),
            (row['CHR2'], row['POS2'], row['POS2']),
            color=get_link_color(row['FUSION'])
        ), axis=1)
        
        # Save the figure
        fig = circos.plotfig(dpi=300)
        fig.savefig(sample_svg, dpi=300, bbox_inches='tight')
        fig.savefig(sample_png, dpi=300, bbox_inches='tight')
        
        # Tag the fusion calls with the sample name
        sample_fusion_result = bnd_only[['CHROM', 'POS', 'CHR2', 'POS2', 'GENE', 'FUSION']].copy()
        sample_fusion_result.insert(0, 'SAMPLE_NAME', sample_name)
        
        # Write the per-sample TSV
        sample_fusion_result.to_csv(sample_tsv, sep='\t', index=False)
        
        # Add to the combined result
        all_fusion_results.append(sample_fusion_result)
        
        print(f"Sample {sample_name} - Circos plot saved to {sample_svg}, {sample_png} and fusion table to {sample_tsv}")
    
    # Merge the fusion calls of all samples into one file
    if all_fusion_results:
        combined_results = pd.concat(all_fusion_results, ignore_index=True)
        combined_results.to_csv(args.output_fusion_tsv, sep='\t', index=False)
        print(f"Combined fusion results saved to {args.output_fusion_tsv}")
    else:
        # Write an empty file
        pd.DataFrame(columns=['SAMPLE_NAME', 'CHROM', 'POS', 'CHR2', 'POS2', 'GENE', 'FUSION']).to_csv(
            args.output_fusion_tsv, sep='\t', index=False
        )
        print("No fusion results found for any sample.")

if __name__ == '__main__':
    main()
