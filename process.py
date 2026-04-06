"""
file for processing MAF/txt files
"""
# HAVEN'T TESTED YET
import pandas as pd
import numpy as np
import gzip

# Function to map CN to mutations
def map_cn(maf_df, cnv_df):
    total_cn_list = []

    for idx, mut in maf_df.iterrows():
        chr_ = mut["Chromosome"]
        pos = mut["Start_Position"] # assuming SNVs, start and end are the same
        # assumption: leukemia mutations are mostly SNVs/splice sites

        segments = cnv_df[cnv_df["Chromosome"] == chr_]

        match = segments[
            (segments["Start"] <= pos) &
            (segments["End"] >= pos)
        ]

        if len(match) > 0:
            total_cn = match.iloc[0]["Copy_Number"]
        else:
            total_cn = 2  # default diploid fallback

        total_cn_list.append(total_cn)

    maf_df["total_cn"] = total_cn_list
    return maf_df

# Create files for MyClone/Pyclone input
program = 'myclone' # 'myclone' or 'pyclone'

# load sample sheet and create new dataframes for maf and cnv files
ss = pd.read_csv('gdc_sample_sheet.2026-03-04.tsv', sep='\t')
maf_key = ss[ss['File Name'].str.contains(r'\.maf.gz$', regex=True)].copy()
cnv_key = ss[ss['File Name'].str.contains(r'\.txt$', regex=True)].copy()
# sort by case ID to ensure maf and cnv files are in same order
maf_key = maf_key.sort_values(by='Case ID') 
cnv_key = cnv_key.sort_values(by='Case ID')

raw_data_path = '/Users/biancakolim/Desktop/Academic/Spring2026/FCBB/final/data/raw/'
output_folder = 'myclone-code/myclone_input/'

for maf, cnv in zip(maf_key.itertuples(index=True, name='Pandas'), cnv_key.itertuples(index=True, name='Pandas')): # uses File ID and File Name to find MAF and CNV file
    folder_maf = maf[1] # File ID
    file_name_maf = maf[2] # File Name
    folder_cnv = cnv[1] # File ID
    file_name_cnv = cnv[2] # File Name

    # Load MAF file
    maf_df = pd.read_csv(raw_data_path+folder_maf+'/'+file_name_maf, sep='\t', compression='gzip', comment='#')
    maf_df = maf_df[[
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "t_ref_count",
        "t_alt_count"
        ]].copy()
    maf_df = maf_df.dropna(subset=["t_ref_count", "t_alt_count"]) # Drop rows with missing counts
    # convert to integers
    maf_df["t_ref_count"] = maf_df["t_ref_count"].astype(int)
    maf_df["t_alt_count"] = maf_df["t_alt_count"].astype(int)
    # compute total depth and VAF
    maf_df["total_depth"] = maf_df["t_ref_count"] + maf_df["t_alt_count"]
    maf_df = maf_df[maf_df["total_depth"] > 0] # avoid division by zero
    maf_df["VAF"] = maf_df["t_alt_count"] / maf_df["total_depth"]
    maf_df["mutation_id"] = ( # create mutation ID in format "chr:pos:ref>alt"
        maf_df["Chromosome"].astype(str) + ":" +
        maf_df["Start_Position"].astype(str) + ":" +
        maf_df["Reference_Allele"] + ">" +
        maf_df["Tumor_Seq_Allele2"])
    
    # Load CNV file
    cnv_df = pd.read_csv(raw_data_path+folder_cnv+'/'+file_name_cnv, sep='\t') # need Copy_Number column

    # Sort for efficient matching
    maf_df = maf_df.sort_values(["Chromosome", "Start_Position"])
    cnv_df = cnv_df.sort_values(["Chromosome", "Start"])
    maf_df = map_cn(maf_df, cnv_df)

    # output file name: ###_caseID.tsv
    case_id = maf[5].split(',')[0]
    idx = f"{maf[0]:03}"
    out_nm = f"{idx}_{case_id}.tsv"

    match program:
        case 'myclone':
            print("MyClone processing not implemented yet")
            myclone_df = maf_df[[
                "mutation_id",
                "t_ref_count",
                "t_alt_count",
                "total_cn",
                "VAF" ]].copy()
            # Rename columns to match MyClone spec
            myclone_df.rename(columns={
                "t_ref_count": "ref_counts",
                "t_alt_count": "var_counts"
                }, inplace=True)
            
            myclone_df.to_csv(output_folder + out_nm, sep="\t", index=False)
        case 'pyclone':
            print("Pyclone processing not implemented yet")


# make dictionary with case ID and relapse or non-relapse status