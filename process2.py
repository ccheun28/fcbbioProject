"""
file for processing MAF/txt files
"""
# HAVEN'T TESTED YET
# copy of process.py

import pandas as pd
import numpy as np
import gzip
import os

# Function to map CN to mutations
def map_cn(maf_df, cnv_df):
    total_cn_list = []
    major_cn_list = []
    minor_cn_list = []
    normal_cn_list = []

    # Infer sex from CNV file
    chromosomes_in_file = set(cnv_df["Chromosome"])
    if "chrY" in chromosomes_in_file:
        sex = "Male"
    else:
        sex = "Female"

    for idx, mut in maf_df.iterrows():
        chr_ = mut["Chromosome"]
        pos = mut["Start_Position"] # assuming SNVs, start and end are the same
        # assumption: leukemia mutations are mostly SNVs/splice sites

        # Determine normal copy number based on chromosome type
        if chr_ == "chrX":
            normal_cn = 1 if sex == "Male" else 2
        elif chr_ == "chrY":
            normal_cn = 1 if sex == "Male" else 0
        else:
            normal_cn = 2
        normal_cn_list.append(normal_cn)

        segments = cnv_df[cnv_df["Chromosome"] == chr_]

        match = segments[
            (segments["Start"] <= pos) &
            (segments["End"] >= pos)
        ]

        if len(match) > 0:
            total_cn = match.iloc[0]["Copy_Number"]
            minor_cn = match.iloc[0]["Minor_Copy_Number"]
            major_cn = match.iloc[0]["Major_Copy_Number"]
        else:
            total_cn = 2  # default diploid fallback
            minor_cn = 1
            major_cn = 1

        total_cn_list.append(total_cn)
        minor_cn_list.append(minor_cn)
        major_cn_list.append(major_cn)

    maf_df["total_cn"] = total_cn_list
    maf_df["major_cn"] = major_cn_list
    maf_df["minor_cn"] = minor_cn_list
    maf_df["normal_cn"] = normal_cn_list
    return maf_df

# Create files for MyClone/Pyclone input
program = 'pyclone' # 'myclone' or 'pyclone'

# load sample sheet and create new dataframes for maf and cnv files
ss = pd.read_csv('gdc_sample_sheet.2026-04-01.tsv', sep='\t') #CHANGED
maf_key = ss[ss['File Name'].str.contains(r'\.maf.gz$', regex=True)].copy()
cnv_key = ss[ss['File Name'].str.contains(r'\.copy_number_variation.seg.txt$', regex=True)].copy()
# sort by case ID to ensure maf and cnv files are in same order
maf_key = maf_key.sort_values(by='Case ID') 
cnv_key = cnv_key.sort_values(by='Case ID')

#print(cnv_key.head(10))
maf_cases = set(maf_key["Case ID"])
cnv_cases = set(cnv_key["Case ID"])


#check for samples without MAF and txt
common_cases = maf_cases & cnv_cases
missing_maf = cnv_cases - maf_cases
missing_cnv = maf_cases - cnv_cases

"""
print("Cases missing MAF files:")
print(missing_maf)

print("\nCases missing CNV files:")
print(missing_cnv)
"""

#only keep samples with both MAF and cnv files
maf_key = maf_key[maf_key["Case ID"].isin(common_cases)].copy()
cnv_key = cnv_key[cnv_key["Case ID"].isin(common_cases)].copy()
maf_key = maf_key.sort_values(by='Case ID')
cnv_key = cnv_key.sort_values(by='Case ID')

raw_data_path = '/Users/charlottecheung/Developer/FCBBioInfo/FCBBiodata/'
output_folder = '/Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_input/' #CHANGED

for maf, cnv in zip(maf_key.itertuples(index=True, name='Pandas'), cnv_key.itertuples(index=True, name='Pandas')): # uses File ID and File Name to find MAF and CNV file
    folder_maf = maf[1] # File ID
    file_name_maf = maf[2] # File Name
    folder_cnv = cnv[1] # File ID
    file_name_cnv = cnv[2] # File Name

    maf_path = os.path.join(raw_data_path, folder_maf, file_name_maf)
    cnv_path = os.path.join(raw_data_path, folder_cnv, file_name_cnv)

    # Check if files exist
    if not os.path.isfile(maf_path):
        print(f"MAF file not found: {maf_path}, skipping sample")
        continue
    if not os.path.isfile(cnv_path):
        print(f"CNV file not found: {cnv_path}, skipping sample")
        continue

    # Load MAF file
    try:
        maf_df = pd.read_csv(raw_data_path+folder_maf+'/'+file_name_maf, sep='\t', compression='gzip', comment='#')
        maf_df["Chromosome"] = maf_df["Chromosome"].astype(str).str.strip()
    except Exception as e:
        print(f"Error reading MAF file {maf_path}: {e}, skipping sample")
        continue

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
    try:
        cnv_df = pd.read_csv(raw_data_path+folder_cnv+'/'+file_name_cnv, sep='\t') # need Copy_Number column
        cnv_df["Chromosome"] = cnv_df["Chromosome"].astype(str).str.strip()
    except Exception as e:
        #print(f"Error reading CNV file {cnv_path}: {e}, skipping sample")
        continue

    # Sort for efficient matching
    maf_df = maf_df.sort_values(["Chromosome", "Start_Position"])
    cnv_df = cnv_df.sort_values(["Chromosome", "Start"])
    maf_df = map_cn(maf_df, cnv_df)

    # output file name: ###_caseID.tsv
    case_id = maf[6].split(',')[0]
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
            pyclone_df = maf_df[[
                "mutation_id",
                "t_ref_count",
                "t_alt_count",
                "total_cn",
                "major_cn",
                "minor_cn",
                "normal_cn"
                 ]].copy()
            # Rename columns to match MyClone spec
            pyclone_df.rename(columns={
                "t_ref_count": "ref_counts",
                "t_alt_count": "var_counts"
                }, inplace=True)
            
            pyclone_df.to_csv(output_folder + out_nm, sep="\t", index=False)
            

# make dictionary with case ID and relapse or non-relapse status