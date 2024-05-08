"""

Description:
------------
A script whichs recopilates all the information obtained from the project. TeloFusVarFinder
All the tables collected from the tables generated in the scripts:
All_protein_class.py or Specific_proteins_class.py (example: zoonomia_proteins_class_specific_mouse.xlsx)
Fusion_analysis.py (example: Fusion_results_coverage.xlsx)
Motif_variant_analysis.py (example: Telomere_variants_catalogue_raw.xlsx)
Extraction_zoonomia_information.py (example: zoonomia_sp_info.xlsx)

Usage:
------

command line:
python FusionVariantProteinInfo.py

Dependencies:
-------------
    - python 3.12.0
    - pandas
    - numpy

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 10-April-2024    

"""

import pandas as pd
import numpy as np

# File paths for data sources
info_file = "zoonomia_sp_info_c.xlsx"
fusion_file = "Final_results_mammals/Fusion_results_coverage.xlsx"
variant_file = "Final_results_mammals/Telomere_variants_catalogue.xlsx"
# Uncomment one of the following lines depending on which protein file you want to use
# protein_file = "zoonomia_all_proteins_class_2.csv"
# protein_file = "zoonomia_all_proteins_class_2_mouse.csv"
# protein_file = "zoonomia_proteins_class_specific.xlsx"
protein_file = "zoonomia_proteins_class_specific_mouse.xlsx"
cancer_file = "SupplementaryData.xls"

# Load species information data from an Excel file
df_info = pd.read_excel(info_file, usecols=["Species", "Order", "Average_bp_read", "Instrument_model", 'Longevity', 'Chrom_num'])
df_protein = pd.read_excel(protein_file)

# Clean up data
df_protein.rename(columns={'Unnamed: 0': 'Species'}, inplace=True)
df_protein['Species'] = df_protein['Species'].str.replace('_', ' ')

# Define a mapping for protein category values as present (true) and absent (False) 
category_mapping = {'I': True, 'PI': True, 'UL': False, 'M': False, 'PM': False, 'L': False, 'PG': False}
df_protein.set_index('Species', inplace=True)
for col in df_protein:
    df_protein[col] = df_protein[col].map(category_mapping)

# Load fusion data from an Excel file and clean it up
df_fusion = pd.read_excel(fusion_file)
df_fusion['Specie_name'] = df_fusion['Specie_name'].str.replace('.csv', '')

# Load and filter variant data from an Excel file 
df_variant = pd.read_excel(variant_file)
df_variant = df_variant[df_variant['Variant'].isin(['Genome_size', 'Total_telomere_reads', 'Porcentage_TTAGGG', 'Porcentage_variability'])]

# Transform variant data to have variants as columns and species as rows
df_variant_t = df_variant.set_index('Variant').T

# Merge fusion data with species information
fusion_info_df = pd.merge(df_fusion, df_info, on='Species', how='left').fillna(np.nan)
print(fusion_info_df)

# Merge fusion data with protein data
fusion_info_protein_df = pd.merge(fusion_info_df, df_protein, on='Species', how='left').fillna(np.nan)
print(fusion_info_protein_df)

# Uncomment the following lines to merge with cancer data (not used in this script)

# df_cancer = pd.read_excel(cancer_file, usecols=["Species", "CMR", "ICM"])
# df_cancer['Species'] = df_cancer['Species'].str.replace('_', ' ')
# fusion_protein_info_cancer_df = pd.merge(fusion_info_protein_df, df_cancer, on='Species', how='left').fillna(np.nan)
# print(fusion_protein_info_cancer_df)

# Set the index of the fusion_info_protein_df to 'Specie_name'
fusion_info_protein_df.set_index('Specie_name', inplace=True)

# Join fusion_info_protein_df with df_variant_t to combine variant data
fusion_protein_info_cancer_variant_df = fusion_info_protein_df.join(df_variant_t, how='outer')

# Fill NaN values with numpy.nan
fusion_protein_info_cancer_variant_df.fillna(np.nan, inplace=True)
fusion_protein_info_cancer_variant_df.reset_index(drop=False, inplace=True)
fusion_protein_info_cancer_variant_df.rename(columns={'index': 'Species_name'}, inplace=True)

# Calculate total telomere size and percentage telomere size
fusion_protein_info_cancer_variant_df['Total_telomere_size'] = (
        fusion_protein_info_cancer_variant_df['Total_telomere_reads'] *
        fusion_protein_info_cancer_variant_df['Average_bp_read'] /
        fusion_protein_info_cancer_variant_df['Genome_coverage']
    )
fusion_protein_info_cancer_variant_df['percentage_telomere_size'] = (
        fusion_protein_info_cancer_variant_df['Total_telomere_size'] /
        fusion_protein_info_cancer_variant_df['Genome_size'] * 100
    )

#merged_excel_filename = 'zoonomia_analysis_data_all_proteins_boolean_coverage_mouse.csv'
merged_excel_filename = 'zoonomia_analysis_data_boolean_coverage_mouse.xlsx'

fusion_protein_info_cancer_variant_df.to_excel(merged_excel_filename, index=False)
print(f"Tabla fusionada guardada en {merged_excel_filename}")
