"""

Description:
------------
A script which analyses the data obtained in the TeloFusVarFinder.
This scripts gets the total frequency of the main motif (TTAGGG or CCCTAA) and its variants.
The script generates 3 tables:

Table_{specie}_1.xlsx for R1, Table_{specie}_2.xlsx for R2
which counts all the variants from TTAGGG and CCCTAA and unifies them in a table with the frequency.

Merged_table_all.xlsx
merge all tables in one where the columns are the file name and the values are their frequency

Telomere_variants_catalogue.xlsx
This table adds CCCTAA motif and its variants to its complementary reversed and the corresponded variants.
Variant: this colum represent the type of the variant included the main motif TTAGGG 
{specie}_1 and {specie}_2 represent the column name from the species and respectively its fastQ file (R1 or R2)

Usage:
------

command line:
python Motif_variant_analysis.py

Dependencies:
-------------
    - python 3.12.0
    - pandas
    - biopython


Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 15-february-2024    

"""


import pandas as pd
import os
from Bio.Seq import Seq
pd.options.mode.chained_assignment = None  

# Define a function to process the csv created in TeloFusVarFinder and extract information on telomere motifs
def process_csv_file(file_path):
    df = pd.read_csv(file_path)  
    total_unique_ids = len(set(df['ID']))  # Calculate the total unique IDs in the file
    ttaggg_dict = {}  # Dictionary to store data for 'TTAGGG' motifs
    ccctaa_dict = {}  # Dictionary to store data for 'CCCTAA' motifs

    # Initialize variables for tracking the current ID and motif data
    current_id = None
    current_motif = None
    current_motif_frequency = None
    current_variants = {}

    for index, row in df.iterrows():
        id_ = row['ID']
        telomere_motif = row.get('TelomereMotif')
        variant = row['Variant']
        frequency = row['Frequency']
        
        # If the ID changes, store the current motif and variants in the dictionary
        if pd.notnull(id_):
            if current_id is not None:
                if current_motif == 'TTAGGG':
                    ttaggg_dict[current_id] = {'Motif': current_motif, 'Motif_Frequency': current_motif_frequency, 'Variant_Frequencies': current_variants}
                elif current_motif == 'CCCTAA':
                    ccctaa_dict[current_id] = {'Motif': current_motif, 'Motif_Frequency': current_motif_frequency, 'Variant_Frequencies': current_variants}
            
            current_id = id_
            current_motif = telomere_motif
            current_motif_frequency = frequency
            current_variants = {}
        
        # Store the variant frequencies in the variants dictionary
        if pd.notnull(variant):
            current_variants[variant] = frequency

    # Store the final motif and variant data in the dictionaries
    if current_id is not None:
        if current_motif == 'TTAGGG':
            ttaggg_dict[current_id] = {'Motif': current_motif, 'Motif_Frequency': current_motif_frequency, 'Variant_Frequencies': current_variants}
        elif current_motif == 'CCCTAA':
            ccctaa_dict[current_id] = {'Motif': current_motif, 'Motif_Frequency': current_motif_frequency, 'Variant_Frequencies': current_variants}
    
    return ttaggg_dict, ccctaa_dict, total_unique_ids  # Return the dictionaries and total telomere reads


folder_path = "Telo_motif_raw/"r
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

# Process all telo_mv_{specie}_{number}.csv files
for file in csv_files:
    # Extract the species name from the file name
    specie = os.path.splitext(file)[0]
    specie = specie.split('_')[2:]
    specie = "_".join(specie)
    

    excel_filename = f"Telo_motif_raw/table_{specie}.xlsx"
    
    
    if not os.path.exists(excel_filename):
        file_path = os.path.join(folder_path, file)
        
        print(f"Processing file: {file_path}")
        
        # Process the CSV file and get the data dictionaries and total telomere reads
        ttaggg_dict, ccctaa_dict, total_unique_ids = process_csv_file(file_path)

        
        ttaggg_data = []
        ttaggg_variant_data = []

        # Process data for 'TTAGGG' motifs and its variants
        for key, value in ttaggg_dict.items():
            ttaggg_data.append({'Motif': value['Motif'], 'Motif_Frequency': value['Motif_Frequency']})
            for variant, frequency in value['Variant_Frequencies'].items():
                ttaggg_variant_data.append({'Variant': variant, 'Frequency': frequency})
                
        # Create dataframes for 'TTAGGG' and variants data and calculate frequencies
        if ttaggg_variant_data:
            motif_df_ttaggg = pd.DataFrame(ttaggg_data)
            motif_ttaggg = motif_df_ttaggg['Motif'].unique()
            motif_ttaggg = ', '.join(motif_ttaggg)
            variant_df_t = pd.DataFrame(ttaggg_variant_data)

            motif_frequency_total_t = motif_df_ttaggg['Motif_Frequency'].sum()
            variant_frequency_total_t = variant_df_t.groupby('Variant').sum().reset_index()

            motif_total_df_t = pd.DataFrame({'Variant': [motif_ttaggg], 'Frequency': [motif_frequency_total_t]})

            variant_frequency_total_t = pd.concat([motif_total_df_t, variant_frequency_total_t]).reset_index(drop=True)
        else:
            variant_frequency_total_t = pd.DataFrame({'Variant': ['TTAGGG'], 'Frequency': [0]})

        # Process data for 'CCCTAA' motifs and its variants
        ccctaa_data = []
        ccctaa_variant_data = []

        for key, value in ccctaa_dict.items():
            ccctaa_data.append({'Motif': value['Motif'], 'Motif_Frequency': value['Motif_Frequency']})
            for variant, frequency in value['Variant_Frequencies'].items():
                ccctaa_variant_data.append({'Variant': variant, 'Frequency': frequency})
        # Create dataframes for 'TTAGGG' and variants data and calculate frequencies
        if ccctaa_variant_data:
            motif_df_ccctaa = pd.DataFrame(ccctaa_data)
            motif_ccctaa = motif_df_ccctaa['Motif'].unique()
            motif_ccctaa = ', '.join(motif_ccctaa)
            variant_df_c = pd.DataFrame(ccctaa_variant_data)

            motif_frequency_total_c = motif_df_ccctaa['Motif_Frequency'].sum()
            variant_frequency_total_c = variant_df_c.groupby('Variant').sum().reset_index()

            motif_total_df_c = pd.DataFrame({'Variant': [motif_ccctaa], 'Frequency': [motif_frequency_total_c]})

            variant_frequency_total_c = pd.concat([motif_total_df_c, variant_frequency_total_c]).reset_index(drop=True)
        else:
            variant_frequency_total_c = pd.DataFrame({'Variant': ['CCCTAA'], 'Frequency': [0]})

        # Combine the frequency totals from 'TTAGGG' and 'CCCTAA' data
        variant_frequency_total = pd.concat([variant_frequency_total_t, variant_frequency_total_c]).reset_index(drop=True)

        # Add the total telomere reads as a final row in the dataframe
        total_unique_ids_df = pd.DataFrame({'Variant': ['Total reads'], 'Frequency': [total_unique_ids]})
        variant_frequency_total = pd.concat([variant_frequency_total, total_unique_ids_df], ignore_index=True)
        
        # Save the processed data 
        variant_frequency_total.to_excel(excel_filename, index=False)
        print(f"Table saved in {excel_filename}")

# Merge all Excel files in the specified folder
folder_path = "Telo_motif_raw"
excel_files = [file for file in os.listdir(folder_path) if file.endswith('.xlsx')]


dfs = []


for file in excel_files:
    file_path = os.path.join(folder_path, file)
    df = pd.read_excel(file_path)
    df['Archivo'] = os.path.splitext(file)[0]  # Add a column with the file name
    dfs.append(df)

# Combine all dataframes into one
merged_df = pd.concat(dfs)

# Remove rows with variants containing 'N' and pivot the table.
merged_df = merged_df[~merged_df['Variant'].str.contains('N')]
merged_df = merged_df.pivot_table(index='Variant', columns='Archivo', values='Frequency', fill_value=0, aggfunc=sum)
merged_df = merged_df.reset_index(drop=False)
merged_excel_filename = 'merged_table_all.xlsx'

# Save the merged dataframe
merged_df.to_excel(merged_excel_filename, index=False)
print(f"Merged table saved in {merged_excel_filename}")

# Function to compare two sequences and return True if they are similar within a specified error margin
def are_similar(sequence1, sequence2, max_error):
    differences = sum(c1 != c2 for c1, c2 in zip(sequence1, sequence2))
    return differences <= max_error


reference_sequence = 'TTAGGG'
df = pd.read_excel('merged_table_all.xlsx')

# Initialize a list for storing similarity scores
similarity_scores = []
for variant in df['Variant']:
    similarity = are_similar(variant, reference_sequence, 1)  # Check similarity within a max error of 1
    similarity_scores.append(similarity)

# Add the similarity scores as a new column in the dataframe
df['Similarity'] = similarity_scores

# Sort the dataframe by similarity and reset the index
df_sorted = df.sort_values(by='Similarity', ascending=False)
df_sorted = df_sorted.reindex(sorted(df_sorted.columns), axis=1)
df_sorted = df_sorted.reset_index(drop=True)

# Remove the similarity column and create a new dataframe without it
df_no_similarity = df_sorted.drop('Similarity', axis=1)
total_reads = {}

# Extract total reads from each column in the dataframe
for column in df_no_similarity.columns:
    value = df_no_similarity[column].iloc[-1]
    
    total_reads[column] = value

# Separate the data into similar and dissimilar rows
true_rows = df_sorted[df_sorted['Similarity'] == True]
true_rows = true_rows.drop('Similarity', axis=1)

false_rows = df_sorted[df_sorted['Similarity'] == False]

# Reverse complement false rows to account for different variants
false_rows_2 = false_rows.copy()

false_rows_2['Variant_r'] = false_rows_2['Variant'].apply(lambda x: str(Seq(x).reverse_complement()))
false_rows_2.rename(columns={'Variant': 'Variant_s'}, inplace=True)
false_rows_2.rename(columns={'Variant_r': 'Variant'}, inplace=True)
false_rows_2 = false_rows_2.drop('Similarity', axis=1)
false_rows_2 = false_rows_2.drop('Variant_s', axis=1)

# Handle specific variants and remove total reads from the false rows
false_rows_2.loc[false_rows_2['Variant'] == 'shtey ltaoA', 'Variant'] = 'Total reads'
false_rows_2 = false_rows_2[false_rows_2['Variant'] != 'Total reads']

# Merge the true rows and reverse complemented false rows on the variant column
merged_df = pd.merge(true_rows, false_rows_2, on='Variant', suffixes=('_TTAGGG', '_CCCTAA'))

# Retrieve list of column names
index_merged_df = list(merged_df.columns)

# Rename columns and set the index to 'Variant'
merged_df = merged_df.rename(columns=lambda x: x.replace('table_', ''))
merged_df = merged_df.rename(columns=lambda x: x.replace('_CCCTAA', ''))
merged_df = merged_df.rename(columns=lambda x: x.replace('_TTAGGG', ''))
print(merged_df['Ochotona_princeps_2'])

merged_df.set_index('Variant', inplace=True)

# Extract species from the column names
species = set()
for index in index_merged_df:
    specie = index.split("_")[1:]
    specie = "_".join(specie)
    specie = specie.split("_")[:-1]
    specie = "_".join(specie)
    species.add(specie)

# Calculate sums of columns for each species
sumas_por_especie = {}
for specie_1 in species:
    columnas_de_especie = [col for col in merged_df.columns if specie_1.lower() in col.lower()]
    
    sumas_por_especie[specie_1] = merged_df[columnas_de_especie].sum(axis=1)

# Create a DataFrame for the sums by specie
df_sumas = pd.DataFrame(sumas_por_especie)
df_sumas = df_sumas.sort_index(axis=1)
df_sumas.reset_index(inplace=True)

# Initialize dictionaries for total spots run, bases run, and genome size
total_spots_run = {}
total_bases_run = {}
genome_size = {}

# Read species information from an Excel file
df_sp_info = pd.read_excel("zoonomia_sp_info_c.xlsx")
for index, row in df_sp_info.iterrows():
    species = row['Species'].replace(" ", "_")
    
    total_spots_run[species] = row['Total_spots_run']
    total_bases_run[species] = row['Total_bases_run']
    genome_size[species] = row['Genome_size']

# Duplicate the data for total spots run, bases run, and genome size
total_spots_run_duplicated = {}
total_bases_run_duplicated = {}
genome_size_duplicated = {}

for species, value in total_spots_run.items():
    total_spots_run_duplicated[f"{species}_1"] = value
    total_spots_run_duplicated[f"{species}_2"] = value

for species, value in total_bases_run.items():
    total_bases_run_duplicated[f"{species}_1"] = value
    total_bases_run_duplicated[f"{species}_2"] = value

for species, value in genome_size.items():
    genome_size_duplicated[f"{species}_1"] = value
    genome_size_duplicated[f"{species}_2"] = value

# Set the display format for float numbers to 1 decimal place
pd.set_option('display.float_format', lambda x: '%.1f' % x)

# Create a DataFrame for total reads
df_total_reads = pd.DataFrame([total_reads])
df_total_reads = df_total_reads.rename(columns=lambda x: x.replace('table_', ''))

# Create DataFrames for total spots run, total bases, and genome size
df_total_spots_run = pd.DataFrame([total_spots_run_duplicated])
df_total_bases = pd.DataFrame([total_bases_run_duplicated])
df_genome_size = pd.DataFrame([genome_size_duplicated])

# Concatenate the sum DataFrame with the total reads DataFrame
df_variants_catalogue = pd.concat([df_sumas, df_total_reads])
df_variants_catalogue.set_index('Variant', inplace=True)
if 'TTAGGG' in df_variants_catalogue.index:
    df_variants_catalogue = df_variants_catalogue.reindex(['TTAGGG'] + 
                                                          [idx for idx in df_variants_catalogue.index if idx != 'TTAGGG'])
df_variants_catalogue = df_variants_catalogue.rename(index={'Total reads': 'Total_telomere_reads'})
df_variants_catalogue = df_variants_catalogue.drop(df_variants_catalogue.columns[0], axis=1)
df_variants_catalogue = df_variants_catalogue.reset_index(drop=False)

# Save the DataFrame
df_variants_catalogue.to_excel("Telomere_variants_catalogue_raw.xlsx", index=False)

df_variants_catalogue.set_index('Variant', inplace=True)
# Concatenate the catalogue DataFrame with total spots run, bases run, and genome size data
df_variants_catalogue = pd.concat([df_variants_catalogue, df_total_spots_run], join='inner')
df_variants_catalogue.index = df_variants_catalogue.index[:-1].tolist() + ['Total_spots_run']
df_variants_catalogue = pd.concat([df_variants_catalogue, df_total_bases], join='inner')
df_variants_catalogue.index = df_variants_catalogue.index[:-1].tolist() + ['Total_bases_run']
df_variants_catalogue = pd.concat([df_variants_catalogue, df_genome_size], join='inner')
df_variants_catalogue.index = df_variants_catalogue.index[:-1].tolist() + ['Genome_size']

# Calculate total frequencies and percentages for TTAGGG and variability
Total_frecuency = df_variants_catalogue.iloc[0:19, :].sum()
Total_frecuency_variant = df_variants_catalogue.iloc[1:19, :].sum()
TTAGGG_frecuency = df_variants_catalogue.iloc[0, :]
porcentage_TTAGGG = TTAGGG_frecuency  / Total_frecuency
porcentage_variability = Total_frecuency_variant  / Total_frecuency

# Add percentages to the DataFrame
df_variants_catalogue.loc['Porcentage_TTAGGG'] = round(porcentage_TTAGGG, 2)
df_variants_catalogue.loc['Porcentage_variability'] = round(porcentage_variability, 2)


df_variants_catalogue = df_variants_catalogue.reset_index(drop=False)
df_variants_catalogue = df_variants_catalogue.rename(columns={'index': 'Variant'})

# Save the final catalogue
df_variants_catalogue.to_excel('Telomere_variants_catalogue.xlsx', index=False)
print(f"Telomere_variants_catalogue.xlsx created")
