"""

Description:
------------
A script which analyses the data obtained in the TeloFusVarFinder.
This script unifies all the fusions from each specie and creates two tables one with the raw data Fusion_results_raw.xlsx 
and the other table with the data normilezed by genome coverage Fusion_results_coverage.xlsx.
This table is devided into:
Inward (distance 0) Fusions of forward + reverse with distance 0 (eg: TTAGGGCCCTAA)
Otward (distance 0) Fusions of reverse + forward with distance 0 (eg: CCCTAATTAGGG)
Inward (other) Fusions of forward + reverse with distance > 0 (eg: TTAGGGxxxCCCTAA)
Otward (other) Fusions of reverse + forward with distance > 0 (eg: CCCTAAxxxTTAGGG)
Total_fusions_0 Total inward and outward fusions with distance 0 
Total_fusions Total of all fusions inward + outward 
Total_inward Total fusions inward 0 + other
Total_outward Total fusions outward 0 + other

Usage:
------

command line:
python Fusion_analysis.py

Dependencies:
-------------
    - python 3.12.0
    - pandas

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 10-february-2024    

"""



import os
import pandas as pd

# Set the directory paths for the data files
directory = '/root/TFM/GCA_all/Final_results_mammals'
directory_fasta = '/root/TFM/GCA_all/Final_results_mammals/Telo_fus_raw'

# Function to count the number of inward and outward fusions
def count_fusions(df):
    inward_distance_0 = df[(df['FusionType'] == 'Inward') & (df['Distance'] == 0)].shape[0]
    outward_distance_0 = df[(df['FusionType'] == 'Outward') & (df['Distance'] == 0)].shape[0]
    inward_other = df[(df['FusionType'] == 'Inward') & (df['Distance'] != 0)].shape[0]
    outward_other = df[(df['FusionType'] == 'Outward') & (df['Distance'] != 0)].shape[0]
    return inward_distance_0, outward_distance_0, inward_other, outward_other

results = []

# Loop through all files in the specified directory
for filename in os.listdir(directory):
    if filename.startswith("combined_distances_") and filename.endswith(".csv"):
        filepath = os.path.join(directory, filename)
        
        # Extract the species name from the file name
        specie_name = filename.split('.')[0]
        specie_name = specie_name.split("_")[2:]
        specie_name = "_".join(specie_name)
        species = filename.split('.')[0]
        species = species.split("_")[:-1]
        species = "_".join(species)
        species = species.split("_")[2:]
        species = "_".join(species)
        species = species.replace("_", " ")

        # Load the data file into a data frame.
        df = pd.read_csv(filepath)
        
        # Count inward and outward fusions at distance 0 and other distances
        inward_distance_0, outward_distance_0, inward_other, outward_other = count_fusions(df)
        
        # Append the results for this species to the list
        results.append((species, specie_name, inward_distance_0, outward_distance_0, inward_other, outward_other))

# Convert the results list into a pandas DataFrame
df_results = pd.DataFrame(results, columns=['Species', 'Specie_name', 'Inward (Distance 0)', 'Outward (Distance 0)', 'Inward (Other)', 'Outward (Other)'])


df_info = pd.read_excel("zoonomia_sp_info_c.xlsx")
df_results_2 = df_results.copy()

# This interates from the data frame to get the values
for index, row in df_results.iterrows():
    species_in_file = row['Species']

    if species_in_file in df_info['Species'].values:

        
        for col in ['Inward (Distance 0)', 'Outward (Distance 0)', 'Inward (Other)', 'Outward (Other)']:
            # Calculate genome coverage for the species
            genome_coverage = df_info[df_info['Species'] == species_in_file]['Total_spots_run'].iloc[0] * df_info[df_info['Species'] == species_in_file]['Average_bp_read'].iloc[0] / df_info[df_info['Species'] == species_in_file]['Genome_size'].iloc[0]
            genome_coverage = round(genome_coverage, 2)
            
            # Calculate the raww fusion count and coverage values
            coverage_value = round(int((df_results.loc[index, col] / genome_coverage) * 100), 2)
            raw_value = df_results.loc[index, col]
            # data frame with raw values

            df_results_2['Total_fusions'] = df_results_2['Inward (Distance 0)'] + df_results_2['Outward (Distance 0)'] + df_results_2['Inward (Other)'] + df_results_2['Outward (Other)']
            df_results_2['Total_fusions_0'] = df_results_2['Inward (Distance 0)'] + df_results_2['Outward (Distance 0)']
            df_results_2['Total_inward'] = df_results_2['Inward (Distance 0)'] + df_results_2['Inward (Other)']
            df_results_2['Total_outward'] = df_results_2['Outward (Other)'] + df_results_2['Outward (Distance 0)']

            # Update df_results with coverage and raw values
            df_results.loc[index, col] = float(coverage_value) if not pd.isnull(coverage_value) else 0
            df_results_2.loc[index, col] = int(raw_value) if not pd.isnull(raw_value) else 0

            # dataframe with coverage values
            df_results['Total_fusions'] = df_results['Inward (Distance 0)'] + df_results['Outward (Distance 0)'] + df_results['Inward (Other)'] + df_results['Outward (Other)']
            df_results['Total_fusions_0'] = df_results['Inward (Distance 0)'] + df_results['Outward (Distance 0)']
            df_results['Total_inward'] = df_results['Inward (Distance 0)'] + df_results['Inward (Other)']
            df_results['Total_outward'] = df_results['Outward (Other)'] + df_results['Outward (Distance 0)']
             
        # Store the genome coverage value for the species in both DataFrames
        df_results_2.loc[index, 'Genome_coverage'] = genome_coverage
        df_results.loc[index, 'Genome_coverage'] = genome_coverage

    else:
        print(f"The species {species_in_file} was not found in df_info")

# Extract directory name from the directory telo_fus_raw where fasta files are stored
# This part obtains the data which has 0 fusions and add it to both dataframes
for filename in os.listdir(directory_fasta):
    # Extract the species name from the file name
    specie_name_fasta = filename.split('.')[0].split("_")[2:]
    specie_name_fasta = "_".join(specie_name_fasta)
    specie_fasta = filename.split('.')[0].split("_")[2:-1]
    specie_fasta = "_".join(specie_fasta).replace("_", " ")

    # If the species name from the fasta file is not in df_results_2, add a new row with zero values
    if specie_name_fasta not in df_results_2['Specie_name'].values:
        genome_coverage = df_info[df_info['Species'] == specie_fasta]['Total_spots_run'].iloc[0] * df_info[df_info['Species'] == specie_fasta]['Average_bp_read'].iloc[0] / df_info[df_info['Species'] == specie_fasta]['Genome_size'].iloc[0]
        genome_coverage = round(genome_coverage, 2)
        new_row = {
            'Species': specie_fasta,
            'Specie_name': specie_name_fasta,
            'Inward (Distance 0)': 0,
            'Outward (Distance 0)': 0,
            'Inward (Other)': 0,
            'Outward (Other)': 0,
            'Total_fusions': 0,
            'Total_fusions_0': 0,
            'Total_inward': 0,
            'Total_outward': 0,
            'Genome_coverage': genome_coverage
        }

        # Add the new row to df_results_2 and df_results
        df_results_2.loc[len(df_results_2)] = new_row
        df_results.loc[len(df_results_2)] = new_row
        
# Define the file names for the output Excel files
output_excel_file_raw = "Fusion_results_raw.xlsx"
output_excel_file_coverage = "Fusion_results_coverage.xlsx"

