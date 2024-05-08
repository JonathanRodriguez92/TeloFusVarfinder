"""

Description:
------------
A script whichs recopilates all the gene information obtained from TOGA at {specie}_orthologsClassification.tsv and {specie}_loss_summ_data.tsv
this recopilates information and creates a table with the information of all genes with the protein name structured like {protein}_g_{t_gene}_{reference}

Usage:
------

command line:
python All_proteins_class.py

Dependencies:
-------------
    - python 3.12.0
    - pandas

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 9-April-2024    

"""


import pandas as pd
import os

# Variables definition
excel_file = "zoonomia_sp_info_c.xlsx"
df = pd.read_excel(excel_file)
#tsv_folders = ["Orthologs_human_hg38/"]
tsv_folders = ['Orthologs_mouse_mm10/']
result_df = pd.DataFrame()
output_excel = "zoonomia_all_proteins_class_2_mouse.csv"
classification_file = "Orthologs_mouse_mm10//Rhinopithecus_roxellana_orthologsClassification.tsv"
classification_data = pd.read_csv(classification_file, sep='\t')

# Dictionary with all the information founde in the classification file which is the t_gene and the protein associated to the gene.
transcript_protein_dict = dict(zip(classification_data['t_gene'], classification_data['t_transcript'].str.split('.').str[1]))

# Process TSV folders
for tsv_folder in tsv_folders:
    for file_name in os.listdir(tsv_folder):
        if file_name.endswith("_loss_summ_data.tsv"):
            
            tsv_file_loss = os.path.join(tsv_folder, file_name)
            # it extracts the specie from the foulder name
            specie = file_name.split('.')[0].split("_")[:-3]
            specie = "_".join(specie)
            print(specie)
            reference = tsv_folder.split("_")[1]
            tsv_data_loss = pd.read_csv(tsv_file_loss, sep='\t')
            for index, row in df.iterrows():
                species_name = str(row['Species']).replace(" ", "_")
            # Filters the rows called 'GENE'
                if species_name == specie:

                    gene_data = tsv_data_loss[tsv_data_loss.iloc[:, 0] == 'GENE']
                    
                    
                    for idx, row in gene_data.iterrows():
                        t_gene = row.iloc[1]  # Gene ID collection: (ENSG0000***)
                        
                        value = row.iloc[2]  # Gets the value next to the gene ID
                        for gene, protein in transcript_protein_dict.items():

                            if gene == t_gene: # verify the gene and t_gene are the same in the classification_file
                        # # it creates a column name from t_gene and the protein founded in the data classification_file
                                column_name_1 = f'{protein}_g_{t_gene}_{reference}'

                                
                                # Creates a column if it doesnt exists
                                if column_name_1 not in result_df.columns:
                                    result_df[column_name_1] = pd.Series(dtype='object')
                                # Assigns the value to the column
                                result_df.loc[species_name, column_name_1] = value
print(result_df)
result_df.to_csv(output_excel)