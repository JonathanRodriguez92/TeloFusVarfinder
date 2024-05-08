"""

Description:
------------
A script whichs recopilates the gene information obtained from TOGA at {specie}_orthologsClassification.tsv and {specie}_loss_summ_data.tsv
this recopilates information and creates a table with the information of specific proteins previously seleceted, in this case:
'ATRX', 'DAXX', 'TERT', 'RTEL1', 'PIF1', "TERF1", "TERF2", "TINF2", "TPP1", 'POT1', "POT1A", 
'POT1B', "DCR1B", "ACD", "WRN", "BLM", "TP53", "SIRT1", "RAD52", "FEN1", "RPA1", 'RPAIN', 
'RPA2', 'ASF1A', 'SETDB1', 'RMI', 'BLM', 'TOP3A', 'PML', 'ZBTB40', 'RAD51'

Usage:
------

command line:
python Specific_protein_class.py

Dependencies:
-------------
    - python 3.12.0
    - pandas

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 29-March-2024    

"""


import pandas as pd
import os

excel_file = "zoonomia_sp_info.xlsx"
df = pd.read_excel(excel_file)
df = df[['Species']]

#tsv_folders = ["Orthologs_human_hg38/"]
tsv_folders = ['Orthologs_mouse_mm10/']

proteins_to_search = ['ATRX', 'DAXX', 'TERT', 'RTEL1', 'PIF1', "TERF1", "TERF2", "TINF2", "TPP1", 'POT1', "POT1A", 'POT1B', "DCR1B", "ACD", "WRN", "BLM", "TP53", "SIRT1", "RAD52", "FEN1", "RPA1", 'RPAIN', 'RPA2', 'ASF1A', 'SETDB1', 'RMI', 'BLM', 'TOP3A', 'PML', 'ZBTB40', 'RAD51']

processed_combinations = set() # This gets unique IDS of the gene information

for tsv_folder in tsv_folders:
    for protein in proteins_to_search:
        for index, row in df.iterrows():
            species = str(row['Species'])
            species = species.replace(" ", "_")
            tsv_file_orth = os.path.join(tsv_folder, f"{species}_orthologsClassification.tsv") 
            tsv_file_loss = os.path.join(tsv_folder, f"{species}_loss_summ_data.tsv")
            
            reference = tsv_folder.split("_")[1]
            if os.path.exists(tsv_file_orth):
                tsv_data_orth = pd.read_csv(tsv_file_orth, sep='\t')
                tsv_data_loss = pd.read_csv(tsv_file_loss, sep='\t')
                # It looks for the protein selected
                filtered_rows_orth = tsv_data_orth[(tsv_data_orth['t_transcript'].str.contains(protein, case=False, regex=True)) |
                                        (tsv_data_orth['q_transcript'].str.contains(protein, case=False, regex=True))]
                # it gets the t_gene id and t_transcript id. 
                if not filtered_rows_orth.empty:
                    for inner_index_1, inner_row_1 in filtered_rows_orth.iterrows():                         
                        t_transcript = inner_row_1['t_transcript'].split(".")[0]
                        protein_1 = inner_row_1['t_transcript'].split(".")[1].upper()
                        if protein == protein_1: # compares if the protein is the one selected.
                            t_gene = inner_row_1['t_gene']
                            print(t_gene)
                            # it creates a column name from t_gene and t_transcript
                            column_name_2 = f'{protein}_g_{t_gene}__{reference}'
                            orthology_class = inner_row_1['orthology_class']
                            column_name_1 = f'{protein}_t_{t_transcript}_{reference}'
                            if column_name_2 not in df.columns:
                                df[column_name_2] = None
                            # t_transcript is avoided
                            #if column_name_1 not in df.columns:
                            #        df[column_name_1] = None
                            #df.loc[index, column_name_1] = orthology_class
                            
                            # it looks for the gene_class found in loss_summ_data.tsv which could be I, PI, UL and so on 
                            # and the added to the value of the column previously named

                            filtered_rows_summ = tsv_data_loss[tsv_data_loss.iloc[:, 1].str.upper().str.contains(t_gene)]
                            
                            if not filtered_rows_summ.empty:
                                for inner_index_2, inner_row_2 in filtered_rows_summ.iterrows(): 
                                    gene_class = inner_row_2.iloc[2]

                            
                            df.loc[index, column_name_2] = gene_class    

                            
                            
                    print(f'{species} {reference} {protein_1} Completed')
print(df) 
output_excel = "zoonomia_proteins_class_specific_mouse.xlsx"
df.to_excel(output_excel, index=False)

print(f"Tabla actualizada guardada en {output_excel}")



