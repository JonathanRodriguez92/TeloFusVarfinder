"""

Description:
------------
Analysing and processing SRR data fastq files
This script looks for the frequency of a motif and its variability in the genome with this parametres:
-error. This represents how many mismatches the motif could have.
-Original motif. This represents the pattern to analyse. In this case TTAGGG or CCCTAA but could be anyone.
-percentage_occurrence. Percentage of frequency of the desired motif in the read. In this case 40 was selected to determine if the sequence is telomeric or not.
-max_coord_diff an extra parameter not currently in use, is a paramater to defind how long a motif can be between one and another.
 for example: TTAGGGCCTTAGGGTTAGGG if the max_coord_diff is 2 the length would be 20 but if is one the max lenght would be 12.
 This is great to use in assembled genomes to find the telomere part.
Usage:
------

command line:
python Telo_motif_process.py

Dependencies:
-------------
    - python 3.12.0
    - pandas
    - biopython
    - numpy
    - csv

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 10-January-2024    

"""


from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import pandas as pd
import os
import csv
import time

class TelomericVariantSearcher:

    def __init__(self, specie, data_dir, srr_number, seq_type, original_motifs=None, error=1, max_coord_diff=3):
        # Variables needed to run the script
        self.original_motifs = original_motifs # A motif is a pattern which can apear in tandem repeats in the genome
        self.seq_type = seq_type # Paired or single
        self.data_dir = data_dir # where fastq file is
        self.srr_number = srr_number # SRR number to analyse
        self.specie = specie # Name of the specie
        self.error = error # max error permited
        self.max_coord_diff = max_coord_diff 
        self.threshold = 40 # This is the minimum porcentage of motif frecuence in the sequence.
    # A function to find the similarity between the original motif and the variant
    def are_similar(self, sequence1, sequence2, max_error):
        differences = sum(c1 != c2 for c1, c2 in zip(sequence1, sequence2))
        return differences <= max_error
        
    # Main fuction which searchs the frequency of the motif in a read and its own variants.
    def search_telomeric_variants(self, fastq_file, min_frequency, number):
        with open(f'telomeric_mv_raw_{self.specie}_{number}.csv', 'w', newline='') as csv_file:
            fieldnames = ['ID', 'Variant', 'Frequency', 'Length', 'MaxTandemRepeatsLength', 'PorcentageTelomereSequence']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            # open the fastq_file
            for records in SeqIO.parse(fastq_file, "fastq"):
                sequence = str(records.seq).upper()
                telomeric_variants = [] 

                for original_motif in self.original_motifs:
                    occurrences = sequence.count(original_motif) # this counts the frequency of the original motif from the read.
                    percentage_occurrence = (occurrences * len(original_motif)) * 100 / len(sequence) # it calculates the percentage depending of the lenght of the read.
                    if percentage_occurrence >= 40:  
                        # A deeper analysis of the sequence.                
                        variant_length = len(original_motif)
                        for i in range(len(sequence) - variant_length + 1):
                            candidate_motif = str(sequence[i:i+variant_length].upper())
                            if self.are_similar(original_motif, candidate_motif, self.error):
                                telomeric_variants.append((candidate_motif, i + 1, i + variant_length)) # this records the composition of the motif/variant, its location and its length.

                            # This part of the analysis looks for extra motifs.
                            #if original_motif == 'TTAGGG':
                                #if self.are_similar(original_motif, candidate_motif, self.error):
                                    #first_three_nucleotides = candidate_motif[:3]
                                    #last_three_nucleotides = candidate_motif[-3:]
                                    #if all(nucleotide in 'ACTG' for nucleotide in first_three_nucleotides) and last_three_nucleotides == 'GGG':
                            #elif original_motif == 'CCCATT':
                                #if self.are_similar(original_motif, candidate_motif, self.error):
                                    #first_three_nucleotides = candidate_motif[:3]
                                    #last_three_nucleotides = candidate_motif[-3:]
                                    #if all(nucleotide in 'ACTG' for nucleotide in last_three_nucleotides) and first_three_nucleotides == 'CCC':
                                        #telomeric_variants.append((candidate_motif, i + 1, i + variant_length))

                        # This part calculates the frequency of the variant/motif
                        variant_frequency = {}
                        for variant, start, end in telomeric_variants:
                            variant_frequency[variant] = variant_frequency.get(variant, {"frequency": 0, "coordinates": []})
                            # This calculates if the variant/motif detected is overlaping with other variant/motif and does not integrate it as an extra frequence
                            is_superposed = any(abs(start - existing_end) == 0 or abs(existing_start - end) == 0 for existing_start, existing_end in variant_frequency[variant]["coordinates"])
                            
                            if not is_superposed:
                                variant_frequency[variant]["frequency"] += 1

                            # this calculates the max lenght of the tamdem repeat. 
                            #    merged_coordinates = []
                            #    for existing_start, existing_end in variant_frequency[variant]["coordinates"]:
                                #    if not (start <= existing_start and end >= existing_end) and not (start >= existing_start and end <= existing_end):
                                #     if abs(start - existing_end) <= self.max_coord_diff or abs(existing_start - end) <= self.max_coord_diff:
                                #         start = min(start, existing_start)
                                #         end = max(end, existing_end)
                                #     else:
                                #         merged_coordinates.append((existing_start, existing_end))
                                # else:
                                #     merged_coordinates.append((existing_start, existing_end))

                                #merged_coordinates.append((start, end))
                                #variant_frequency[variant]["coordinates"] = merged_coordinates
                        # The data is stored in the csv.
                        if variant_frequency:
                            for variant, variant_info in variant_frequency.items():
                                frequency = variant_info["frequency"]
                                # max_tandem_length = max(end - start + 1 for start, end in variant_info["coordinates"])
                                porcentage_seq = round((frequency * len(variant)) * 100 / len(sequence), 2) # This calculates the percentage of the variant/motif in the current read
                                
                                if frequency >= min_frequency:
                                    writer.writerow({'ID': records.id, 
                                                    'Variant': variant, 
                                                    'Frequency': frequency,
                                                    'Length': len(sequence), # 'MaxTandemRepeatsLength': max_tandem_length,
                                                    'PorcentageTelomereSequence': porcentage_seq})
        time.sleep(60)                          
        print(f'telomeric_mv_raw_{self.specie}_{number}.csv created')
        
        # A second part to get a better and more visible table with the data obtained.
        df = pd.read_csv(f'telomeric_mv_raw_{self.specie}_{number}.csv')
        
        df = df.sort_values(by='PorcentageTelomereSequence', ascending=False) # data is sort by the percentage telomere sequence to get only the most frequenced (motif)
        df = df[df['PorcentageTelomereSequence'] >= self.threshold] # data is stored with a maximun threshold previously selected
        
        df = df.sort_values(by=['ID', 'Frequency'], ascending=[True, False]) 
        result = df.groupby('ID', sort=False).first().reset_index() # All Ids are stored to futher analysis
        df2 = pd.read_csv(f'telomeric_mv_raw_{self.specie}_{number}.csv')
        result2 = pd.merge(df2, result[['ID', 'Variant']], on='ID', how='inner') # Both dfs are merged to get one with variant and motif separated
        result2['Variant_x'] = np.where(result2['Variant_x'] != result2['Variant_y'], result2['Variant_x'], '') 
        result2 = result2.rename(columns={'Variant_x': 'Variant', 'Variant_y': 'TelomereMotif'}) # Names are changed to get a clearer perspective
        result2.sort_values(by='ID', ascending=True)
        result2 = result2[['ID', 'TelomereMotif', 'Variant', 'Frequency','PorcentageTelomereSequence', 'Length']] 
        result2['ID_read'] = result2['ID'].str.split('.').str[1].astype(int) # Organisation of ID
        result2 = result2.sort_values(by=['ID_read', 'Frequency'], ascending=[True, False])
        result2 = result2.drop('ID_read', axis = 1)

        duplicates_mask = result2.duplicated(subset=['ID'])
        result2.loc[duplicates_mask, ['ID','TelomereMotif']] = np.nan
        os.makedirs('Telomeric_motifs_variants_raw', exist_ok=True)
        result2.to_csv(f'Telomeric_motifs_variants_raw/telo_mv_{self.specie}_{number}.csv')
        
        print(f'telo_mv_{self.specie}_{number}.csv created')
        os.remove(f'telomeric_mv_raw_{self.specie}_{number}.csv') 