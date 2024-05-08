"""

Description:
------------
Analysing and processing SRR data fastq files (SINGLE OR PAIRED).
This script looks for telomere fusions in the genome with this parametres:
-Mismatches. This represents how many mismatches depending of the score
-Pattern. This represents the pattern to analyse in this case TTAGGG or CCCTAA but could be anyone
-Max distance. The distance between the two parts of the fusion

Usage:
------

command line:
python Telo_fusion_process.py

Dependencies:
-------------
    - python 3.12.0
    - pandas
    - biopython
    - numpy

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 12-February-2024    

"""


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align 
import pandas as pd
import numpy as np
import os


class TeloFusProcessor:
    def __init__(self, srr_number, seq_type, data_dir, specie):
        # Variables needed to work
        self.srr_number = srr_number
        self.seq_type = seq_type
        self.specie = specie
        self.data_dir = data_dir
    # Function to generate circular rotations or permutations to the pattern selected. for example TTAGGG, TAGGGT, AGGGTT and so on.
    def circular_rotations(self, pattern):
        return set(pattern[i:] + pattern[:i] for i in range(len(pattern)))
    # Function to generate a bed file where the pattern is founded in any fastq file.      
    def process_pattern(self, circ_rotation, pattern ,direction, bed_file, fastq_file):
        #This small function was added to reduce time as not all reads has the pattern TTAGGG or CCCTAA, so this increase the speed considerably.
        def search_pattern(seq, patterns):
            for pattern in patterns:
                if pattern in seq:
                    return True
            return False

        # Pairwise align to detect all possible options where the pattern could appear.
        mismatch_penalty = 0
        gap_open_penalty = -1
        gap_extend_penalty = -1
        min_score = 17 # This represents one mismatch in the pattern
        # Aligns the sequence previously detected and the pattern which could be a list of patterns.
        aligner = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=mismatch_penalty,
                                    open_gap_score=gap_open_penalty, extend_gap_score=gap_extend_penalty)
        # It opens a bed file to store the data
        with open(bed_file, "a") as bed:
            for record in SeqIO.parse(fastq_file, "fastq"): 
                ADN_seq = record.seq.upper()
                record_id = record.id
                if search_pattern(ADN_seq, pattern):
                    for rotation in circ_rotation:
                        alignments = aligner.align(rotation, ADN_seq)
                        for alignment in alignments:
                            if alignment.score >= min_score:
                                coordinates = np.array(alignment.coordinates)
                                start, end = min(coordinates[-1:, 0]), max(coordinates[-1:, -1])
                                score = alignment.score
                                bed.write(f"{record_id}\t{start}\t{end}\t{direction}\t{score}\n")
    
    # Function to analyse the coordinates of the fusion found which can be inward or otward.
    def merge_coord(self, coord_dict, data, max_distance):

        for record_id, coordinates in coord_dict.items():
            forward_coords = coordinates['forward']
            reverse_coords = coordinates['reverse']
            
            for forward_coord in forward_coords:
                for reverse_coord in reverse_coords:
                    # Evaluates the distance of each convination.
                    overlap_start = max(forward_coord[0], reverse_coord[0])
                    overlap_end = min(forward_coord[1], reverse_coord[1])

                    if overlap_start == overlap_end:
                        # If there is a overlap the distance would be 0
                        distance = 0
                    else:
                        # if there is not overlap the distance would be calculate depending of the coordinates
                        distance = min(
                            abs(forward_coord[1] - reverse_coord[0]),
                            abs(reverse_coord[1] - forward_coord[0])
                        )
                    # After calculating the distance a condition is used to determine if the distance is enough or not.
                    if distance <= max_distance:
                        total_score = forward_coord[2] + reverse_coord[2]
                        if total_score >= 35:
                            fusion_start = min(forward_coord[0], reverse_coord[0])
                            fusion_end = max(forward_coord[1], reverse_coord[1])
                            # inward is the pattern where the fusion is forward + reverse
                            # outward is the pattern where the fusion is reverse + forward
                            fusion_type = "Inward" if forward_coord[0] <= reverse_coord[0] else "Outward" 
                            data.append({
                                'ID': record_id,
                                'Distance': distance,
                                'FusionStart': fusion_start,
                                'FusionEnd': fusion_end,
                                'FusionType': fusion_type,
                                'TotalScore': total_score
                            })
    # A function to create a fasta file with the sequence and the header.
    def fusion_coord(self, fusion_coordinates, fastq_file, output_fasta_file):
        with open(output_fasta_file, "a") as output_fasta:
            for seq_record in SeqIO.parse(fastq_file, "fastq"):
                record_id = seq_record.id
                if record_id in fusion_coordinates:
                    coordinates = fusion_coordinates[record_id]
                    DNA_seq = seq_record.seq
                    total_length = len(DNA_seq)
                    for fusion_start, fusion_end, orientation, total_score, distance in coordinates:
                        fusion = DNA_seq[fusion_start:fusion_end]
                        # Header has the  orientation, th fusion, the length of the read, the total score, the fusion sequence and the distance within the fusion between the two patterns.
                        header = f"{orientation} {fusion_start}-{fusion_end} TotalLength={total_length} bp TotalScore={total_score} fusion {fusion} distance = {distance}"
                        seq = DNA_seq
                        seq_record = SeqRecord(seq, id=record_id, description=header)
                        SeqIO.write(seq_record, output_fasta, "fasta")

    # Main analysis function.                    
    def candidate_fus(self, bed_forward, bed_reverse, fastq_file, number):
        print('Loading bed_files')
        # To free space the data, in the bed files generated are analysed in chunks.
        # This unifies the bed files into one data set to get all the fusions founded.
        # This unifies the IDs of the bed files which are the same and records all the information provided in one dictionary
        coord_dict = {}
        for df_forward_chunk in pd.read_csv(bed_forward, sep='\t', header=None, names=['ID', 'Start', 'End', 'Direction', 'Score'], chunksize=chunksize):
            for index, row in df_forward_chunk.iterrows():
                record_id = row['ID'] 
                start = row['Start']
                end = row['End']
                score = row['Score']
                if record_id in coord_dict:
                    coord_dict[record_id]['forward'].append((start, end, score),)
                else:
                    coord_dict[record_id] = {'forward': [(start, end, score)], 'reverse': []}

        for df_reverse_chunk in pd.read_csv(bed_reverse, sep='\t', header=None, names=['ID', 'Start', 'End', 'Direction', 'Score'], chunksize=chunksize):
            for index, row in df_reverse_chunk.iterrows():
                record_id = row['ID']
                start = row['Start']
                end = row['End']
                score = row['Score']
                if record_id in coord_dict:
                    coord_dict[record_id]['reverse'].append((start, end, score))

        coord_dict = {record_id: coordinates for record_id, coordinates in coord_dict.items() if coordinates['reverse']}

        data = []
        max_distance = 20
        # After the data is unified in one dictionary "coord_dict" the distance is calculated and all information is stored in "data".
        self.merge_coord(coord_dict, data, max_distance)

        combined_df = pd.DataFrame(data)
        # "data" is orgined and prepared to remove redundances and get the less distance of the fusion aswell as the better scored one.
        if not combined_df.empty:

            combined_df.sort_values(by=['ID', 'FusionStart', 'Distance', 'TotalScore'], ascending=[True, True, True, False], inplace=True) # getting same ID fusions
            combined_df.reset_index(drop=True, inplace=True)
            combined_df = combined_df.drop_duplicates(subset=['ID', 'FusionStart'], keep='first') # Duplicates of the same ID and the same fusion start are removed
            combined_df.sort_values(by=['ID', 'FusionEnd', 'FusionStart'], ascending=[True, True, False], inplace=True) # Looking for the closest fusion, start coordinate is descending and end coordinate is ascending
            combined_df = combined_df.drop_duplicates(subset=['ID', 'FusionEnd'], keep='first') # Duplicates of the same ID and same fusion end are removed
            combined_df = combined_df.drop_duplicates(subset=['ID', 'Distance'], keep='first') # Duplicates with the same distance aswell are removed
            combined_df['ID_read'] = combined_df['ID'].str.split('.').str[1].astype(int) 
            combined_df = combined_df.sort_values(by=('ID_read'), ascending=True) # The ID is devided by a "." so to sort values, the ID was splitted to get a better organization
            combined_df = combined_df.drop('ID_read', axis = 1)
            combined_df = combined_df.drop_duplicates(subset='ID', keep='first')
            combined_df.to_csv(f"combined_distances_{self.specie}_{number}.csv", index=False) # Data is stored in CSV file

            csv_file = f"combined_distances_{self.specie}_{number}.csv"
            os.makedirs('Telo_fus_raw', exist_ok=True)
            output_fasta_file = f"Telo_fus_raw/Telo_fus_{self.specie}_{number}.fasta"

            fusion_coordinates = {}

            
            df = pd.read_csv(csv_file)
            # data is stored in a dictionary
            for _, row in df.iterrows():
                record_id = row['ID']
                fusion_start = row['FusionStart']
                fusion_end = row['FusionEnd']
                orientation = row['FusionType']
                total_score = row['TotalScore']
                distance = row['Distance']
                if record_id not in fusion_coordinates:
                    fusion_coordinates[record_id] = []
                fusion_coordinates[record_id].append((fusion_start, fusion_end, orientation, total_score, distance))
            # Fuction is called to create the fasta file.
            self.fusion_coord(fusion_coordinates, fastq_file, output_fasta_file)
        # If there are no fusions in the file a fasta file is created with no sequences.
        else:
            sequence = Seq("None")
            record = SeqRecord(sequence, id="No sequences", description="")
            output_fasta_file = f"Telo_fus_raw/Telo_fus_{self.specie}_{number}.fasta"

            with open(output_fasta_file, "w") as output_fasta:
                SeqIO.write(record, output_fasta, "fasta")                

