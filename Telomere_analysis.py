
"""

Description:
------------
This script unifies Telo_motif_process.py, Telo_fusion_process.py and SRR_download_process working in a multiprocessing environment.
Analysing and processing SRR data fastq files (SINGLE OR PAIRED).

Usage:
------

command line:
python Telomere_analysis.py {number} / this represents the row of the table used with a SRR number.


example of usage with multiprocessing in a regular environment:

def process_row_wrapper(args):
    index, row = args
    process_species(row)

#if __name__ == "__main__":
    df = pd.read_excel('zoonomia_subgroup.xlsx', engine='openpyxl')
    num_processes = 2

    # Process rows in batches of 3
    for i in range(0, len(df), num_processes):
        batch = df.iloc[i:i+num_processes]
        processes = []
        for index, rows in batch.iterrows():
            process = Process(target=process_row_wrapper, args=((index, rows),))
            process.daemon = False  # Set daemon attribute explicitly
            processes.append(process)
            process.start()

        for process in processes:
            process.join()

Dependencies:
-------------
    - python 3.12.0
    - pandas
    - multiprocessing

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 10-april-2024    

"""

from Telo_fusion_process import TeloFusProcessor
from Telo_motif_process import TelomericVariantSearcher
from SRR_download_process import SrrDownloader
from multiprocessing import Process
import shutil
import pandas as pd
import os


# funtion to analyse the SRR fastq files
def process_species(row):
    # Variables needed to work with the program
    srr_number = row['SRR_number'] 
    specie = row['Species'].replace(" ", "_")
    seq_type = row['Library_layout']

    # Paths where the data is stored
    dir = f'{specie}'
    path_SRA = f'/data3/genomes/ZOONOMIA/SRA/'
    path_1 = os.path.join(path_SRA, dir) 
    os.makedirs(path_1, exist_ok=True)
    download_file = f'/data3/genomes/ZOONOMIA/SRA/{specie}'
    descompresor_file = f'/data3/genomes/ZOONOMIA/SRA/{specie}/{srr_number}'
    complete_file = f'/data3/genomes/ZOONOMIA/SRA/Complete_{specie}'
    path_compresed = os.path.join(descompresor_file, f'{srr_number}.sra')
    path_paired = os.path.join(download_file, f'{srr_number}_1.fastq')
    path_paired_f_1 = os.path.join(download_file, f'{srr_number}_filtered_1.fastq')
    path_paired_f_2 = os.path.join(download_file, f'{srr_number}_filtered_2.fastq')
    path_single = os.path.join(download_file, f'{srr_number}.fastq')
    path_single_f = os.path.join(download_file, f'{srr_number}_filtered.fastq')

    # Data for telo_motif_process
    original_motifs = ['TTAGGG', 'CCCTAA']
    print(f'Processing {specie.replace("_", " ")}')

    # If the SRR is paired ended
    if seq_type == 'PAIRED' and not pd.isnull(srr_number):

        # Using SRR_download to download the SRR runs of fastq files.
        if not os.path.exists(complete_file):
            if not os.path.exists(path_paired) and not os.path.exists(path_compresed):
                print(f'downloading {specie}')
                srr_down = SrrDownloader(srr_number, download_file)     # Downloading process             
                srr_down.srr_downloader()
                
            if os.path.exists(path_compresed) and not os.path.exists(path_paired):
                print(f'decompressing {specie}')
                srr_down = SrrDownloader(srr_number, download_file)     # Decompressing process              
                srr_down.srr_descompresor()
                shutil.rmtree(descompresor_file)

            if not os.path.exists(path_paired_f_1) or not os.path.exists(path_paired_f_2):
                print('filtering fastq')
                srr_down = SrrDownloader(srr_number, download_file)     # Filtering process 
                srr_down.srr_filter(path_paired_f_1)
                print('All filtered')

            if os.path.exists(path_paired_f_1) and os.path.exists(path_paired_f_2):
                
                # Calling Telo_fusion_process and Telo_motif_process
                telo_fus_processor = TeloFusProcessor(srr_number, seq_type, download_file, specie)
                telomeric_searcher = TelomericVariantSearcher(specie, download_file, srr_number, seq_type, original_motifs=original_motifs)
                # Variables needed to call the functions
                circ_forward = telo_fus_processor.circular_rotations("TTAGGGTTAGGGTTAGGG")
                circ_reverse = telo_fus_processor.circular_rotations("CCCTAACCCTAACCCTAA")

                pattern_forward = telo_fus_processor.circular_rotations("TTAGGG")
                pattern_reverse = telo_fus_processor.circular_rotations("CCCTAA")

                fastq_file_r1 = os.path.join(download_file, f"{srr_number}_filtered_1.fastq")
                fastq_file_r2 = os.path.join(download_file, f"{srr_number}_filtered_2.fastq")
                
                # BED files to store data
                bed_forward_1 = f"forward{specie}_1.bed"
                bed_reverse_1 = f"reverse{specie}_1.bed"
                bed_forward_2 = f"forward{specie}_2.bed"
                bed_reverse_2 = f"reverse{specie}_2.bed"
                    
                print(f'Processing pattern')
                print(f'Searching for telomeric variants and motifs')

                # Multiprocessing analysis where 6 process work simultaneously

                # Loking for R1 and R2 fusions, 4 processes for processing the patherns TTAGGGTTAGGTTAGG and CCCTAACCCTAACCCTAA in each fastQ file.
                process1 = Process(target=telo_fus_processor.process_pattern, args=(circ_forward, pattern_forward, 'forward', bed_forward_1, fastq_file_r1))
                process2 = Process(target=telo_fus_processor.process_pattern, args=(circ_reverse, pattern_reverse, 'reverse', bed_reverse_1, fastq_file_r1))
                process3 = Process(target=telo_fus_processor.process_pattern, args=(circ_forward, pattern_forward, 'forward', bed_forward_2, fastq_file_r2))
                process4 = Process(target=telo_fus_processor.process_pattern, args=(circ_reverse, pattern_reverse, 'reverse', bed_reverse_2, fastq_file_r2))
                # Calculating frequency of TTAGGG repeats and its variants with one mismatch in R1 and R2.
                process5 = Process(target=telomeric_searcher.search_telomeric_variants, args=(fastq_file_r1, 1, '1'))
                process6 = Process(target=telomeric_searcher.search_telomeric_variants, args=(fastq_file_r2, 1, '2'))

                process1.start()
                process2.start()
                process3.start()
                process4.start()
                process5.start()
                process6.start()

                process1.join()
                process2.join()
                process3.join()
                process4.join()
                process5.join()
                process6.join()

                print('Processing fusions')

                # Multiprocessing analysis where 2 process work simultaneously

                # The unification of R1 and R2 forward and reverse data.
                processx = Process(target=telo_fus_processor.candidate_fus(bed_forward_1, bed_reverse_1, fastq_file_r1, '1'))
                processy = Process(target=telo_fus_processor.candidate_fus(bed_forward_2, bed_reverse_2, fastq_file_r2, '2'))

                processx.start()
                processy.start()

                processx.join()
                processy.join()

            else:
                print('Not filtered properly')
            # Files to delete after analysis to free space            
            files_to_delete = [f"reverse{specie}_1.bed", f"forward{specie}_1.bed",f"reverse{specie}_2.bed", f"forward{specie}_2.bed" ] 
            for file_to_delete in files_to_delete:
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)

            print(f"Deleted files:: {', '.join(files_to_delete)}")
            dir_2 = f'Complete_{specie}'
            os.makedirs(os.path.join(path_SRA, dir_2), exist_ok=True)
            shutil.rmtree(download_file)        
    # If the SRR is Single ended
    elif seq_type == 'SINGLE' and not pd.isnull(srr_number):
        if not os.path.exists(complete_file):
            # Using SRR_download to download the SRR runs of fastq files.
            if not os.path.exists(path_single) and not os.path.exists(path_compresed):
                print(f'downloading {specie}')
                srr_down = SrrDownloader(srr_number, download_file)         # Downloading process          
                srr_down.srr_downloader()
                
            if os.path.exists(path_compresed) and not os.path.exists(path_single):
                print(f'decompressing {specie}')
                srr_down = SrrDownloader(srr_number, download_file)         # Decompressing process   
                srr_down.srr_descompresor()
                shutil.rmtree(descompresor_file)

            if not os.path.exists(path_single_f):
                print('filtering fastq')
                srr_down = SrrDownloader(srr_number, download_file)         # Filtering process
                srr_down.srr_filter()

            if os.path.exists(path_single_f):
                # Calling Telo_fusion_process and Telo_motif_process
                telo_fus_processor = TeloFusProcessor(srr_number, seq_type, download_file, specie)
                telomeric_searcher = TelomericVariantSearcher(specie, download_file, srr_number, seq_type, original_motifs=original_motifs)
                # Variables needed to call the functions
                circ_forward = telo_fus_processor.circular_rotations("TTAGGGTTAGGGTTAGGG")
                circ_reverse = telo_fus_processor.circular_rotations("CCCTAACCCTAACCCTAA")

                pattern_forward = telo_fus_processor.circular_rotations("TTAGGG")
                pattern_reverse = telo_fus_processor.circular_rotations("CCCTAA")

                fastq_file = os.path.join(download_file, f"{srr_number}_filtered.fastq")

                # BED files to store data
                bed_forward = f"forward{specie}.bed"
                bed_reverse = f"reverse{specie}.bed"

                print("Looking for coordinates")

                print(f'Processing pattern')
                # Multiprocessing analysis where 3 process work simultaneously

                # Loking for fusions, 2 processes for processing the patherns TTAGGGTTAGGTTAGG and CCCTAACCCTAACCCTAA in each fastQ file.
                process1 = Process(target=telo_fus_processor.process_pattern, args=(circ_forward, pattern_forward, 'forward', bed_forward, fastq_file))
                process2 = Process(target=telo_fus_processor.process_pattern, args=(circ_reverse, pattern_reverse, 'reverse', bed_reverse, fastq_file))
                # Calculating frequency of TTAGGG repeats and its variants with one mismatch.
                process3 = Process(target=telomeric_searcher.search_telomeric_variants, args=(fastq_file, 1, '0'))

                process1.start()
                process2.start()
                process3.start()

                process1.join()
                process2.join()
                process3.join()

                telo_fus_processor.candidate_fus(bed_forward, bed_reverse, fastq_file, '0')
                # Files to delete to free space
                files_to_delete = [f"reverse{specie}.bed", f"forward{specie}.bed" ] 
                for file_to_delete in files_to_delete:
                    if os.path.exists(file_to_delete):
                        os.remove(file_to_delete)

            dir_2 = f'Complete_{specie}'
            os.makedirs(os.path.join(path_SRA, dir_2), exist_ok=True)
            shutil.rmtree(download_file)

    else:
        print('No SRR available')

# This function creates a variable which has to be assigned before, in the command line, where we can decide which row to analyse.
def process_row(row_number):
    df = pd.read_excel('zoonomia_sp_info_c.xlsx', engine='openpyxl')

    if row_number < 1 or row_number > len(df):
        print("row number invalid")
        return
    row = df.iloc[row_number - 1]

    process_species(row)
    print("Processing row", row_number) 
# Here is assigned to a task ID in the Cluster
if __name__ == "__main__":
    row_number = int(os.environ.get('SGE_TASK_ID', 0))
    if row_number < 1 or row_number > 250:
        print("Número de fila inválido")
    else:
        process_row(row_number)