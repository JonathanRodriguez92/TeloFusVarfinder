"""

Description:
------------
Download, decompresing and filtering fastQ data

Usage:
------

command line:
python SRR_download_process.py srr_number download_file

Dependencies:
-------------
    - python 3.12.0
    - prefetch
    - fastq-dump
    - fastq_quality_filter
    - sra-tools

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 20-December-2024    

"""
import subprocess
import os
import time
import sys


class SrrDownloader:
    def __init__(self, srr_number, download_file):
        self.srr_number = srr_number 
        self.download_file = download_file
        self.max_attempts = 5
    # Function to download the fastq data from https://trace.ncbi.nlm.nih.gov/ using prefetch
    def srr_downloader(self):
        current_attempt = 1
        # Attemps if internet does not get a good signal
        while current_attempt <= self.max_attempts:
            try:
                download = f'prefetch {self.srr_number} -X 500g -p -O {self.download_file}'
                print(f'Downloading {self.srr_number} (Attempt {current_attempt}/{self.max_attempts})')
                subprocess.run(download, shell=True, check=True)
                break  
            except subprocess.CalledProcessError as e:
                print(f"Error during download (Attempt {current_attempt}/{self.max_attempts}): {e}")
                current_attempt += 1
                time.sleep(10)  

        if current_attempt > self.max_attempts:
            print(f"Failed to download {self.srr_number} after {self.max_attempts} attempts.")
            sys.exit()
    # Function to decompres into fastq_R1 and fastq_R2 data using fastq-dump
    def srr_descompresor(self):
        current_attempt = 1
        # Attemps if internet does not get a good signal
        while current_attempt <= self.max_attempts:
            try:
                print(f'Decompressing {self.srr_number}.sra')
                command = f'fastq-dump --outdir {self.download_file} --split-3 --clip {self.download_file}/{self.srr_number}'
                print(f'Decompressing {self.srr_number} (Attempt {current_attempt}/{self.max_attempts})')
                subprocess.run(command, shell=True, check=True)
                break  
            except subprocess.CalledProcessError as e:
                print(f"Error during decompression (Attempt {current_attempt}/{self.max_attempts}): {e}")
                current_attempt += 1
                time.sleep(10)  

        if current_attempt > self.max_attempts:
            print(f"Failed to decompress {self.srr_number} after {self.max_attempts} attempts.")
    # Function to filters fastq_R1 and fastq_R2 data using fastq_quality_filter
    def srr_filter(self, path_paired=None):
        # Diferenciates if data is SINGLE.
        if os.path.exists(f'{self.download_file}/{self.srr_number}.fastq'):
            current_attempt = 1
            # Attemps if internet does not get a good signal
            while current_attempt <= self.max_attempts:
                try:
                    input_fastq = f'{self.srr_number}.fastq'
                
                    print(f'Filtering {input_fastq}')
                    # Quality score higher than 30 and stores the sequences with more than 50 % of the bases with that score.
                    filtered_fastq = f'fastq_quality_filter -q 30 -p 50 -i {self.download_file}/{input_fastq} -o {self.download_file}/{self.srr_number}_filtered.fastq'
                    subprocess.run(filtered_fastq, shell=True, check=True) 
                    print(f'Filtered {input_fastq} successfully.')
                    break
                except subprocess.CalledProcessError as e:
                    print(f"Error during filtering (Attempt {current_attempt}/{self.max_attempts}): {e}")
                    current_attempt += 1
                    time.sleep(10)
            if current_attempt > self.max_attempts:
                print(f"Failed to filtering {input_fastq} after {self.max_attempts} attempts.")    

        # Diferenciates if data is PAIRED.
        elif os.path.exists(f'{self.download_file}/{self.srr_number}_1.fastq'):
            current_attempt = 1
            # Attemps if internet does not get a good signal

            while current_attempt <= self.max_attempts:
                try:
                    input_fastq_1 = f'{self.srr_number}_1.fastq'
                    input_fastq_2 = f'{self.srr_number}_2.fastq'
                    
                    print(f'Filtering {input_fastq_1}')

                    # if the filtered process fails and only one fastQ file was filtered this does not repeat the process of filtration.
                    if os.path.exists(path_paired):
                        # Quality score higher than 30 and stores the sequences with more than 50 % of the bases with that score.
                        filtered_fastq_2 = f'fastq_quality_filter -q 30 -p 50 -i {self.download_file}/{input_fastq_2} -o {self.download_file}/{self.srr_number}_filtered_2.fastq'
                        subprocess.run(filtered_fastq_2, shell = True, check=True)
                    else:
                        # Quality score higher than 30 and stores the sequences with more than 50 % of the bases with that score.
                        filtered_fastq_1 = f'fastq_quality_filter -q 30 -p 50 -i {self.download_file}/{input_fastq_1} -o {self.download_file}/{self.srr_number}_filtered_1.fastq'
                        filtered_fastq_2 = f'fastq_quality_filter -q 30 -p 50 -i {self.download_file}/{input_fastq_2} -o {self.download_file}/{self.srr_number}_filtered_2.fastq'

                        subprocess.run(filtered_fastq_1, shell = True, check=True)
                        subprocess.run(filtered_fastq_2, shell = True, check=True)

                    print(f'Filtered {input_fastq_1} and {input_fastq_2}  successfully.')
                    break
                except subprocess.CalledProcessError as e:
                    print(f"Error during filtering (Attempt {current_attempt}/{self.max_attempts}): {e}")
                    current_attempt += 1
                    time.sleep(10)
            if current_attempt > self.max_attempts:
                print(f"Failed to filtering {input_fastq_1} and {input_fastq_2} after {self.max_attempts} attempts.")    
        
                
