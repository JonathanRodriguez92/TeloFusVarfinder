# TeloFusVarfinder
A tool to search telomere variability and fusions in the genome

## Overview

TeloFusVarfinder is a powerful tool designed for the analysis of telomere variability and potential telomere fusions within genomes. This tool facilitates the examination of telomere motifs, their variants, and potential fusion events by providing streamlined data processing and analysis functionalities.

## Features

- **Telomere Fusion Detection:** Identifies potential telomere fusion events and provides relevant data for further analysis.
- **Telomere Variability Analysis:** The tool allows for the investigation of different telomere motifs and their frequencies in the genome.
- **Data Management:** Efficiently handles raw data, processes CSV files, and merges data across different files and species.
  
### Telomere Fusion Detection `TeloFusProcessor` 

The process of telomere fusion detection in TeloFusVarfinder involves several key functions and methodologies:

- **`Circular_rotations` Function:** This function is crucial for detecting telomere fusions as it enables the circular rotation of sequences. This allows for the identification of potential fusion, even thought the principal motif is TTAGGG or CCCTAA the fusion point of the telomere could be any permutation as the shortening of the telomere could ocurr in any base.
  - Inputs:
    - Pattern: In the zoonomia project was TTAGGGTTAGGGTTAGGG and CCCTAACCCTAACCCTAA
  - Outputs:
    - A list of patterns with their posible circular permutations   

- **`Process_pattern` Function:** This function leverages the `PairwiseAligner` class from the `Bio.Align` module to perform local alignment on patterns and stores the data in 2 bed files forward and reverse. The alignment is configured with the following parameters:
    - `mismatch_penalty = 0`
    - `gap_open_penalty = -1`
    - `gap_extend_penalty = -1`
    - `min_score = 17` 
    
    These parameters allow for the identification of potential fusions with high sensitivity while providing a slightly flexibility in matching (1 mismatch).
  - Inputs:
    - Circular rotations: A pattern created in the circular_rotation function
    - Pattern: A simple pattern to detect potencial reads
    - Direction: Forward or reverse depending of what you analyse
    - Bed file: The name of the bed file (e.g.: {specie}{direction}{number}.bed)
    - Fastq file: The file to analyse
  - Outputs:
    - Bed file: Project case was {Specie}{direction}{number}.bed as it could be for R1 or R2.  
      
- **`Merge_coord` Function:** The `merge_coord` function processes 2 BED files (forward and reverse) to identify potential fusion events by analysing the distance between patterns. It merges coordinates with the same ID from different files to provide a comprehensive view of potential fusions..
  - Inputs:
    - Coord_dict: A dictionary that involves the coordinates with the same ID previously calculated.
    - Data: An empty list
    - Max_distance: The maximum distance permitted between patterns (e.g.: if 20 is the maximum distance permitted and the distance between pattern is 10 the fusion would look like TTAGGGTTAGGGxxxxxxxxxxCCCTAACCTAA where x is a random base)
      
  - Outputs:
    -  A CSV table with fusion details: ID, coordinates, type of fusion Inward(forward-reverse) or Outward(reverse-forward) and total score of the fusion

- **`Fusion_coord` Function:** This function is responsible for searching through FASTA sequence data and creating a new FASTA file containing the sequences of interest. It uses the coordinates identified in the `merge_coord` function to extract the sequences related to potential fusions.
  - Inputs:
    - Fusion coordinates: A dictionary which reads the csv table and coordinates the data previously calculated.
    - Fastq_file: To look for the reads and store them in a fasta file
  - Outputs:
    - Fasta file
      
- **`candidate_fus` Function:** The main function called ties all of these processes together to facilitate the search for telomere fusions. It manages the overall workflow, ensuring that each function is called in the correct sequence to efficiently identify potential fusions within the genome. At the end Bed files are delected to free space.
  
### Telomere Variability Analysis `TelomericVariantSearcher`

The process of telomere variability analysis in TeloFusVarfinder involves two key functions:

- **`are_similar` Function:** This function compares two sequences to determine their similarity. It accepts two sequences as input and a maximum allowed difference (`max_error`). The function computes the number of differing characters between the sequences and returns a boolean indicating whether the difference is within the specified error threshold. This function helps identify telomeric variants that closely match a reference sequence.

- **`search_telomeric_variants` Function:** This function is responsible for searching telomeric variants within the genome data. It scans through sequence data, utilizing the `are_similar` function to identify sequences that closely resemble known telomeric motifs. The function can handle large datasets and generate insights into the distribution and frequency of different telomeric motifs across various species.
  - Input:
    - min_frequency: A threshold which stores the data with a minimun frequency (e.g: 1)
    - fastq_file: Fastq_file to analyse
    - number: A number to diferenciates the output file name with others.
  - Output:
    - A CSV file with
      - ID: the ID of the read or scaffold analysed
      - TelomereMotif: The original motif analysed (e.g: TTAGGG).
      - Variant: The variant analysed from the original motif with max_errors permitted.
      - Frequency: How many times appear the motif in the read and do not overlaps.
      - PorcentageTelomereSequence: How big in porcentage is the motif/variant appearing in the read (e.g: 30% of TTAGGG).
      - Length: Total length of the read.

Together, these functions enable the tool to provide comprehensive insights into telomere variability across different genomes and species. By identifying telomeric variants and their frequencies, users can gain valuable information about telomere structure and potential areas for further investigation.

### Data Management `SrrDownloader`

The process to download SRR files from sra-toolkits prefetch, fastq-dump and filtration with fastq_quality_filter.

## Scripts

Here are the main scripts in this repository and their descriptions:

- **[Telo_fusion_process.py](Telo_fusion_process.py)**: 
    - Description: Handles telomere fusion detection, including circular rotation function, local alignment using the `PairwiseAligner` with specified parameters, merging coordinates to find fusions from BED files, and searching FASTA sequences to create a new FASTA file.
    - [View the script](https://github.com/JonathanRodriguez92/TeloFusVarfinder/blob/4dd5a032f68f10269ed0521ab7bc86b334ee53a8/Telo_fusion_process.py)
  
- **[Telo_motif_process.py](Telo_motif_process.py)**:
    - Description: Manages telomere variability analysis, including functions for comparing sequences for similarity and searching for telomeric variants.
    - [View the script](https://github.com/JonathanRodriguez92/TeloFusVarfinder/blob/3b199b1d46f3c0bce1ce718c241652ab3d02ce50/Telo_motif_process.py)
  
- **[SRR_download_process.py](SRR_download_process.py)**:
    - Description: Responsible for downloading and filtering data.
    - [View the script](https://github.com/JonathanRodriguez92/TeloFusVarfinder/blob/a8b36d2f228cee7212a8d6f4bfeb888119aee90c/SRR_download_process.py)

## Getting Started

To use TeloFusVarfinder, follow these steps:

1. **Set Up Your Environment:**
    - Make sure you have Python 3.12 installed on your machine.
    - Install necessary dependencies using `pip` or your preferred package manager.
    - **All depencdices are here**
    - [see the script](https://github.com/JonathanRodriguez92/TeloFusVarfinder/blob/e79e59c390402ca8c2273e304a654428488712fa/miniconda_env_lib.sh)

2. **Download the Scripts:**
    - Clone the repository or download the scripts directly to your project directory.

3. **Prepare Your Data:**
    - Organize your data files.
      - A table which has columns called
        - Species
        - srr_number
          
      - Your own motifs and patterns to analyse 

4. **Run the Scripts:**
    - Execute the scripts Telomere_analysis.py and whichs parallelizes all the scripts in one (`Telo_fusion_process.py`, `Telo_motif_process.py`, and `SRR_download_process.py`).
    - Important: The script was runned in a cluster, so the end of it should be modified to ensure you can use it other computers
    - [see the script](https://github.com/JonathanRodriguez92/TeloFusVarfinder/blob/892fd542f1347c6d50d0b80235fc15faf33c8970/Telomere_analysis.py)
      

## Contact

For questions, feedback, or support, please contact jonatan.rodriguez@estudiante.uam.es.

Enjoy using TeloFusVarfinder for your genome analysis!

