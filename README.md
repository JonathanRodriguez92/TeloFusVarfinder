# TeloFusVarfinder
A tool to search telomere variability and fusions in the genome

## Overview

TeloFusVarfinder is a powerful tool designed for the analysis of telomere variability and potential telomere fusions within genomes. This tool facilitates the examination of telomere motifs, their variants, and potential fusion events by providing streamlined data processing and analysis functionalities.

## Features

- **Telomere Fusion Detection:** Identifies potential telomere fusion events and provides relevant data for further analysis.
- **Telomere Variability Analysis:** The tool allows for the investigation of different telomere motifs and their frequencies in the genome.
- **Data Management:** Efficiently handles raw data, processes CSV files, and merges data across different files and species.
  
### Telomere Fusion Detection

The process of telomere fusion detection in TeloFusVarfinder involves several key functions and methodologies:

- **Circular Rotation Function:** This function is crucial for detecting telomere fusions as it enables the circular rotation of sequences. This allows for the identification of potential fusion, even thought the principal motif is TTAGGG or CCCTAA the fusion point of the telomere could be any permutation as the shortening of the telomere could ocurr in any base.
  - Inputs:
    - Pattern: In the zoonomia project was TTAGGGTTAGGGTTAGGG and CCCTAACCCTAACCCTAA
  - Outputs:
    - A list of patterns with their posible circular permutations   

- **Process Pattern with Local Alignment:** This function leverages the `PairwiseAligner` class from the `Bio.Align` module to perform local alignment on patterns and stores the data in 2 bed files forward and reverse. The alignment is configured with the following parameters:
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
      
- **Merge Coordinates Function:** The `merge_coord` function processes 2 BED files (forward and reverse) to identify potential fusion events by analysing the distance between patterns. It merges coordinates with the same ID from different files to provide a comprehensive view of potential fusions..
  - Inputs:
    - Coord_dict: A dictionary that involves the coordinates with the same ID previously calculated.
    - Data: An empty list
    - Max_distance: The maximum distance permitted between patterns (e.g.: if 20 is the maximum distance permitted and the distance between pattern is 10 the fusion would look like TTAGGGTTAGGGxxxxxxxxxxCCCTAACCTAA where x is a random base)
      
  - Outputs:
    -  A table with fusion details: ID, coordinates, type of fusion Inward(forward-reverse) or Outward(reverse-forward), total score selected

- **Fusion Coordinates Function:** This function is responsible for searching through FASTA sequence data and creating a new FASTA file containing the sequences of interest. It uses the coordinates identified in the `merge_coord` function to extract the sequences related to potential fusions.

- **Main Function:** The main function ties all of these processes together to facilitate the search for telomere fusions. It manages the overall workflow, ensuring that each function is called in the correct sequence to efficiently identify potential fusions within the genome.

## File Descriptions

The project consists of several Python scripts that provide different aspects of telomere variability and fusion analysis:

- **Telo_fusion_process.py:** This script is responsible for detecting potential telomere fusions in the genome. It processes data to identify and analyze possible fusion events.
  
- **Telo_motif_process.py:** This script analyzes the variability of telomere motifs in the genome. It processes CSV files to extract information about telomere motifs and variants, then creates summary Excel files with the relevant data.
  
- **SRR_download_process.py:** This script handles the downloading of raw data from the Sequence Read Archive (SRA) for analysis. It ensures that the required data is available for further processing.

## Getting Started

To use TeloFusVarfinder, follow these steps:

1. **Set Up Your Environment:**
    - Make sure you have Python 3.x installed on your machine.
    - Install necessary dependencies using `pip` or your preferred package manager.

2. **Download the Scripts:**
    - Clone the repository or download the scripts directly to your project directory.

3. **Prepare Your Data:**
    - Organize your data files (CSV, Excel) in the appropriate directory as per the script requirements.

4. **Run the Scripts:**
    - Execute the scripts (`Telo_fusion_process.py`, `Telo_motif_process.py`, and `SRR_download_process.py`) in the desired sequence depending on your analysis needs.

5. **Review the Results:**
    - The scripts will produce Excel files with the analyzed data. Review these files to gain insights into telomere variability and potential fusion events.

## Contact

For questions, feedback, or support, please contact [your email address].

Enjoy using TeloFusVarfinder for your genome analysis!

