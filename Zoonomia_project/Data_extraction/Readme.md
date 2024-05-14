# Extraction of information
This is how I stracted all data from different pathways
## **Extraction_zoonomia_information.py**
### Description

This is a Python script designed for text mining to extract information from various websites. It specifically targets:

- [Zoonomia Project](https://zoonomiaproject.org/the-mammal-tree-list-view/): Extracts project information.
- [NCBI Datasets API](https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/): Fetches genome information.
- [NCBI Assembly](https://ncbi.nlm.nih.gov/assembly): Retrieves SRR number and desired information.
- [Animal Diversity Web](https://animaldiversity.org): Gathers longevity data.
  
The script utilizes web scraping techniques along with API requests to extract specific information from the provided websites. It then processes the extracted data for further analysis or storage.

### Usage

To use the script, run the following command: **python Extraction_zoonomia_information.py**

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **jsonlines**
- **requests**
- **beautifulSoup**
- **pandas**
- **numpy**
- **zipfile**

You can install these dependencies using pip or conda:
```
pip install jsonlines requests beautifulsoup4 pandas numpy zipfile

```
## **Fusion_analysis.py**

### Description

This script analyzes the data obtained in the TeloFusVarFinder project. It unifies all the fusions from each species and creates two tables: one with the raw data (`Fusion_results_raw.xlsx`) and the other table with the data normalized by genome coverage (`Fusion_results_coverage.xlsx`). The tables are divided into:

- **Inward (distance 0) Fusions**: Fusions of forward + reverse with distance 0 (e.g., TTAGGGCCCTAA)
- **Outward (distance 0) Fusions**: Fusions of reverse + forward with distance 0 (e.g., CCCTAATTAGGG)
- **Inward (other) Fusions**: Fusions of forward + reverse with distance > 0 (e.g., TTAGGGxxxCCCTAA)
- **Outward (other) Fusions**: Fusions of reverse + forward with distance > 0 (e.g., CCCTAAxxxTTAGGG)
- **Total_fusions_0**: Total inward and outward fusions with distance 0
- **Total_fusions**: Total of all fusions inward + outward
- **Total_inward**: Total fusions inward 0 + other
- **Total_outward**: Total fusions outward 0 + other

### Usage

To use the script, run the following command: `python Fusion_analysis.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**

You can install these dependencies using pip:

```
pip install pandas
```

## **Motif_variant_analysis.py**

### Description

This script analyzes the data obtained in the TeloFusVarFinder project. It calculates the total frequency of the main motif (TTAGGG or CCCTAA) and its variants. The script generates three tables:

- **Table_{specie}_1.xlsx** for R1
- **Table_{specie}_2.xlsx** for R2

These tables count all the variants from TTAGGG and CCCTAA and unify them in a table with their frequency.

- **Merged_table_all.xlsx**: Merges all tables into one where the columns are the file name and the values are their frequency.

- **Telomere_variants_catalogue.xlsx**: Adds CCCTAA motif and its variants to its complementary reversed and the corresponding variants. The "Variant" column represents the type of the variant including the main motif TTAGGG. `{specie}_1` and `{specie}_2` represent the column name from the species and respectively its fastQ file (R1 or R2).

### Usage

To use the script, run the following command: `python Motif_variant_analysis.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**
- **biopython**

You can install these dependencies using pip:

```
pip install pandas biopython
```
## **Fasta_fusion_analysis.py**

### Description

This script analyzes the data obtained in the TeloFusVarFinder to examine fusions with distance 0 and attempt to identify any correlation of the breakpoint.

### Usage

To use the script, run the following command: `python Fasta_fusion_analysis.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**

You can install these dependencies using pip:
```
pip install pandas 
```

## **Extraction_gene_information.py**

### Description

This script utilizes BeautifulSoup and a URL to download TOGA (Transcript and Gene Ortholog Assignment) information, specifically transcript and gene information from the assembled data obtained in the Zoonomia project. It downloads and organizes the data into `{species}_orthologsClassification.tsv` files, representing the orthology classification given by TOGA, which includes the following classes:

- **one2one**: One-to-one, a gene from the reference genome to one region of the query genome.
- **one2many**: One-to-many, a gene from the reference genome to many regions of the query genome.
- **many2many**: Many-to-many, many genes from the reference genome to many regions of the query genome.
- **one2zero**: One-to-zero, no orthologous genes found in the assembled genome.

### Usage

To use the script, run the following command: `python Extraction_gene_information.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **gzip**
- **requests**
- **BeautifulSoup**
- **pandas**

You can install these dependencies using pip:
```
pip install gzip requests beautifulsoup4 pandas 

```

## **Extraction_gene_to_gene_information.py**

### Description

This script utilizes BeautifulSoup and a URL to download TOGA (Transcript and Gene Ortholog Assignment) information, specifically transcript and gene information from the assembled data obtained in the Zoonomia project. It downloads and organizes the data into `{species}_loss_summ_data.tsv` files, representing the gene, transcript, and projection data with different classes:

- **I (intact)**: Represents intact genes/transcripts.
- **PI (partially intact)**: Represents partially intact genes/transcripts.
- **PG (Paralogous Projection)**: Represents paralogous projections where orthologous chains are not identified.
- **UL (Uncertain loss)**: Represents uncertain loss where there is not enough data for resolution.
- **PM (Partially Missing)**: Represents partially missing genes/transcripts.
- **M (Missing)**: Represents missing genes/transcripts.
- **L (Lost)**: Represents lost genes/transcripts.

### Usage

To use the script, run the following command: `python Extraction_gene_to_gene_information.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **gzip**
- **requests**
- **BeautifulSoup**
- **pandas**

You can install these dependencies using pip:
```
pip install gzip requests beautifulsoup4 pandas 

```

## **Specific_proteins_class.py**

### Description

This script compiles gene information obtained from TOGA (`{species}_orthologsClassification.tsv` and `{species}_loss_summ_data.tsv`) and creates a table with information about specific proteins previously selected. The selected proteins include: 'ATRX', 'DAXX', 'TERT', 'RTEL1', 'PIF1', "TERF1", "TERF2", "TINF2", "TPP1", 'POT1', "POT1A", 'POT1B', "DCR1B", "ACD", "WRN", "BLM", "TP53", "SIRT1", "RAD52", "FEN1", "RPA1", 'RPAIN', 'RPA2', 'ASF1A', 'SETDB1', 'RMI', 'BLM', 'TOP3A', 'PML', 'ZBTB40', and 'RAD51'.

### Usage

To use the script, run the following command: `python Specific_proteins_class.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**

You can install these dependencies using pip:

```
pip install pandas 

```

## **All_proteins_class.py**

### Description

This script compiles all the gene information obtained from TOGA (`{species}_orthologsClassification.tsv` and `{species}_loss_summ_data.tsv`) and creates a table with the information of all genes with the protein name structured like `{protein}_g_{t_gene}_{reference}`.

### Usage

To use the script, run the following command: `python All_proteins_class.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**

You can install these dependencies using pip:

```
pip install pandas 

```

## **FusionVariantProteinInfo.py**

### Description

This script compiles all the information obtained from the project TeloFusVarFinder. It gathers tables generated in the following scripts:

- **All_protein_class.py** or **Specific_proteins_class.py** (example: `zoonomia_proteins_class_specific_mouse.xlsx`)
- **Fusion_analysis.py** (example: `Fusion_results_coverage.xlsx`)
- **Motif_variant_analysis.py** (example: `Telomere_variants_catalogue_raw.xlsx`)
- **Extraction_zoonomia_information.py** (example: `zoonomia_sp_info.xlsx`)

### Usage

To use the script, run the following command: `python FusionVariantProteinInfo.py`

### Dependencies

Ensure you have the following dependencies installed:

- **Python 3.12.0**
- **pandas**
- **numpy**

You can install these dependencies using pip:

```
pip install pandas numpy

```

## Acces to the scripts:




