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
### Acces to the script




