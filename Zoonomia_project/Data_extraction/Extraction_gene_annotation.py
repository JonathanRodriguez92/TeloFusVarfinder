"""

Description:
------------
A script that uses beatifulsoup and a url to download TOGA information specially transcript and gene information from the assembled data obtained in the project named Zoonomia.
It downloads and creates a folder from the reference selected and organises the data into {specie}_orthologsClassification.tsv
This data has the representation the orthology classification given by TOGA with this different classes:
one2one -one to one, a gene from the reference genome to one regions of the query genome
one2many - one to many, a gene from the reference genome to many regions of the query genome
many2many - many to many, many genes from the reference genome to many regions of the query genome
one2zero - one 2 zero, no orthologous genes found in the assembled genome.

Usage:
------

command line:
python Extraction_gene_to_gene_information.py

Dependencies:
-------------
    - python 3.12.0
    - gzip
    - requests
    - beautifulSoup
    - pandas

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 15-March-2024    

"""

import os
import requests
import gzip
import pandas as pd
import urllib3
from io import BytesIO
from bs4 import BeautifulSoup

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

reference_list = ['human_hg38', 'mouse_mm10']

df = pd.read_excel('zoonomia_sp_info.xlsx', engine='openpyxl')

# Gets all data from human and mouse
for reference in reference_list:    
    for index, row in df.iterrows():
        output_directory = f'Orthologs_{reference}'
        os.makedirs(output_directory, exist_ok=True)
        specie = row['Species'].replace(" ", "_")
        url = f'https://genome.senckenberg.de/download/TOGA/{reference}_reference/'
        response = requests.get(url, verify=False)
        soup = BeautifulSoup(response.text, 'html.parser')
        directories = soup.find_all('a', href=True)
        # It looks through out all the folders and get all names
        for directory in directories:
            directory_url = directory['href']

            if directory_url.endswith('/'):
                directory_response = requests.get(url + directory_url, verify=False)

                directory_soup = BeautifulSoup(directory_response.text, 'html.parser')

                # it gets all the links to start the download
                links = directory_soup.find_all('a', href=True)
                for link in links:
                    link_url =link['href']
                    if specie in link_url:
                        if link_url.endswith('/'):
                            # Recopilates the data orthologsClassification.tsv from each folder
                            orthologos_classification = "orthologsClassification.tsv.gz"
                            orthologos_url = url + directory_url + link_url + orthologos_classification
                            orthologos_response = requests.get(orthologos_url, verify=False)
                            if orthologos_response.status_code == 200:
                                # it opens a directorey and places the data with a new name.
                                with open(os.path.join(output_directory, f'{specie}_orthologsClassification.tsv.gz'), 'wb') as f_out:
                                    f_out.write(orthologos_response.content)
                                    print(f'{specie}_orthologsClassification.tsv.gz downloaded')
                                    # it decompresses the data
                                    with gzip.open(BytesIO(orthologos_response.content), 'rb') as f_in:
                                        with open(os.path.join(output_directory, f'{specie}_orthologsClassification.tsv'), 'wb') as f_out_uncompressed:
                                            f_out_uncompressed.write(f_in.read())
                                            print(f'{specie}_orthologsClassification.tsv decompressed')

                                    os.remove(os.path.join(output_directory, f'{specie}_orthologsClassification.tsv.gz'))
                                    print(f'{specie}_orthologsClassification.tsv.gz removed')
                                    break  
        print("Proceso completado.")




