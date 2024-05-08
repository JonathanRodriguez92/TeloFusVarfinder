"""

Description:
------------
A script that uses beatifulsoup and a url to download TOGA information specially transcript and gene information from the assembled data obtained in the project named Zoonomia.
It downloads and creates a folder from the reference selected and organises the data into {specie}_loss_summ_data.tsv
This data has the representation of Gene, Transcript and Projection with different classes:
I (intact), PI (partially intact), PG(Paralogous Projection, not orthologos chains identified), 
UL (Uncertain loss, not enough data to have a resolution), PM (Partially Missing), M(Missing), L(lost)

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
        specie = row['Species'].replace(" ", "_")
        output_directory = f'Orthologs_{reference}'
        directory_out = os.path.join(output_directory, f'{specie}_loss_summ_data.tsv')
        if not os.path.exists(directory_out):
            os.makedirs(output_directory, exist_ok=True)
            
            print(specie)
            url = f'https://genome.senckenberg.de/download/TOGA/{reference}_reference/'
            response = requests.get(url, verify=False)
            soup = BeautifulSoup(response.text, 'html.parser')
            # It looks through out all the folders and get all names.
            directories = soup.find_all('a', href=True) 
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
                                # Recopilates the data loss_summ_data.tsv from each folder
                                loss_summ = "loss_summ_data.tsv.gz"
                                loss_summ_url = url + directory_url + link_url + loss_summ
                                loss_summ_response = requests.get(loss_summ_url, verify=False)
                                if loss_summ_response.status_code == 200:
                                    # it opens a directorey and places the data with a new name.
                                    with open(os.path.join(output_directory, f'{specie}_loss_summ_data.tsv.gz'), 'wb') as f_out:
                                        f_out.write(loss_summ_response.content)
                                        print(f'{specie}_loss_summ_data.tsv.gz downloaded')
                                        # it decompresses the data
                                        with gzip.open(BytesIO(loss_summ_response.content), 'rb') as f_in:
                                            with open(os.path.join(output_directory, f'{specie}_loss_summ_data.tsv'), 'wb') as f_out_uncompressed:
                                                f_out_uncompressed.write(f_in.read())
                                                print(f'{specie}_loss_summ_data.tsv decompressed')

                                        os.remove(os.path.join(output_directory, f'{specie}_loss_summ_data.tsv.gz'))
                                        break  
    print("Process completed.")
