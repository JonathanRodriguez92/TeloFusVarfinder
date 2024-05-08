"""

Description:
------------
A text mining script designed to extract information from the websites: https://zoonomiaproject.org/the-mammal-tree-list-view/ project information, 
https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/ for genome information, https://ncbi.nlm.nih.gov/assembly for SRR number and desired information 
and https://animaldiversity.org for longevity.

Usage:
------

command line:
python Extraction_zoonomia_information.py

Dependencies:
-------------
    - python 3.12.0
    - jsonlines
    - requests
    - beautifulSoup
    - pandas
    - numpy
    - zipfile

Author:
------
    - Author: Jonatan Rodriguez del Valle
    - Date: 10-December-2024    

"""


import jsonlines
import requests
import zipfile
import codecs
import os
import shutil
import re
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

# Extraction from https://zoonomiaproject.org/the-mammal-tree-list-view/
def parse_html(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    species_divs = soup.find_all('div', class_='tree')
    data = []
    # Data extracted
    for species_div in species_divs:
        species_name = species_div.find('h3', class_='species').text
        ncbi_link = species_div.find('a', href=lambda href: href and 'ncbi.nlm.nih.gov/assembly/GC' in href)
        gc_number = ncbi_link['href'].split('/')[-1] if ncbi_link else ""
        order_elem = species_div.find('p', class_='order')
        family_elem = species_div.find('p', class_='family')
        order = order_elem.text.replace("ORDER: ", "") if order_elem else "NA"
        order = order.capitalize()
        family = family_elem.text.replace("FAMILY: ", "") if family_elem else "NA"
        data.append((species_name, order, family, gc_number))

    return pd.DataFrame(data, columns=['Species', 'Order', 'Family', 'NCBI_AssemblyAccNo'])
# Downloading and extrar from a url  
def download_and_extract_zip(url, temp_folder):
    response = requests.get(url)
    response.raise_for_status()
    # downloading all content from the url
    with open(os.path.join(temp_folder, "temp.zip"), "wb") as temp_file:
        temp_file.write(response.content)
    # Unziping file
    with zipfile.ZipFile(os.path.join(temp_folder, "temp.zip"), "r") as zip_ref:
        zip_ref.extractall(temp_folder)

# Extraction of information from the url https://www.ncbi.nlm.nih.gov/sra?term={sra_value}&cmd=DetailsSearch&report=FullXml 
# Extraction of information from the jsonl file extracted in the url https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/...
def extract_values(jsonl_file):
    with codecs.open(jsonl_file, "rb", encoding="utf-8", errors="ignore") as file:
        records = jsonlines.Reader(file)
        for record in records:

            accession = re.sub(r'\.\d+', '', record.get("accession", "")) # accession number (GC number)
            assembly_stat = record.get('assemblyStats', {}) 
            genome_coverage = assembly_stat.get('genomeCoverage',"").replace("x", "") # Genome coverage
            genome_size = assembly_stat.get("totalSequenceLength", "") # Total sequence Length
            genome_size = int(genome_size) # Genome size
            organism = record.get("organism", {}) # Organisim
            common_name = organism.get("commonName", "") # Common name
            assembly_info = record.get("assemblyInfo", {})
            assembly_method = assembly_info.get("assemblyMethod", "") # Assembly method 
            biosample = assembly_info.get("biosample", {})
            sample_ids = biosample.get("sampleIds", [{}])
            sra_value = next((sample.get("value", "NA") for sample in sample_ids if sample.get("db") == "SRA"), "N/A") # SRA value
            url2 = f"https://www.ncbi.nlm.nih.gov/sra?term={sra_value}&cmd=DetailsSearch&report=FullXml"
            response = requests.get(url2)
            
            if sra_value != np.nan:
                if response.status_code == 200:
                    data = []
                    max_gb_multiple_runs = False
                    max_gb_one_run = False
                    # Extraction of SRR and all information needed for the project
                    # Important: This is a strict selection of which srr number to select.
                    # A project (SRA) could have multiples experiments (SRX) and each experiment could have multiple runs (SRR). 
                    # regresion expresion used to find all the data
                    match_experiment = re.findall(r'<b>&lt;EXPERIMENT_PACKAGE&gt(.*?)<b>&lt;/EXPERIMENT_PACKAGE&gt;', response.text)
                    for match in match_experiment:
                        match_bytes = re.search(r'bytes="<span>(.*?)</span>"', match)
                        num_bytes = int(match_bytes.group(1))
                        num_gb = round(num_bytes / (1024 ** 3), 3)
                        match_runs = re.search(r'runs="<span>(.*?)</span>', match)
                        num_runs = int(match_runs.group(1))
                        data.append((num_gb, num_runs))
                        # Calculations of experments with more than 1 run (SRR).
                        experiments_with_multiple_runs = [exp for exp in data if exp[1] > 1]
                        
                        if experiments_with_multiple_runs:
                            # Calculating how big the run is and get the higer one
                            max_gb_multiple_runs = max(experiments_with_multiple_runs, key=lambda x: x[0])
                            index_max_gb_mult_run = data.index(max_gb_multiple_runs)
                            print("Experiment whith higher weight and has multiple runs:", max_gb_multiple_runs, "position:",  index_max_gb_mult_run )
                            
                        else:
                            print("There is not more than one run.")
                        # Calculating experiments with 1 run (SRR).
                        experiments_with_one_run = [exp for exp in data if exp[1] == 1]

                        if experiments_with_one_run:
                            # Calculating how big the run is and get the higher one
                            max_gb_one_run = max(experiments_with_one_run, key=lambda x: x[0])
                            index_max_gb_one_run = data.index(max_gb_one_run)
                            print("Experiment whith higher weight and has 1 runs:", max_gb_one_run, "positionn:" , index_max_gb_one_run)
                            
                        else:
                            print("There is not more than one run.")
                    
                    # Prioritising experiments with only one run and bigger than 20 GB
                    if max_gb_one_run and max_gb_one_run[0] > 20:
                        position = 0
                        for match in match_experiment:
                            if position == index_max_gb_one_run:
                                
                                match_instrument = re.search(r'<b>&lt;INSTRUMENT_MODEL&gt;</b>(.*?)<b>&lt;/INSTRUMENT_MODEL&gt;</b>', match)
                                # Gets the instrument model (eg. illumina)
                                instrument_model = match_instrument.group(1) if match_instrument else np.nan 
                                # Gets the SRR value
                                match_srr = re.search(r'accession="<span>(.?RR.*?)</span>"', match)
                                srr_value = match_srr.group(1) 

                                match_total_bases = re.search( r'total_bases="<span>(.*?)</span>"', match)
                                total_bases = match_total_bases.group(1) if match_total_bases else np.nan
                                # Gets the total number of bases of the run
                                total_bases = int(total_bases) 
                                match_spots = re.search(r'total_spots="(.*?)"', match)
                                if match_spots:
                                    total_spots_1 = match_spots.group(1)
                                    total_spots = re.sub(r'</?span>', '', total_spots_1)
                                    # Gets total number of reads in the run
                                    total_spots = float(total_spots)
                                else:
                                    total_spots = np.nan 
                                match_average = re.search(r'average="(.*?)"', match)
                                if match_average:
                                    average_1 = match_average.group(1)
                                    average = re.sub(r'</?span>', '', average_1)
                                    # Gets the total number of bases (average) per read
                                    average = float(average) 
                                else:
                                    average = np.nan
                                match_library = re.search(r'SINGLE|PAIRED', match)
                                # Gets if is Single or Paired
                                library_layout = match_library.group() if match_library else np.nan 
                                break

                            position += 1
                        
                    # Prioritising experiments which has multiple experiments and the biggest run is taken only the bigger one.         
                    elif max_gb_multiple_runs and max_gb_multiple_runs[0] > max_gb_multiple_runs[1] * 10:
                        position=0
                        for match in match_experiment:
                            if position == index_max_gb_mult_run:
                                matchs_srr_all = re.findall(r'RUN <(.*?)&lt;IDENTIFIERS&gt', match)
                                data_2 = []
                                count = 0
                                max_gb_per_run = False
                                for match_srr_all in matchs_srr_all:
                                    match_srr_bytes = re.search(r'size="<span>(.*?)</span>', match_srr_all)
                                    num_bytes = int(match_srr_bytes.group(1))
                                    num_gb = round(num_bytes / (1024 ** 3), 3)
                                    data_2.append((num_gb, count))
                                    print(num_gb)
                                    count += 1
                            
                                max_gb_per_run = max(data_2, key=lambda x: x[0])
                                index_max_gb_per_run = data_2.index(max_gb_per_run)
                                print("Bigger run:", max_gb_per_run, "position:" , index_max_gb_per_run)

                                if max_gb_per_run and max_gb_per_run[0] > 20:
                                    position = 0
                                    for match_srr_all in matchs_srr_all:
                                        if position == index_max_gb_per_run:

                                            match_srr = re.search(r'accession="<span>(.?RR.*?)</span>"',  match_srr_all)
                                            # Gets the SRR value
                                            srr_value = match_srr.group(1) 
                                            print('SRR selected:', srr_value)                                  
                                            match_instrument = re.search(r'<b>&lt;INSTRUMENT_MODEL&gt;</b>(.*?)<b>&lt;/INSTRUMENT_MODEL&gt;</b>', match_srr_all)
                                            # Gets the instrument model (eg. illumina)
                                            instrument_model = match_instrument.group(1) if match_instrument else np.nan 

                                            match_total_bases = re.search( r'total_bases="<span>(.*?)</span>"', match_srr_all)
                                            total_bases = match_total_bases.group(1) if match_total_bases else np.nan
                                            # Gets the total number of bases of the run
                                            total_bases = int(total_bases) 
                                            print('total bp:', total_bases)
                                            match_spots = re.search(r'total_spots="(.*?)"', match_srr_all)
                                            if match_spots:
                                                total_spots_1 = match_spots.group(1)
                                                total_spots = re.sub(r'</?span>', '', total_spots_1)
                                                # Gets total number of reads in the run
                                                total_spots = int(total_spots) 
                                            else: 
                                                total_spots = np.nan
                                            match_average = re.search(r'average="(.*?)"', match_srr_all)
                                            if match_average:
                                                average_1 = match_average.group(1)
                                                average = re.sub(r'</?span>', '', average_1)
                                                # Gets the total number of bases (average) per read
                                                average = int(average) 
                                            else:
                                                average = np.nan
                                            match_library = re.search(r'SINGLE|PAIRED', match_srr_all) 
                                            # Gets if is Single or Paired
                                            library_layout = match_library.group() if match_library else np.nan 
                                            break

                                        position+=1
                                
                                else:
                                    instrument_model = np.nan 
                                    srr_value = np.nan 
                                    total_spots =np.nan 
                                    library_layout =np.nan 
                                    average = np.nan 
                                    total_bases = np.nan 
                                break
                            position+=1
                    else:
                        instrument_model = np.nan
                        srr_value = np.nan
                        total_spots =np.nan
                        library_layout =np.nan
                        average = np.nan
                        total_bases = np.nan
                else:
                    instrument_model = np.nan
                    srr_value = np.nan
                    total_spots =np.nan
                    library_layout =np.nan
                    average = np.nan
                    total_bases = np.nan
            else:
                instrument_model =np.nan
                srr_value = np.nan
                total_spots =np.nan
                library_layout =np.nan
                average = np.nan
                total_bases = np.nan
        
            
            yield assembly_method, common_name, accession,  genome_coverage, total_bases, genome_size,  library_layout, sra_value,  srr_value,  instrument_model, total_spots, average

# Extraction of longevity form the url: https://animaldiversity.org

def extractionLongevity(url):

    response = requests.get(url)
    if response.status_code == 200: 
        soup = BeautifulSoup(response.content, "html.parser")
        try:                
            lifespan_element = soup.find("h3", id="lifespan_longevity").find_next("dd")
            lifespan_text = lifespan_element.get_text(strip=True)
            words = lifespan_text.split()
            for word in words:
                try:
                    years = float(word.strip())
                    print(years)
                    return years
                except ValueError:
                    pass
        except AttributeError:
            pass

        
    else:
        print("Webpage not accessible") 

# This unifies all the function and stores all the data in a df
def main():
    with open('zoonomia.html', 'r', encoding='utf-8') as file:
        html_content = file.read()
    print('Searching for data')
    df = parse_html(html_content)
    
    df['Common_name'] = ""  
    df['Accession_number'] = ""
    df['Genome_coverage'] = ""
    df['Genome_size'] = ""
    df['AssemblyMethod'] = ""    
    df['SRA_number'] = ""
    df['SRR_number'] = ""
    df['Library_layout'] = "" 
    df['Instrument_model'] = ""
    df['Total_spots_run'] = ""
    df['Total_bases_run'] = ""
    df['Average_bp_read'] = ""
    
    temp_folder = "temp_folder"
    os.makedirs(temp_folder, exist_ok=True)
    # https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/
    print('Searching SRA and SRR numbers and data')
    for index, row in df.iterrows():
        gc_number = row['NCBI_AssemblyAccNo']
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{gc_number}/download?&filename={gc_number}.zip"
        
        download_and_extract_zip(url, temp_folder)

        jsonl_file_path = os.path.join(temp_folder, "ncbi_dataset", "data", "assembly_data_report.jsonl")

        for assembly_method, common_name, accession,  genome_coverage, total_bases, genome_size,  library_layout, sra_value,  srr_value,  instrument_model, total_spots, average in extract_values(jsonl_file_path):
            
            df.at[index, 'Common_name'] = common_name
            df.at[index, 'Accession_number'] = accession
            df.at[index, 'Genome_coverage'] = genome_coverage
            df.at[index, 'Genome_size'] = genome_size
            df.at[index, 'AssemblyMethod'] = assembly_method
            df.at[index, 'SRA_number'] = sra_value
            df.at[index, 'SRR_number'] = srr_value
            df.at[index, 'Library_layout'] = library_layout
            df.at[index, 'Instrument_model'] = instrument_model           
            df.at[index, 'Total_spots_run'] = total_spots
            df.at[index, 'Total_bases_run'] = total_bases
            df.at[index, 'Average_bp_read'] = average

    print('Searching Longevity')
    longevities = []
    for index, row in df.iterrows():
        species = str(row['Species'])
        species = species.replace(" ", "_")
        url = f"https://animaldiversity.org/accounts/{species}"
        longevity = extractionLongevity(url)
        longevities.append(longevity)

    df['Longevity'] = longevities
    df.to_excel('zoonomia_sp_info.xlsx', index=False, engine='openpyxl')
    shutil.rmtree(temp_folder)
    

if __name__ == "__main__":
    main()