#!/bin/bash



mkdir -p /home/jrodriguezdv/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/jrodriguezdv/miniconda3/miniconda.sh
bash /home/jrodriguezdv/miniconda3/miniconda.sh -b -u -p /home/jrodriguezdv/miniconda3
rm -rf /home/jrodriguezdv/miniconda3/miniconda.sh
miniconda3/bin/conda init bash
# Creación del environment
conda create -n telomere_env python=3.12
conda activate telomere_env
conda install -c bioconda sra-tools
conda install biopython
conda install jsonlines
conda install openpyxl
conda install pandas
conda install requests
conda install BeautifulSoup
conda install seaborn
conda install gprofiler-official
conda install xlrd
conda install sklearn

