#!/bin/bash
# Nombre del trabajo
#$ -t 1-250
#$ -tc 6
#$ -N RAW_2_TELO_zoo
# Salida estándar y de error
#$ -e error_telo_zoo.$TASK_ID.log
#$ -o output_telo_zoo.$TASK_ID.log
# Cola de trabajos 
#$ -P AG
# cores necesarios.
#$ -pe pthreads 6
# memoria RAM requerida
#$ -l h_vmem=6G
#$ -A bioinformatica

export PYTHONPATH=
unset PKG_CONFIG_PATH
unset LD_LIBRARY_PATH
unset R_HOME
unset JAR_HOME
unset PERL5LIB
unset SAMTOOLS
.  /home/jrodriguezdv/miniconda3/bin/activate
conda activate telomere_env
python Telomere_analysis.py

