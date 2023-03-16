# AI BLAST genotyping

command line tool for running the BLAST genotyping method 

# Requirements
To run this tool you will need to have BLAST locally installed. 
Recommend installing with conda: 

``
conda install -c bioconda blast
``


You will also need to have a local copy of the BLAST database that has been created based on the reference strains for the UK AI genotypes. 
This can be copied from the HPC:

[sftp://hpc1mgmt.unix.phe.gov.uk/phengs/hpc_projects/nicc80_ai/geno_blastdb](url)

You will need the files all files labelled "ai_geno_ref.*"

If you do not have access to the folder then please discuss access with the project lead. 

# Command line tool 

```python ai_genotyping_tool.py -h 

usage: ai_genotyping_tool.py [-h] --output_dir OUTPUT_DIR [--testing {yes,no}] 
[--extension EXTENSION] 
[--input_file INPUT_FILE | --input_folder INPUT_FOLDER] --blastdb BLASTDB

AI UK Genotyping command line tool

`optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder. Default: CWD.
  --testing {yes,no}, -t {yes,no}
                        Debugging mode. Specify by either "yes" or "no"
  --extension EXTENSION, -e EXTENSION
                        FASTA file extension if not default
  --input_file INPUT_FILE, -i INPUT_FILE
                        Input FASTA file
  --input_folder INPUT_FOLDER, -f INPUT_FOLDER
                        Input FASTA file
  --blastdb BLASTDB, -b BLASTDB
                        Reference BLAST database
                        `

```


## Example usage

****Process all FASTA files (based on extension) in a folder****

python ai_genotyping_tool_v01.py -e fasta -f ../gisaid -b ai_geno_groups_ref -o ../gisaid/genotyping_testing


****Process all samples in single FASTA file****

python ai_genotyping_tool_v01.py -i ../gisaid/test_seq_20230307.fasta -b ai_geno_groups_ref -o ../gisaid/genotyping_testing

## Assumptions:

Assumes that the FASTA header ends with the segment: "|HA" or "_HA"

## Features to add still

1. Testing data
2. Historic table for joint genotyping calls


# Details of how the BLAST database was set up & initial testing

[https://phecloud.sharepoint.com/:p:/r/teams/GenomicsCell/Shared%20Documents/Avian_Flu/Genotyping/AI_genotyping_20230209%20-%20Copy.pptx?d=wb923b858ed304b7480dbe572050270ab&csf=1&web=1&e=u6aFJT](url)