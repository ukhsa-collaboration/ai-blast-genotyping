# AI BLAST genotyping

command line tool for running the BLAST genotyping method 

# Requirements & Dependencies
To run this tool you will need to numpy, pandas and BLAST installed. The versions are present in the requirements.txt file. 
Installation can be performed with either pip or conda. 

# Running on the HPC

To run this script on the HPC you will need to load pre-existing modules, and have copied over and un-zip the gitlab code for both the AI BLAST genotyping gitlab repo. 
You do not have to install further python modules or use the requirements.txt due to the limitations of updating python code on the HPC. 

```
# standard setup to get off the head node
module load sge
module load qrsh 
# specific requirements
module load sge blast+/2.2.27
```


## Files required
1. **BLAST database** - You will also need to have a local copy of the BLAST database that has been created based on the reference strains for the UK AI genotypes. 
This can be copied from the HPC:

### THE INPUT FILES HAVE BEEN UPDATED SO PLEASE COPY ACROSS THE LATEST VERSION FROM THE HPC

[sftp://hpc1mgmt.unix.phe.gov.uk/phengs/hpc_projects/nicc80_ai/geno_blastdb](url)

You will need the files all files labelled "ai_geno_groups_ref.*"

2. **Genotypes reference table** - "genotype_groups_examples.csv" copy also available in the same HPC folder. 

The reference table should be placed in the main repo folder to be accessed by the genotyping script

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

****Process all FASTA files (based on extension) in a folder or all samples in single FASTA file****
```

python /path/ai-blast-genotyping/ai_genotyping_tool.py -o OUTPUT_DIR -t {yes,no} -e FASTA_EXTENSION [--input_file INPUT_FILE OR --input_folder INPUT_FOLDER] -b BLASTDB –n FILETAG 
```
## Assumptions:

Assumes that the FASTA header ends with the segment: "|HA" or "_HA"

## Test data from the UK reference sequences

A test data set for the UK reference genomes is available on the HPC folder in the location below
sftp://hpc1mgmt.unix.phe.gov.uk/phengs/hpc_projects/nicc80_ai/geno_blastdb/test_data




# Details of how the BLAST database was set up & initial testing

[https://phecloud.sharepoint.com/:p:/r/teams/GenomicsCell/Shared%20Documents/Avian_Flu/Genotyping/AI_genotyping_20230209%20-%20Copy.pptx?d=wb923b858ed304b7480dbe572050270ab&csf=1&web=1&e=u6aFJT](url)