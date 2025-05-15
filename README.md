# AI BLAST genotyping

This repo contains the command line tool for running the universal BLAST genotyping method as well as helper scripts to process new genotype references and generate FASTA files for genotyping. 

# Requirements & Dependencies
To run this tool you will need to numpy, pandas and BLAST installed. The versions are present in the requirements.txt file. 
Installation can be performed with conda:
```
conda install --file requirements.txt
```


## Files required
There are a series of required files for either adding new genotypes or running the BLAST genotyping. These are included in the reference_files folder within the repo and include:
- The latest BLAST db - all_genotype_references.db.*
- The FASTA file from which the BLAST db is built (all_genotype_references.fasta)
- the constellation table - blast_geno_threshold_table98.csv
- a genotyping key for the meta-data - genotype.key.csv


### CHANGE LOG

- The universal genotyping tool no longer requires the "genotype_groups_examples.csv" table and instead uses the constellations table blast_geno_threshold_table98.csv. 
- The method itself has changed to incorporate groups per segment which are matched up to a genotype constellation (combination of groups across the segments). Some genotype constellations are known duplicates. 
- The genotype prioritisation based on prevalence is no longer part of the universal genotyping process as is not feasible with the current data availability. The top BLAST hit is instead reported in addition to the genotype constellation, which is matched back to the genotype. 

# Adding new genotypes

### Basic process
1. If adding GenoFlu references there is a script to prepare the data for that: 

```
python genoflu_prep.py -h
usage: genoflu_prep.py [-h] --genoflu GENOFLU --previous_geno_tab PREVIOUS_GENO_TAB --output_dir OUTPUT_DIR

New US genotype prep tool

optional arguments:
  -h, --help            show this help message and exit
  --genoflu GENOFLU, -g GENOFLU
                        Genoflu folder
  --previous_geno_tab PREVIOUS_GENO_TAB, -p PREVIOUS_GENO_TAB
                        Previous genotype key from ai-blast-genotyping prep
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder location

```
The FASTA file can then be input into step 3. 

2. If adding other references you will need you can create the FASTA manually but the FASTA header must be in the following format: Sequence_identifier*|Genotype|Subtype|Scheme|Segment 
If the isolate data is already in aiseqdb skip to step 3. 

3. Run the new_genotyping_prep.py script using a BLAST identity threshold of **98%**

```
python new_genotype_prep.py -h
usage: new_genotype_prep.py [-h] --geno_key GENO_KEY --reffasta REFFASTA --output_dir OUTPUT_DIR
                            [--epiid EPIID | --fasta FASTA] [--username USERNAME] --tagname TAGNAME --threshold
                            THRESHOLD

Script to update BLAST genotyping schemas

optional arguments:
  -h, --help            show this help message and exit
  --geno_key GENO_KEY, -g GENO_KEY
                        existing genotyping key
  --reffasta REFFASTA, -r REFFASTA
                        existing genotype reference fasta
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder. Default: CWD.
  --epiid EPIID, -e EPIID
                        CSV file of new genotype sequences with meta-data
  --fasta FASTA, -f FASTA
                        Input FASTA file for new genotype references (SPECIFIC INPUT FORMAT REQUIRED)
  --username USERNAME, -u USERNAME
                        username for aiseqdb query
  --tagname TAGNAME, -n TAGNAME
                        Analysis name for output files
  --threshold THRESHOLD, -t THRESHOLD
                        Threshold for BLAST filtering
```

4. Once you have reviewed the outputs and are happy you need to add the latest version of the BLASTdb and the blast_genotyping_threshold_table98.csv to this git repo.

# BLAST genotyping command line tool 

### Overall process

1. *NOTE: This is an internal script to retrieve sequences from bespoke database.* First you will need to create a FASTA file to run the BLAST genotyping on using the prep_for_genotyping.py script. 

```
python prep_for_genotyping.py -h
usage: prep_for_genotyping.py [-h] --input_file INPUT_FILE --output_dir OUTPUT_DIR --tagname TAGNAME

Script to create AI segment alignments for mutation scan

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        CSV file retrieved from AI database.
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder. Default: CWD.
  --tagname TAGNAME, -n TAGNAME
                        Analysis name for output files

```
2. Run the BLAST genotyping tool on your new FASTA file: 

``` 
python ai_genotyping_tool.py -h
usage: ai_genotyping_tool.py [-h] --output_dir OUTPUT_DIR [--testing {yes,no}] --tagname TAGNAME
                             [--extension EXTENSION] [--input_file INPUT_FILE | --input_folder INPUT_FOLDER]
                             --blastdb BLASTDB --strict STRICT --identity IDENTITY

AI UK Genotyping command line tool

optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder. Default: CWD.
  --testing {yes,no}, -t {yes,no}
                        Debugging mode. Specify by either "yes" or "no"
  --tagname TAGNAME, -n TAGNAME
                        Analysis name for output files
  --extension EXTENSION, -e EXTENSION
                        FASTA file extension if not default
  --input_file INPUT_FILE, -i INPUT_FILE
                        Input FASTA file
  --input_folder INPUT_FOLDER, -f INPUT_FOLDER
                        Input FASTA file
  --blastdb BLASTDB, -b BLASTDB
                        Reference BLAST database
  --strict STRICT, -s STRICT
                        Strict version, all 8 segments for genotype
  --identity IDENTITY, -d IDENTITY
                        Percentage identity threshold. Default = 98.
```

```
## Example usage

python /path/ai-blast-genotyping/ai_genotyping_tool.py -o OUTPUT_DIR -t {yes,no} -e FASTA_EXTENSION [--input_file INPUT_FILE OR --input_folder INPUT_FOLDER] -b BLASTDB –n FILETAG 

```

## Assumptions:

Assumes that the FASTA header ends with the segment: "|HA" or "_HA" and will check for these delimiters & otherwise will fail. 

