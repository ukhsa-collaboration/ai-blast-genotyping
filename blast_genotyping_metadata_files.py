# Concatenated tree prep for mutation scan 
#
# outputs: per cds alignment, table of any excluded sequences
# use mafft 



from pathlib import Path
import argparse
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys
import logging
from datetime import datetime, date
import pandas as pd
from collections import Counter
from Bio import SeqIO
import subprocess


now = date.today()

__version__ = 0.1
__author__ = 'Kate Howell'

def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"Script to create AI segment alignments for mutation scan")
    parser.add_argument('--input_file', '-i', required=True, help='CSV file retrieved from AI database.')
    parser.add_argument('--blast_summary', '-b', required=True, help='BLAST summary file')
    parser.add_argument('--output_dir', '-o', required=True, default=os.getcwd(), help='Output folder. Default: CWD.')
    parser.add_argument('--tagname', '-n', required=True, help='Analysis name for output files')

    # Add threads
    args = parser.parse_args()

    # Need to handle output dir before setting up logging files.
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)

    return args

def logging_file_setup(output_folder):
    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    datetime_stamp = now.strftime("%Y%m%d-%H%M%S")
    logging_file_output = os.path.join(output_folder, str(str(now.strftime("%Y%m%d")) + str('_logging_file.log')))
    print("Logging output in the file: {}".format(logging_file_output))
    logging.basicConfig(
            filename=logging_file_output,  # filename=logging_file_output,
            filemode='w',
            format='%(asctime)s:%(levelname)s:%(message)s',
            level=logging.INFO)
    print("Logging Level: INFO")
    logging.info("Logging Started.")

def check_arguments(args):
	if not os.path.exists(args.input_file):
		print(f'File does not appear to exists: {args.input_file}. Please check')
		return 1
	if not os.path.isdir(args.output_dir):
		os.mkdir(args.output_dir)
	return 0


path = os.path.dirname(sys.argv[0])

def subset_metadata(g, newtab,filetag):
     subblasttab = newtab[newtab['Top_Hit']==g]
     subblasttab.to_csv(os.path.join(args.output_dir,str(f'{filetag}_{now}_dbextract_metadata_{g}.csv')))
     

def main(args):
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    output_dir = os.path.abspath(args.output_dir)
    print(output_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # Set up logging
    logging_file_setup(output_dir)

    #main steps of script
    # 1. Input metadata file and create fasta file for each segment
    filetag = args.tagname
    metatab = pd.read_csv(args.input_file)
    metatab['segment_name'] = metatab['segment_name'].fillna('NA')

    blasttab = pd.read_csv(args.blast_summary)

    newtab = pd.merge(metatab,blasttab, left_on='isolate_epi_id', right_on='isolate_epi_id',how='left')
    newtab.to_csv(os.path.join(args.output_dir,str(f'{filetag}_{now}_dbextract_metadata_blast.csv')))
    print("Genotyping results summary:")
    print(Counter(blasttab['Top_Hit']))
    genotypes = list(set(blasttab['Top_Hit']))
    for g in genotypes:
         subset_metadata(g, newtab,filetag)


    print("Genotype group metadata files now produced!")
    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f'Pipeline time to completion: {start_time - end_time}')


if __name__ == '__main__':
    args = read_commandline()
    # Setting up testing and logging
    logging_file_setup(args.output_dir)
    check = check_arguments(args)
    if check == 1:
        sys.exit(logging.error("Arguments provided were not expected. Please check log."))
    sys.exit(main(args))
