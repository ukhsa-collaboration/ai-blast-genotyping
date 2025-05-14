import argparse
import logging
import os
import sys

import pandas as pd

from datetime import date, datetime

import subprocess
import src.utilities as util 
import src.fasta_processing as fap
import src.blast_wrangle as bl
import src.genotype_groups as geno

now = date.today()

__version__ = 1.0
__author__ = "Kate Howell"



def read_commandline():
    """
    BLAST genotyping setup tool - command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(
        description=f"Script to update BLAST genotyping schemas"
    )
    parser.add_argument(
        "--geno_key", "-g", required=True, help="existing genotyping key",default = "reference_files/genotype_key.csv"
    )
    parser.add_argument(
        "--reffasta", "-r", required=True, help="existing genotype reference fasta",default = "reference_files/all_genotype_references.fasta"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        required=True,
        default=os.getcwd(),
        help="Output folder. Default: CWD.",
    )
    input_group = (
        parser.add_mutually_exclusive_group()
    )  # input group requires one or the other
    input_group.add_argument(
        "--epiid", "-e", required=False, help="CSV file of new genotype sequences with meta-data"
    )
    input_group.add_argument(
        "--fasta", "-f", required=False, help="Input FASTA file for new genotype references (SPECIFIC INPUT FORMAT REQUIRED)"
    )
    parser.add_argument(
        "--username", "-u", required=False, help="username for aiseqdb query"
    )
    parser.add_argument(
        "--tagname", "-n", required=True, help="Analysis name for output files"
    )
    parser.add_argument(
        "--threshold", "-t", required=True, help="Threshold for BLAST filtering"

    )
    # Add threads
    args = parser.parse_args()

    # Need to handle output dir before setting up logging files.
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)

    return args


repopath = os.path.dirname(sys.argv[0])
        

def main(args):
    """
    Main running of the script to run the BLAST query and wrangle the results to provide a per segment and sample summary of the genotyping results.

    :return: N/A
    """
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    global output_dir
    output_dir = os.path.abspath(args.output_dir)

    # Set up logging
    util.logging_file_setup(output_dir)
    # Read in required files depending if CSV or FASTA provided
    outputfasta = os.path.join(args.output_dir,f'all_genotype_references_{now}.fasta')
    if args.reffasta == "reference_files/all_genotype_references.fasta":
        subprocess.call(f"cp {os.path.join(os.path.dirname(sys.argv[0]),args.reffasta)} {outputfasta}",shell=True)
    else:
        subprocess.call(f"cp {args.reffasta} {outputfasta}",shell=True)
    
    if args.epiid is not None:
        epitab = pd.read_csv(args.epiid)
        ids = "','".join(epitab['sequence'])
        print(f"{len(epitab['sequence'])} new genotypes to be added to the database")
        db_extract = util.query_aiseqdb(args.username, ids)
        for i in epitab['sequence']:
            subepitab = db_extract[db_extract['isolate_epi_id'] == i]
            print(f'{i} has {subepitab.shape[0]} segments available')
            if subepitab.shape[0] < 6:
                print(f'{i} has less than 6 segments available. Constellation will be skipped. Genotype NOT added.')
                logging.info(f'{i} has less than 6 segments available. Constellation will be skipped. Genotype NOT added.')
                indexseq = epitab[epitab['sequence'] == i].index
                epitab.drop(indexseq, inplace=True)
            elif subepitab.shape[0] < 8:
                print(f'{i} does not have all 8 segments available. Constellation will be incomplete')
        db_extract['isolate_name']=db_extract['isolate_name'].str.replace(" ","_")
        db_extract2 = pd.merge(epitab,db_extract,right_on='isolate_epi_id',left_on='sequence',how="left")
        fastafile, filtereddf = fap.create_fasta_seqs(args,now,db_extract2)
        fap.per_segment_fasta(now,util.segments, fastafile, output_dir)
        epitab = epitab[epitab['sequence'].isin(list(filtereddf['isolate_epi_id']))]
        subepitab = epitab[['Genotype', 'sequence', 'subtype', 'schema']]
        if args.geno_key == "":
            genotab = geno.update_geno_key(args,now,os.path.join(os.path.dirname(sys.argv[0]),args.geno_key), subepitab)
        else:
            genotab = geno.update_geno_key(args,now,args.geno_key,subepitab)

    if args.fasta is not None:
        newgenotab = fap.new_fasta_parsing(args,now,args.fasta)
        print(f"{len(newgenotab['sequence'])} new genotypes to be added to the database")
        if args.geno_key == "":
            genotab = geno.update_geno_key(args,now,os.path.join(os.path.dirname(sys.argv[0]),args.geno_key), newgenotab)
        else:
            genotab = geno.update_geno_key(args,now,args.geno_key, newgenotab)
        epitab = newgenotab
    try:
        blasttab =bl.blast_processing(args,now)
    except:
        print("Unable to run the BLAST processing, exiting...")
        sys.exit()
    blast,queries = bl.filter_blast_results(args,now,blasttab)
    geno.geno_groups(queries, blast,int(args.threshold),epitab['Genotype'])
    geno.files_for_gitrepo(now,output_dir)
    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: {start_time - end_time}")


if __name__ == "__main__":
    args = read_commandline()

    # Setting up testing and logging
    util.logging_file_setup(args.output_dir)
    check = util.check_arguments(args)
    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )
    sys.exit(main(args))