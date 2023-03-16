


'''
### AI Genotyping command line tool

## Overview of the tool:
# For each sample we need to blast the sequence of the 8 segments against the reference blast database and get the top genotype hit for each segment
# BLAST tabular output format means that we can automate the results processing
# Wrangle the results of the blast search results together and make a call about the genotype
# Compare the genotypes to the results that we have from the existing sequence database 

'''

import os, glob,sys
import pandas as pd
import numpy as np
import logging
import argparse
from datetime import datetime
import shutil
import subprocess
from collections import Counter

now = datetime.now()
x_int = int(now.strftime("%Y%m%d%H%M%S"))

__version__ = 0.1
__author__ = 'Kate Howell'


def logging_file_setup(output_folder, testing):
    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    datetime_stamp = now.strftime("%Y%m%d-%H%M%S")
    logging_file_output = os.path.join(output_folder, str(str(now.strftime("%Y%m%d")) + str('_logging_file.log')))
    print("Logging output in the file: {}".format(logging_file_output))

    if testing == 'yes':
        logging.basicConfig(
            filename=logging_file_output,  # filename=logging_file_output,
            filemode='w',
            format='%(asctime)s:%(levelname)s:%(message)s',
            level=logging.DEBUG)
        print("Logging Level: DEBUG")
    else:
        logging.basicConfig(
            filename=logging_file_output,  # filename=logging_file_output,
            filemode='w',
            format='%(asctime)s:%(levelname)s:%(message)s',
            level=logging.INFO)
        print("Logging Level: INFO")
    logging.info("Logging Started.")


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"AI UK Genotyping command line tool")
    parser.add_argument('--output_dir', '-o', required=True, default=os.getcwd(),
                        help='Output folder. Default: CWD.')
    parser.add_argument('--testing', '-t', required=False, default=False, choices=['yes', 'no'],
                        help='Debugging mode. Specify by either "yes" or "no"')  # Change this to a list of options
    parser.add_argument('--extension', '-e', required=False, default='.fasta', help='FASTA file extension if not default')
    input_group = parser.add_mutually_exclusive_group() # input group requires one or the other
    input_group.add_argument('--input_file', '-i', required=False, help='Input FASTA file')
    input_group.add_argument('--input_folder', '-f', required=False, help='Input FASTA file')
    parser.add_argument('--blastdb', '-b', required=True, help='Reference BLAST database')
    args = parser.parse_args()
    print(args)
    # Need to handle output dir before setting up logging files.
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)

    return args


def check_arguments(args):
    """
    Check that paths provided exist and create output folder if it doesn't exist already.

    :param args: output from arg parse
    :return: status
    """
    # Check if input file or input folder are present
    if args.input_folder is not None:
        if not os.path.exists(args.input_folder):
            logging.error(f'File does not appear to exists: {args.input_folder}. Please check')
            return 1

    if not os.path.exists(args.blastdb+".nin"):
        logging.error(f'BLAST database does not appear to exists: {args.blastdb}. Please check')
        return 1
    return 0


def testing_functions(testing_functions_parameters):
    # Generate output folder
    if testing_functions_parameters[1] == str('output_folder'):
        if testing_functions_parameters[0] == True:
            if os.path.exists(testing_functions_parameters[2]):
                print('Deleting output_folder {}'.format(testing_functions_parameters[2]))
                try:  ## Try to remove tree; if failed show an error using try...except on screen
                    shutil.rmtree(testing_functions_parameters[2])
                except OSError as e:
                    print("Error: %s - %s." % (e.filename, e.strerror))
                os.mkdir(testing_functions_parameters[2])
            else:
                os.mkdir(testing_functions_parameters[2])
        elif testing_functions_parameters[0] == False:
            # Create output folder
            if not os.path.exists(testing_functions_parameters[2]):
                print(
                    "Output folder does not exist. Attempting to create folder: %s" % (testing_functions_parameters[2]))
                os.mkdir(testing_functions_parameters[2])


# Reference data & input parameters - move to Lib folder
reftabdict = {'PB2_Group1': 'A/chicken/England/053052/2021', 'PB2_Group2': 'A/chicken/Wales/053969/2021', 
'PB2_Group3': 'A/chicken/Scotland/054477/2021', 'PB1_Group1': 'A/chicken/England/053052/2021', 'PB1_Group2': 'A/chicken/England/152082/2022',
 'PA_Group1': 'A/chicken/England/053052/2021', 'PA_Group2': 'A/chicken/Scotland/054477/2021', 
 'PA_Group3': 'A/herring_gull/England/324803/2022', 'PA_Group4': 'A/Humboldt_penguin/England/161651/2022', 
 'HA_Group1': 'A/chicken/England/053052/2021', 'HA_Group2': 'A/Greylag_goose/England/054503/2021',
  'NP_Group1': 'A/chicken/England/053052/2021', 'NP_Group2': 'A/turkey/England/016515/2022', 'NP_Group3': 'A/chicken/England/069816/2021',
   'NP_Group4': 'A/chicken/England/085598/2022', 'NP_Group5': 'A/chicken/England/118935/2022', 'NA_Group1': 'A/chicken/England/053052/2021', 
   'NA_Group2': 'A/herring_gull/England/324803/2022', 'M_Group1': 'A/chicken/England/053052/2021', 
   'NS_Group1': 'A/chicken/England/053052/2021', 'NS_Group2': 'A/chicken/England/085598/2022', 'NS_Group3': 'A/pheasant/England/251536/2022'}

genogroups = pd.read_csv("~/Documents/avian_influenza/apha/genotype_groups_examples.csv",dtype=str)
genogroups['Segment'] = genogroups['Segment'].replace(np.nan, "NA")

groups = genogroups['Group'].tolist()
genotypes = genogroups['Genotypes'].tolist()
genodict = dict(zip(genogroups.Labels, genogroups.Genotypes))
blast_cols = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
threshold = 90

def tidy_fasta_files(sample):
    logging.info("Tidying FASTA headers for the following non-standard characters.... () ,-\n")
    os.system("perl -pi -e 's/\(/_/g' "+sample)
    os.system("perl -pi -e 's/\)/_/g' "+sample)
    os.system("perl -pi -e 's/ /_/g' "+sample)
    os.system("perl -pi -e 's/,/_/g' "+sample)
    os.system("perl -pi -e 's/--//g' "+sample)
    os.system("perl -pi -e 's/-\n-//g' "+sample)
    
# function to run the blast searches to temp file

def run_blast(db,folder,sample):
    tidy_fasta_files(sample)
    logging.info("Performing the BLAST searches per FASTA file")
    subprocess.call("blastn -db "+db+" -query "+os.path.join(folder,sample)+" -out "+os.path.join(args.output_dir,sample)+".blast.out -outfmt 6 -max_target_seqs 1",shell=True)

def match_genotype_dict(blast_pass):
    segdict = dict(zip(genogroups["Segment"]+"_"+genogroups.Example_sequence, genogroups.Labels))
    #map the query sequence back to the genotyping table
    blast_pass["group_match"] = blast_pass['ref_match'].apply(lambda x: segdict.get(x)).fillna('')
    blast_pass["genotype_match"] = blast_pass['group_match'].apply(lambda x: genodict.get(x)).fillna('')
    return blast_pass

def tidy_sample_name(table,sample_name):
    test = table[sample_name].str.split(pat = "|",expand = True)
    if test.shape[1]==3:
        table[sample_name] = test[test.columns[0]]+"|"+test[test.columns[1]]
    if test.shape[1]==2:
        table[sample_name] = test[test.columns[0]]
    else:
        logging.info("Unsure of sample name format so not performing tidying step")
    return table

def tidy_blast_table(folder,sample):
    logging.info("Reading in BLAST output:")
    logging.info(os.path.join(args.output_dir,sample+".blast.out"))
    blasttab = pd.read_csv(os.path.join(args.output_dir,sample+".blast.out"),sep="\t",header = None)
    blasttab.columns = blast_cols
    blast_pass = blasttab[blasttab["pident"] >= threshold]
    blast_fail = blasttab[blasttab["pident"] < threshold]
    if blast_fail.shape[0]>= 1:
        logging.info(blast_fail.shape[0]," sequences do not meet a minimum percentage sequence identity to the ref seq!!")
    logging.info("BLAST pass table dimensions:")
    logging.info(blast_pass.shape)
    if "|" in blast_pass['qseqid'][0]:
        test = blast_pass['qseqid'].str.split(pat = "|",expand = True)
        blast_pass['segment'] = test[test.columns[-1]]
    else: 
        test = blast_pass['qseqid'].str.split(pat = "_",expand = True)
        blast_pass['segment'] = test[test.columns[-1]]
    blast_pass['ref_match'] = blast_pass['segment']+"_"+blast_pass['sseqid'].map(lambda x: x.split('|')[0]).str.replace("_","")
    blast_pass = match_genotype_dict(blast_pass)
    blast_pass.to_csv(os.path.join(folder,sample+".blast.out2"))
    results_df = pd.DataFrame(blast_pass['qseqid'].to_list(), columns=['sample',])
    results_df['genotype_match'] = blast_pass['genotype_match']
    results_df['segment'] = blast_pass['segment']
    if blast_fail.shape[0]>=1:
        fail = blast_fail['qseqid'].to_list()
        for f in fail:
            missing =[[f,"NA","NA"]]
            results_df = results_df.append(pd.DataFrame(missing, columns=results_df.columns),ignore_index=True)
    results_df.sort_values(by=['sample'],inplace=True)
    results_df = tidy_sample_name(results_df,'sample')

    return(results_df)


def run_full_blasts(folder,mode,extension):
    results_tabs = []
    logging.info('Input folder provided:')
    logging.info(os.path.join(folder,"*"+extension))
    fastas = glob.glob(os.path.join(folder,"*"+extension))
    logging.info("Fasta files found:\n")
    logging.info(fastas)
    if mode == "all_in_folder":   
        for f in fastas:
            logging.info("Running BLAST")
            run_blast(args.blastdb,folder,f)
            logging.info("Reviewing BLAST results")
            sresults = tidy_blast_table(folder,f)
            logging.info("Combining results across FASTA files...")
            results_tabs.append(sresults)
            newdf = pd.concat(results_tabs)
            newdf.to_csv(os.path.join(args.output_dir,str(x_int)+"_summary.csv"))
            logging.info("Output file written to:")
            logging.info(os.path.join(args.output_dir,str(x_int)+"_summary.csv"))
        return(newdf)
    
    elif mode =="single":
        logging.info("Running BLAST on single FASTA file")
        head_tail = os.path.split(folder)
        logging.info(head_tail)
        run_blast(args.blastdb,head_tail[0],os.path.join(head_tail[0],head_tail[1]))
        logging.info("Reviewing BLAST results")
        sresults = tidy_blast_table(head_tail[0],head_tail[1])
        sresults.to_csv(os.path.join(args.output_dir,str(x_int)+"_summary.csv")) 
        logging.info("Output file written to:")
        logging.info(os.path.join(args.output_dir,str(x_int)+"_summary.csv"))
        return sresults
    
def create_persample_summary(summarytab):
    pivot = pd.pivot_table(summarytab, values='genotype_match', 
                                index='sample', 
                                columns='segment', fill_value = "N/A", aggfunc='first')
    consensus = []
    freq=[]
    for index, row in pivot.iterrows():
            newlist = []
            topgeno = []
            for s in pivot.columns:
                try:
                    test = row[s].split("|")
                    for t in test:
                        newlist.append(t)
                except:
                    newlist.append("")
            genfreq = Counter(newlist)
            counts = list(genfreq.values())
            maxgeno = max(counts)
            keys = list(genfreq.keys())
            for n,c in enumerate(counts):
                if c==maxgeno:
                    topgeno.append(keys[n])
            freq.append("|".join(topgeno))
            consensus.append(genfreq)
    pivot["consensus"] = consensus
    pivot["Top_Hit"] = freq
    pivot.to_csv(os.path.join(args.output_dir,str(x_int)+"_ukgenotyping_summary.csv"))
    return pivot
# Need a function to summarise the tables / files processed

def overall_summary(pivot):
    # Number of samples included
    # Summary of genotypes called
    print(f"Samples processed: ",pivot.shape[0])
    print("Genotypes called: ")
    print(Counter(pivot['Top_Hit']))
# Need a function to use historic genotype rules when multiple genotype calls are returned
### e.g. AIV07-B1|AIV07-B2 = AIV07-B2

# Run script

def main(args):
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    output_dir = os.path.abspath(args.output_dir)

    # Generate output folder
    testing_functions([args.testing, str('output_folder'), output_dir])

    # Set up logging
    logging_file_setup(output_dir, args.testing)

    if args.input_folder is not None:
        logging.info("Processing folder full of Genbank input files")
        input_folder = os.path.abspath(args.input_folder)
        summarytab = run_full_blasts(input_folder,"all_in_folder",args.extension)
        pivot_out = create_persample_summary(summarytab)
        overall_summary(pivot_out)
    if args.input_file is not None:
        logging.info("Processing single Genbank input files")
        input_file= os.path.abspath(args.input_file)
        summarytab = run_full_blasts(input_file,"single",args.extension)
        pivot_out = create_persample_summary(summarytab)
        overall_summary(pivot_out)

    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f'Pipeline time to completion: {start_time - end_time}')


if __name__ == '__main__':
    args = read_commandline()

    # Setting up testing and logging
    logging_file_setup(args.output_dir, args.testing)
    check = check_arguments(args)
    if check == 1:
        sys.exit(logging.error("Arguments provided were not expected. Please check log."))
    sys.exit(main(args))
