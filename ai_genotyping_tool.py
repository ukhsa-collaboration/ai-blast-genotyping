"""
### AI Genotyping command line tool

## Overview of the tool:
# For each sample we need to blast the sequence of the 8 segments against the reference blast database and get the top genotype hit for each segment
# BLAST tabular output format means that we can automate the results processing
# Wrangle the results of the blast search results together and make a call about the genotype
# Compare the genotypes to the results that we have from the existing sequence database 

"""

import os
import glob
import sys
import pandas as pd
import numpy as np
import logging
import argparse
from datetime import datetime, date
import shutil
import subprocess
from collections import Counter
from Bio import SeqIO

now = date.today()


__version__ = 2.0
__author__ = "Kate Howell"


def logging_file_setup(output_folder, testing):
    '''
    Return the sum of two decimal numbers in binary digits.

            Parameters:
                    a (int): A decimal integer
                    b (int): Another decimal integer

            Returns:
                    binary_sum (str): Binary string of the sum of a and b
    '''
    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    logging_file_output = os.path.join(
        output_folder, str(str(now.strftime("%Y%m%d")) + str("_geno_logging_file.log"))
    )
    print("Logging output in the file: {}".format(logging_file_output))

    if testing == "yes":
        logging.basicConfig(
            filename=logging_file_output,
            filemode="w",
            format="%(asctime)s:%(levelname)s:%(message)s",
            level=logging.DEBUG,
        )
        print("Logging Level: DEBUG")
    else:
        logging.basicConfig(
            filename=logging_file_output,
            filemode="w",
            format="%(asctime)s:%(levelname)s:%(message)s",
            level=logging.INFO,
        )
        print("Logging Level: INFO")
    logging.info("Logging Started.")


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"AI UK Genotyping command line tool")
    parser.add_argument(
        "--output_dir",
        "-o",
        required=True,
        default=os.getcwd(),
        help="Output folder. Default: CWD.",
    )
    parser.add_argument(
        "--testing",
        "-t",
        required=False,
        default=False,
        choices=["yes", "no"],
        help='Debugging mode. Specify by either "yes" or "no"',
    )  # Change this to a list of options
    parser.add_argument(
        "--tagname", "-n", required=True, help="Analysis name for output files"
    )
    parser.add_argument(
        "--extension",
        "-e",
        required=False,
        default=".fasta",
        help="FASTA file extension if not default",
    )
    input_group = (
        parser.add_mutually_exclusive_group()
    )  # input group requires one or the other
    input_group.add_argument(
        "--input_file", "-i", required=False, help="Input FASTA file"
    )
    input_group.add_argument(
        "--input_folder", "-f", required=False, help="Input FASTA file"
    )
    parser.add_argument(
        "--blastdb", "-b", required=True, help="Reference BLAST database"
    )
    parser.add_argument(
        "--strict",
        "-s",
        required=True,
        help="Strict version, all 8 segments for genotype",
        default="no",
    )
    parser.add_argument(
        "--identity",
        "-d",
        required=True,
        help="Percentage identity threshold. Default = 97",
        default=97,
    )
    args = parser.parse_args()
    # Need to handle output dir before setting up logging files.hat's the timeline for the analysis
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        logging.info(f"Creating output folder:\n {args.output_dir}")
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
            logging.error(
                f"File does not appear to exists: {args.input_folder}. Please check"
            )
            return 1

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    if not os.path.exists(args.blastdb + ".nin"):
        logging.error(
            f"BLAST database does not appear to exists: {args.blastdb}. Please check"
        )
        return 1

    return 0


def testing_functions(testing_functions_parameters):
    """
    Basic testing function

    :return: N/A
    """
    # Generate output folder
    if testing_functions_parameters[1] == str("output_folder"):
        if testing_functions_parameters[0] == True:
            if os.path.exists(testing_functions_parameters[2]):
                print(
                    "Deleting output_folder {}".format(testing_functions_parameters[2])
                )
                logging.info(
                    "Deleting output_folder {}".format(testing_functions_parameters[2])
                )
                try:  ## Try to remove tree; if failed show an error using try...except on screen
                    shutil.rmtree(testing_functions_parameters[2])
                except OSError as e:
                    print("Error: %s - %s." % (e.filename, e.strerror))
                    logging.info("Error: %s - %s." % (e.filename, e.strerror))
                os.mkdir(testing_functions_parameters[2])
            else:
                os.mkdir(testing_functions_parameters[2])
        elif testing_functions_parameters[0] == False:
            # Create output folder
            if not os.path.exists(testing_functions_parameters[2]):
                print(
                    "Output folder does not exist. Attempting to create folder: %s"
                    % (testing_functions_parameters[2])
                )
                logging.info(
                    "Output folder does not exist. Attempting to create folder: %s"
                    % (testing_functions_parameters[2])
                )
                os.mkdir(testing_functions_parameters[2])


path = os.path.dirname(sys.argv[0])
blast_cols = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]

segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

if os.path.exists(os.path.join(path, "blast_geno_threshold_table98.csv")):
    genoblast = pd.read_csv(
        os.path.join(path, "blast_geno_threshold_table98.csv"), dtype=str
    )
    genogroups = pd.melt(
        genoblast,
        id_vars=["sequence", "genotype", "subtype", "constellation"],
        value_vars=segments,
    )
    print(genogroups.head)
    genogroups = genogroups.rename(columns={"variable": "Segment", "value": "Group"})
    genogroups["Segment"] = genogroups["Segment"].replace(np.nan, "NA")
    groups = genogroups["Group"].tolist()
    genotypes = genogroups["genotype"].tolist()
    genodict = genogroups.groupby("Group")["genotype"].apply(list).to_dict()
    genodict2 = {}
    for key in genodict:
        genodict2[key] = "|".join(genodict[key])
    genodict = genodict2
else:
    logging.error(
        f'Reference genotypes table does not exists:"genotype_groups_examples.csv". Please check if the file is in the correct location'
    )


def tidy_fasta_files(sample):
    """
    Run regex on FASTA files to ensure BLAST searches can run

    :return: None
    """
    logging.info(
        "Tidying FASTA headers for the following non-standard characters.... () ,-\n"
    )
    os.system("perl -pi -e 's/\(/_/g' " + sample)
    os.system("perl -pi -e 's/\)/_/g' " + sample)
    os.system("perl -pi -e 's/ /_/g' " + sample)
    os.system("perl -pi -e 's/,/_/g' " + sample)
    os.system("perl -pi -e 's/--//g' " + sample)
    os.system("perl -pi -e 'spath/-\n-//g' " + sample)


def duplicate_fasta_check(fasta):
    """
    Read in FASTA file, check if sequence in dictionary, report any duplicates

    :return: duplicates_list
    """
    seqdict = {}
    duplist = []
    fasta_sequences = SeqIO.parse(open(fasta), "fasta")
    for seq in fasta_sequences:
        # check if sequence already in dictionary
        if seq.id in seqdict:
            print(f"Duplicate sequence found for {seq.id}")
            logging.info(f"Duplicate sequence found for {seq.id}")
            duplist.append(seq.id)
        else:
            seqdict[seq.id] = str(seq.seq)

   # print(
  #      f"{len(duplist)} duplicates found in FASTA files. Any duplicates should be reviewed prior to interpreting results!"
   # )
    logging.info(
        print(
            f"{len(duplist)} duplicates found in FASTA files. Any duplicates should be reviewed prior to interpreting results!"
        )
    )
    return duplist


def missing_fasta_check(fasta, segdict):
    """
    Read in FASTA file, check if segments are present for each sample, report any missing segments

    :return: segmissing
    """
    print("Check which FASTAs available...")
    delim = "None"
    fasta_sequences = SeqIO.parse(open(fasta), "fasta")

    for seq in fasta_sequences:
        if "|" in seq.id:
            delim = "|"
        elif "." in seq.id:
            delim = "."
        elif "_" in seq.id:
            delim = "."
        identifier = seq.id.split(delim)
        newlist = []
        try:
            sofar = segdict.get(identifier[0])
            for s in sofar:
                newlist.append(s)
        except:
            pass
        newlist.append(identifier[-1])
        segdict[identifier[0]] = newlist

    return segdict


def create_segment_tab(segdict):
    """
    Create a segment table that checks the FASTA file for presence of all 8 segments and writes to a csv table

    :return: segment dataframe,fasta_count
    """
    print("Creating segment tab....")
    segdf = {}
    samples = segdict.keys()
    fasta_count = 0
    for s in samples:
        segfound = segdict.get(s)
        notfound = set(segments) - set(segfound)
        if len(notfound) >= 1:
            print(f"Missing segments identified for {s}: {notfound}")
            logging.info(f"Missing segments identified for {s}: {notfound}")
        row = []
        for seg in segments:
            if seg in segfound:
                row.append(seg)
                fasta_count += 1
            elif seg in notfound:
                row.append("missing")

        segdf[s] = row
    return segdf, fasta_count


def run_blast(db, folder, sample, output_dir):
    """
    Run the BLAST query using sub-process

    :return: N/A
    """
    tidy_fasta_files(os.path.join(folder, sample))
    logging.info("Performing the BLAST searches per FASTA file")
    subprocess.call(
        f"blastn -db {db} -query {os.path.join(folder,sample)} -out {os.path.join(output_dir,sample)}.blast.out -outfmt 6 -max_target_seqs 1",
        shell=True,
    )


def match_genotype_dict(blast_pass):
    """
    Match the BLAST hit to a genotype group and possible genotype match per segment

    :return: BLAST hit dataframe
    """
    print("Matching genotype to dictionary")
    segdict = dict(
        zip(
            genogroups["Segment"] + "_" + genogroups.sequence.str.replace("_", ""),
            genogroups.Group,
        )
    )
    # map the query sequence back to the genotyping table
    # need to look up which group the reference is in
    blast_pass["group_match"] = (
        blast_pass["ref_match"].apply(lambda x: segdict.get(x)).fillna("")
    )
    blast_pass["genotype_match"] = (
        blast_pass["group_match"].apply(lambda x: genodict.get(x)).fillna("")
    )
    return blast_pass


def create_hit_dict(blast_pass):
    hitdict = {}
    for q in list(set(blast_pass['isolate_epi_id'])):
        subblast = blast_pass[blast_pass['isolate_epi_id'] == q]
     #   print(subblast.shape)
       # for s in segments:
       #     subblastseg = subblast[subblast['segment']==s]
        if subblast.shape[0] >= 1:
         #  print(Counter(list(subblast['hit_geno'])))
            hitdict[subblast["isolate_epi_id"].iloc[0]] = [dict(Counter(list(subblast['hit_isolate_name']))),Counter(list(subblast['hit_geno']))]
    return hitdict

def tidy_blast_table(folder, sample, segmissing, fasta_count):
    """
    Wrange the BLAST results table, exclude results that do not meet threshold, separate no sequence vs. no hit, split the query and hit ids into isolate_epi_id and segment.

    :return: newresultstab dataframe
    """
    print("Tidying BLAST tables....break")
    logging.info("Reading in BLAST output:")
    logging.info(os.path.join(args.output_dir, f"{sample}.blast.out"))
    blasttab = pd.read_csv(
        os.path.join(args.output_dir, f"{sample}.blast.out"), sep="\t", header=None
    )
    if blasttab.shape[0] < fasta_count:
        print(
            f"{fasta_count - blasttab.shape[0]} FASTAs do not meet minimum BLAST thresholds"
        )
        logging.info(
            f"{fasta_count - blasttab.shape[0]} FASTAs do not meet minimum BLAST thresholds"
        )
    blasttab.columns = blast_cols
    threshold = int(args.identity)
    blast_pass = blasttab[blasttab["pident"] >= threshold]
    blast_pass = blast_pass.sort_values(['qseqid','pident'],ascending=[True,False])
    blast_fail = blasttab[blasttab["pident"] < threshold]
    if blast_fail.shape[0] >= 1:
        logging.info(
            f"{blast_fail.shape[0]} sequences do not meet a minimum percentage sequence identity to the ref seq!!"
        )
        print(
            f"{blast_fail.shape[0]} sequences do not meet a minimum percentage sequence identity to the ref seq!!"
        )
    logging.info("BLAST pass table dimensions:")
    logging.info(blast_pass.shape)
    if blast_pass.shape[0] == 0:
        print(
            "No sequences found to meet the minimum threshold. Tool now exiting, please check that you are using the same subtype of sequences as in the database!"
        )
        sys.exit(
            logging.error(
                "No sequences found to meet the minimum threshold. Tool now exiting, please check that you are using the same subtype of sequences as in the database!"
            )
        )
  #  print("Checking sequence identifier for delimiters.... ")
    logging.info("Checking sequence identifier for delimiters.... ")
    if "|" in list(blast_pass["qseqid"])[0]:
        delim = "|"

    elif "." in list(blast_pass["qseqid"])[0]:
        delim = "."
    elif "_" in list(blast_pass["qseqid"])[0]:
        delim = "_"
    else:
        logging.error(
            f'Sample header delimiter unknown! Please check, "." or "|" expected!'
        )
        print(f'Sample header delimiter unknown! Please check, "." or "|" expected!')
 #   print(f"Delimiter found: {delim}")
    logging.info(f"Delimiter found: {delim}")
    logging.info("Preparing BLAST table for summary")
 #   print("Preparing BLAST table for summary")
    try:
        test = blast_pass["qseqid"].str.split(pat=delim, expand=True)
    except:
        print(f"Issue with delimiter! {delim}")
    blast_pass.insert(len(blast_pass.columns), "segment", test[test.columns[-1]])
    blast_pass.insert(len(blast_pass.columns), "isolate_epi_id", test[test.columns[0]])
    blast_pass.insert(
        len(blast_pass.columns),
        "ref_match",
        blast_pass["segment"]
        + "_"
        + blast_pass["sseqid"].map(lambda x: x.split("|")[0]).str.replace("_", ""),
    )
    blast_pass[['hit_isolate_name', 'hit_geno','hit_subtype','hit_segment']] = blast_pass['sseqid'].str.split('|', expand=True)
    blast_pass = match_genotype_dict(blast_pass)
    blast_pass.to_csv(os.path.join(folder, f"{sample}.blast.out2"))
   
    results_df = pd.DataFrame()
    results_df.insert(
        len(results_df.columns), "genotype_match", blast_pass["genotype_match"]
    )
    results_df.insert(len(results_df.columns), "segment", blast_pass["segment"])
    results_df.insert(
        len(results_df.columns), "isolate_epi_id", blast_pass["isolate_epi_id"]
    )
    results_df.insert(len(results_df.columns), "top_hit", blast_pass["ref_match"])
    newresults_df = [results_df]
     
    hitdict = create_hit_dict(blast_pass)
    print("hit dictionary created....")
    if blast_fail.shape[0] >= 1:
        test = blast_fail["qseqid"].str.split(pat=delim, expand=True)
        blast_fail.insert(len(blast_fail.columns), "segment", test[test.columns[-1]])
        blast_fail.insert(
            len(blast_fail.columns), "isolate_epi_id", test[test.columns[0]]
        )
        missing = blast_fail[["segment", "isolate_epi_id"]]
        missing.insert(
            len(missing.columns),
            "genotype_match",
            ["No match via BLAST"] * missing.shape[0],
        )
        newresults_df.append(missing)
    newresultstab = pd.concat(newresults_df)
    samples = list(set(newresultstab["isolate_epi_id"]))

    for s in samples:
        subresults = newresultstab[newresultstab["isolate_epi_id"] == s]
        if subresults.shape[0] == 1:
            pass
        elif subresults.shape[0] == 0:
            submissing = segmissing[segmissing["sample"] == s]
            segfile = list(set(newresultstab["segment"]))[0]
            if submissing[segfile] == "missing":
                newresultstab.loc[len(newresultstab)] = ["No sequence", segfile, s]
            else:
                newresultstab.loc[len(newresultstab)] = ["BLAST FAIL", segfile, s]
        elif subresults.shape[0] < 8:
            submissing = segmissing[segmissing["sample"] == s]

            segfound = list(set(subresults["segment"]))
            segcheck = list(set(segments) - set(segfound))
            for seg in segcheck:

                if submissing[seg].iloc[0] == "missing":
                    newresultstab.loc[len(newresultstab)] = ["No sequence", seg, s,"Null"]
                else:
                    newresultstab.loc[len(newresultstab)] = ["BLAST FAIL", seg, s,"Null"]

    return newresultstab,hitdict


def run_full_blasts(folder, mode, extension, output_dir):
    """

    Perform several checks on the FASTA file(s) and run the BLAST queries. For the multiple file option, concatenate the results. Match the BLAST hits to the genotyping reference table.

    :return: BLAST results dataframe
    """
    print("Running BLASTs")
    results_tabs = []
    logging.info("Input folder provided:")
    logging.info(os.path.join(folder, "*" + extension))
    fastas = glob.glob(os.path.join(folder, "*" + extension))

    segdict = {}
    if mode == "all_in_folder":
        logging.info("Fasta files found:")
        logging.info(fastas)
        segtabs = []
        logging.info("Running BLAST")
        print("Running BLAST on input folder")
        for f in fastas:

            
            segdict = missing_fasta_check(f, segdict)
            duplist = duplicate_fasta_check(f)
            #if len(duplist) >= 1:
           #     print(f"Duplicates identified: {duplist}")
           #     logging.info(f"Duplicates identified: {duplist}")
            run_blast(args.blastdb, folder, f, output_dir)
           # logging.info("Reviewing BLAST results")
          #  print("Reviewing BLAST results")
            segmissing, fasta_count = create_segment_tab(segdict)

            segdicttab = pd.DataFrame.from_dict(segmissing, orient="index")
            segdicttab = segdicttab.reset_index()
            segtabs.append(segdicttab)
            segcols = ["sample"]
            for s in segments:
                segcols.append(s)
            segdicttab.columns = segcols

            sresults,hitdict = tidy_blast_table(args.output_dir, f, segdicttab, fasta_count)
           # logging.info("Combining results across FASTA files...")
           # print("Combining results across FASTA files...")
            results_tabs.append(sresults)
        newdf = pd.concat(results_tabs)
        segtab = pd.concat(segtabs)
        segtab.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_segment_table.csv")
        )
        logging.info(
            f"Segment table written to : {os.path.join(output_dir,f'{now}_{args.tagname}_segment_table.csv')}"
        )
        # add step to sort out missing segments table

        newdf.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        logging.info("Output file written to:")
        logging.info(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )

        return newdf

    elif mode == "single":
        logging.info("Running BLAST on single FASTA file")
        print("Running BLAST on single FASTA file")
        head_tail = os.path.split(folder)
        logging.info("Fasta file found:")
        logging.info(folder)
        logging.info(head_tail)
        duplist = duplicate_fasta_check(folder)
        segdict = missing_fasta_check(folder, segdict)
        segmissing, fasta_count = create_segment_tab(segdict)
        segdicttab = pd.DataFrame.from_dict(segmissing, orient="index")
        segdicttab = segdicttab.reset_index()
        segdicttab.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_segment_table.csv")
        )
        logging.info(
            f"Segment table written to : {os.path.join(output_dir,f'{now}_{args.tagname}_segment_table.csv')}"
        )
        segdicttab.columns = [
            "sample",
            "PB2",
            "PB1",
            "PA",
            "HA",
            "NP",
            "NA",
            "MP",
            "NS",
        ]
        if len(duplist) >= 1:
            print(f"Duplicates identified: {duplist}")
            logging.info(f"Duplicates identified: {duplist}")
        run_blast(args.blastdb, head_tail[0], head_tail[1], output_dir)
        logging.info("Reviewing BLAST results")
        print("Reviewing BLAST results")
        sresults,hitdict = tidy_blast_table(output_dir, head_tail[1], segdicttab, fasta_count)

        sresults.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        logging.info("Output file written to:")
        logging.info(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        print("Output file written to:")
        print(os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv"))
        return sresults, hitdict


#def details_lookup(q, output_dir):
#  #  print("look up details")
 #   blastout = glob.glob(os.path.join(output_dir, "*.blast.out2"))
  #  if len(blastout)==1:
   #     blasttab = pd.read_csv(os.path.join(output_dir, blastoutdetails_lookup[0]))
   # else:
    #    blasts = []
    #    for b in blastout:
    #        blasts.append(pd.read_csv(b))
    #    blasttab = pd.concat(blasts)
    #subblast = blasttab[blasttab['isolate_epi_id'] == q]
   # match = []
   # geno = []
   # for index, row in subblast.iterrows():
   #     hit = row['sseqid'].split("|")
   #     match.append(hit[0])
   #     geno.append(hit[1])
   # print(dict(Counter(match)),dict(Counter(geno)))
    #return dict(Counter(list(subblast['hit_isolate_name']))),Counter(list(subblast['hit_geno']))


def tidydict(pivot, collist):
    print("tidy up step")
    tidyup = ['{','}',"'"]
    for c in collist:
        for t in tidyup:
            pivot[c] = pivot[c].replace(t,"")
    return pivot

#def geno_prevalence():
#    


def create_persample_summary(summarytab, segthreshold,hitdict):
    """

    Wrangle the results from the BLAST table (matched to genotypes) into a summary and likely top result.
    :return: results pivot table


    """
    print("Per sample summary step...")
    pivot = pd.pivot_table(
        summarytab,
        values="genotype_match",
        index="isolate_epi_id",
        columns="segment",
        fill_value="No sequence",
        aggfunc="first",
    )
   # path = os.path.dirname(sys.argv[0])
    consensus = []
   # prioritytab = pd.read_csv(os.path.join(path, "genotype_prioritisation.csv"))
    freq = []
    results = []
    details = []
    genohits = []
    genoblast = []
   # print(hitdict)
  #  print(len(set(hitdict.keys() - set(pivot.index))))
    
    for index, row in pivot.iterrows():
        result = ""
        newlist = []
        topgeno = []
        
        newlist = split_geno_matches(pivot, row, newlist)
        genfreq = Counter(newlist)
        counts = list(genfreq.values())
        maxgeno = max(counts)
        keys = list(genfreq.keys())
        tophits = count_geno_hits(
           segthreshold, topgeno, counts, maxgeno, keys
        )
      #  print(result)
        freq.append(tophits)
        consensus.append(genfreq)
        results.append(result)
     #   print(hitdict)
      #  print(index)
        # try:
        keylist = hitdict.keys()
        if index in keylist:
            topmatch, genodict = hitdict[index]
            gencounts = list(genodict.values())
        # except:
        #     pass
        #   print(index,max(gencounts))
            if len(gencounts)>=1:
                    maxgenohit = max(gencounts)
                    genkeys = list(genodict.keys())
                #  topgeno = []
                    if maxgenohit >=6:
                        for n, c in enumerate(gencounts):
                            if c == maxgenohit:
                                genoblast.append(genkeys[n])
                    else:
                        genoblast.append("No genotype had >=6 segments, review top hit")
            else:
                genoblast.append("")
            details.append(topmatch)
        #    genohits.append(tophits)
        else:
            print(index, " not in blast dictionary")   
            details.append("No top hits found, BLAST fail??")  
            genoblast.append("")
    print(genoblast)
   # print(len(details),len(genohits))
    pivot["consensus"] = consensus
    pivot["Top_Hit"] = freq
    pivot['Top_hit_details'] = details
    
  #  pivot['Top_hit_genotype'] = dict(genodict)
    #pivot = tidydict(pivot,['Top_hit_details','Top_hit_genotype'])
    pivot["Best_result"] = genoblast
    pivot["Best_result"] = pivot["Best_result"].str.replace(
        "No sequence", "Insufficient sequence data"
    )

    pivot.to_csv(
        os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv")
    )
    logging.info("Summary output file written to:")
    logging.info(
        os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv")
    )
    print("Summary output file written to:")
    print(os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv"))
    return pivot


def count_geno_hits(segthreshold, topgeno, counts, maxgeno, keys):
    """

    Review genotype matches across the 8 segments, return the top results and where necessary choose a top hit out of multiple matches
    :return: top result per sample

    #try to simplify method....
    """
    for n, c in enumerate(counts):
        if c == maxgeno:
            topgeno.append(keys[n])
            if int(maxgeno) >= segthreshold:
                tophit = "|".join(topgeno)
                #if "|" in tophit:
                   # freqs = []
                    #for t in topgeno:
                    #    hittab = prioritytab[prioritytab["Genotype"] == t]
                    #    if hittab.shape[0] > 0:
                    ##        freqs.append(hittab["Frequency"].iloc[0])
                    #    else:
                    #        logging.info(f"Genotype not found in frequency table: {t}")
                    #      #  print(f"Genotype not found in frequency table: {t}")
                    #        freqs.append(0)

                   # maxhit = max(freqs)
                   # result = topgeno[freqs.index(maxhit)]

               # else:
                #    result = keys[n]
            else:
                tophit = "No known genotype: Please review individual segments results"
    return tophit


def split_geno_matches(pivot, row, newlist):
    """

    Split the genotype matches into a list
    :return: update newlist list


    """
    for s in pivot.columns:
        # try:
           # print(row[s])
            test = row[s].split("|")
           # print(test)
            for t in test:
                newlist.append(t)
        #except:
        #    newlist.append("")
    return newlist


# Need a function to summarise the tables / files processed


def overall_summary(pivot):
    """
    Produce a summary of the per sample results to the command line and log

    :return: N/A
    """
    # Number of samples included
    # Summary of genotypes called
    print(f"Samples processed: {pivot.shape[0]}")
    print("Genotypes called: ")
    print(Counter(pivot["Best_result"]))
    logging.info(f"Samples processed: {pivot.shape[0]}")
    logging.info("Genotypes called: ")
    logging.info(Counter(pivot["Best_result"]))


# Need a function to use historic genotype rules when multiple genotype calls are returned
### e.g. AIV07-B1|AIV07-B2 = AIV07-B2


# Run script


def main(args):
    """
    Main running of the script to run the BLAST query and wrangle the results to provide a per segment and sample summary of the genotyping results.

    :return: N/A
    """
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    global output_dir
    output_dir = os.path.abspath(args.output_dir)
    
    if args.strict == "yes":
        segthreshold = 8
    else:
        segthreshold = 7
    # Generate output folder
    testing_functions([args.testing, str("output_folder"), output_dir])
    logging.info(args)

    # Set up logging
    logging_file_setup(output_dir, args.testing)
    # Read in required files
    if args.input_folder is not None:
        logging.info("Processing folder full of input files")
        print("Processing folder full of input files")
        input_folder = os.path.abspath(args.input_folder)

        summarytab,hitdict = run_full_blasts(
            input_folder, "all_in_folder", args.extension, output_dir
        )
        pivot_out = create_persample_summary(summarytab, segthreshold,hitdict)
        overall_summary(pivot_out)
    if args.input_file is not None:
        print("Processing single FASTA files")
        logging.info("Processing single FASTA files")
        input_file = os.path.abspath(args.input_file)
        summarytab,hitdict = run_full_blasts(input_file, "single", args.extension, output_dir)
        pivot_out = create_persample_summary(summarytab, segthreshold,hitdict)
        overall_summary(pivot_out)

    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: {start_time - end_time}")


if __name__ == "__main__":
    args = read_commandline()

    # Setting up testing and logging
    logging_file_setup(args.output_dir, args.testing)
    check = check_arguments(args)
    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )
    sys.exit(main(args))
