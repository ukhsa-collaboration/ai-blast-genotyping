"""
### AI Genotyping command line tool

## Overview of the tool:
# For each sample we need to blast the sequence of the 8 segments against the reference blast database and get the top genotype hit for each segment
# BLAST tabular output format means that we can automate the results processing
# Wrangle the results of the blast search results together and make a call about the genotype
# Compare the genotypes to the results that we have from the existing sequence database 

"""

import os
import sys

import logging
import argparse
from datetime import datetime, date

import src.utilities as util
import src.blast_wrangle as bl


__version__ = 2.0
__author__ = "Kate Howell"


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
        "--input_folder", "-f", required=False, help="Input FASTA folder"
    )
    parser.add_argument(
        "--blastdb",
        "-b",
        required=True,
        help="Reference BLAST database",
        default="reference_files/all_genotype_references.db",
    )
    parser.add_argument(
        "--constellation",
        "-c",
        required=True,
        help="Constellation table",
        default="reference_files/blast_geno_threshold_table98.csv",
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
        help="Percentage identity threshold. Default = 98",
        default=98,
    )
    args = parser.parse_args()
    # Need to handle output dir before setting up logging files.hat's the timeline for the analysis
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        logging.info(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)
    return args


# Run script


def main(args):
    """
    Main running of the script to run the BLAST query and wrangle the results to provide a per segment and sample summary of the genotyping results.

    :return: N/A
    """
    start_time = datetime.now()  # Start time for calculating performance improvements

    ### global declarations
    global output_dir
    global repopath
    global now
    global constellation_path

    now = date.today()
    output_dir = os.path.abspath(args.output_dir)

    if args.blastdb == "reference_files/all_genotype_references.db":
        repopath = os.path.dirname(sys.argv[0])
    else:
        repopath = ""

    if args.strict == "yes":
        segthreshold = 8
    else:
        segthreshold = 7
    # Generate output folder
    util.testing_functions([args.testing, str("output_folder"), output_dir])
    logging.info(args)

    # Set up logging
    util.logging_file_setup(output_dir, args.testing)
    # Read in required files
    if args.input_folder is not None:
        logging.info("Processing folder full of input files")
        print("Processing folder full of input files")
        input_folder = os.path.abspath(args.input_folder)

        summarytab, hitdict = bl.run_full_blasts(
            repopath,
            args,
            now,
            input_folder,
            "all_in_folder",
            args.extension,
            output_dir,
        )
        pivot_out = bl.create_persample_summary(
            args, now, summarytab, segthreshold, hitdict
        )
        bl.overall_summary(pivot_out)
    if args.input_file is not None:
        print("Processing single FASTA files")
        logging.info("Processing single FASTA files")
        input_file = os.path.abspath(args.input_file)
        summarytab, hitdict = bl.run_full_blasts(
            repopath, args, now, input_file, "single", args.extension, output_dir
        )
        pivot_out = bl.create_persample_summary(
            args, now, summarytab, segthreshold, hitdict
        )
        bl.overall_summary(pivot_out)

    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: {start_time - end_time}")


if __name__ == "__main__":
    args = read_commandline()

    # Setting up testing and logging
    util.logging_file_setup(args.output_dir, args.testing)

    check = util.check_arguments(args)

    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )
    sys.exit(main(args))
