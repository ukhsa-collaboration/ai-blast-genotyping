# Concatenated tree prep for mutation scan
#
# outputs: per cds alignment, table of any excluded sequences
# use mafft


import argparse
import os
import sys
import logging
from datetime import datetime, date
import pandas as pd
import src.utilities as util
import src.fasta_processing as fap

now = date.today()

__version__ = 1.1
__author__ = "Kate Howell"


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(
        description=f"Script to create AI segment alignments for mutation scan"
    )
    parser.add_argument(
        "--input_file", "-i", required=True, help="CSV file retrieved from AI database."
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        required=True,
        default=os.getcwd(),
        help="Output folder. Default: CWD.",
    )
    parser.add_argument(
        "--tagname", "-n", required=True, help="Analysis name for output files"
    )

    # Add threads
    args = parser.parse_args()

    # Need to handle output dir before setting up logging files.
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)

    return args





path = os.path.dirname(sys.argv[0])




def write_metadata_table(df, filetag, output_dir):
    """
    Create a reduced metadata dataframe that removes the na_sequence columns
    Write the dataframe to a new file

    :param args: df, file name tag, output directory
    :return: None
    """
    df_edited = df.drop(columns=["na_sequence"])
    df_output_file = os.path.join(
        output_dir, str(f"{filetag}_{now}_edited_metadata.csv")
    )
    df_edited.to_csv(df_output_file, index=False)
    logging.info(f"Meta-data for segment file created {df_output_file}.")


path = os.path.dirname(sys.argv[0])


def main(args):
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    output_dir = os.path.abspath(args.output_dir)
    print(output_dir)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # Set up logging
    util.logging_file_setup(output_dir)
    # main steps of script
    # 1. Input metadata file and create fasta file for each segment
    filetag = args.tagname
    metatab = pd.read_csv(args.input_file)
    metacheck = util.metatab_checks(args.input_file,metatab)
    if metacheck == 1:
        sys.exit(
            logging.error("Metadata file is not as not expected. Please check log.")
        )
    metatab["segment_name"] = metatab["segment_name"].fillna("NA")
    metatab_edited = fap.create_fasta_headers(metatab)
    print(f"Generating full fasta files based on given arguments. Filetag {filetag}")
    logging.info(
        f"Generating full fasta files based on given arguments. Filetag {filetag}"
    )
    if metatab_edited.shape[0] == 0:
        print(f"No sequence data in input table!")
        sys.exit(logging.info(f"No sequence data in input table !"))
    else:
        fastapath = fap.write_fasta_file(now,metatab_edited, filetag, output_dir)

    write_metadata_table(metatab_edited, filetag, output_dir)
    print(f"Fasta file created: {fastapath}")
    logging.info(f"Fasta file created: {fastapath}")
    print("Ready for genotyping")
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
