# Concatenated tree prep for mutation scan
#
# outputs: per cds alignment, table of any excluded sequences
# use mafft


from pathlib import Path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys
import logging
from datetime import datetime, date
import pandas as pd
import subprocess


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


def logging_file_setup(output_folder):
    """
    Log file setup in output folder

    :return: None
    """

    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    datetime_stamp = now.strftime("%Y%m%d-%H%M%S")
    logging_file_output = os.path.join(
        output_folder,
        str(str(now.strftime("%Y%m%d")) + str("_fileprep_logging_file.log")),
    )
    print("Logging output in the file: {}".format(logging_file_output))
    logging.basicConfig(
        filename=logging_file_output,  # filename=logging_file_output,
        filemode="w",
        format="%(asctime)s:%(levelname)s:%(message)s",
        level=logging.INFO,
    )
    print("Logging Level: INFO")
    logging.info("Logging Started.")
    logging.info(f"Version: {__version__}")


def check_arguments(args):
    """
    Check that paths provided exist and create output folder if it doesn't exist already.

    :param args: output from arg parse
    :return: status
    """
    if not os.path.exists(args.input_file):
        print(f"File does not appear to exists: {args.input_file}. Please check")
        return 1
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    return 0


def metatab_checks(metatab):
    """
    Check that the metadata file is as expected, with > 1 row, correct columns and na_sequence with data.

    :param args: metatab
    :return: status
    """
    if metatab.shape[0] == 0:
        logging.error(
            f"No rows found in the input meta-data csv file {args.input_file}. Please check that you have the correct file and try again!"
        )
        return 1
    metacols = [
        "na_sequence",
        "isolate_epi_id",
        "segment_name",
        "isolate_name",
        "isolate_id",
    ]
    for m in metacols:
        if m in metatab.columns:
            pass
        else:
            logging.error(
                f"Required columns are not present {m}. Please check that you have the correct file and try again!"
            )
            return 1
    if metatab["na_sequence"].isnull().sum():
        logging.error(
            f"No sequence information is present in na_sequence column. Please check that you have the correct file and try again!"
        )
        return 1
    return 0


path = os.path.dirname(sys.argv[0])

segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def create_fasta_headers(df):
    """
    Produce a FASTA header to accomanpy the nucleotide sequence from the database table
    Includes the isolate_epi_id, isolate_id and segment_name
    Excludes any sequences where the segment is not specified

    :param args: df
    :return: filtered_df
    """
    filtered_df = df[df["segment_name"].notnull()]
    filtered_df["fasta_header"] = (
        filtered_df["isolate_epi_id"] + str("|") + filtered_df['id_1'].astype(str)+str("|") + filtered_df["segment_name"]
    )
    logging.info("Fasta headers created and added to meta-data....")
    return filtered_df


def write_fasta_file(df, filetag, output_dir):
    """
    Create fasta files per segment using the fasta_header and
    na_sequence columns from the database extract dataframe


    :param args: df, segment, file name tag, output directory
    :return: output_file fasta file name
    """
    logging.info(f"Processing fasta file for input file....")
    sequence_df = df[["fasta_header", "na_sequence", "segment_name", "isolate_epi_id"]]
    # change method to remove warnings
    sequence_df.insert(
        len(sequence_df.columns), "length", len(sequence_df["na_sequence"])
    )
    sequence_df = sequence_df.sort_values(
        ["fasta_header", "length"], ascending=[False, False]
    )
    # look at duplicates
    uniqueids = set(sequence_df["fasta_header"])
    #   print(f'{len(uniqueids)} unique ids in FASTA headers')
    logging.info(f"{len(uniqueids)} unique ids in FASTA headers")
    sequence_df = sequence_df.reset_index()
    sequence_df = sequence_df.drop_duplicates(subset=["fasta_header"], keep="first")
    sequence_dictionary_list = sequence_df.to_dict("records")
    fasta_filename = f"{filetag}_{now}.fasta"
    output_file = os.path.join(output_dir, fasta_filename)
    seq_file = open(output_file, "w")
    logging.info(f"Creating fasta file for segment: {output_file}....")
    for sequence in sequence_dictionary_list:
        fasta_header = sequence["fasta_header"]
        fasta_sequence = sequence["na_sequence"]
        file_content = SeqRecord(Seq(fasta_sequence), id=fasta_header, description="")
        SeqIO.write(file_content, seq_file, "fasta")
    seq_file.close()
    logging.info(f"Fasta file for segment complete....")
    logging.info(f"Meta-data for segment contains {sequence_df.shape[0]} records.")
    return output_file


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
    logging_file_setup(output_dir)
    # main steps of script
    # 1. Input metadata file and create fasta file for each segment
    filetag = args.tagname
    metatab = pd.read_csv(args.input_file)
    metacheck = metatab_checks(metatab)
    if metacheck == 1:
        sys.exit(
            logging.error("Metadata file is not as not expected. Please check log.")
        )
    metatab["segment_name"] = metatab["segment_name"].fillna("NA")
    metatab_edited = create_fasta_headers(metatab)
    print(f"Generating full fasta files based on given arguments. Filetag {filetag}")
    logging.info(
        f"Generating full fasta files based on given arguments. Filetag {filetag}"
    )
    if metatab_edited.shape[0] == 0:
        print(f"No sequence data in input table!")
        sys.exit(logging.info(f"No sequence data in input table !"))
    else:
        fastapath = write_fasta_file(metatab_edited, filetag, output_dir)

    write_metadata_table(metatab_edited, filetag, output_dir)
    print(f"Fasta file created: {fastapath}")
    logging.info(f"Fasta file created: {fastapath}")
    print("Ready for genotyping")
    end_time = datetime.now()  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: {start_time - end_time}")


if __name__ == "__main__":
    args = read_commandline()
    # Setting up testing and logging
    logging_file_setup(args.output_dir)
    check = check_arguments(args)
    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )
    sys.exit(main(args))
