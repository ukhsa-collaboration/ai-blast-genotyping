"""Module providing a range of utilities functions for the mutation scan wrapper."""

import os
import sys
import logging
from datetime import date
import shutil
import sqlalchemy
from sqlalchemy import text, inspect, create_engine, inspect, URL
import pandas as pd



segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


now = date.today()


def logging_file_setup(output_folder, testing):
    """
    Return the sum of two decimal numbers in binary digits.

            Parameters:
                    a (int): A decimal integer
                    b (int): Another decimal integer

            Returns:
                    binary_sum (str): Binary string of the sum of a and b
    """
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
    if args.blastdb == "reference_files/all_genotype_references.db":
        repopath = os.path.dirname(sys.argv[0])
    else:
        repopath = ""
    #   print(os.path.join(repopath,args.blastdb + ".nin"))
    if not os.path.exists(os.path.join(repopath, args.blastdb + ".nin")):
        logging.error(
            f"BLAST database does not appear to exist: {os.path.join(repopath,args.blastdb)}. Please check"
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
        if testing_functions_parameters[0] is True:
            if os.path.exists(testing_functions_parameters[2]):
                print(
                    "Deleting output_folder {}".format(testing_functions_parameters[2])
                )
                logging.info(
                    "Deleting output_folder {}".format(testing_functions_parameters[2])
                )
                try:  
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


def tidydict(pivot, collist):
    print("tidy up step")
    tidyup = ["{", "}", "'"]
    for c in collist:
        for t in tidyup:
            pivot[c] = pivot[c].replace(t, "")
    return pivot



def query_aiseqdb(username, epiids):
    """
    Run query of aiseqdb to retrieve records for specific isolate_epi_ids.

    :param username: user name for the database (e.g. kate.howell)
    :param epiids: list of isolate_epi_ids from input table


    :return: database extract from aiseqdb
    """
    print("Preparing to query aiseqdb....")
    SQLquery = f"""select *
    from isolates  
    left join segment_sequences on (segment_sequences.isolate_id = isolates.id) 
    left join files as fasta_files 
    on (segment_sequences.file_id = fasta_files.id) where isolates.isolate_epi_id in ('{epiids}')"""
    
    hostname = "SqlpgEpiDevCol01.unix.phe.gov.uk"
    database_name = "ai_seq_db_preprod"
    username = f'{username}@phe.gov.uk'
    db_conn = create_engine(f'postgresql+psycopg2://{username}@{hostname}/{database_name}')
    insp = inspect(db_conn)

    print(insp.get_table_names())
    
    db_extract = pd.read_sql(
            sql=text(SQLquery), 
            con=db_conn.connect())
    print(f'{db_extract.shape[0]} segments retrieved from aiseqdb')

    return db_extract

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def get_aiseqdb_conn():
    return sqlalchemy.create_engine(
        URL.create(
            "postgresql+psycopg2",
            username=os.environ["AISEQDB_PROD_USER"],
            password=os.environ["AISEQDB_PROD_PASSWORD"],
            host=os.environ["AISEQDB_PROD_HOST"],
            database=os.environ["AISEQDB_PROD_DB"]
        )
    )


def tidy_files(output_dir,subtype):
    # remove fasta files in alns complete
    # create folder of alns and meta and zip up 

    shutil.make_archive(f'genotyping_run_{now}', format='zip', root_dir=output_dir)



def metatab_checks(input,metatab):
    """
    Check that the metadata file is as expected, with > 1 row, correct columns and na_sequence with data.

    :param args: metatab
    :return: status
    """
    if metatab.shape[0] == 0:
        logging.error(
            f"No rows found in the input meta-data csv file {input}. Please check that you have the correct file and try again!"
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
