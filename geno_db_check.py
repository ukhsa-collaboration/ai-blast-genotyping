import os
import glob
import sys
import pandas as pd
import subprocess
from pathlib import Path
import uuid
import shutil
import logging
import sqlalchemy
from sqlalchemy import (
    select,
    MetaData,
    inspect,
    and_,
    URL,
    func
)

import argparse

from datetime import date, datetime
from dateutil.relativedelta import relativedelta

now = date.today()
x_int = int(now.strftime("%Y%m%d%H%M%S"))

__version__ = 0.1
__author__ = 'Kate Howell'

def logging_file_setup(output_folder):
    """
    Log file setup in output folder

    :return: None
    """
    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    datetime_stamp = now.strftime("%Y%m%d-%H%M%S")
    logging_file_output = os.path.join(output_folder, str(str(now.strftime("%Y%m%d")) + str('_geno_id_check.log')))
    print("Logging output in the file: {}".format(logging_file_output))
    logging.basicConfig(
            filename=logging_file_output,  # filename=logging_file_output,
            filemode='w',
            format='%(asctime)s:%(levelname)s:%(message)s',
            level=logging.INFO)
    print("Logging Level: INFO")
    logging.info("Logging Started.")
    logging.info(f"Version: {__version__}")


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


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"AI DB PDAG / RVU Query Extract")
    parser.add_argument('--username', '-u', required=True, help='PHE email address')
    parser.add_argument('--password_file', '-p', required=False, help='Password file path if not running with kerberos authentication')
    parser.add_argument('--no-rewrite-username', required=False, action='store_true', help='Skips rewriting usernames to include @phe.gov.uk')
    parser.add_argument('--use-aiseqdb-environment-variables', required=False, action='store_true', help='Read aiseqdb credentials from bash environment variables. Overrides username and password flags')
    parser.add_argument('--output_dir','-o', required=True,default=os.getcwd(),help='Output directory path')
        
    args = parser.parse_args()
    print(args)
    logging.info(f'Arguments provided: {args}')
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
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    return 0



def main(args):
    start_time = datetime.now()  # Start time for calculating performance improvements

    
    global output_dir
    output_dir = args.output_dir
    logging.info(f"Ouput directory: {output_dir}")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # Set up logging
    logging_file_setup(output_dir)

    db_conn = get_aiseqdb_conn()

    #Set up aiseqdb connection and tables
    insp = inspect(db_conn)
    global aiseqdb_metadata
    aiseqdb_metadata = MetaData()
    aiseqdb_metadata.reflect(bind=db_conn)
    global isolates_table
    global segment_sequences_table
    global files_table
    global genotype

    isolates_table = aiseqdb_metadata.tables["isolates"]
    segment_sequences_table = aiseqdb_metadata.tables["segment_sequences"]
    files_table = aiseqdb_metadata.tables["files"]
    genotype = aiseqdb_metadata.tables["genotype"]

    # query genotyping table to see what data already exists
    lastingest = func.max(genotype.c.ingest_date)
    
    datecheck = pd.read_sql(
            sql=lastingest, 
            con=db_conn.connect()
        )
    querydate = datecheck.iloc[0][0] - relativedelta(weeks=int(2))
    existing_id = select(func.distinct(genotype.c.isolate_epi_id))
    processed = datecheck = pd.read_sql(
        sql=existing_id, 
        con=db_conn.connect()
    )
    
    # query the main tables for new data
    dbquery = select(
        isolates_table.join(
            segment_sequences_table,
            isolates_table.c.id == segment_sequences_table.c.isolate_id,
            isouter=True,
        ).join(
            files_table,
            segment_sequences_table.c.file_id == files_table.c.id,
            isouter=True,
        ))
    dbquery = dbquery.where(files_table.c.processed_date >=(querydate))
    dbextract = pd.read_sql(
            sql=dbquery, 
            con=db_conn.connect()
        )

    #check which isolate_epi_ids still need genotyping
    need_genotyping = set(dbextract.isolate_epi_id) - set(processed.distinct_1)


    subdbextract = dbextract[dbextract['isolate_epi_id'].isin(need_genotyping)]
    subdbextract.to_csv(os.path.join(output_dir,f"{now}_db_extract_newseqs_forgenotyping.csv"))

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