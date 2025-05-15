import os
import sys
import pandas as pd
import logging
from sqlalchemy import inspect, MetaData, func,select
import argparse

from datetime import date, datetime
from dateutil.relativedelta import relativedelta
import src.utilities as util 


now = date.today()
x_int = int(now.strftime("%Y%m%d%H%M%S"))

__version__ = 0.1
__author__ = 'Kate Howell'




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


def main(args):
    start_time = datetime.now()  # Start time for calculating performance improvements

    
    global output_dir
    output_dir = args.output_dir
    logging.info(f"Ouput directory: {output_dir}")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # Set up logging
    util.logging_file_setup(output_dir)

    db_conn = util.get_aiseqdb_conn()

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
    util.logging_file_setup(args.output_dir)
    check = util.check_arguments(args)
    if check == 1:
        sys.exit(logging.error("Arguments provided were not expected. Please check log."))
    
    sys.exit(main(args))