# run pdag export once a week

# steps required
# copy genotyping linelist.csv from the HPC
# run main script
# delete temp data if zip file with today's date exists
# copy data to the HPC 
# remove zip folders 

import logging
import os
import sys
import glob
from datetime import date, datetime
import subprocess

import shutil


from sqlalchemy import (
    select,
    MetaData,
    inspect,
    and_,
    URL,
    func
)


from datetime import date, datetime
from dateutil.relativedelta import relativedelta
now = date.today()
x_int = int(now.strftime("%Y%m%d%H%M%S"))

__version__ = 0.1
__author__ = 'Kate Howell'


 
def main():
    start_time = datetime.now()  # Start time for calculating performance improvements
    #for f in glob.glob("deletable_temp_folder*"):
    #    rmtree(f)
    output_dir = os.path.join("/data/genotyping",str(now))
    dbpath = "/data/PDAG_exports/aiseqdb_genotypes_latest.db"
    #Path(f"deletable_temp_folder_{uuid.uuid4()}").resolve()
    if not os.path.isdir(output_dir):  # Set up output folder
        print(f"Creating output folder:\n {output_dir}")
        os.mkdir(output_dir)
    logging.info(f"Ouput directory: {output_dir}")
    p = subprocess.Popen(["scp", "-i", "~/.ssh/id_rsa",f"kate.howell@158.119.147.128://phengs/hpc_projects/nicc80_ai/genoflugelbinder_runs/aiseqdb_genotypes_latest.db", "/data/PDAG_exports/"])
    sts = os.waitpid(p.pid, 0)
    logging.info("Genotyping db copied from HPC")

    p = subprocess.Popen(["sqlite3", "-header", "-csv", dbpath, "'select * from genotypes;'", ">", "/data/PDAG_exports/linelist.csv"])
    sts = os.waitpid(p.pid, 0)
    logging.info("Genotyping db queried for line list")
    print("Genotyping db queried for line list")

   # print("Genotyping db copied from HPC")
    p = subprocess.Popen(["conda","run","-n","blastgeno","python","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/geno_db_check.py","-u","kate.howell","--no-rewrite-username","--use-aiseqdb-environment-variables","-o",output_dir])
    sts = os.waitpid(p.pid, 0)
    logging.info("BLAST genotyping run list created")


    p = subprocess.Popen(["conda","run","-n","blastgeno","python","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/prep_for_genotyping.py","-i",os.path.join(output_dir,f"{now}_db_extract_newseqs_forgenotyping.csv"),"-o",output_dir,"-n","GenoRun"])
    sts = os.waitpid(p.pid, 0)
    logging.info("FASTA file created")

    p = subprocess.Popen(["conda","run","-n","blastgeno","python","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/ai_genotyping_tool.py","-o",output_dir,"-n","GenoRun","-i",os.path.join(output_dir,f"GenoRun_{now}.fasta"),"-b","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/reference_files/all_genotype_references.db","-c","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/reference_files/blast_geno_threshold_table98.csv","-s","yes","-d","98"])
    sts = os.waitpid(p.pid, 0)
    logging.info("Genotyping run complete")

    p = subprocess.Popen(["conda","run","-n","blastgeno","python","/home/phe.gov.uk/kate.howell/ai-blast-genotyping/genotype_ingest_fileprep.py","-o",output_dir,"-f",os.path.join(output_dir,f"GenoRun_{now}.fasta"),"-g","/home/phe.gov.uk/kate.howell/data/PDAG_exports/linelist.csv","-b",os.path.join(output_dir,f"{now}_GenoRun_BLAST_summary.csv"),"-o",output_dir])
    sts = os.waitpid(p.pid, 0)
    logging.info("Genotyping run complete")

    shutil.make_archive(f'genotyping_run_{now}', format='zip', root_dir=output_dir)
    # copy over folders to the HPC
    p = subprocess.Popen(["scp", "-i", "~/.ssh/id_rsa",f"genotyping_run_{now}.zip","kate.howell@158.119.147.128://phengs/hpc_projects/nicc80_ai/geno_blastdb/analysis_results" ])
    sts = os.waitpid(p.pid, 0)
    for f in glob.glob("*zip"):
        os.remove(f)
    end_time = datetime.now()  # end time for calculating performance improvements
    logging.info(f'Pipeline time to completion: {start_time - end_time}')

  

if __name__ == '__main__':
    
    # Setting up testing and logging
   # logging_file_setup(args.output_dir)
   # check = check_arguments(args)
 #   if check == 1:
   #     sys.exit(logging.error("Arguments provided were not expected. Please check log."))
    
    sys.exit(main())