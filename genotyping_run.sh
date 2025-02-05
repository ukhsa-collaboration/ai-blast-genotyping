#!/bin/sh
source ~/.bashrc

## change to routine dir
cd ~/scheduled_tasks

## log the run starting
echo "`date` Starting run" >> genotyping_run.log

## activate conda and run script
conda activate blastgeno
/home/phe.gov.uk/kate.howell/ai-blast-genotyping/run_genotyping.py &>> genotyping_run.log

date=$(date '+%Y-%m-%d')
source /home/phe.gov.uk/kate.howell/ai_seq_db_env2/bin/activate

ai_seq_db import_genotypes --data $(date '+%Y-%m-%d')_genotype_ingest.csv --user kate.CRON --notes "db ingest CRON job $(date '+%Y-%m-%d')" --data_dir /home/phe.gov.uk/kate.howell/data/genotyping/$(date '+%Y-%m-%d')


## log the run completing
echo "`date` Run complete" >> genotyping_run.log
