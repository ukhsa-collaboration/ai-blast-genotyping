import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from sqlalchemy import text, inspect, schema,create_engine
from datetime import date, datetime
from collections import Counter
import json
import subprocess
from subprocess import PIPE, CalledProcessError, check_call, Popen, STDOUT

now = date.today()

__version__ = 1.0
__author__ = "Kate Howell"



def read_commandline():
    """
    BLAST genotyping setup tool - command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(
        description=f"Script to update BLAST genotyping schemas"
    )
    parser.add_argument(
        "--geno_key", "-g", required=True, help="existing genotyping key"
    )
    parser.add_argument(
        "--reffasta", "-r", required=True, help="existing genotype reference fasta"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        required=True,
        default=os.getcwd(),
        help="Output folder. Default: CWD.",
    )
    input_group = (
        parser.add_mutually_exclusive_group()
    )  # input group requires one or the other
    input_group.add_argument(
        "--epiid", "-e", required=False, help="CSV file of new genotype sequences with meta-data"
    )
    input_group.add_argument(
        "--fasta", "-f", required=False, help="Input FASTA file for new genotype references (SPECIFIC INPUT FORMAT REQUIRED)"
    )
    parser.add_argument(
        "--username", "-u", required=False, help="username for aiseqdb query"
    )
    parser.add_argument(
        "--tagname", "-n", required=True, help="Analysis name for output files"
    )
    parser.add_argument(
        "--threshold", "-t", required=True, help="Threshold for BLAST filtering"

    )
    # Add threads
    args = parser.parse_args()

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


segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

repopath = os.path.dirname(sys.argv[0])

def logging_file_setup(output_folder):
    """
    Set up log

    :return: N/A
    """
    print("Setting up logging information for debugging.")
    # Need to clean up log files, currently generates multiple log files.
    logging_file_output = os.path.join(
        output_folder, str(str(now.strftime("%Y%m%d")) + str("_geno_logging_file.log"))
    )
    print("Logging output in the file: {}".format(logging_file_output))


    logging.basicConfig(
            filename=logging_file_output,
            filemode="w",
            format="%(asctime)s:%(levelname)s:%(message)s",
            level=logging.INFO,
    )
    print("Logging Level: INFO")
    logging.info("Logging Started.")


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


def per_segment_fasta(segments, fastafile, output_dir):
    """
    Create FASTA file of the genotype sequences per segment, run alignments & trees

    :param segments: segment list

    :param fastafile: all genotypes references FASTA file

    :param output_dir: output directory

    :return: N/A
    """
    for s in segments:
        seq_file = open(os.path.join(output_dir,f"all_genotype_references_{now}_{s}.fasta"), "w")
        fasta_sequences = SeqIO.parse(open(fastafile), 'fasta')
        for seq in fasta_sequences:
            if s in seq.id:
                SeqIO.write(seq, seq_file, "fasta")
        #try:
            #print('mafft', os.path.join(output_dir, f"all_genotype_references_{now}_{s}.fasta"), '>',os.path.join(output_dir,f"all_genotype_references_{now}_{s}.aln"))
        #    aln_seqs = subprocess.Popen(['mafft', os.path.join(output_dir, f"all_genotype_references_{now}_{s}.fasta"), '>',os.path.join(output_dir,f"all_genotype_references_{now}_{s}.aln")],stdout=PIPE, stderr=PIPE)
         #   out, err = aln_seqs.communicate()
         #   result = out.decode()
         #   print("Error : ",err )
         #   print("Result : ",result )
         #   aln_seqs.wait()
        #except CalledProcessError as exc:
        #    print(exc.output)  


def check_geno_exists(dbextract):
    """
    Check if the exact sequence already exists in the FASTA file and remove any duplicates

    :param db_extract: aiseqdb table extract

    :return: reduced dbextract
    """
   # dbextract = check_fasta_forids(dbextract,'isolate_epi_id')
    dbextract = check_fasta_forids(dbextract,'isolate_name')
    return dbextract

def check_fasta_forids(dbextract,columnname):
    for i in list(set(dbextract[columnname])):
        idfound = False
        fasta_sequences = SeqIO.parse(open(os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta")), 'fasta')
        for seq in fasta_sequences:
            if i in seq.id:
                idfound = True
        if idfound == True:        
            print(f"{columnname} already in genotyping base files, skipping, {i}")
            indexseq = dbextract[dbextract[columnname] == i].index
            dbextract.drop(indexseq, inplace=True)
    return dbextract

def create_fasta_seqs(db_extract):
    """
    Create FASTA file of the new genotype sequences and update the existing file

    :param db_extract: aiseqdb table extract

    :return: All gentoype fasta file
    """
    db_extract["segment_name"] = db_extract["segment_name"].fillna("NA")
    db_extract['header'] = db_extract.isolate_name.astype(str).replace(" ","_") + '|'  + \
    db_extract.Genotype.astype(str) + '|' +db_extract.subtype_x.astype(str) + '|' +\
    db_extract.segment_name.astype(str) 
    filtered_df = db_extract[db_extract["na_sequence"].notnull()]
    filtered_df = check_geno_exists(filtered_df)
    seq_file = open(os.path.join(args.output_dir,f"new_genotype_references_{now}.fasta"), "w")
    updated_seq_file = open(os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta"), "a")
    for index, row in filtered_df.iterrows():
        fasta_header = row["header"]
        fasta_sequence = row["na_sequence"]
        file_content = SeqRecord(Seq(fasta_sequence), id=fasta_header, description="")
        SeqIO.write(file_content, seq_file, "fasta")
        SeqIO.write(file_content, updated_seq_file, "fasta")

    return os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta"), filtered_df
    
def update_geno_key(genokey, subepitab):
    """
    Update the genotyping key table with new genotypes

    :param genokey: existing genotyping key table 
    :param subepitab: new rows for genotyping key table


    :return: combined genotyping key table
    """
    genotab = pd.read_csv(genokey)
    newdf = pd.concat([genotab, subepitab])
    newname = genokey.replace(".csv",f"_{now}.csv")
    newdf.to_csv(os.path.join(args.output_dir,f'{newname}'))
    return newdf


def new_fasta_parsing(seqfile):
    """
    Parse provided FASTA file to retrieve new genotype meta-data and sequences

    :param seqfile: FASTA file for new genotype sequences

    :return: new rows for genotyping key table
    """
    newgenodict = {}
    fasta_sequences = SeqIO.parse(open(seqfile), 'fasta')
    updated_seq_file = open(os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta"), "a")
    for seq in fasta_sequences:
        info = seq.id.split("|")
        newgenodict[info[1]] = [info[0], info[2], info[3], info[4]] 
        seq.id = f'{info[0].replace(" ","_")}|{info[1]}|{info[2]}|{info[4]}'
        seq.description = ""
        SeqIO.write(seq, updated_seq_file, "fasta")
    newgenotab = pd.DataFrame.from_dict(newgenodict, orient="index")
    newgenotab = newgenotab.reset_index()
    newgenotab.columns = ['Genotype', 'sequence', 'subtype', 'schema', 'segment']

    #remove column for segment and de-duplicate
    newgenotab = newgenotab.drop('segment', axis=1)
    newgenotab = newgenotab.drop_duplicates()

    return newgenotab

def blast_processing():
    """
    Run commands to create new BLASTdb and run the query 

    :return: BLAST query table output
    """
  #  subprocess create blast db
    logging.info("Creating the new BLAST database")

    subprocess.call(
        f"makeblastdb -in {os.path.join(args.output_dir,f'all_genotype_references_{now}.fasta')} -out {os.path.join(args.output_dir,f'all_genotype_references_{now}.db')} -dbtype nucl",
        shell=True,
    )
    logging.info("Performing the BLAST searches per FASTA file")
    blasttab = os.path.join(args.output_dir, f'all_genotype_references_{now}.blast.out')
    #subprocess query blast db with file
    subprocess.call(
        f"blastn -db {os.path.join(args.output_dir,f'all_genotype_references_{now}.db')} -query {os.path.join(args.output_dir,f'all_genotype_references_{now}.fasta')} -out {blasttab} -outfmt 6",
        shell=True,
    )
    return blasttab


def filter_blast_results(blasttab):
    """
    Tidy up the BLAST results table and output csv version of table

    :param blasttab: BLAST output table
    :return: tidied BLAST query table output and list of isolates
    """
    blast = pd.read_csv(blasttab,sep="\t")
    blast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    blast[['query_isolate_name', 'query_genotype','query_subtype','query_segment']] = blast['qseqid'].str.split('|', expand=True)
    blast[['hit_isolate_name', 'hit_genotype','hit_subtype','hit_segment']] = blast['sseqid'].str.split('|', expand=True)
    blast.to_csv(os.path.join(args.output_dir,f"blast_h5n1_geno_refs_{now}.csv"))
    queries = list(set(blast['query_isolate_name']))
    return blast,queries

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def check_overlap(q,grouplist,groupdict,ingroup):
    """
    Review groups per segment and merge highly overlapping groups
    

    :return: N/A
    """
    maxoverlap = 0
    groupfound = ""
    for key, val in groupdict.items():
        overlap = intersection(list(set(grouplist)), val)   
        if len(overlap) > maxoverlap:
            maxoverlap = len(overlap)
            groupfound = key

    groupdict, groupdecision, ingroup = decide_group_merge(q, maxoverlap, groupfound, grouplist, groupdict, ingroup)
    return groupdict, groupdecision, ingroup

def decide_group_merge(q,maxoverlap,groupfound, grouplist, groupdict, ingroup):
    """
    Review groups per segment and merge highly overlapping groups
    

    :return: N/A
    """
    if float(maxoverlap) / float(len(set(grouplist))) >= 0.75:
        groupmembers = groupdict[groupfound]
        print(f"Group found with > 75% matches to existing group, {q} {groupfound}")
        logging.info(f"Group found with > 75% matches to existing group, {q} {groupfound}")
        logging.info(f"Adding {q} to {groupfound}")
        newgrouplist = []
        #add query and other members of group to new list
        for g in grouplist:
            newgrouplist.append(g)
        #add in the existing group members
        for g in groupmembers:
            newgrouplist.append(g)
        groupdict[groupfound] = list(set(newgrouplist))
        groupdecision = "merge"
        ingroup = True
        return groupdict, groupdecision, ingroup
    else:
        groupdecision = "separate"
        return groupdict, groupdecision, ingroup
    
def create_group_json(queries, t, subblast, s,subtype):
    """
    Create the json table 

    :param queries: list of isolate_names in BLAST query
    :param t: filter threshold for percentage identity (default: 98%)
    :param subblast: Filtered BLAST table per segment
    :param s: Segment

    :return: N/A
    """
    segblast = subblast[subblast['query_segment'] == s]
    groupdict = {}
    groups = 1
    #remove duplicate sequences using set
    for n, q in enumerate(list(set(queries))):
       
        hitblast = segblast[segblast['query_isolate_name'] == q]
       
       
        #if no hits then no segment available, otherwise would match to itself
        if hitblast.shape[0] == 0:
            print(f"{q} no {s} found")
            pass

        else:
            if subtype == True:
                
                h5type = hitblast['query_subtype'].iloc[0]
                h5type = h5type.replace("A/", "")
                hitblast['hit_subtype'] = hitblast['hit_subtype'].str.replace("A/", "")
                hitblast = hitblast[hitblast['hit_subtype'] == h5type]
            
            ingroup = False
                #check if query is already in a group in dictionary 
            for key, val in groupdict.items():
                    
                    if q in list(val):
                        ingroup = True
                        pass
            while not ingroup:
                    #if not in the group
                    grouplist = list(hitblast['hit_isolate_name'])
                    grouplist.append(q) 
                    #check if it overlaps with an existing group, and by how much
                    groupdict, groupdecision,ingroup = check_overlap(q, grouplist, groupdict, ingroup)
                    if groupdecision == "separate":
                        #new group needed
                        groupdict, groups, ingroup = new_group(s, groupdict, groups, grouplist)
                        pass
                    elif groupdecision=="merge":
                        #already added to a group
                        pass
                    else:
                        print("Not merge or separate??")
            qfound = False
            count = 0
            for key, val in groupdict.items():
                    if q in list(val):
                        #   print(q, key)
                        qfound = True
                        count = count+1
            if not qfound:
                    print(q,"Not added to dictionary", ingroup)
            if count>1:
                    print(q,s,"duplicated!!!")
                    remove = 0
                    for key, val in groupdict.items():
                        if q in list(val):
                            remove = remove + 1
                            if remove > 1:
                                val = val.remove(q)
                                groupdict[key] = val

    print(f'Segment {s}:{len(groupdict.keys())} groups')
    with open(os.path.join(args.output_dir,f'blast_geno_threshold_{t}_{s}.json'), 'w') as fp:
        json.dump(groupdict, fp)

def new_group(s, groupdict, groups, grouplist):
    newgrouplist = []
    for g in grouplist:
        ng = False
        for key, val in groupdict.items():
            if g in list(val):
                ng = True
                exit   
        if not ng:
             newgrouplist.append(g)
    # print(ne)
    logging.info(f"New group defined: {s}_group{groups}, {list(set(newgrouplist))}")
    groupdict[f'{s}_group{groups}'] = list(set(newgrouplist))
    ingroup = True
    groups = groups+1
    return groupdict, groups, ingroup

def create_constellation(queries, blast, tabcols, t):
    """
    Create the constellation table per genotype

    :param queries: list of isolate_names in BLAST query
    :param t: filter threshold for percentage identity (default: 98%)
    :param blast: BLAST table 
    :param tabcols: columns to include in table

    :return: constellations dataframe for specific threshold
    """
    testdf = []
    for q in queries:
        hitblast = blast[blast['query_isolate_name']==q]
        
        seqinfo = [q,hitblast['query_genotype'].iloc[0],hitblast['query_subtype'].iloc[0]]
        for s in segments:
            with open(os.path.join(args.output_dir,f'blast_geno_threshold_{t}_{s}.json'), 'r') as file:
                data = file.read()
            tsdict = json.loads(data)
            found = False
            if not found:
                for key, val in tsdict.items():
                        if q in list(val):
                            seqinfo.append(key)
                            found = True

        if len(seqinfo)>11:
            print(seqinfo)
        testdf.append(seqinfo)
    thresholddf = pd.DataFrame(testdf, columns = tabcols)
    for s in segments:
        thresholddf[s] = thresholddf[s].fillna("Absent") 
    thresholddf["constellation"] = thresholddf[['PB2','PB1','PA','HA','NP','NA','MP','NS']].agg("|".join, axis=1)
    thresholddf.to_csv(os.path.join(args.output_dir,f'blast_geno_threshold_table{t}_{now}.csv'))
    return thresholddf



def describe_geno(genotype,constellation,thresholddf,segments):
    """
    Describe the new genotypes and relatives in groups

    :param constellation: constellation for new genotype
    :param thresholddf:constellation table

    :return: N/A
    """
    groups = constellation.split("|")
    genotype = genotype.replace("\/", "_")
    try:
        with open(os.path.join(args.output_dir,f"genotype_description_{genotype}.txt"),"w") as filex:
            filex.write(f'{genotype} has the following relatives per segment\n')
      #  print(f'{genotype} has the following relatives per segment')
            for n,g in enumerate(groups):
                seg_threshold = thresholddf[thresholddf[segments[n]] == g]
        #     print(f'{segments[n]}: {"|".join(list(seg_threshold["genotype"]))}')
                filex.write(f'{g}: {"|".join(list(seg_threshold["genotype"]))}\n')
    except:
        pass


def report_new_genotypes(thresholddf, genolist):
    """
    Report on results for the new genotypes and flag duplicates

    :param thresholddf: constellation table
    :param genolist: new genotype list


    :return: N/A
    """
    #print(genolist)
    subthresholds = thresholddf[thresholddf['genotype'].isin(genolist)]
    constellations = []
    for index, row in subthresholds.iterrows():
        print(f"{row['genotype']} has the following constellation: {row['constellation']}")
        constellations.append(row['constellation'])
        describe_geno(row['genotype'],row['constellation'],thresholddf,segments)
    dup_check = thresholddf[thresholddf['constellation'].isin(constellations)]
    if dup_check.shape[0] > subthresholds.shape[0]:
        print("Duplicate constellation found for at least one genotype....")
        thresholddf.to_csv(f'new_genotype_duplicate_constellations_tab_{now}.csv')
    for c in constellations:
        dup_check = thresholddf[thresholddf['constellation'] == c]
        if dup_check.shape[0] > 1:
            print(f"{c} has {len(dup_check['genotype'])} genotypes:{', '.join(list(dup_check['genotype']))}")
            print("Additional review required to determine if extra layer of SNP typing is required!!")

def geno_groups(queries,blast,t,newgeno):
    """
    Wrapper for running the grouping process per segment an creating the constellation

    :param queries: list of sequences
    :param blast: BLAST output table
    :param t: percentage identity threshold for filtering
    :param genolist: list of genotypes

    :return: N/A
    """
    tabcols = ['sequence','genotype','subtype']
    tabcols.extend(segments)
    # we want the different combinations / constellations in a table


    subblast = blast[blast['pident']>=t]

    for s in segments:
        if s == "NA":
            create_group_json(queries, t, subblast, s,subtype=True)
        else:
            create_group_json(queries, t, subblast, s,subtype=False)
    thresholddf = create_constellation(queries, blast, tabcols, t)
    report_new_genotypes(thresholddf,newgeno)



def main(args):
    """
    Main running of the script to run the BLAST query and wrangle the results to provide a per segment and sample summary of the genotyping results.

    :return: N/A
    """
    start_time = datetime.now()  # Start time for calculating performance improvements
    # Create paths for qsub submission
    output_dir = os.path.abspath(args.output_dir)


    # Set up logging
    logging_file_setup(output_dir)
    # Read in required files depending if CSV or FASTA provided
    outputfasta = os.path.join(args.output_dir,f'all_genotype_references_{now}.fasta')
   
    subprocess.call(f"cp {args.reffasta} {outputfasta}",shell=True)
    
    if args.epiid is not None:
        epitab = pd.read_csv(args.epiid)
        ids = "','".join(epitab['sequence'])
        print(f"{len(epitab['sequence'])} new genotypes to be added to the database")
        db_extract = query_aiseqdb(args.username, ids)
        for i in epitab['sequence']:
            subepitab = db_extract[db_extract['isolate_epi_id'] == i]
            print(f'{i} has {subepitab.shape[0]} segments available')
            if subepitab.shape[0] < 6:
                print(f'{i} has less than 6 segments available. Constellation will be skipped. Genotype NOT added.')
                logging.info(f'{i} has less than 6 segments available. Constellation will be skipped. Genotype NOT added.')
                indexseq = epitab[epitab['sequence'] == i].index
                epitab.drop(indexseq, inplace=True)
            elif subepitab.shape[0] < 8:
                print(f'{i} does not have all 8 segments available. Constellation will be incomplete')
        db_extract['isolate_name']=db_extract['isolate_name'].str.replace(" ","_")
        db_extract2 = pd.merge(epitab,db_extract,right_on='isolate_epi_id',left_on='sequence',how="left")
        fastafile, filtereddf = create_fasta_seqs(db_extract2)
        per_segment_fasta(segments, fastafile, output_dir)
        epitab = epitab[epitab['sequence'].isin(list(filtereddf['isolate_epi_id']))]
        subepitab = epitab[['Genotype', 'sequence', 'subtype', 'schema']]
        genotab = update_geno_key(args.geno_key,subepitab)

    if args.fasta is not None:
        newgenotab = new_fasta_parsing(args.fasta)
        print(f"{len(newgenotab['sequence'])} new genotypes to be added to the database")
        genotab = update_geno_key(args.geno_key, newgenotab)
        epitab = newgenotab
    try:
        blasttab = blast_processing()
    except:
        print("Unable to run the BLAST processing, exiting...")
        sys.exit()
    blast,queries = filter_blast_results(blasttab)
    geno_groups(queries, blast,int(args.threshold),epitab['Genotype'])
    

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