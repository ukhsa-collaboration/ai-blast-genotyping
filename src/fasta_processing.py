import logging
import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

    logging.info(
        print(
            f"{len(duplist)} duplicates found in FASTA files. Any duplicates should be reviewed prior to interpreting results!"
        )
    )
    return duplist



def get_isolate_id_fasta(fasta):
    seqdict = {}
    fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
    for seq in fasta_sequences:
        if "|" in seq.id:
            seporator = "|"
        elif "." in seq.id:
             seporator = "."
        else:
            logging.warning("no known seporator found, exiting!")
            sys.exit()
        info = seq.id.split(seporator)
        seqdict[info[1]] = [info[0],info[2]]
    seqdf = pd.DataFrame.from_dict(seqdict, orient='index')
    seqdf =seqdf.reset_index()
    seqdf.columns = ['segment_sequence_id','isolate_epi_id','segment_name']
    return seqdf



def per_segment_fasta(now,segments, fastafile, output_dir):
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


def check_geno_exists(args,now,dbextract):
    """
    Check if the exact sequence already exists in the FASTA file and remove any duplicates

    :param db_extract: aiseqdb table extract

    :return: reduced dbextract
    """

    dbextract = check_fasta_forids(args,now,dbextract,'isolate_name')
    return dbextract

def check_fasta_forids(args,now,dbextract,columnname):
    """
    Look up if the isolate_name is already part of the FASTA file, skip if present. 

    :param db_extract: aiseqdb table extract
    :param columnname: aiseqdb table extract
    :return: reduced db_extract
    """
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

def create_fasta_seqs(args,now,db_extract):
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
    filtered_df = check_geno_exists(args,now,filtered_df)
    seq_file = open(os.path.join(args.output_dir,f"new_genotype_references_{now}.fasta"), "w")
    updated_seq_file = open(os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta"), "a")
    for index, row in filtered_df.iterrows():
        fasta_header = row["header"]
        fasta_sequence = row["na_sequence"]
        file_content = SeqRecord(Seq(fasta_sequence), id=fasta_header, description="")
        SeqIO.write(file_content, seq_file, "fasta")
        SeqIO.write(file_content, updated_seq_file, "fasta")

    return os.path.join(args.output_dir,f"all_genotype_references_{now}.fasta"), filtered_df


def new_fasta_parsing(args,now,seqfile):
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


def write_fasta_file(now,df, filetag, output_dir):
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
