import logging
import os
from Bio import SeqIO


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
