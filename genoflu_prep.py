import pandas as pd
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import logging
import argparse
import datetime
from Bio import SeqIO
from src.utilities import util 

now = datetime.date.today()



__version__ = 1.0
__author__ = "Kate Howell"


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"New US genotype prep tool")


    parser.add_argument(
        "--genoflu", "-g", required=True, help="Genoflu folder"
    )

    parser.add_argument(
        "--previous_geno_tab", "-p", required=True, help="Previous genotype key from ai-blast-genotyping prep"
    )
    parser.add_argument(
        "--output_dir", "-o", required=True, help="Output folder location"
    )    
    args = parser.parse_args()
    # Need to handle output dir before setting up logging files.hat's the timeline for the analysis
    if not os.path.isdir(args.output_dir):  # Set up output folder
        print(f"Creating output folder:\n {args.output_dir}")
        logging.info(f"Creating output folder:\n {args.output_dir}")
        os.mkdir(args.output_dir)
    return args




def total_keys(test_dict):
        count = 0
        for key in test_dict:
            count += 1
            if isinstance(test_dict[key], dict):
                count += total_keys(test_dict[key])
        return count


def identify_new_genotypes(args, genotab):
    oldgenotab = pd.read_csv(args.previous_geno_tab)
    suboldgenotab = oldgenotab[oldgenotab['schema']=="US"]
    newgenolist = set(list(genotab['Genotype'])) - set(list(suboldgenotab['Genotype']))
    changedgeno = set(list(suboldgenotab['Genotype'])) - set(list(genotab['Genotype']))
    print(f'New genotyes found in Genoflu are: {", ".join(newgenolist)}')
    print(f'The following genotypes are no longer found (likely minor to major genotype): {", ".join(changedgeno)}')
    logging.info(f'New genotyes found in Genoflu are: {", ".join(newgenolist)}')
    logging.info(f'The following genotypes are no longer found (likely minor to major genotype): {", ".join(changedgeno)}')

def create_genotype_fastas(output_dir, genotab, subtype_list, seqdict, segments):
   # seqdict = {}
    for index, row in genotab.iterrows():
        with open(os.path.join(output_dir,f"{row['Genotype']}_ref.fasta"),"w") as seq_file:
            for s in segments:
                key = f'{row[s]}_{s}'
               
                seqdescript = seqdict[key][0].split(" ")
                seqdict[f'{seqdescript[1]}_{s}'] = [row[s],f"{row['Genotype']}_ref",row['Genotype'],f'H5N{subtype_list[index]}']
                newid = f'{row["Genotype"]}_ref|{row["Genotype"]}|H5N{subtype_list[index]}|US|{s}'
                file_content = SeqRecord(Seq(seqdict[key][1]), id=newid, description="")
                SeqIO.write(file_content, seq_file, "fasta")
    seqdf = pd.DataFrame.from_dict(seqdict, orient='index')
    seqdf.to_csv(os.path.join(output_dir,f"genoflu_reference_table_record_{now}.csv"))

def create_sequence_dict(fastas):
    seqdict = {}
    seqcount = 0
    for f in fastas:
        fasta_sequences = SeqIO.parse(open(f), 'fasta')
        for seq in fasta_sequences:
            info = seq.description.split(" ")

            seqcount += 1
            if f'{info[0]}_{info[2]}' in seqdict:
                if seq.description == seqdict[f'{info[0]}_{info[2]}'][0]:
                    pass
                else:
                    print(seq.description,seqdict[f'{info[0]}_{info[2]}'][0])
            seqdict[f'{info[0]}_{info[2]}'] = [seq.description,seq.seq]
            
    print(total_keys(seqdict))
    print(seqcount)
    return seqdict

def na_subtype(genotab):
    natype = genotab['NA']
    subtype_list = []
    for n in natype:
        if "N" in n:
            character = n.find('N')
            subtype_list.append(n[character+1])
        else:
            subtype_list.append("1")
    
    return subtype_list


def main(args):
    start_time = datetime.datetime.now()
      # Start time for calculating performance improvements
    folder = args.genoflu
    print(folder)
    genotab = pd.read_excel(os.path.join(folder, "genotype_key.xlsx"))
    logging.info(f'Genoflu table read in, {len(genotab["Genotype"])} in table....')
    #print(genotab.head)
    fastas = glob.glob(os.path.join(folder,"fastas/*fasta"))
    logging.info(f'{len(fastas)} FASTA files found....')
    logging.info('Parsing subtypes....')
    subtype_list = na_subtype(genotab)
    seqdict = create_sequence_dict(fastas)
    segments = list(genotab.columns[1:-1])
    segments.append("NS")
    create_genotype_fastas(args.output_dir, genotab, subtype_list, seqdict, segments)
    #which are in the new genotypes?   
    identify_new_genotypes(args, genotab)

    end_time = datetime.datetime.now()  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: { end_time - start_time}")



if __name__ == "__main__":
    args = read_commandline()
    # Setting up testing and logging
    util.logging_file_setup(args.output_dir,"no")
    check = util.check_arguments(args)

    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )

    sys.exit(main(args))

