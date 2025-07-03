import pandas as pd
import os
import sys

import logging
import argparse
import datetime
import subprocess
import src.utilities as util
import src.fasta_processing as fap

now = datetime.date.today()


__version__ = 1.0
__author__ = "Kate Howell"


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(
        description=f"Genotype ingest preparation for aiseqdb"
    )
    

    parser.add_argument("--genoflu", "-g", required=True, help="Genoflu results")
    # parser.add_argument(
    #     "--gdir", "-d", required=True, help="Genoflu directory"
    # )
    parser.add_argument("--fasta", "-f", required=True, help="Genotyping FASTA file")

    parser.add_argument("--blast", "-b", required=True, help="BLAST genotyping results")
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


def get_git_tag(folderdir):
    dir_path = os.getcwd()
    os.chdir(folderdir)
    cmd2 = ["git", "describe", "--tags"]
    proc = subprocess.Popen(cmd2, stdout=subprocess.PIPE)
    version = proc.stdout.read()
    os.chdir(dir_path)
    version = version.decode("utf-8")
    version = version.replace("\n", "")
    print(f"Git tag of {folderdir} is {version}")
    return version


def main(args):
    start_time = datetime.datetime.now()
    # Start time for calculating performance improvements
    blasttab = pd.read_csv(args.blast)

    genotab = pd.read_csv(args.genoflu)
    newtab = pd.merge(
        blasttab, genotab, right_on="sample", left_on="isolate_epi_id", how="left"
    )
    #
    subtab, subblast2 = wrangle_blast_tab(newtab)
    subgeno2, seqdf = genoflu_table_wrangling(subtab)

    # combine all the tables & rename columns

    comb1 = pd.merge(
        subgeno2,
        seqdf,
        left_on=["isolate_epi_id", "segment_name"],
        right_on=["isolate_epi_id", "segment_name"],
        how="left",
    )

    subblast2["segment"] = subblast2["segment"].fillna("NA")
    comb1["segment_name"] = comb1["segment_name"].fillna("NA")
    comb2 = pd.merge(
        comb1,
        subblast2,
        left_on=["isolate_epi_id", "segment_name"],
        right_on=["isolate_epi_id", "segment"],
        how="left",
    )
    comb2["segment_sequence_id"] = comb2["segment_sequence_id"].fillna("")
    comb2["final_blast_genotype"] = comb2["final_blast_genotype"].str.replace(
        "No known genotype: Please review individual segments results", "Not assigned"
    )
    comb2 = comb2.drop(columns=["segment"])
    comb2["final_genoflu_genotype"] = comb2["final_genoflu_genotype"].fillna(
        "Not assigned"
    )
    comb2.to_csv(
        os.path.join(args.output_dir, f"{now}_genotype_ingest.csv"), index=False
    )
    end_time = (
        datetime.datetime.now()
    )  # end time for calculating performance improvements

    logging.info(f"Pipeline time to completion: { end_time - start_time}")


def genoflu_table_wrangling(subtab):
    genoflucols = ["isolate_epi_id", "Genotype List Used, >=98%", "Genotype_simplified"]
    # melt the dataframe, how to get the count of segments, how we get the version of the BLAST genotyping
    gftab = subtab[genoflucols]
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

    for s in segments:
        info = []
        for i in list(gftab["Genotype List Used, >=98%"]):
            found = False
            try:
                segs = i.split(",")

                for g in segs:
                    if s in g:
                        details = g.split(":")
                        info.append(details[1])
                        found = True
            except:
                pass
            if not found:
                info.append("missing or below threshold")
        # gftab[s] = info
        gftab.insert(len(gftab.columns), s, info)
    genoflucols2 = [
        "isolate_epi_id",
        "Genotype_simplified",
        "PB2",
        "PB1",
        "PA",
        "HA",
        "NP",
        "NA",
        "MP",
        "NS",
    ]
    subbgeno = pd.melt(
        gftab[genoflucols2], id_vars=["isolate_epi_id", "Genotype_simplified"]
    )

    # rename the columns
    subbgeno_nomissing = subbgeno[subbgeno["value"] != "missing or below threshold"]
    df2 = (
        subbgeno.groupby(["isolate_epi_id"]).size().reset_index(name="number_segments")
    )
    subgeno2 = pd.merge(
        subbgeno, df2, left_on="isolate_epi_id", right_on="isolate_epi_id", how="left"
    )
    subgeno2["genoflu_version"] = get_git_tag(
        "/home/phe.gov.uk/kate.howell/Documents/GenoFLU"
    )
    subgeno2["ingest_date"] = now
    subgeno2.columns = [
        "isolate_epi_id",
        "final_genoflu_genotype",
        "segment_name",
        "genoflu_group",
        "number_segments",
        "genoflu_version",
        "ingest_date",
    ]
    seqdf = fap.get_isolate_id_fasta(args.fasta)
    return subgeno2, seqdf


def wrangle_blast_tab(newtab):
    droplist = [
        "consensus",
        "sample",
        "Genotype",
        "Genotype Sample Title List",
        "Genotype Percent Match List",
        "Genotype Mismatch List",
        "Genotype Average Depth of Coverage List",
    ]

    keepcols = set(newtab.columns) - set(droplist)

    keepcols = list(set(newtab.columns)-set(droplist))
    #print(keepcols)
    #print(set(keepcols) - set(newtab.columns))
    subtab = newtab[keepcols]

    blastcols = [
        "isolate_epi_id",
        "PB2",
        "PB1",
        "PA",
        "HA",
        "NP",
        "NA",
        "MP",
        "NS",
        "Top_Hit",
    ]
    # melt the dataframe, how to get the count of segments, how we get the version of the BLAST genotyping
    subblast = pd.melt(subtab[blastcols], id_vars=["isolate_epi_id", "Top_Hit"])

    # rename the columns
    subblast.columns = [
        "isolate_epi_id",
        "final_blast_genotype",
        "segment",
        "blast_group",
    ]

    return subtab, subblast


if __name__ == "__main__":
    args = read_commandline()
    # Setting up testing and logging
    util.logging_file_setup(args.output_dir, "no")
    check = util.check_arguments(args)

    if check == 1:
        sys.exit(
            logging.error("Arguments provided were not expected. Please check log.")
        )

    sys.exit(main(args))
