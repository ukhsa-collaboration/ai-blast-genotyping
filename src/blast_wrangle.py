import os
import logging
import src.utilities as util
import src.fasta_processing as fap
from collections import Counter
import pandas as pd
import sys
import glob
import subprocess
import numpy as np

global genogroups
global genodict
global constellation_path


def process_ref_tables(constellation_path):
    genoblast = pd.read_csv(os.path.join(constellation_path), dtype=str)

    genogroups = pd.melt(
        genoblast,
        id_vars=["sequence", "genotype", "subtype", "constellation"],
        value_vars=util.segments,
    )

    genogroups = genogroups.rename(columns={"variable": "Segment", "value": "Group"})
    genogroups["Segment"] = genogroups["Segment"].replace(np.nan, "NA")

    genodict = genogroups.groupby("Group")["genotype"].apply(list).to_dict()
    genodict2 = {}
    for key in genodict:
        genodict2[key] = "|".join(genodict[key])
    genodict = genodict2
    return genogroups, genodict


def match_genotype_dict(blast_pass, genogroups, genodict):
    """
    Match the BLAST hit to a genotype group and possible genotype match per segment

    :return: BLAST hit dataframe
    """
    print("Matching genotype to dictionary")
    segdict = dict(
        zip(
            genogroups["Segment"] + "_" + genogroups.sequence.str.replace("_", ""),
            genogroups.Group,
        )
    )
    # map the query sequence back to the genotyping table
    # need to look up which group the reference is in
    blast_pass["group_match"] = (
        blast_pass["ref_match"].apply(lambda x: segdict.get(x)).fillna("")
    )
    blast_pass["genotype_match"] = (
        blast_pass["group_match"].apply(lambda x: genodict.get(x)).fillna("")
    )
    return blast_pass


#


def create_hit_dict(blast_pass):
    hitdict = {}
    for q in list(set(blast_pass["isolate_epi_id"])):
        subblast = blast_pass[blast_pass["isolate_epi_id"] == q]
        if subblast.shape[0] >= 1:
            hitdict[subblast["isolate_epi_id"].iloc[0]] = [
                dict(Counter(list(subblast["hit_isolate_name"]))),
                Counter(list(subblast["hit_geno"])),
            ]
    return hitdict


def check_constellations(args):
    if args.constellation == "reference_files/blast_geno_threshold_table98.csv":
        repopath = os.path.dirname(sys.argv[0])
        constellation_path = os.path.join(repopath, args.constellation)
    else:
        constellation_path = args.constellation

    if os.path.exists(constellation_path):
        genogroups, genodict = process_ref_tables(constellation_path)
    else:
        logging.error(
            'Reference genotypes table does not exist:"genotype_groups_examples.csv". Please check if the file is in the correct location'
        )
    return constellation_path


def tidy_blast_table(args, folder, sample, segmissing, fasta_count):
    """
    Wrange the BLAST results table, exclude results that do not meet threshold, separate no sequence vs. no hit, split the query and hit ids into isolate_epi_id and segment.

    :return: newresultstab dataframe
    """
    print("Tidying BLAST tables....break")
    logging.info("Reading in BLAST output:")
    logging.info(os.path.join(args.output_dir, f"{sample}.blast.out"))
    blasttab = pd.read_csv(
        os.path.join(args.output_dir, f"{sample}.blast.out"), sep="\t", header=None
    )
    if blasttab.shape[0] < fasta_count:
        print(
            f"{fasta_count - blasttab.shape[0]} FASTAs do not meet minimum BLAST thresholds"
        )
        logging.info(
            f"{fasta_count - blasttab.shape[0]} FASTAs do not meet minimum BLAST thresholds"
        )
    blasttab.columns = util.blast_cols
    threshold = int(args.identity)
    blast_pass = blasttab[blasttab["pident"] >= threshold]
    blast_pass = blast_pass.sort_values(["qseqid", "pident"], ascending=[True, False])
    blast_fail = blasttab[blasttab["pident"] < threshold]
    if blast_fail.shape[0] >= 1:
        logging.info(
            f"{blast_fail.shape[0]} sequences do not meet a minimum percentage sequence identity to the ref seq!!"
        )
        print(
            f"{blast_fail.shape[0]} sequences do not meet a minimum percentage sequence identity to the ref seq!!"
        )
    logging.info("BLAST pass table dimensions:")
    logging.info(blast_pass.shape)
    if blast_pass.shape[0] == 0:
        print(
            "No sequences found to meet the minimum threshold. Tool now exiting, please check that you are using the same subtype of sequences as in the database!"
        )
        sys.exit(
            logging.error(
                "No sequences found to meet the minimum threshold. Tool now exiting, please check that you are using the same subtype of sequences as in the database!"
            )
        )
    logging.info("Checking sequence identifier for delimiters.... ")
    if "|" in list(blast_pass["qseqid"])[0]:
        delim = "|"

    elif "." in list(blast_pass["qseqid"])[0]:
        delim = "."
    elif "_" in list(blast_pass["qseqid"])[0]:
        delim = "_"
    else:
        logging.error(
            'Sample header delimiter unknown! Please check, "." or "|" expected!'
        )
        print('Sample header delimiter unknown! Please check, "." or "|" expected!')
    logging.info(f"Delimiter found: {delim}")
    logging.info("Preparing BLAST table for summary")
    try:
        test = blast_pass["qseqid"].str.split(pat=delim, expand=True)
    except:
        print(f"Issue with delimiter! {delim}")
    blast_pass.insert(len(blast_pass.columns), "segment", test[test.columns[-1]])
    blast_pass.insert(len(blast_pass.columns), "isolate_epi_id", test[test.columns[0]])
    blast_pass.insert(
        len(blast_pass.columns),
        "ref_match",
        blast_pass["segment"]
        + "_"
        + blast_pass["sseqid"].map(lambda x: x.split("|")[0]).str.replace("_", ""),
    )
    blast_pass[["hit_isolate_name", "hit_geno", "hit_subtype", "hit_segment"]] = (
        blast_pass["sseqid"].str.split("|", expand=True)
    )
    constellation_path = check_constellations(args)
    genogroups, genodict = process_ref_tables(constellation_path)
    blast_pass = match_genotype_dict(blast_pass, genogroups, genodict)
    blast_pass.to_csv(os.path.join(folder, f"{sample}.blast.out2"))

    results_df = pd.DataFrame()
    results_df.insert(
        len(results_df.columns), "genotype_match", blast_pass["genotype_match"]
    )
    results_df.insert(len(results_df.columns), "segment", blast_pass["segment"])
    results_df.insert(
        len(results_df.columns), "isolate_epi_id", blast_pass["isolate_epi_id"]
    )
    results_df.insert(len(results_df.columns), "top_hit", blast_pass["ref_match"])
    newresults_df = [results_df]

    hitdict = create_hit_dict(blast_pass)
    print("hit dictionary created....")
    if blast_fail.shape[0] >= 1:
        test = blast_fail["qseqid"].str.split(pat=delim, expand=True)
        blast_fail.insert(len(blast_fail.columns), "segment", test[test.columns[-1]])
        blast_fail.insert(
            len(blast_fail.columns), "isolate_epi_id", test[test.columns[0]]
        )
        missing = blast_fail[["segment", "isolate_epi_id"]]
        missing.insert(
            len(missing.columns),
            "genotype_match",
            ["No match via BLAST"] * missing.shape[0],
        )
        newresults_df.append(missing)
    newresultstab = pd.concat(newresults_df)
    samples = list(set(newresultstab["isolate_epi_id"]))

    for s in samples:
        subresults = newresultstab[newresultstab["isolate_epi_id"] == s]
        if subresults.shape[0] == 1:
            pass
        elif subresults.shape[0] == 0:
            submissing = segmissing[segmissing["sample"] == s]
            segfile = list(set(newresultstab["segment"]))[0]
            if submissing[segfile] == "missing":
                newresultstab.loc[len(newresultstab)] = ["No sequence", segfile, s]
            else:
                newresultstab.loc[len(newresultstab)] = ["BLAST FAIL", segfile, s]
        elif subresults.shape[0] < 8:
            submissing = segmissing[segmissing["sample"] == s]

            segfound = list(set(subresults["segment"]))
            segcheck = list(set(util.segments) - set(segfound))
            for seg in segcheck:

                if submissing[seg].iloc[0] == "missing":
                    newresultstab.loc[len(newresultstab)] = [
                        "No sequence",
                        seg,
                        s,
                        "Null",
                    ]
                else:
                    newresultstab.loc[len(newresultstab)] = [
                        "BLAST FAIL",
                        seg,
                        s,
                        "Null",
                    ]

    return newresultstab, hitdict


def run_blast(db, folder, sample, output_dir):
    """
    Run the BLAST query using sub-process

    :return: N/A
    """
    fap.tidy_fasta_files(os.path.join(folder, sample))
    logging.info("Performing the BLAST searches per FASTA file")
    print(db)
    subprocess.call(
        f"blastn -db {db} -query {os.path.join(folder,sample)} -out {os.path.join(output_dir,sample)}.blast.out -outfmt 6 -max_target_seqs 1",
        shell=True,
    )


def run_full_blasts(repopath, args, now, folder, mode, extension, output_dir):
    """

    Perform several checks on the FASTA file(s) and run the BLAST queries. For the multiple file option, concatenate the results. Match the BLAST hits to the genotyping reference table.

    :return: BLAST results dataframe
    """
    print("Running BLASTs")
    results_tabs = []
    logging.info("Input folder provided:")
    logging.info(os.path.join(folder, "*" + extension))
    fastas = glob.glob(os.path.join(folder, "*" + extension))

    segdict = {}
    if mode == "all_in_folder":
        logging.info("Fasta files found:")
        logging.info(fastas)
        segtabs = []
        logging.info("Running BLAST")
        print("Running BLAST on input folder")
        hitdict = {}
        for f in fastas:

            segdict = fap.missing_fasta_check(f, segdict)
            duplist = fap.duplicate_fasta_check(f)

            run_blast(os.path.join(repopath, args.blastdb), folder, f, output_dir)

            segmissing, fasta_count = util.create_segment_tab(segdict)

            segdicttab = pd.DataFrame.from_dict(segmissing, orient="index")
            segdicttab = segdicttab.reset_index()
            segtabs.append(segdicttab)
            segcols = ["sample"]
            for s in util.segments:
                segcols.append(s)
            segdicttab.columns = segcols

            sresults, hitdict_perfile = tidy_blast_table(
                args, args.output_dir, f, segdicttab, fasta_count
            )
            hitdict.update(hitdict_perfile)
            results_tabs.append(sresults)
        newdf = pd.concat(results_tabs)
        segtab = pd.concat(segtabs)
        segtab.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_segment_table.csv")
        )
        logging.info(
            f"Segment table written to : {os.path.join(output_dir,f'{now}_{args.tagname}_segment_table.csv')}"
        )
        # add step to sort out missing segments table

        newdf.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        logging.info("Output file written to:")
        logging.info(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )

        return newdf, hitdict

    elif mode == "single":
        logging.info("Running BLAST on single FASTA file")
        print("Running BLAST on single FASTA file")
        head_tail = os.path.split(folder)
        logging.info("Fasta file found:")
        logging.info(folder)
        logging.info(head_tail)
        duplist = fap.duplicate_fasta_check(folder)
        segdict = fap.missing_fasta_check(folder, segdict)
        segmissing, fasta_count = util.create_segment_tab(segdict)
        segdicttab = pd.DataFrame.from_dict(segmissing, orient="index")
        segdicttab = segdicttab.reset_index()
        segdicttab.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_segment_table.csv")
        )
        logging.info(
            f"Segment table written to : {os.path.join(output_dir,f'{now}_{args.tagname}_segment_table.csv')}"
        )
        segdicttab.columns = [
            "sample",
            "PB2",
            "PB1",
            "PA",
            "HA",
            "NP",
            "NA",
            "MP",
            "NS",
        ]
        if len(duplist) >= 1:
            print(f"Duplicates identified: {duplist}")
            logging.info(f"Duplicates identified: {duplist}")
        run_blast(
            os.path.join(repopath, args.blastdb), head_tail[0], head_tail[1], output_dir
        )
        logging.info("Reviewing BLAST results")
        print("Reviewing BLAST results")
        sresults, hitdict = tidy_blast_table(
            args, output_dir, head_tail[1], segdicttab, fasta_count
        )

        sresults.to_csv(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        logging.info("Output file written to:")
        logging.info(
            os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv")
        )
        print("Output file written to:")
        print(os.path.join(output_dir, f"{now}_{args.tagname}_BLAST_summary.csv"))
        return sresults, hitdict


def create_persample_summary(args, now, summarytab, segthreshold, hitdict):
    """

    Wrangle the results from the BLAST table (matched to genotypes) into a summary and likely top result.
    :return: results pivot table


    """
    print("Per sample summary step...")
    pivot = pd.pivot_table(
        summarytab,
        values="genotype_match",
        index="isolate_epi_id",
        columns="segment",
        fill_value="No sequence",
        aggfunc="first",
    )

    consensus = []
    freq = []
    results = []
    details = []
    genoblast = []

    for index, row in pivot.iterrows():
        process_pivot(
            segthreshold,
            hitdict,
            pivot,
            consensus,
            freq,
            results,
            details,
            genoblast,
            index,
            row,
        )
    pivot["consensus"] = consensus
    pivot["Top_Hit"] = freq
    pivot["Top_hit_details"] = details
    pivot["Best_result"] = genoblast
    pivot["Best_result"] = pivot["Best_result"].str.replace(
        "No sequence", "Insufficient sequence data"
    )

    pivot.to_csv(
        os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv")
    )
    logging.info("Summary output file written to:")
    logging.info(
        os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv")
    )
    print("Summary output file written to:")
    print(os.path.join(args.output_dir, f"{now}_{args.tagname}_genotyping_summary.csv"))
    return pivot


def process_pivot(
    segthreshold,
    hitdict,
    pivot,
    consensus,
    freq,
    results,
    details,
    genoblast,
    index,
    row,
):
    """

    Review pivot counts across geno groups
    :return: info to update pivot


    """
    result = ""
    newlist = []
    topgeno = []

    newlist = split_geno_matches(pivot, row, newlist)
    genfreq = Counter(newlist)
    counts = list(genfreq.values())
    maxgeno = max(counts)
    keys = list(genfreq.keys())
    tophits = count_geno_hits(segthreshold, topgeno, counts, maxgeno, keys)

    freq.append(tophits)
    consensus.append(genfreq)
    results.append(result)

    keylist = hitdict.keys()
    if index in keylist:
        topmatch, genodict = hitdict[index]
        gencounts = list(genodict.values())
        if len(gencounts) >= 1:
            maxgenohit = max(gencounts)
            genkeys = list(genodict.keys())
            if maxgenohit >= 6:
                for n, c in enumerate(gencounts):
                    if c == maxgenohit:
                        genoblast.append(genkeys[n])
            else:
                genoblast.append("No genotype had >=6 segments, review top hit")
        else:
            genoblast.append("")
        details.append(topmatch)
    else:
        print(index, " not in blast dictionary")
        details.append("No top hits found, BLAST fail??")
        genoblast.append("BLAST fail")


def count_geno_hits(segthreshold, topgeno, counts, maxgeno, keys):
    """

    Review genotype matches across the 8 segments, return the top results and where necessary choose a top hit out of multiple matches
    :return: top result per sample

    #try to simplify method....
    """
    for n, c in enumerate(counts):
        if c == maxgeno:
            topgeno.append(keys[n])
            if int(maxgeno) >= segthreshold:
                tophit = "|".join(topgeno)
            else:
                tophit = "No known genotype: Please review individual segments results"
    return tophit


def split_geno_matches(pivot, row, newlist):
    """

    Split the genotype matches into a list
    :return: update newlist list

    """
    for s in pivot.columns:
        test = row[s].split("|")
        for t in test:
            newlist.append(t)

    return newlist


def overall_summary(pivot):
    """
    Produce a summary of the per sample results to the command line and log

    :return: N/A
    """
    # Number of samples included
    # Summary of genotypes called
    print(f"Samples processed: {pivot.shape[0]}")
    print("Genotypes called: ")
    print(Counter(pivot["Top_Hit"]))
    logging.info(f"Samples processed: {pivot.shape[0]}")
    logging.info("Genotypes called: ")
    logging.info(Counter(pivot["Top_Hit"]))



def filter_blast_results(args,now,blasttab):
    """
    Tidy up the BLAST results table and output csv version of table

    :param blasttab: BLAST output table
    :return: tidied BLAST query table output and list of isolates
    """
    blast = pd.read_csv(blasttab,sep="\t")
    blast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    blast[['query_isolate_name', 'query_genotype','query_subtype','query_segment']] = blast['qseqid'].str.split('|', expand=True)
    blast[['hit_isolate_name', 'hit_genotype','hit_subtype','hit_segment']] = blast['sseqid'].str.split('|', expand=True)
    blast.to_csv(os.path.join(args.output_dir,f"blast_h5n1_geno_refs_{now}.csv"),index=False)
    queries = list(set(blast['query_isolate_name']))
    return blast,queries



def blast_processing(args,now):
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

