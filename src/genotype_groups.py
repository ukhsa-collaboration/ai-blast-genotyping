import src.utilities as util 
import src.fasta_processing as fap
import src.blast_wrangle as bl
import pandas as pd
import os
import logging
import subprocess
import json
import glob

def update_geno_key(args,now,genokey, subepitab):
    """
    Update the genotyping key table with new genotypes

    :param genokey: existing genotyping key table 
    :param subepitab: new rows for genotyping key table


    :return: combined genotyping key table
    """
    genotab = pd.read_csv(genokey)
    newdf = pd.concat([genotab, subepitab])
    head_tail = os.path.split(genokey)
    newname = head_tail[1].replace(".csv",f"_{now}.csv")
    newdf.to_csv(os.path.join(args.output_dir,f'{newname}'),index=False)
    return newdf


def check_overlap(q,grouplist,groupdict,ingroup):
    """
    Review groups per segment and merge highly overlapping groups
    
    :param q: query list to check overlap
    :param grouplist: query relatives list
    :param groupdict: existing group dictionary
    :param ingroup: decision variable (merge / separate)
    :return: groupdict, groupdecision, ingroup

    """
    maxoverlap = 0
    groupfound = ""
    for key, val in groupdict.items():
    
        if val is None:
            pass
        else:
            overlap = util.intersection(list(set(grouplist)), val)   
            if len(overlap) > maxoverlap:
                maxoverlap = len(overlap)
                groupfound = key

    groupdict, groupdecision, ingroup = decide_group_merge(q, maxoverlap, groupfound, grouplist, groupdict, ingroup)
    return groupdict, groupdecision, ingroup

def decide_group_merge(q,maxoverlap,groupfound, grouplist, groupdict, ingroup):
    """
    Review groups per segment and merge highly overlapping groups
    :param q: query list to check overlap
    :param maxoverlap: max overlap with any existing groups
    :param groupfound: query found in group list
    :param grouplist: query relatives list
    :param groupdict: existing group dictionary
    :param ingroup: decision variable (merge / separate)

    :return: groupdict, groupdecision, ingroup
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
    
def create_group_json(args,queries, t, subblast, s,subtype):
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
            if subtype is True:
                
                h5type = hitblast['query_subtype'].iloc[0]
                h5type = h5type.replace("A/", "")
                hitblast['hit_subtype'] = hitblast['hit_subtype'].str.replace("A/", "")
                hitblast = hitblast[hitblast['hit_subtype'] == h5type]
            
            ingroup = False
                #check if query is already in a group in dictionary 
            for key, val in groupdict.items():
                    if val is None:
                        pass
                    else:
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
    """
    Create new group

    :param s: segment
    :param groupdict: group dictionary
    :param groups: running group count
    :param grouplist: list of relatives for group

    :return: groupdict, groups, ingroup (T/F variable)
    """
    newgrouplist = []
    for g in grouplist:
        ng = False
        for key, val in groupdict.items():
            if g in list(val):
                ng = True
                exit   
        if not ng:
             newgrouplist.append(g)

    logging.info(f"New group defined: {s}_group{groups}, {list(set(newgrouplist))}")
    groupdict[f'{s}_group{groups}'] = list(set(newgrouplist))
    ingroup = True
    groups = groups+1
    return groupdict, groups, ingroup

def create_constellation(args,now,queries, blast, tabcols, t):
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
        for s in util.segments:
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
    for s in util.segments:
        thresholddf[s] = thresholddf[s].fillna("Absent") 
    thresholddf["constellation"] = thresholddf[['PB2','PB1','PA','HA','NP','NA','MP','NS']].agg("|".join, axis=1)
    thresholddf.to_csv(os.path.join(args.output_dir,f'blast_geno_threshold_table{t}_{now}.csv'),index=False)
    return thresholddf



def describe_geno(args,genotype,constellation,thresholddf,segments):
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


def report_new_genotypes(now,thresholddf, genolist):
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
        describe_geno(row['genotype'],row['constellation'],thresholddf,util.segments)
    dup_check = thresholddf[thresholddf['constellation'].isin(constellations)]
    if dup_check.shape[0] > subthresholds.shape[0]:
        print("Duplicate constellation found for at least one genotype....")
        thresholddf.to_csv(f'new_genotype_duplicate_constellations_tab_{now}.csv',index=False)
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
    tabcols.extend(util.segments)
    # we want the different combinations / constellations in a table


    subblast = blast[blast['pident']>=t]

    for s in util.segments:
        if s == "NA":
            create_group_json(queries, t, subblast, s,subtype=True)
        else:
            create_group_json(queries, t, subblast, s,subtype=False)
    thresholddf = create_constellation(queries, blast, tabcols, t)
    report_new_genotypes(thresholddf,newgeno)

def files_for_gitrepo(now,output_dir):
    try:
        subprocess.call(f"mkdir {os.path.join(output_dir,'for_git')}",shell=True)
    except:
        pass
    constpath = os.path.join(output_dir,f"blast_geno_threshold_table98_{now}.csv")
    newconstpath = os.path.join(output_dir,"for_git","blast_geno_threshold_table98.csv")
    subprocess.call(f"cp {constpath} {newconstpath}",shell=True)
    
    dbfiles = glob.glob(os.path.join(output_dir,'all_genotype_references_*'))
    for d in dbfiles:
        if "db" in d:
            newdbloc = d.replace(f"_{now}","")
            newdbloc = newdbloc.replace(output_dir,os.path.join(output_dir,"for_git"))
            subprocess.call(f'cp {d} {newdbloc}',shell=True)