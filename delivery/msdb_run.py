
import os
import re
import sys
import json
import time
import random
import sqlite3

import masstool

from msdb import RTFR



closeT = 5


# loc = ":memory:"
loc = "test2.db"
rtfr = RTFR(loc)

if False:
    rtfr.add_csv("pred.csv")
    rtfr.create_index()
    rtfr.create_index_mw()

protein_count = {}
frag2prid = {}
if os.path.exists("frag2prid.json"):
    frag2prid = json.load(open("frag2prid.json"))
    for uniprotids in frag2prid.values():
        for prid in uniprotids:
            protein_count[prid] = 0
else:
    f = open("fragments_lookup.tsv")
    for l in f:
        prid, fr = l.strip().split("\t")
        if fr not in frag2prid:
            frag2prid[fr] = set()
        frag2prid[fr].add(prid)
    f.close()
    for k, v in frag2prid.items():
        frag2prid[k] = list(v)
    json.dump(frag2prid, open("frag2prid.json", "w"))



"""
0 Raw file
1 Scan number
2 Retention time
3 Ion injection time
4 Total ion current
5 Collision energy
6 Summations
7 Base peak intensity
8 Elapsed time
9 Identified
10 Matched
11 Reverse
12 MS/MS IDs
13 Sequence
14 Length
15 Filtered peaks
16 m/z
17 Mass
18 Charge
19 Type
20 Fragmentation
21 Mass analyzer
22 Parent intensity fraction
23 Fraction of total spectrum
24 Base peak fraction
25 Precursor full scan number
26 Precursor intensity
27 Precursor apex fraction
28 Precursor apex offset
29 Precursor apex offset time
30 Scan event number
31 Modifications
32 Modified sequence
33 Proteins
34 Score
35 PEP
36 Experiment
37 Reporter PIF
38 Reporter fraction
39 Intens Comp Factor
40 CTCD Comp
41 RawOvFtT
42 AGC Fill
43 Scan index
44 MS scan index
45 MS scan number
"""
x=True

last_fr = 0
last_sn = 0

missedcount = 0
savedcount1 = 0
savedcount2 = 0
unindentifiedcount = 0
ms2counts = 0


logf = open("ided_missed.txt", "w")
def missed(detail):
    tmp = list(map(str, detail))
    logf.write("\t".join(tmp) + "\n")
    return

percent_tmp = [x/100 for x in range(100)]


start_ts = time.time()
query_time = 0
query_count = 0
for l in open("msmsScans96.txt"):
    l = l.strip().split("\t")
    if x:
        x = False
        continue

    ms2counts += 1
    running_percent = ms2counts / 3793412


    if len(percent_tmp) > 0:
        if abs(running_percent-percent_tmp[0]) < pow(0.1, 6):
            log = "%0.2f\t%0.2f (%0.2f+%0.2f)\t%0.4f\t%0.2fs" % (
                running_percent * 100,
                (savedcount1 + savedcount2) / ms2counts * 100,
                savedcount1 / ms2counts * 100,
                savedcount2 / ms2counts * 100,
                missedcount / ms2counts * 100,
                time.time() - start_ts
            )
            print(log)
            percent_tmp.pop(0)



    f_num = int(re.compile(r"Fraction(\d{1,2})").findall(l[0])[0])
    scan_num = int(l[1])
    rt = float(l[2])
    seq = l[13].strip()

    if l[17] == "NaN":
        # No charge
        continue
    mw = float(l[17])


    # Make sure it is in the right order print(f_num, scan_num)
    assert f_num >= last_fr
    if f_num > last_fr:
        last_sn = 0
    last_fr = f_num
    assert scan_num >= last_sn
    last_sn = scan_num

    mods = l[31].strip()
    protein = l[33].strip()

    mz = float(l[16])
    charge = int(l[18])
    assert charge >= 0
    if charge != 0:
        diff = mz*charge - charge - float(mw)
        if diff > 1:
            print(mz, charge, mw, diff)

    assert l[9] in "+-"
    id_flag = l[9] == "+"

    if seq != "":
        doublecheckmw = masstool.peptide_mw(seq) #+ mod_mass["Carbamidomethyl"]
        if "Oxidation" in mods:
            continue
        assert mods in ["Acetyl (Protein N-term)", "Unmodified", ""]
        for aa in seq:
            if aa == "K":
                doublecheckmw += masstool.mod_mass["TMT"]
            if aa == "C":
                doublecheckmw += masstool.mod_mass["Carbamidomethyl"]
        if mods == "Acetyl (Protein N-term)":
            doublecheckmw += masstool.mod_mass["Acetyl"]
        else:
            doublecheckmw += masstool.mod_mass["TMT"]

        diff = abs(mw - doublecheckmw)
        if diff > 1:
            # Issue
            print("ISSUE")

    query_start_ts = time.time()
    #dbres = list(rtfr.query(mw, rt, f_num))
    dbres = list(rtfr.query_mass(mw, rt, f_num))
    query_time += time.time()-query_start_ts
    query_count += 1
    if len(dbres) == 0:
        """
        if id_flag:
            missedcount += 1
            missed([seq, mods, mw, charge, mz, rt, f_num, protein])
        """
        savedcount1 += 1
        if len(frag2prid.get(seq, [])) > 0:
            missedcount += 1
            missed([seq, mods, mw, charge, mz, rt, f_num, protein])

        continue

    closeout = True
    allpossibleprids = set()
    for r0 in dbres:
        prids = frag2prid[r0[0]]
        allpossibleprids = allpossibleprids.union(prids)
        for prid in prids:
            if protein_count[prid] < closeT:
                closeout = False
                break

    if closeout:
        savedcount2 += 1
    else:
        if seq in frag2prid:
            for prid in frag2prid[seq]:
                protein_count[prid] += 1
        else:
            unindentifiedcount += 1

    #uniprot_ids = list(filter(lambda x: not x.startswith("REV_") and not x.startswith("CON_") and x!="", protein.split(";")))
    #if set(uniprot_ids) != set(frag2prid.get(seq, [])):
    #    print(seq, uniprot_ids, frag2prid.get(seq, []))

    """
    uniprot_ids = list(filter(lambda x: x != "", protein.split(";")))
    for pr in uniprot_ids:
        if pr.startswith("REV__"):
            pr = pr[5:]

        if pr not in allpossibleprids:
            pass
            #print(pr, allpossibleprids, rev, seq)
        protein_count[prid] += 1
    """



    '''
    closedout = True
    uniprot_ids = list(filter(lambda x:x != "", protein.split(";")))
    for pr in uniprot_ids:
        if pr not in protein_count:
            print(pr)
            protein_count[pr] = 0
        if protein_count[pr] < closeT:
            closedout = False
    if len(uniprot_ids) == 0:
        unindentifiedcount += 1
        closedout = False

    if closedout:
        savedcount2 += 1
        continue

    for pr in uniprot_ids:
        protein_count[pr] += 1
    '''


total_time_took = time.time() - start_ts


stat = f"""
Total MS2 scans: {ms2counts}
{savedcount1} {round(100 * savedcount1 / ms2counts, 2)}% could be saved due to wrong mw-rt-fr
{savedcount2} {round(100 * savedcount2 / ms2counts, 2)}% could be saved due to closed out
{missedcount} {round(100 * missedcount / ms2counts, 2)}% are missed due to incorrect rt/fr
{unindentifiedcount} {round(100 * unindentifiedcount / ms2counts, 2)}% found peptide sequences are not in any proteins 
Took {round(total_time_took)} seconds
Took {round(query_time)} seconds on database query, average for {round(query_time/query_count*1000, 2)} ms each query
"""

print(stat)



