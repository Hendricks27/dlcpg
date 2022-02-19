
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

issue = 0
savedcount1 = 0
savedcount2 = 0
ms2counts = 0

protein_count = {}
for l in open("msmsScans96.txt"):
    l = l.strip().split("\t")
    if x:
        x = False
        continue

    ms2counts += 1
    f_num = int(re.compile(r"Fraction(\d{1,2})").findall(l[0])[0])
    scan_num = int(l[1])
    rt = float(l[2])
    seq = l[13].strip()
    if l[17] == "NaN":
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
            print(seq, mods, mw, doublecheckmw, diff, "%0.2f" % (diff/57.021463735))

    dbres = list(rtfr.query(mw, rt, f_num))
    if len(dbres) == 0:
        if id_flag:
            print(seq, mods, protein)
        savedcount1 += 1


    closedout = True
    uniprot_ids = list(filter(lambda x:x != "", protein.split(";")))
    for pr in uniprot_ids:
        if pr not in protein_count.items():
            protein_count[pr] = 0
        if protein_count[pr] < closeT:
            closedout = False

    if closedout:
        savedcount2 += 1
    else:
        for pr in uniprot_ids:
            protein_count[pr] += 1







print(savedcount2, ms2counts, "%0.2f" % (100 * savedcount2 / ms2counts))



