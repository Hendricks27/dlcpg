import os
import sys

dname = ["f96", "f25"][1]
input_file = open("%s.csv" % dname)

lines_rt = []
lines_fr = []
for l in input_file:

    seq, modtmp, mw, oxinum, rt, fr = l.strip().split("\t")
    # print(seq, modtmp, oxinum, rt, fr)

    newmod = ""
    oxi = False
    if modtmp == "Unmodified":
        pass
    elif modtmp == "Acetyl (Protein N-term)":
        newmod = "0|Acetyl|"
    elif "Oxidation" in modtmp:
        if oxinum != "1":
            continue
        if seq.count("M") > 1:
            continue
        assert modtmp in ["Oxidation (M)", "Acetyl (Protein N-term);Oxidation (M)"]
        if modtmp == "Acetyl (Protein N-term);Oxidation (M)":
            newmod = "0|Acetyl|"
        oxi = True
    else:
        print(modtmp)

    for i, aa in enumerate(seq):
        pos = i+1

        if i == 0 and newmod != "0|Acetyl|":
            newmod = "0|TMT|"

        if aa == "K":
            newmod += "%s|TMT|" % pos

        if aa == "C":
            newmod += "%s|Carbamidomethyl|" % pos

        if aa == "M" and oxi:
            newmod += "%s|Oxidation|" % pos

    newmod = newmod[:-1]

    lines_rt.append(",".join([seq, newmod, rt]))
    lines_fr.append(",".join([seq, newmod, fr]))





lines_rt.sort()
lines_fr.sort()

output_file = open("../%s_fraction.csv" % dname, "w")
output_file.write("seq,modifications,tr\n")
for l in lines_fr:
    output_file.write(l + "\n")

output_file = open("../%s_retention.csv" % dname, "w")
output_file.write("seq,modifications,tr\n")
for l in lines_rt:
    output_file.write(l + "\n")
