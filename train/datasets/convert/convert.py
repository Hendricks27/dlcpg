import os
import sys

input_file = open("fraction_df_0196.csv")


init=True
lines = []
for l in input_file:
    l = l.strip()
    if init:
        init = False
        continue

    seq, modtmp, fr, x = l.split(",")
    # fr = float(fr)

    newmod = ""
    if modtmp == "Unmodified":
        pass
    elif modtmp == "Acetyl (Protein N-term)":
        newmod = "0|Acetyl|"
    elif "Oxidation" in modtmp:
        continue
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

    newmod = newmod[:-1]
    newl = ",".join([seq, newmod, fr])
    lines.append(newl)



lines.sort()

output_file = open("../fraction_df_0196_woo.csv", "w")
output_file.write("seq,modifications,tr\n")
for l in lines:
    output_file.write(l + "\n")


