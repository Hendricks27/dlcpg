import os
import sys
import json


f = open("human_library.csv")
f2 = open("human_library_w_mass.csv").read().strip().split("\n")
rt = json.load(open("pred_f96_retention_rt.json"))
fr = json.load(open("pred_f96_fraction_fr.json"))

out = open("../../delivery/pred.csv", "w")

i = 0
for l in f:
    line = l.strip().split(",")
    s, m, mw = f2[i].split(",")
    assert line[0] == s
    assert line[1] == m

    line += [mw, rt[i], fr[i]]
    out.write("\t".join(map(str, line)) + "\n")

    i += 1




