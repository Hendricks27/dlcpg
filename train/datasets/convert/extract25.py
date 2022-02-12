import os
import sys


"""
0 Sequence
1 Modifications
2 Mass
3 Mass Fractional Part
4 Protein Groups
5 Proteins
6 Gene Names
7 Protein Names
8 Unique (Groups)
9 Unique (Proteins)
10 Acetyl (Protein N-term)
11 Oxidation (M)
12 Missed cleavages
13 Experiment 1
14 Experiment 10
15 Experiment 11
16 Experiment 12
17 Experiment 13
18 Experiment 14
19 Experiment 15
20 Experiment 16
21 Experiment 17
22 Experiment 18
23 Experiment 19
24 Experiment 2
25 Experiment 20
26 Experiment 21
27 Experiment 22
28 Experiment 23
29 Experiment 24
30 Experiment 3
31 Experiment 4
32 Experiment 5
33 Experiment 6
34 Experiment 7
35 Experiment 8
36 Experiment 9
37 Experiment A
38 Retention time
39 Calibrated retention time
40 Charges
41 PEP
42 MS/MS scan number
43 Raw file
44 Score
45 Delta score
46 Reverse
47 Potential contaminant
48 Intensity
49 Intensity 1
50 Intensity 10
51 Intensity 11
52 Intensity 12
53 Intensity 13
54 Intensity 14
55 Intensity 15
56 Intensity 16
57 Intensity 17
58 Intensity 18
59 Intensity 19
60 Intensity 2
61 Intensity 20
62 Intensity 21
63 Intensity 22
64 Intensity 23
65 Intensity 24
66 Intensity 3
67 Intensity 4
68 Intensity 5
69 Intensity 6
70 Intensity 7
71 Intensity 8
72 Intensity 9
73 Intensity A
74 id
75 Protein group IDs
76 Peptide ID
77 Evidence IDs
78 MS/MS IDs
79 Best MS/MS
80 Oxidation (M) site IDs
81 MS/MS Count
82 Taxonomy IDs
"""



input_file_path = "./20211203-combined25/txt/modificationSpecificPeptides.txt"
out_file_path = "dfc25.csv"


outputxxxx = open(out_file_path, "w")

def tocorrectorder(l):
    res = [None] * 25
    convert_index = [0, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 19, 20, 21, 22, 23, 2, 3, 4, 5, 6, 7, 8, 24]
    for i, entry in enumerate(l):
        res[convert_index[i]] = entry
    assert None not in res
    return res

with open(input_file_path) as f:
    i = 0
    for l in f:
        i += 1
        l = l.split("\t")
        if i == 1:
            x = 0
            tmp = []
            for e in l:

                if "Experiment" not in e:
                    continue

                if e == "Experiment A":
                    e = "Experiment 25"

                x += 1
                tmp.append(int(e[10:])-1)

            # convert_index in tocorrectorder
            # print(tmp)

            continue

        seq = l[0]
        mod = l[1]
        oxi = l[11]
        mw = l[2]

        rt = float(l[38])

        pep = float(l[41])
        score = float(l[44])

        if pep > 0.01:
            continue
        if score < 90:
            continue

        frs = l[13:38]
        frsc = tocorrectorder(frs)
        frsd = {}
        for j, fr in enumerate(frsc):
            if fr == "":
                continue
            fr = int(fr)
            frsd[j] = fr

        intensity = tocorrectorder(l[49:74])
        intensityd = {}
        for j, ints in enumerate(intensity):
            if ints in "0":
                continue
            ints = int(ints)
            intensityd[j] = ints

        fraction = None
        if len(frsd) == 0:
            raise RuntimeError
        elif len(frsd) == 1:
            fraction = list(frsd.keys())[0]
        else:
            keys = set(frsd.keys()).union(set(intensityd.keys()))
            assert set(frsd.keys()).issuperset(intensityd.keys())

            wfr1 = 0
            s1 = sum(frsd.values())
            for k in keys:
                wfr1 += frsd[k] / s1 * k
            fraction = wfr1

            wfr2 = 0
            s2 = sum(intensityd.values())
            if s2 == 0:
                wfr2 = wfr1
            else:
                for k in keys:
                    wfr2 += intensityd.get(k, 0) / s2 * k

            fraction = wfr2

        line = [seq, mod, mw, oxi, rt, fraction]

        outputxxxx.write("\t".join(map(str, line)) + "\n")













