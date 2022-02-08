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
30 Experiment 25
31 Experiment 26
32 Experiment 27
33 Experiment 28
34 Experiment 29
35 Experiment 3
36 Experiment 30
37 Experiment 31
38 Experiment 32
39 Experiment 33
40 Experiment 34
41 Experiment 35
42 Experiment 36
43 Experiment 37
44 Experiment 38
45 Experiment 39
46 Experiment 4
47 Experiment 40
48 Experiment 41
49 Experiment 42
50 Experiment 43
51 Experiment 44
52 Experiment 45
53 Experiment 46
54 Experiment 47
55 Experiment 48
56 Experiment 49
57 Experiment 5
58 Experiment 50
59 Experiment 51
60 Experiment 52
61 Experiment 53
62 Experiment 54
63 Experiment 55
64 Experiment 56
65 Experiment 57
66 Experiment 58
67 Experiment 59
68 Experiment 6
69 Experiment 60
70 Experiment 61
71 Experiment 62
72 Experiment 63
73 Experiment 64
74 Experiment 65
75 Experiment 66
76 Experiment 67
77 Experiment 68
78 Experiment 69
79 Experiment 7
80 Experiment 70
81 Experiment 71
82 Experiment 72
83 Experiment 73
84 Experiment 74
85 Experiment 75
86 Experiment 76
87 Experiment 77
88 Experiment 78
89 Experiment 79
90 Experiment 8
91 Experiment 80
92 Experiment 81
93 Experiment 82
94 Experiment 83
95 Experiment 84
96 Experiment 85
97 Experiment 86
98 Experiment 87
99 Experiment 88
100 Experiment 89
101 Experiment 9
102 Experiment 90
103 Experiment 91
104 Experiment 92
105 Experiment 93
106 Experiment 94
107 Experiment 95
108 Experiment 96
109 Retention time
110 Calibrated retention time
111 Charges
112 PEP
113 MS/MS scan number
114 Raw file
115 Score
116 Delta score
117 Reverse
118 Potential contaminant
119 Intensity
120 Intensity 1
121 Intensity 10
122 Intensity 11
123 Intensity 12
124 Intensity 13
125 Intensity 14
126 Intensity 15
127 Intensity 16
128 Intensity 17
129 Intensity 18
130 Intensity 19
131 Intensity 2
132 Intensity 20
133 Intensity 21
134 Intensity 22
135 Intensity 23
136 Intensity 24
137 Intensity 25
138 Intensity 26
139 Intensity 27
140 Intensity 28
141 Intensity 29
142 Intensity 3
143 Intensity 30
144 Intensity 31
145 Intensity 32
146 Intensity 33
147 Intensity 34
148 Intensity 35
149 Intensity 36
150 Intensity 37
151 Intensity 38
152 Intensity 39
153 Intensity 4
154 Intensity 40
155 Intensity 41
156 Intensity 42
157 Intensity 43
158 Intensity 44
159 Intensity 45
160 Intensity 46
161 Intensity 47
162 Intensity 48
163 Intensity 49
164 Intensity 5
165 Intensity 50
166 Intensity 51
167 Intensity 52
168 Intensity 53
169 Intensity 54
170 Intensity 55
171 Intensity 56
172 Intensity 57
173 Intensity 58
174 Intensity 59
175 Intensity 6
176 Intensity 60
177 Intensity 61
178 Intensity 62
179 Intensity 63
180 Intensity 64
181 Intensity 65
182 Intensity 66
183 Intensity 67
184 Intensity 68
185 Intensity 69
186 Intensity 7
187 Intensity 70
188 Intensity 71
189 Intensity 72
190 Intensity 73
191 Intensity 74
192 Intensity 75
193 Intensity 76
194 Intensity 77
195 Intensity 78
196 Intensity 79
197 Intensity 8
198 Intensity 80
199 Intensity 81
200 Intensity 82
201 Intensity 83
202 Intensity 84
203 Intensity 85
204 Intensity 86
205 Intensity 87
206 Intensity 88
207 Intensity 89
208 Intensity 9
209 Intensity 90
210 Intensity 91
211 Intensity 92
212 Intensity 93
213 Intensity 94
214 Intensity 95
215 Intensity 96
216 id
217 Protein group IDs
218 Peptide ID
219 Evidence IDs
220 MS/MS IDs
221 Best MS/MS
222 Oxidation (M) site IDs
223 MS/MS Count
224 Taxonomy IDs"""

outputxxxx = open("df0196.csv", "w")

def tocorrectorder(l):
    res = [None] * 96
    convert_index = [0, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 2, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 3, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 4, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 5, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 6, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 7, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 8, 89, 90, 91, 92, 93, 94, 95]
    for i, entry in enumerate(l):
        res[convert_index[i]] = entry
    assert None not in res
    return res

with open("./20211130_fraction1-96/txt/modificationSpecificPeptides.txt") as f:
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
                x += 1
                tmp.append(int(e[10:])-1)

            continue


        seq = l[0]
        mod = l[1]
        oxi = l[11]

        rt = float(l[109])

        pep = float(l[112])
        score = float(l[115])

        if pep > 0.01:
            continue
        if score < 90:
            continue

        frs = l[13:109]
        frsc = tocorrectorder(frs)
        frsd = {}
        for j, fr in enumerate(frsc):
            if fr == "":
                continue
            fr = int(fr)
            frsd[j] = fr

        intensity = tocorrectorder(l[120:215+1])
        intensityd = {}
        for j, ints in enumerate(intensity):
            if ints in "0":
                continue
            ints = int(ints)
            intensityd[j] = ints

        fraction = None
        if len(frsd) == 0 :
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

        line = [seq, mod, oxi, rt, fraction]

        outputxxxx.write("\t".join(map(str, line)) + "\n")













