import os
import random
import sys
import json
import sqlite3

input1 = open("human_library.csv")

rt_output = open("lib_rt.csv", "w")
fr_output = open("lib_fr.csv", "w")

i = 0
for l in input1:

    if i == 0:
        rtl = l.strip() + ",rt\n"
        frl = l.strip() + ",fraction\n"
    else:
        rtl = l.strip() + ",%0.2f\n" % (random.randint(15, 125) * 1.0)
        frl = l.strip() + ",%0.2f\n" % (random.randint(1, 970) / 10.)

    i += 1

    rt_output.write(rtl)
    fr_output.write(frl)



input1.close()
rt_output.close()
fr_output.close()



