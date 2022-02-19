
import os
import sys
import json
import time
import random
import sqlite3

from msdb import RTFR


# loc = ":memory:"
loc = "test2.db"
rtfr = RTFR(loc)

if False:
    rtfr.add_csv("pred.csv")
    rtfr.create_index()



for i in range(1000):
    ts = time.time()

    alllines = 3717909

    rt = rtfr.random_real(20, 115)
    mw = rtfr.random_real(400, 4200)
    fr = rtfr.random_real(3, 96)
    res = list(rtfr.query(mw, rt, fr))

    duration = time.time() - ts

    line = "\t".join(map(str, [len(res), alllines, "%0.3f P" % (len(res)/alllines*100), "%0.2f ms" % (duration*1000)]))

    print(line)




