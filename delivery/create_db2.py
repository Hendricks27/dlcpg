import os
import random
import sys
import json
import time
import sqlite3



class RTFR():


    def __init__(self, dbp):
        self.connection = None
        self.cursor = None
        self.connect(dbp)

    def connect(self, dbp):
        self.connection = sqlite3.connect(dbp)
        self.cursor = self.connection.cursor()

    def initialize(self):
        self.cursor.execute('''CREATE TABLE info (seq text, mods text, mw real, rt real, fr real)''')
        self.connection.commit()

    def random_real(self, a, b):
        return random.randint(a*1000, b*1000) / 1000


    def add_csv(self, file_path):
        i = 0
        csv = open(file_path)
        for l in csv:
            if i == 0:
                i += 1
                continue
            seq, mod, d = l.strip().split(",")
            self.add_data_single((seq, mod, float(d), self.random_real(15, 125), self.random_real(2, 97)))
        self.commit()
        csv.close()
        return


    def add_data_single(self, d):
        self.cursor.execute("INSERT INTO info VALUES %s" % (d, ))

    def commit(self):
        self.connection.commit()




# loc = ":memory:"
loc = "test2.db"
rtfr = RTFR(loc)


if loc == ":memory:":
    rtfr.initialize()
    rtfr.add_csv("human_library.csv")



#a = list(rtfr.query("fraction", (50, 62), limit=100000))
# b = list(rtfr.query("retention_time", (3000, 3500), limit=100000))


qtest = """
SELECT seq, mods, mw, rt, fr
FROM info
WHERE mw BETWEEN 2000 AND 2005
AND rt BETWEEN 40 AND 58
AND fr BETWEEN 50 AND 62
"""




ts = time.time()
i = 0
alllines = 3717909
for row in rtfr.cursor.execute(qtest):
    i += 1
    # print(row)

duration = time.time() - ts

print(i, alllines, "%0.5f" % (i/alllines*100), duration)


