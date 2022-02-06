import os
import sys
import json
import time
import sqlite3



class RTFR():

    allowed_table = [
        "retention_time",
        "fraction"
    ]

    def __init__(self, dbp):
        self.connection = None
        self.cursor = None
        self.connect(dbp)

    def connect(self, dbp):
        self.connection = sqlite3.connect(dbp)
        self.cursor = self.connection.cursor()

    def initialize(self):
        self.cursor.execute(
            '''
            CREATE TABLE retention_time
            (seq text, mods text, rt real)
            ''')

        self.cursor.execute(
            '''
            CREATE TABLE fraction
            (seq text, mods text, fr real)
            ''')
        self.connection.commit()


    def add_csv(self, table, file_path):
        assert table in self.allowed_table
        i = 0
        csv = open(file_path)
        for l in csv:
            if i == 0:
                i += 1
                continue
            seq, mod, d = l.strip().split(",")
            self.add_data_single(table, (seq, mod, float(d)))
        self.commit()
        csv.close()
        return


    def add_data_single(self, table, triplet):
        assert table in self.allowed_table
        self.cursor.execute("INSERT INTO %s VALUES %s" % (table, triplet))

    def commit(self):
        self.connection.commit()

    def query(self, table, drange, limit=500):
        assert table in self.allowed_table

        data_type = None
        if table == "retention_time":
            data_type = "rt"
        elif table == "fraction":
            data_type = "fr"
        else:
            raise RuntimeError

        q = """
        SELECT * FROM %s  
        WHERE %s BETWEEN %s AND %s
        ORDER BY seq
        LIMIT %s;
        """ % (table, data_type, drange[0], drange[1], limit)
        for row in self.cursor.execute(q):
            yield row


loc = ":memory:" # test.db
rtfr = RTFR(loc)


if loc == ":memory:":
    rtfr.initialize()
    rtfr.add_csv("retention_time", "lib_rt.csv")
    rtfr.add_csv("fraction", "lib_fr.csv")


#a = list(rtfr.query("fraction", (50, 62), limit=100000))
# b = list(rtfr.query("retention_time", (3000, 3500), limit=100000))


qtest = """
SELECT fraction.seq, fraction.mods, fraction.fr
FROM fraction
WHERE fraction.fr BETWEEN 50 AND 62
"""

qtest = """
SELECT retention_time.seq, retention_time.mods, retention_time.rt
FROM retention_time
WHERE retention_time.rt BETWEEN 1000 AND 1300
"""

qtest = """
SELECT fraction.seq, fraction.mods, fraction.fr, retention_time.rt
FROM fraction
INNER JOIN retention_time ON fraction.seq = retention_time.seq AND fraction.mods = retention_time.mods
WHERE fraction.fr BETWEEN 50 AND 62
AND retention_time.rt BETWEEN 70 AND 87
"""


ts = time.time()
i = 0
alllines = 3717909
for row in rtfr.cursor.execute(qtest):
    i += 1
    # print(row)

duration = time.time() - ts

print(i, alllines, "%0.2f" % (i/alllines*100), duration)

time.sleep(20)