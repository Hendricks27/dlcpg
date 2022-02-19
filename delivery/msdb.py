import os
import sys
import json
import time
import random
import sqlite3



class RTFR():

    def __init__(self, dbp):
        # Confidence interval 99%
        # mw range 10ppm
        self._mw_range = 4600 * 10 / pow(10, 6)
        self._rt_range = 11
        self._fr_range = 5

        self.connection = None
        self.cursor = None
        self.connect(dbp)

    def connect(self, dbp):
        init = os.path.exists(dbp)
        self.connection = sqlite3.connect(dbp)
        self.cursor = self.connection.cursor()

        if not init:
            self.cursor.execute('CREATE TABLE info (seq text, mods text, mw real, rt real, fr real)')
            self.commit()

    def random_real(self, a, b):
        return random.randint(a*1000, b*1000) / 1000

    def add_csv(self, file_path):
        i = 0
        csv = open(file_path)
        for l in csv:
            seq, mod, mw, rt, fr = l.strip().split("\t")
            self.add_data_single((seq, mod, float(mw), float(rt), float(fr)))
        self.commit()
        csv.close()
        return

    def create_index(self):
        self.cursor.execute("CREATE INDEX ind ON info(mw, rt, fr)")
        self.commit()

    def add_data_single(self, d):
        self.cursor.execute("INSERT INTO info VALUES %s" % (d, ))

    def commit(self):
        self.connection.commit()

    def set_mw_range(self, d):
        self._mw_range = d

    def set_rt_range(self, d):
        self._rt_range = d

    def set_fr_range(self, d):
        self._fr_range = d

    def query(self, mw, rt, fr):
        mw_start = mw - self._mw_range
        mw_end   = mw + self._mw_range

        rt_start = rt - self._rt_range
        rt_end   = rt + self._rt_range

        fr_start = fr - self._fr_range
        fr_end   = fr + self._fr_range

        q = """
        SELECT seq, mods, mw, rt, fr
        FROM info
        WHERE mw BETWEEN %0.2f AND %0.2f
        AND rt BETWEEN %0.2f AND %0.2f
        AND fr BETWEEN %0.2f AND %0.2f
        """ % (mw_start, mw_end, rt_start, rt_end, fr_start, fr_end)

        for row in self.cursor.execute(q):
            yield row


if __name__ == "__main__":
    print("hello")



