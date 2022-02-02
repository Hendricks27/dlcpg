"""
Main code used to evaluate DeepLC performance
"""

__author__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__credits__ = ["Robbin Bouwmeester", "Ralf Gabriels", "Prof. Lennart Martens", "Sven Degroeve"]
__license__ = "Apache License, Version 2.0"
__maintainer__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__email__ = ["Robbin.Bouwmeester@ugent.be", "Ralf.Gabriels@ugent.be"]

import os
import sys
import time


import hashlib
import itertools
import matplotlib
matplotlib.use('Agg')

import cnn_functions

import numpy as np
import pandas as pd
import tensorflow as tf

from matplotlib import pyplot as plt

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model

from sklearn.model_selection import KFold


config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)


pd.set_option("display.max_rows", None, "display.max_columns", None)


dataset_name = "tinytestmod"



if __name__ == "__main__":
    df = cnn_functions.read_infile("datasets/%s.csv" % (dataset_name))
    aa_comp = cnn_functions.read_aa_lib("aa_comp_rel.csv")
    print(df)

    df = cnn_functions.get_feat_df(df, aa_comp=aa_comp)

    matrix_index =["matrix", "matrix_sum", "pos_matrix", "matrix_all", "matrix_hc"]

    for pep in ["Pep_0", "Pep_1"]:
        print(pep)
        for m in matrix_index:
            print(m)
            print(df.loc[pep, m])
            print("\n-------------------------\n")
        print("\n=========================\n")








