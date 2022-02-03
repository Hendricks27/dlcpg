"""
Code used to run the retention time predictor
"""

__author__ = "Robbin Bouwmeester"
__copyright__ = "Copyright 2019"
__credits__ = ["Robbin Bouwmeester","Prof. Lennart Martens","Sven Degroeve"]
__license__ = "Apache License, Version 2.0"
__version__ = "1.0"
__maintainer__ = "Robbin Bouwmeester"
__email__ = "Robbin.Bouwmeester@ugent.be"

from deeplc import DeepLC
from feat_extractor import FeatExtractor

# Native imports
import pickle
import sys
import os
import random
import itertools
import argparse
from collections import Counter
import re

# Pandas
import pandas as pd

# Matplotlib
from matplotlib import pyplot as plt

# Numpy
import numpy as np

def parse_arguments():
    """
    Read arguments from the command line

    Parameters
    ----------
        
    Returns
    -------
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--file_pred", type=str, dest="file_pred",default="",
                        help="Specify the file and path to make predictions for")

    parser.add_argument("--file_cal", type=str, dest="file_cal", default="",
                        help="Specify the file and path for calibrating the predictions (leave empty for no calibration)")

    parser.add_argument("--file_pred_out", type=str, dest="file_pred_out", default="",
                        help="Specify the outputfile for the (calibrated) predictions")

    parser.add_argument("--file_model",
                        help="Specify the model(s) to use to make the (calibrated) predictions seperate with spaces",
                        nargs="+",
                        default=["mods/full_dia_fixed_mods.hdf5","mods/full_integrated_dataset_v3.hdf5", "mods/full_seqs_21_ptm_alltype_fixed_mods.hdf5"])

    parser.add_argument("--n_threads", type=int, dest="n_threads", default=16,
                        help="Number of peaks to extract and consider for combinations in a spectrum")

    parser.add_argument("--split_cal", type=int, dest="split_cal", default=50,
                        help="Number of peaks to extract and consider for combinations in a spectrum")

    parser.add_argument("--dict_divider", type=int, dest="dict_divider", default=50,
                        help="Number of peaks to extract and consider for combinations in a spectrum")

    parser.add_argument("--batch_num", type=int, dest="batch_num", default=250000,
                        help="Batch size (of peptides) to use for predicting the retention time. Lower to decrease memory footprint")

    parser.add_argument("--version", action="version", version="%(prog)s 1.0")

    results = parser.parse_args()

    return results

def main():
    """
    Main function for the CLI interface
    
    Parameters
    ----------
        
    Returns
    -------
    """
    argu = parse_arguments()

    run(file_pred=argu.file_pred,
        file_cal=argu.file_cal,
        file_pred_out=argu.file_pred_out,
        file_model=argu.file_model,
        n_threads=argu.n_threads,
        split_cal=argu.split_cal,
        dict_divider=argu.dict_divider,
        batch_num=argu.batch_num)

def run(file_pred="",
        file_cal="",
        file_pred_out="",
        file_model="",
        n_threads=32,
        split_cal=50,
        dict_divider=50,
        batch_num=50000):
    """
    Main function to run the DeepLC code
    
    Parameters
    ----------
    file_pred : str
        the file in peprec format that we need to make predictions for
        this file is not required to contain a tr column
    file_cal : str
        the file in peprec format that we use for calibrating the prediction
        model. This file is required to contain a tr column
    file_pred_out : str
        outfile for predictions, the file is in peprec format and predictions
        are added in the column !!!
    file_model : str | list
        the model(s) to try for retention time prediction can be a single location
        or several locations for multiple models to try
    n_threads : int
        number of threads to run mainly the feature extraction on
    split_cal : int
        number of splits or divisions to use for the calibration
    dict_divider : int
        !!!
        
    Returns
    -------
    """
    df_pred = pd.read_csv(file_pred)
    df_pred = df_pred.fillna("")

    if len(file_cal) > 1:
        df_cal = pd.read_csv(file_cal)
        df_cal = df_cal.fillna("")

    # Make a feature extraction object; you can skip this if you do not want to use the default settings
    # for DeepLC. Here we want to use a model that does not use RDKit features so we skip the chemical
    # descriptor making procedure.
    f_extractor = FeatExtractor(add_sum_feat=False,
                                ptm_add_feat=False,
                                ptm_subtract_feat=False,
                                standard_feat = False,
                                chem_descr_feat = False,
                                add_comp_feat = False,
                                cnn_feats = True,
                                verbose = True)
    
    # Make the DeepLC object that will handle making predictions and calibration
    pepper = DeepLC(path_model=file_model,
                f_extractor=f_extractor,
                cnn_model=True,
                n_jobs=n_threads,
                verbose=False,
                batch_num=batch_num)

    # Calibrate the original model based on the new retention times
    if len(file_cal) > 1:
        pepper.calibrate_preds(seq_df=df_cal)
    
    # Make predictions; calibrated and uncalibrated
    if len(file_cal) > 1:
        preds = pepper.make_preds(seq_df=df_pred)
    else:
        preds = pepper.make_preds(seq_df=df_pred,calibrate=False)

    df_pred["Predicted tR"] = preds
    df_pred.to_csv(file_pred_out)

    if len(file_cal) > 1 and "tr" in df_pred.columns:
        plt.figure(figsize=(11.5,9))
        plt.scatter(df_pred["tr"],df_pred["Predicted tR"],s=3)
        plt.title("Predicted retention times")
        plt.xlabel("Observed tr")
        plt.ylabel("Predicted tr")
        plt.savefig("preds_%s.png" % (os.path.basename(file_pred).split(".")[0]), dpi=150)

if __name__ == "__main__":
    main()