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
import json


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


def run():

    a_blocks=[2, 3]
    a_kernel=[2, 4, 8]
    a_max_pool=[2]
    a_filters_start=[128, 256, 1024]
    a_stride=[1]
    b_blocks=[2, 3]
    b_kernel=[2, 4, 8]
    b_max_pool=[2]
    b_filters_start=[128, 256, 1024]
    b_stride=[1]
    global_neurons=[8, 16, 32]
    global_num_dens=[3]
    regularizer_val=[0.0000025]

    params = list(itertools.product(*[a_blocks,
                                        a_kernel,
                                        a_max_pool,
                                        a_filters_start,
                                        a_stride,
                                        b_blocks,
                                        b_kernel,
                                        b_max_pool,
                                        b_filters_start,
                                        b_stride,
                                        global_neurons,
                                        global_num_dens,
                                        regularizer_val]))

    """
    Use this cell to train setup specific models ready for production environments

        PRODUCTION MODEL FIT

    """

    from tensorflow.keras.layers import ReLU
    import importlib

    lrelu = lambda x: tf.keras.activations.relu(x, alpha=0.1, max_value=20.0)
    importlib.reload(cnn_functions)

    n_epochs = 100
    batch_size = 256
    ratio_test = 0.9
    ratio_valid = 0.95
    sel_cols = ["seq", "modifications", "tr", "index_name", "predictions"]
    fit_hc = True
    use_correction_factor = True
    hc_str = "_hc"

    datasets = [
        "fraction_df_0196_woo",
        # "tinytest"
    ]
    for dataset_name in datasets:
        df = cnn_functions.read_infile("datasets/%s.csv" % (dataset_name))
        aa_comp = cnn_functions.read_aa_lib("aa_comp_rel.csv")

        df = cnn_functions.get_feat_df(df,aa_comp=aa_comp)

        if use_correction_factor:
            correction_factor = (df["tr"].max())/10.0 #*0.05
        else:
            correction_factor = 1.0

        # TODO K-FORD cross validation

        df_train,df_test = cnn_functions.train_test(df,ratio_train=ratio_test)
        df_train,df_valid = cnn_functions.train_test(df_train,ratio_train=ratio_valid)
        #df_valid = df_test

        X_train, X_train_sum, X_train_global, X_train_hc, y_train = cnn_functions.get_feat_matrix(df_train)
        X_valid, X_valid_sum, X_valid_global, X_valid_hc, y_valid = cnn_functions.get_feat_matrix(df_valid)
        X_test, X_test_sum, X_test_global, X_test_hc, y_test = cnn_functions.get_feat_matrix(df_test)

        y_train = y_train/correction_factor
        y_valid = y_valid/correction_factor
        y_test = y_test/correction_factor

        mods_optimized = []

        for p in params:
            a_blocks, a_kernel, a_max_pool, a_filters_start, a_stride, b_blocks, b_kernel, b_max_pool, b_filters_start, b_stride, global_neurons, global_num_dens, regularizer_val = p
            param_hash = hashlib.md5(",".join(map(str,p)).encode()).hexdigest()

            param_history_path = "./mod_para/%s.json" % param_hash
            if os.path.exists(param_history_path):
                continue

            json.dump(p, open(param_history_path, "w"))

            mod_name = "mods/full%s_%s_%s.hdf5" % (hc_str,dataset_name,param_hash)

            model = cnn_functions.init_model(X_train, X_train_sum, X_train_global,X_test_hc,
                                            a_blocks=a_blocks,
                                            a_kernel=a_kernel,
                                            a_max_pool=a_max_pool,
                                            a_filters_start=a_filters_start,
                                            a_stride=a_stride,
                                            b_blocks=b_blocks,
                                            b_kernel=b_kernel,
                                            b_max_pool=b_max_pool,
                                            b_filters_start=b_filters_start,
                                            b_stride=b_stride,
                                            global_neurons=global_neurons,
                                            global_num_dens=global_num_dens,
                                            regularizer_val=regularizer_val,
                                            fit_hc=fit_hc)

            mcp_save = ModelCheckpoint(mod_name,
                                    save_best_only=True,
                                    monitor='val_mean_absolute_error',
                                    mode='min')
            if fit_hc:
                history = model.fit([X_train,X_train_sum,X_train_global,X_train_hc], y_train,
                                    validation_data=([X_valid,X_valid_sum,X_valid_global,X_valid_hc],y_valid),
                                    epochs=n_epochs,
                                    verbose=1,
                                    batch_size=batch_size,
                                    callbacks=[mcp_save],
                                    shuffle=True)
            else:
                history = model.fit([X_train,X_train_sum,X_train_global], y_train,
                                    validation_data=([X_valid,X_valid_sum,X_valid_global],y_valid),
                                    epochs=n_epochs,
                                    verbose=1,
                                    batch_size=batch_size,
                                    callbacks=[mcp_save],
                                    shuffle=True)

            with open("full_hc_%s_params.txt" % (dataset_name), "w") as file:
                file.write("correction_factor,%s\n" % (correction_factor))
                file.write("hc,%s\n" % (fit_hc))

            mods_optimized.append(load_model(mod_name, custom_objects={'<lambda>': lrelu}))


            # Evaluate

            cnn_functions.write_preds(df_train.loc[:, sel_cols],
                                        X_train,
                                        X_train_sum,
                                        X_train_global,
                                        X_train_hc,
                                        mods_optimized,
                                        fit_hc=fit_hc,
                                        correction_factor=correction_factor,
                                        outfile_name="predictions/%s_full%s_train_%s.csv" % (dataset_name, hc_str, param_hash))

            cnn_functions.write_preds(df_valid.loc[:, sel_cols],
                                    X_valid,
                                    X_valid_sum,
                                    X_valid_global,
                                    X_valid_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_full%s_valid_%s.csv" % (dataset_name, hc_str, param_hash))

            cnn_functions.write_preds(df_test.loc[:, sel_cols],
                                    X_test,
                                    X_test_sum,
                                    X_test_global,
                                    X_test_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_full%s_test_%s.csv" % (dataset_name, hc_str, param_hash))




            cnn_functions.plot_preds(X_test,
                    X_test_sum,
                    X_test_global,
                    X_test_hc,
                    y_test,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_full_hc_test_%s.png" % (dataset_name, param_hash),
                    plot_title="%s_full%s_test" % (dataset_name, hc_str))

            cnn_functions.plot_preds(X_train,
                    X_train_sum,
                    X_train_global,
                    X_train_hc,
                    y_train,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_full_hc_train_%s.png" % (dataset_name, param_hash),
                    plot_title="%s_full%s_train" % (dataset_name, hc_str))

            cnn_functions.plot_preds(X_valid,
                    X_valid_sum,
                    X_valid_global,
                    X_valid_hc,
                    y_valid,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_full_hc_valid_%s.png" % (dataset_name, param_hash),
                    plot_title="%s_full%s_valid" % (dataset_name, hc_str))



if __name__ == "__main__":

    start_ts = time.time()

    run()
    duration = time.time() - start_ts
    print("Time took: " + str(duration))
    print()


