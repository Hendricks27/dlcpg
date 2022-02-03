"""
Main code used to run AA evaluation
"""

__author__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__credits__ = ["Robbin Bouwmeester", "Ralf Gabriels", "Prof. Lennart Martens", "Sven Degroeve"]
__license__ = "Apache License, Version 2.0"
__maintainer__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__email__ = ["Robbin.Bouwmeester@ugent.be", "Ralf.Gabriels@ugent.be"]

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
from tensorflow.keras.callbacks import ModelCheckpoint

from sklearn.model_selection import KFold

import cnn_functions

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd

import itertools

import hashlib

# Keras
import tensorflow as tf

from tensorflow.keras.models import load_model

try: from tensorflow.keras.backend import set_session
except ImportError: from tensorflow.compat.v1.keras.backend import set_session
try: from tensorflow.keras.backend import clear_session
except ImportError: from tensorflow.compat.v1.keras.backend import clear_session
try: from tensorflow.keras.backend import get_session
except ImportError: from tensorflow.compat.v1.keras.backend import get_session

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

tf.__version__

import time


def run():
    a_blocks=[3]
    a_kernel=[2,4,8]
    a_max_pool=[2]
    a_filters_start=[256]
    a_stride=[1]
    b_blocks=[2]
    b_kernel=[2]
    b_max_pool=[2]
    b_filters_start=[128]
    b_stride=[1]
    global_neurons=[16]
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
    from tensorflow.keras.layers import ReLU
    import importlib

    lrelu = lambda x: tf.keras.activations.relu(x, alpha=0.1, max_value=20.0)
    importlib.reload(cnn_functions)

    n_epochs = 100
    batch_size = 256
    ratio_test = 0.9
    ratio_valid = 0.95
    sel_cols = ["seq","modifications","tr","index_name"]
    fit_hc = True
    use_correction_factor = True
    hc_str = "_hc"

    datasets = [
                "PXD005573_mcp",
                "mod_fixed_mods",
    ]

    aas = ['D','I','V','H','C','P','Q','S','N','Y','K','W', 'Y', 'E', 'A', 'F','L','M','R','T']

    for dataset_name in datasets:
        for aa in aas:
            df = cnn_functions.read_infile("datasets/%s.csv" % (dataset_name))

            aa_comp = cnn_functions.read_aa_lib("aa_comp_rel.csv",reset_to_glycine=aa)

            df = cnn_functions.get_feat_df(df,aa_comp=aa_comp)

            if use_correction_factor: correction_factor = (df["tr"].max())/10.0
            else: correction_factor = 1.0

            df_aa_count = cnn_functions.count_aa(df)
            df_aa_count.fillna(0.0,inplace=True)

            df_aa = df.loc[df_aa_count[df_aa_count[aa] > 0.0].index]
            df_train,df_test = cnn_functions.train_test(df.loc[df_aa_count[df_aa_count[aa] == 0.0].index],
                                                        seed=42,
                                                        ratio_train=0.9)

            df_train["aa_left_out"] = len(df_train.index)*[aa]
            df_train["dataset"] = len(df_train.index)*[dataset_name]

            df_test["aa_left_out"] = len(df_test.index)*[aa]
            df_test["dataset"] = len(df_test.index)*[dataset_name]

            df_aa["aa_left_out"] = len(df_aa.index)*[aa]
            df_aa["dataset"] = len(df_aa.index)*[dataset_name]
            df_aa["reset_to_glycine"] = len(df_aa.index)*["yes"]

            X_train, X_train_sum, X_train_global, X_train_hc, y_train = cnn_functions.get_feat_matrix(df_train)
            X_test, X_test_sum, X_test_global, X_test_hc, y_test = cnn_functions.get_feat_matrix(df_test)
            X_aa, X_aa_sum, X_aa_global, X_aa_hc, y_aa = cnn_functions.get_feat_matrix(df_aa)


            df_t = cnn_functions.read_infile("datasets/%s.csv" % (dataset_name))

            aa_comp_t = cnn_functions.read_aa_lib("aa_comp_rel.csv")
            df_t = cnn_functions.get_feat_df(df_t,aa_comp=aa_comp_t)
            df_aa_t = df_t.loc[df_aa_count[df_aa_count[aa] > 0.0].index]
            df_aa_t["aa_left_out"] = len(df_aa_t.index)*[aa]
            df_aa_t["dataset"] = len(df_aa_t.index)*[dataset_name]
            df_aa_t["reset_to_glycine"] = len(df_aa_t.index)*["no"]
            X_aa_t, X_aa_sum_t, X_aa_global_t, X_aa_hc_t, y_aa_t = cnn_functions.get_feat_matrix(df_aa_t)

            y_train = y_train/correction_factor
            y_test = y_test/correction_factor
            y_aa = y_aa/correction_factor
            y_aa_t = y_aa_t/correction_factor

            mods_optimized = []

            for p in params:
                a_blocks, a_kernel,a_max_pool,a_filters_start,a_stride,b_blocks,b_kernel,b_max_pool,b_filters_start,b_stride,global_neurons,global_num_dens,regularizer_val = p
                param_hash = hashlib.md5(",".join(map(str,p)).encode()).hexdigest()
                mod_name = "mods/%s_aa_cv%s_%s_train_%s.hdf5" % (dataset_name, hc_str, aa, param_hash)            

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
                
                time.sleep(2)

                mcp_save = ModelCheckpoint(mod_name,
                                        save_best_only=True,
                                        monitor='val_mean_absolute_error',
                                        mode='min')

                time.sleep(2)

                if fit_hc:
                    history = model.fit([X_train,X_train_sum,X_train_global,X_train_hc], y_train,
                                        validation_data=([X_test,X_test_sum,X_test_global,X_test_hc],y_test),
                                        epochs=n_epochs, 
                                        verbose=1, 
                                        batch_size=batch_size,
                                        callbacks=[mcp_save],
                                        shuffle=True)
                else:
                    history = model.fit([X_train,X_train_sum,X_train_global], y_train,
                                        validation_data=([X_test,X_test_sum,X_test_global],y_test),
                                        epochs=n_epochs, 
                                        verbose=1, 
                                        batch_size=batch_size,
                                        callbacks=[mcp_save],
                                        shuffle=True)

                with open("full_hc_%s_params.txt" % (dataset_name), "w") as file:
                    file.write("correction_factor,%s\n" % (correction_factor))
                    file.write("hc,%s\n" % (fit_hc))

                mods_optimized.append(load_model(mod_name,
                                                custom_objects = {'<lambda>': lrelu}))
            
            cnn_functions.write_preds(df_train.loc[:,sel_cols], 
                                X_train,
                                X_train_sum,
                                X_train_global,
                                X_train_hc,
                                mods_optimized,
                                fit_hc=fit_hc,
                                correction_factor=correction_factor,
                                outfile_name="predictions/%s_%s_aa_cv%s_train.csv" % (dataset_name,aa,hc_str))

            cnn_functions.write_preds(df_test.loc[:,sel_cols], 
                                    X_test,
                                    X_test_sum,
                                    X_test_global,
                                    X_test_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_aa_cv%s_valid.csv" % (dataset_name,aa,hc_str))

            cnn_functions.write_preds(df_aa.loc[:,sel_cols], 
                                    X_aa,
                                    X_aa_sum,
                                    X_aa_global,
                                    X_aa_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_aa_cv%s_test_glycine.csv" % (dataset_name,aa,hc_str))

            cnn_functions.plot_preds(X_train, 
                    X_train_sum,
                    X_train_global,
                    X_train_hc,
                    y_train,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_%s_aa_cv%s_train.png" % (dataset_name,aa,hc_str),
                    plot_title="%s_%s_aa_cv%s_train" % (dataset_name,aa,hc_str))

            cnn_functions.plot_preds(X_test, 
                    X_test_sum,
                    X_test_global,
                    X_test_hc,
                    y_test,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_%s_aa_cv%s_valid.png" % (dataset_name,aa,hc_str),
                    plot_title="%s_%s_aa_cv%s_valid" % (dataset_name,aa,hc_str))

            cnn_functions.plot_preds(X_aa, 
                    X_aa_sum,
                    X_aa_global,
                    X_aa_hc,
                    y_aa,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_%s_aa_cv%s_test_glycine.png" % (dataset_name,aa,hc_str),
                    plot_title="%s_%s_aa_cv%s_test_glycine" % (dataset_name,aa,hc_str))

            cnn_functions.write_preds(df_aa_t.loc[:,sel_cols], 
                                    X_aa_t,
                                    X_aa_sum_t,
                                    X_aa_global_t, 
                                    X_aa_hc_t,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_aa_cv%s_test.csv" % (dataset_name,aa,hc_str))

            cnn_functions.plot_preds(X_aa_t,
                    X_aa_sum_t,
                    X_aa_global_t,
                    X_aa_hc_t,
                    y_aa_t,
                    mods_optimized,
                    fit_hc=fit_hc,
                    correction_factor=correction_factor,
                    file_save="figures/%s_%s_aa_cv%s_test.png" % (dataset_name,aa,hc_str),
                    plot_title="%s_%s_aa_cv%s_test" % (dataset_name,aa,hc_str))
 	
            time.sleep(5)

if __name__ == "__main__":
    run()