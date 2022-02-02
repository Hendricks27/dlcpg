"""
Main code used to generate learning curves
"""

__author__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__credits__ = ["Robbin Bouwmeester", "Ralf Gabriels", "Prof. Lennart Martens", "Sven Degroeve"]
__license__ = "Apache License, Version 2.0"
__maintainer__ = ["Robbin Bouwmeester", "Ralf Gabriels"]
__email__ = ["Robbin.Bouwmeester@ugent.be", "Ralf.Gabriels@ugent.be"]

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import load_model

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

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

import tensorflow as tf
tf.__version__


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
    sel_cols = ["seq","modifications","tr","index_name","predictions"]
    fit_hc = True
    use_correction_factor = True
    hc_str = "_hc"

    sel_cols = ["seq","modifications","tr","index_name","num_train","dataset","iter","predictions"]

    div_size = 11
    num_iters = 5
    n_epochs = 100
    batch_size = 256
    fit_hc = True
    use_correction_factor = True
    hc_str = "_hc"

    datasets = [
                "dia_fixed_mods",
                "PXD005573_mcp",
                "hela_hf_psms_aligned",
                "mod_fixed_mods",
                "LUNA_HILIC_fixed_mods",
                "unmod_fixed_mods",
                "ATLANTIS_SILICA_fixed_mods",
                "LUNA_SILICA_fixed_mods",
                "SCX_fixed_mods",
                "prosit_ptm_2020",
                "Xbridge_fixed_mods",
                "PXD008783_median_calibrate",
                "hela_lumos_2h_psms_aligned",
                "pancreas_psms_aligned",
                "yeast_120min_psms_aligned",
                "hela_lumos_1h_psms_aligned",
                "yeast_60min_psms_aligned",
                "arabidopsis_psms_aligned",
                "plasma_lumos_2h_psms_aligned",
                "plasma_lumos_1h_psms_aligned",
                "proteometools_library"
                "PXD005573_mcp"
    ]
    for dataset_name in datasets:
        df = cnn_functions.read_infile("datasets/%s.csv" % (dataset_name))
        aa_comp = cnn_functions.read_aa_lib("aa_comp_rel.csv")
        df = cnn_functions.get_feat_df(df,aa_comp=aa_comp)
        
        if use_correction_factor: correction_factor = (df["tr"].max())/10.0
        else: correction_factor = 1.0
            
        df["tr"] = df["tr"]/correction_factor

        for i in range(num_iters):
            df_train,df_test_temp = cnn_functions.train_test(df,seed=i,ratio_train=0.9)
            df_test,df_test_final = cnn_functions.train_test(df_test_temp,seed=i,ratio_train=0.5)
            model = False
            for max_val in np.linspace(0,len(df_train.index),div_size):
                max_val = int(max_val)
                if max_val == 0: continue

                df_train["num_train"] = len(df_train.index)*[max_val]
                df_train["dataset"] = len(df_train.index)*[dataset_name]
                df_train["iter"] = len(df_train.index)*[i]
                df_train_iter = df_train.loc[df_train.index[0:max_val]]

                df_test["num_train"] = len(df_test.index)*[max_val]
                df_test["dataset"] = len(df_test.index)*[dataset_name]
                df_test["iter"] = len(df_test.index)*[i]

                df_test_final["num_train"] = len(df_test_final.index)*[max_val]
                df_test_final["dataset"] = len(df_test_final.index)*[dataset_name]
                df_test_final["iter"] = len(df_test_final.index)*[i]
                
                X_train, X_train_sum, X_train_global, X_train_hc, y_train = cnn_functions.get_feat_matrix(df_train_iter)
                X_test, X_test_sum, X_test_global, X_test_hc, y_test = cnn_functions.get_feat_matrix(df_test)
                X_test_final, X_test_sum_final, X_test_global_final, X_test_hc_final, y_test_final = cnn_functions.get_feat_matrix(df_test_final)

                mods_optimized = []
            
                for p in params:
                    a_blocks, a_kernel,a_max_pool,a_filters_start,a_stride,b_blocks,b_kernel,b_max_pool,b_filters_start,b_stride,global_neurons,global_num_dens,regularizer_val = p
                    param_hash = hashlib.md5(",".join(map(str,p)).encode()).hexdigest()
                    mod_name = "mods/%s_learning_curve_iter_%s_train_%s_%s.hdf5" % (dataset_name, i, max_val, param_hash)

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
                                            validation_data=([X_test, X_test_sum, X_test_global, X_test_hc],y_test),
                                            epochs=n_epochs, 
                                            verbose=1, 
                                            batch_size=batch_size,
                                            callbacks=[mcp_save],
                                            shuffle=True)
                    else:
                        history = model.fit([X_train,X_train_sum,X_train_global], y_train, 
                                            validation_data=([X_test, X_test_sum, X_test_global],y_test),
                                            epochs=n_epochs, 
                                            verbose=1, 
                                            batch_size=batch_size,
                                            callbacks=[mcp_save],
                                            shuffle=True)

                    with open("full_hc_%s_params.txt" % (dataset_name), "w") as file:
                        file.write("correction_factor,%s\n" % (correction_factor))
                        file.write("hc,%s\n" % (fit_hc))

                    mods_optimized.append(load_model(mod_name,
                                                    custom_objects = {'<lambda>': lrelu},
                                                    compile=True))

                cnn_functions.write_preds(df_train_iter.loc[:,sel_cols], 
                                        X_train,
                                        X_train_sum,
                                        X_train_global,
                                        X_train_hc,
                                        mods_optimized,
                                        fit_hc=fit_hc,
                                        correction_factor=correction_factor,
                                        outfile_name="predictions/%s_learning_curve_train_iter_%s_train_%s.csv" % (dataset_name,i,max_val))

                cnn_functions.write_preds(df_test.loc[:,sel_cols], 
                                        X_test, X_test_sum, X_test_global, X_test_hc,
                                        mods_optimized,
                                        fit_hc=fit_hc,
                                        correction_factor=correction_factor,
                                        outfile_name="predictions/%s_learning_curve_valid_iter_%s_train_%s.csv" % (dataset_name,i,max_val))
                                        
                perf = cnn_functions.write_preds(df_test_final.loc[:,sel_cols], 
                                        X_test_final, X_test_sum_final, X_test_global_final, X_test_hc_final,
                                        mods_optimized,
                                        fit_hc=fit_hc,
                                        correction_factor=correction_factor,
                                        outfile_name="predictions/%s_learning_curve_test_iter_%s_train_%s.csv" % (dataset_name,i,max_val))
                
                cnn_functions.plot_preds(X_train,
                        X_train_sum,
                        X_train_global,
                        X_train_hc,
                        y_train,
                        mods_optimized,
                        fit_hc=fit_hc,
                        correction_factor=correction_factor,
                        file_save="figures/%s_learning_curve_train_iter_%s_%s.png" % (dataset_name,i,max_val),
                        plot_title="%s_learning_curve_train_iter_%s_train_%s" % (dataset_name,i,max_val))

                cnn_functions.plot_preds(X_test, X_test_sum, X_test_global, X_test_hc,
                        y_test,
                        mods_optimized,
                        fit_hc=fit_hc,
                        correction_factor=correction_factor,
                        file_save="figures/%s_learning_curve_valid_iter_%s_%s.png" % (dataset_name,i,max_val),
                        plot_title="%s_learning_curve_valid_iter_%s_train_%s" % (dataset_name,i,max_val))


                cnn_functions.plot_preds(
                        X_test_final, X_test_sum_final, X_test_global_final, X_test_hc_final,
                        y_test_final,
                        mods_optimized,
                        fit_hc=fit_hc,
                        correction_factor=correction_factor,
                        file_save="figures/%s_learning_curve_test_iter_%s_%s.png" % (dataset_name,i,max_val),
                        plot_title="%s_learning_curve_test_iter_%s_train_%s" % (dataset_name,i,max_val))

if __name__ == "__main__":
    run()