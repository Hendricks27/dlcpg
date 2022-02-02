"""
Main code used to evalaute DeepLC on modifications
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
from collections import Counter

import itertools

import hashlib

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)

import tensorflow as tf
tf.__version__

def contains_mod(mod,contain_str):
    if contain_str in mod: return True
    else: return False
    
def get_native_version(index_name):
    return index_name.split("|")[0]+"|"

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

    def get_all_mods(mods):
        split_mods = mods.split("|")
        return [split_mods[m] for m in range(1,len(split_mods),2)]
        
    df = pd.read_csv("datasets/prosit_ptm_2020.csv")
    df.fillna("",inplace=True)

    all_pos_mods = list(itertools.chain(*list(df["modifications"].apply(get_all_mods))))
    all_pos_mods = Counter(all_pos_mods)

    n_epochs = 100
    batch_size = 256
    ratio_test = 0.9
    ratio_valid = 0.95
    sel_cols = ["seq","modifications","tr","index_name","predictions"]
    fit_hc = True
    use_correction_factor = True
    hc_str = "_hc"

    datasets = [
         "PXD008783_median_calibrate"
    ]
    sel_cols = ["seq","modifications","tr","index_name","num_train","dataset","iter","predictions"]
    sel_cols_native_mod_comparison = ["idents","seq","modifications","tr","idents","mod_id","mod_pred","mod_str","mod_tr"]

    for dataset_name in datasets:
        
        

        for m,num in all_pos_mods.most_common()[0:20]:
            print(m)
            
            aa_comp = cnn_functions.read_aa_lib("aa_comp_rel.csv")
            
            df = pd.read_csv("datasets/%s.csv" % (dataset_name))
            df.fillna("",inplace=True)
            df["idents"] = df["seq"]+"|"+df["modifications"]
            
            if use_correction_factor: correction_factor = (df["tr"].max())/10.0
            else: use_correction_factor = 1.0
            df["tr"] = df["tr"]/correction_factor
            
            df = cnn_functions.get_feat_df(df,aa_comp=aa_comp)
            
            mod_presence = df["modifications"].apply(contains_mod,args=(m,))
            df_test = df[mod_presence]
            df_train = df[~mod_presence]

            df_train,df_valid = cnn_functions.train_test(df_train,ratio_train=0.95)

            X_train, X_train_sum, X_train_global, X_train_hc, y_train = cnn_functions.get_feat_matrix(df_train)
            X_valid, X_valid_sum, X_valid_global, X_valid_hc, y_valid = cnn_functions.get_feat_matrix(df_valid)
            X_test, X_test_sum, X_test_global, X_test_hc, y_test = cnn_functions.get_feat_matrix(df_test)

            mods_optimized = []
        
            for p in params:
                a_blocks, a_kernel,a_max_pool,a_filters_start,a_stride,b_blocks,b_kernel,b_max_pool,b_filters_start,b_stride,global_neurons,global_num_dens,regularizer_val = p
                param_hash = hashlib.md5(",".join(map(str,p)).encode()).hexdigest()
                mod_name = "mods/pfind_%s_%s_%s.hdf5" % (dataset_name, m.replace("-","").replace(">",""),param_hash)

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

                mods_optimized.append(load_model(mod_name,
                                                custom_objects = {'<lambda>': lrelu}))

            df.index = df["idents"]

            m = m.replace("-","").replace(">","")

            cnn_functions.write_preds(df_train.loc[:,sel_cols], 
                                    X_train,
                                    X_train_sum,
                                    X_train_global, 
                                    X_train_hc, 
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_modanalysis_train.csv" % (dataset_name,m))

            cnn_functions.write_preds(df_valid.loc[:,sel_cols], 
                                    X_valid,
                                    X_valid_sum,
                                    X_valid_global,
                                    X_valid_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_modanalysis_valid.csv" % (dataset_name,m))
            
            cnn_functions.write_preds(df_test.loc[:,sel_cols], 
                                    X_test,
                                    X_test_sum,
                                    X_test_global,
                                    X_test_hc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_modanalysis_test.csv" % (dataset_name,m))

            cnn_functions.plot_preds(X_train,
                                    X_train_sum,
                                    X_train_global,
                                    X_train_hc,
                                    y_train,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    file_save="figures/%s_%s_modanalysis_train.png" % (dataset_name,m),
                                    plot_title="%s_%s_modanalysis_train" % (dataset_name,m))
            
            cnn_functions.plot_preds(X_valid,
                                    X_valid_sum,
                                    X_valid_global,
                                    X_valid_hc,
                                    y_valid,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    file_save="figures/%s_%s_modanalysis_valid.png" % (dataset_name,m),
                                    plot_title="%s_%s_modanalysis_valid" % (dataset_name,m))

            cnn_functions.plot_preds(X_test,
                                    X_test_sum,
                                    X_test_global,
                                    X_test_hc,
                                    y_test,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    file_save="figures/%s_%s_modanalysis_test.png" % (dataset_name,m),
                                    plot_title="%s_%s_modanalysis_test" % (dataset_name,m))
    
            # Ignore the mods for predictions below
            df_test.drop(["matrix","pos_matrix","matrix_sum","matrix_all","matrix_hc"],axis=1,inplace=True)
            

            df_test_glyc = cnn_functions.get_feat_df(df_test,aa_comp=aa_comp,ignore_mods=m.lower())
            X_test_glyc, X_test_sum_glyc, X_test_global_glyc, X_test_hc_glyc, y_test_glyc = cnn_functions.get_feat_matrix(df_test_glyc)
            
            cnn_functions.write_preds(df_test_glyc.loc[:,sel_cols], 
                                    X_test_glyc,
                                    X_test_sum_glyc,
                                    X_test_global_glyc,
                                    X_test_hc_glyc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    outfile_name="predictions/%s_%s_no_mod_modanalysis_train.csv" % (dataset_name,m))
            
            cnn_functions.plot_preds(X_test_glyc,
                                    X_test_sum_glyc,
                                    X_test_global_glyc,
                                    X_test_hc_glyc,
                                    y_test_glyc,
                                    mods_optimized,
                                    fit_hc=fit_hc,
                                    correction_factor=correction_factor,
                                    file_save="figures/%s_%s_no_mod_modanalysis_train.png" % (dataset_name,m),
                                    plot_title="%s_%s_no_mod_modanalysis_train" % (dataset_name,m))


if __name__ == "__main__":
    run()