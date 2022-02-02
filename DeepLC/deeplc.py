"""
Main code used to generate LC retention time predictions. This provides the main interface.

For the library versions see the .yml file
"""

__author__ = "Robbin Bouwmeester"
__copyright__ = "Copyright 2019"
__credits__ = ["Robbin Bouwmeester","Prof. Lennart Martens","Sven Degroeve"]
__license__ = "Apache License, Version 2.0"
__version__ = "1.0"
__maintainer__ = "Robbin Bouwmeester"
__email__ = "Robbin.Bouwmeester@ugent.be"

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_DIR)

# Native imports
import time
import pickle
from operator import itemgetter
import sys
from configparser import ConfigParser
import time
import gc

# Pandas
import pandas as pd

# Numpy
import numpy as np

# XGBoost
import xgboost as xgb

# Keras
import tensorflow as tf

from tensorflow.keras.models import load_model

try: from tensorflow.keras.backend import set_session
except ImportError: from tensorflow.compat.v1.keras.backend import set_session
try: from tensorflow.keras.backend import clear_session
except ImportError: from tensorflow.compat.v1.keras.backend import clear_session
try: from tensorflow.keras.backend import get_session
except ImportError: from tensorflow.compat.v1.keras.backend import get_session


# Set to force CPU calculations
#os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

# Set for TF V1.0 (counters some memory problems of nvidia 20 series GPUs)
#config = tf.ConfigProto()
#config.gpu_options.allow_growth = True
#session = tf.Session(config=config)

# Set for TF V2.0 (counters some memory problems of nvidia 20 series GPUs)
#config = tf.compat.v1.ConfigProto()
#config.gpu_options.allow_growth = True
#session = tf.compat.v1.Session(config=config)

# Feature extraction
from feat_extractor import FeatExtractor

# Multiproc
from multiprocessing import Pool

import copy

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

# Reset Keras Session
def reset_keras():
    sess = get_session()
    clear_session()
    sess.close()
    sess = get_session()
    gc.collect()

class DeepLC():
    """
    Place holder, fill later

    Parameters
    ----------

    Returns
    -------

    """
    def __init__(self,
                 main_path=os.path.dirname(os.path.realpath(__file__)),
                 path_model=None,
                 verbose=True,
                 bin_dist=1,
                 dict_cal_divider = 100,
                 split_cal = 25,
                 n_jobs=32,
                 config_file=None,
                 f_extractor=None,
                 cnn_model=False,
                 batch_num=350000):        
        # if a config file is defined overwrite standard parameters
        if config_file:
            cparser = ConfigParser()
            cparser.read(config_file)
            dict_cal_divider = cparser.getint("lcPep","dict_cal_divider")
            split_cal = cparser.getint("lcPep","split_cal")
            n_jobs = cparser.getint("lcPep","n_jobs")

        self.main_path = main_path
        self.verbose = verbose
        self.bin_dist = bin_dist
        self.calibrate_dict = {}
        self.calibrate_min = float("inf")
        self.calibrate_max = 0
        self.cnn_model = cnn_model

        self.batch_num = batch_num
        self.dict_cal_divider = dict_cal_divider
        self.split_cal = split_cal
        self.n_jobs = n_jobs

        if path_model:
            if self.cnn_model:
                self.model = path_model
            else:
                with open(path_model, "rb") as handle:
                    self.model = pickle.load(handle)
        
        if f_extractor:
            self.f_extractor = f_extractor
        else:
            self.f_extractor = FeatExtractor()
    
    def __str__(self):
        return("""
  _____                  _      _____ 
 |  __ \                | |    / ____|
 | |  | | ___  ___ _ __ | |   | |     
 | |  | |/ _ \/ _ \ '_ \| |   | |     
 | |__| |  __/  __/ |_) | |___| |____ 
 |_____/ \___|\___| .__/|______\_____|
                  | |                 
                  |_|                   
              """)
        

    def do_f_extraction(self,
                        seqs,
                        mods,
                        identifiers):
        """
        Extract all features we can extract; without parallelization; use if you want to run feature extraction
        with a single core

        Parameters
        ----------
        seqs : list
            peptide sequence list; should correspond to mods and identifiers
        mods : list
            naming of the mods; should correspond to seqs and identifiers
        identifiers : list
            identifiers of the peptides; should correspond to seqs and mods

        Returns
        -------
        pd.DataFrame
            feature matrix
        """
        return self.f_extractor.full_feat_extract(seqs,mods,identifiers)

    def do_f_extraction_pd(self,
                        df_instances):
        """
        Extract all features we can extract; without parallelization; use if you want to run feature extraction
        with a single thread; and use a defined dataframe

        Parameters
        ----------
        df_instances : object :: pd.DataFrame
            dataframe containing the sequences (column:seq), modifications (column:modifications) and naming (column:index)

        Returns
        -------
        pd.DataFrame
            feature matrix
        """
        return self.f_extractor.full_feat_extract(df_instances["seq"],df_instances["modifications"],df_instances.index)
    
    def do_f_extraction_pd_parallel(self,
                        df_instances):
        """
        Extract all features we can extract; with parallelization; use if you want to run feature extraction
        with multiple threads; and use a defined dataframe

        Parameters
        ----------
        df_instances : object :: pd.DataFrame
            dataframe containing the sequences (column:seq), modifications (column:modifications) and naming (column:index)
        
        Returns
        -------
        pd.DataFrame
            feature matrix
        """
        df_instances_split = np.array_split(df_instances, self.n_jobs)
        pool = Pool(self.n_jobs)
        if self.n_jobs == 1: df = self.do_f_extraction_pd(df_instances)
        else: df = pd.concat(pool.map(self.do_f_extraction_pd, df_instances_split))
        pool.close()
        pool.join()
        return df

    def make_preds_core(self,
                seqs=[],
                mods=[],
                identifiers=[],
                calibrate=True,
                seq_df=None,
                correction_factor=1.0,
                mod_name=False):
        """
        Make predictions for sequences

        Parameters
        ----------
        seq_df : object :: pd.DataFrame
            dataframe containing the sequences (column:seq), modifications (column:modifications) and naming (column:index);
            will use parallel by default!
        seqs : list
            peptide sequence list; should correspond to mods and identifiers
        mods : list
            naming of the mods; should correspond to seqs and identifiers
        identifiers : list
            identifiers of the peptides; should correspond to seqs and mods
        calibrate : boolean
            calibrate predictions or just return the predictions

        Returns
        -------
        np.array
            predictions
        """

        # See if we got a list; if not assume we got a df
        if len(seqs) == 0:
            # Make a copy, because we do not want to change to original df
            seq_df = seq_df.copy()
        else:
            # Make a df out of provided lists
            seq_df = pd.DataFrame([seqs,mods]).T
            seq_df.columns = ["seq","modifications"]
            seq_df.index = identifiers
        
        # Only run on unique peptides, defined by seq+mod
        # TODO sort the mods in the peprec on both position and alphabet mod; to not let duplicates through!
        seq_df["idents"] = seq_df["seq"]+"|"+seq_df["modifications"]
        identifiers = list(seq_df.index)

        # Save a row identifier to seq+mod mapper so output has expected return shapes
        identifiers_to_seqmod = dict(zip(seq_df.index,seq_df["idents"]))
        
        # Drop duplicated seq+mod
        seq_df.drop_duplicates(subset=["idents"],inplace=True)            

        if self.verbose:
            cnn_verbose = 1
        else:
            cnn_verbose = 0

        # If we need to apply deep NN
        if self.cnn_model:
            if self.verbose: print("Extracting features for the CNN model ...")
            X = self.do_f_extraction_pd_parallel(seq_df)
            X = X.loc[seq_df.index]
            
            X_sum = np.stack(X["matrix_sum"])
            X_global = np.concatenate((np.stack(X["matrix_all"]),
                                    np.stack(X["pos_matrix"])),
                                    axis=1)

            X = np.stack(X["matrix"])
        else:
            if self.verbose: print("Extracting features for the predictive model ...")
            seq_df.index
            X = self.do_f_extraction_pd_parallel(seq_df)
            X = X.loc[seq_df.index]

            X = X[self.model.feature_names]
        
        ret_preds = []

        # If we need to calibrate
        if calibrate:
            if self.verbose: print("Predicting with calibration ...")

            cal_preds = []

            # Load the model differently if we are going to use a CNN
            if self.cnn_model:
                if mod_name == False:
                    mod = load_model(self.model)
                else:
                    mod = load_model(mod_name)
                uncal_preds = mod.predict([X,X_sum,X_global],batch_size=5120).flatten()/correction_factor
            else:
                # first get uncalibrated prediction
                uncal_preds = self.model.predict(X)/correction_factor

            for uncal_pred in uncal_preds:
                try:
                    slope,intercept,x_correction = self.calibrate_dict[str(round(uncal_pred,self.bin_dist))]
                    cal_preds.append(slope * (uncal_pred-x_correction) + intercept)
                except KeyError:
                    # outside of the prediction range ... use the last calibration curve
                    if uncal_pred <= self.calibrate_min:
                        slope,intercept,x_correction = self.calibrate_dict[str(round(self.calibrate_min,self.bin_dist))]
                        cal_preds.append(slope * (uncal_pred-x_correction) + intercept)
                    elif uncal_pred >= self.calibrate_max:
                        slope,intercept,x_correction = self.calibrate_dict[str(round(self.calibrate_max,self.bin_dist))]
                        cal_preds.append(slope * (uncal_pred-x_correction) + intercept)
                    else:
                        slope,intercept,x_correction = self.calibrate_dict[str(round(self.calibrate_max,self.bin_dist))]
                        cal_preds.append(slope * (uncal_pred-x_correction) + intercept)
            ret_preds = np.array(cal_preds)
        else:
            if self.verbose: print("Predicting values ...")

            # Load the model different if we use CNN
            if self.cnn_model:
                if mod_name == False:
                    mod = load_model(self.model)
                else:
                    mod = load_model(mod_name)
                ret_preds = mod.predict([X,X_sum,X_global],batch_size=5120,verbose=cnn_verbose).flatten()/correction_factor
            else:
                ret_preds = self.model.predict(X)/correction_factor

        pred_dict = dict(zip(seq_df["idents"],ret_preds))

        # Map from unique peptide identifiers to the original dataframe
        ret_preds_shape = []
        for ident in identifiers:
            ret_preds_shape.append(pred_dict[identifiers_to_seqmod[ident]])

        if self.verbose: print("Predictions done ...")

        # Below can cause freezing on some systems
        # It is meant to clear any remaining vars in memory
        reset_keras()
        del mod

        return ret_preds_shape

    def make_preds(self,
            seqs=[],
            mods=[],
            identifiers=[],
            calibrate=True,
            seq_df=None,
            correction_factor=1.0,
            mod_name=False):
        if self.batch_num == 0:
            return self.make_preds_core(seqs=seqs,
                                        mods=mods,
                                        identifiers=identifiers,
                                        calibrate=calibrate,
                                        seq_df=seq_df,
                                        correction_factor=correction_factor,
                                        mod_name=mod_name)
        else:
            ret_preds = []
            if len(seqs) > 0:
                seq_df = pd.DataFrame({"seq": seqs,
                                       "modifications": mods},
                                       index=identifiers)
            for g, seq_df_t in seq_df.groupby(np.arange(len(seq_df)) // self.batch_num):
                temp_preds = self.make_preds_core(identifiers=identifiers,
                                                  calibrate=calibrate,
                                                  seq_df=seq_df_t,
                                                  correction_factor=correction_factor,
                                                  mod_name=mod_name)
                ret_preds.extend(temp_preds)

                #if self.verbose:
                print("Finished predicting retention time for: %s/%s" % (len(ret_preds),len(seq_df)))
            return ret_preds


    def calibrate_preds_func(self,
                             seqs=[],
                             mods=[],
                             identifiers=[],
                             measured_tr=[],
                             correction_factor=1.0,
                             seq_df=None,
                             use_median=True,
                             mod_name=False):
        """
        Make calibration curve for predictions

        Parameters
        ----------
        seqs : list
            peptide sequence list; should correspond to mods and identifiers
        mods : list
            naming of the mods; should correspond to seqs and identifiers
        identifiers : list
            identifiers of the peptides; should correspond to seqs and mods
        measured_tr : list
            measured tr of the peptides; should correspond to seqs, identifiers, and mods
        correction_factor : float
            correction factor that needs to be applied to the supplied measured trs
        seq_df : object :: pd.DataFrame
            a pd.DataFrame that contains the sequences, modifications and observed
            retention times to fit a calibration curve
        use_median : boolean
            flag to indicate we need to use the median valuein a window to perform calibration
        mod_name
            specify a model to use instead of the model assigned originally to this instance
            of the object

        Returns
        -------
        float
            the minimum value where a calibration curve was fitted, lower values
            will be extrapolated from the minimum fit of the calibration curve
        float
            the maximum value where a calibration curve was fitted, higher values
            will be extrapolated from the maximum fit of the calibration curve
        dict
            dictionary with keys for rounded tr, and the values concern a linear model
            that should be applied to do calibration (!!! what is the shape of this?)        
        """
        if len(seqs) == 0:
            seq_df.index
            predicted_tr = self.make_preds(seq_df=seq_df,calibrate=False,correction_factor=correction_factor,mod_name=mod_name)
            measured_tr = seq_df["tr"]
        else:
            predicted_tr = self.make_preds(seqs=seqs,mods=mods,identifiers=identifiers,calibrate=False,correction_factor=correction_factor,mod_name=mod_name)
        
        # sort two lists, predicted and observed based on measured tr
        tr_sort = [(mtr,ptr) for mtr,ptr in sorted(zip(measured_tr,predicted_tr), key=lambda pair: pair[0])]
        measured_tr = [mtr for mtr,ptr in tr_sort]
        predicted_tr = [ptr for mtr,ptr in tr_sort]

        mtr_mean = []
        ptr_mean = []

        calibrate_dict = {}
        calibrate_min = float('inf')
        calibrate_max = 0

        if self.verbose: print("Selecting the data points for calibration (used to fit the linear models between)")

        # smooth between observed and predicted
        for mtr,ptr in zip(self.split_seq(measured_tr,self.split_cal),self.split_seq(predicted_tr,self.split_cal)):
            if use_median:
                mtr_mean.append(np.median(mtr))
                ptr_mean.append(np.median(ptr))
            else:
                mtr_mean.append(sum(mtr)/len(mtr))
                ptr_mean.append(sum(ptr)/len(ptr))

        if self.verbose: print("Fitting the linear models between the points")

        # calculate calibration curves
        for i in range(0,len(ptr_mean)):
            if i >= len(ptr_mean)-1: continue
            delta_ptr = ptr_mean[i+1]-ptr_mean[i]
            delta_mtr = mtr_mean[i+1]-mtr_mean[i]

            slope = delta_mtr/delta_ptr
            intercept = mtr_mean[i]
            x_correction = ptr_mean[i]

            # optimized predictions using a dict to find calibration curve very fast
            for v in np.arange(round(ptr_mean[i],self.bin_dist),round(ptr_mean[i+1],self.bin_dist),1/((self.bin_dist)*self.dict_cal_divider)):
                if v < calibrate_min:
                    calibrate_min = v
                if v > calibrate_max:
                    calibrate_max = v
                calibrate_dict[str(round(v,1))] = [slope,intercept,x_correction]

        if self.verbose: print("Time to calibrate: %s seconds" % (time.time() - t0))

        return calibrate_min, calibrate_max, calibrate_dict

    def calibrate_preds(self,
                        seqs=[],
                        mods=[],
                        identifiers=[],
                        measured_tr=[],
                        correction_factor=1.0,
                        seq_df=None,
                        use_median=True):
        """
        Make calibration curve for predictions

        Parameters
        ----------
        seqs : list
            peptide sequence list; should correspond to mods and identifiers
        mods : list
            naming of the mods; should correspond to seqs and identifiers
        identifiers : list
            identifiers of the peptides; should correspond to seqs and mods
        measured_tr : list
            measured tr of the peptides; should correspond to seqs, identifiers, and mods
        correction_factor : float
            correction factor that needs to be applied to the supplied measured trs
        seq_df : object :: pd.DataFrame
            a pd.DataFrame that contains the sequences, modifications and observed
            retention times to fit a calibration curve
        use_median : boolean
            flag to indicate we need to use the median valuein a window to perform calibration

        Returns
        -------

        """
        if type(self.model) == str:
            self.model = [self.model]
        
        if self.verbose: print("Start to calibrate predictions ...")
        if self.verbose: print("Ready to find the best model out of: %s" % (self.model))

        best_perf = float("inf")
        best_calibrate_min = 0.0
        best_calibrate_max = 0.0
        best_calibrate_dict = {}
        best_model = ""
        
        for m in self.model:
            if self.verbose: print("Trying out the following model: %s" % (m))
            calibrate_output = self.calibrate_preds_func(seqs=seqs,
                                                        mods=mods,
                                                        identifiers=identifiers,
                                                        measured_tr=measured_tr,
                                                        correction_factor=correction_factor,
                                                        seq_df=seq_df,
                                                        use_median=use_median,
                                                        mod_name=m)

            self.calibrate_min, self.calibrate_max, self.calibrate_dict = calibrate_output

            if len(self.calibrate_dict.keys()) == 0: continue

            preds = self.make_preds(seqs=seqs,
                                    mods=mods,
                                    identifiers=identifiers,
                                    calibrate=True,
                                    seq_df=seq_df,
                                    correction_factor=correction_factor,
                                    mod_name=m)

            if len(measured_tr) == 0:
                perf = sum(abs(seq_df["tr"]-preds))
            else:
                perf = sum(abs(measured_tr-preds))

            if self.verbose: print("For current model got a performance of: %s" % (perf/len(preds)))

            if perf < best_perf:
                # TODO is deepcopy really required?
                best_calibrate_dict = copy.deepcopy(self.calibrate_dict)
                best_calibrate_min = copy.deepcopy(self.calibrate_min)
                best_calibrate_max = copy.deepcopy(self.calibrate_max)

                best_model = copy.deepcopy(m)                
                best_perf = perf
                
                if self.verbose: print("Model with the best performance got selected: %s" % (best_model))
        
        self.calibrate_dict = best_calibrate_dict
        self.calibrate_min = best_calibrate_min
        self.calibrate_max = best_calibrate_max
        self.model = best_model

        if self.verbose: print("Model with the best performance got selected: %s" % (best_model))

    def split_seq(self,
                a,
                n):
        """
        Split a list (a) into multiple chunks (n)

        Parameters
        ----------
        a : list
            list to split
        n : list
            number of chunks

        Returns
        -------
        list
            chunked list
        """

        # since chunking is not alway possible do the modulo of residues
        k, m = divmod(len(a), n)
        return(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
