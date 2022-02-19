import os
import sys
import time
import json
import pandas as pd
from deeplc import DeepLC, FeatExtractor




data_name = "f96_retention"
calibration_file_path = "./dataset/f96_retention.csv"
model_path = [
    "full_hc_d0196_retention_1fd8363d9af9dcad3be7553c39396960.hdf5",
    "full_hc_d0196_retention_8c22d89667368f2f02ad996469ba157e.hdf5",
    "full_hc_d0196_retention_cb975cfdd4105f97efa0b3afffe075cc.hdf5",
]



if __name__ == "__main__":

    peptide_file = "human_library.csv"
    pep_df = pd.read_csv(peptide_file, sep=",")
    pep_df['modifications'] = pep_df['modifications'].fillna("")

    calibration_file = calibration_file_path
    cal_df = pd.read_csv(calibration_file, sep=",")
    cal_df['modifications'] = cal_df['modifications'].fillna("")

    print(pep_df)

    # Initiate a DeepLC instance that will perform the calibration and predictions
    dlc = DeepLC(
        path_model=model_path,
        cnn_model=True,
        f_extractor=FeatExtractor(
            add_sum_feat=False,
            ptm_add_feat=False,
            ptm_subtract_feat=False,
            standard_feat=False,
            chem_descr_feat=False,
            add_comp_feat=False,
            cnn_feats=True,
            verbose=True),
        verbose=True,
    )
    dlc.calibrate_preds(seq_df=cal_df)


    preds = dlc.make_preds(seq_df=pep_df)
    json.dump(preds, open("pred_%s_rt.json" % data_name, "w"))




