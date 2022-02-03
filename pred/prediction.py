import os
import sys
import time
import pandas as pd
from deeplc import DeepLC, FeatExtractor



if __name__ == "__main__":
    peptide_file = "test.csv"

    pep_df = pd.read_csv(peptide_file, sep=",")
    pep_df['modifications'] = pep_df['modifications'].fillna("")

    calibration_file = "fraction_df_0196_woo.csv"
    cal_df = pd.read_csv(calibration_file, sep=",")
    cal_df['modifications'] = cal_df['modifications'].fillna("")

    print(pep_df)

    # Initiate a DeepLC instance that will perform the calibration and predictions
    dlc = DeepLC(
        path_model=["full_hc_fraction_df_0196_woo_1fd8363d9af9dcad3be7553c39396960.hdf5"],
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
    print(preds)


