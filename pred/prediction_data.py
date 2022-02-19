import os
import sys
import time
import json

data_name = ""
calibration_file_path = ""
model_path = []

data_name = "f25_fraction"
calibration_file_path = "./dataset/f25_fraction.csv"
model_path = []

data_name = "f25_retention"
calibration_file_path = "./dataset/f25_retention.csv"
model_path = []

data_name = "f96_fraction"
calibration_file_path = "./dataset/f96_fraction.csv"
model_path = [
    "full_hc_d0196_fraction_1fd8363d9af9dcad3be7553c39396960.hdf5",
    "full_hc_d0196_fraction_8c22d89667368f2f02ad996469ba157e.hdf5",
    "full_hc_d0196_fraction_cb975cfdd4105f97efa0b3afffe075cc.hdf5",
]

data_name = "f96_retention"
calibration_file_path = "./dataset/f96_retention.csv"
model_path = [
    "full_hc_d0196_retention_1fd8363d9af9dcad3be7553c39396960.hdf5",
    "full_hc_d0196_retention_8c22d89667368f2f02ad996469ba157e.hdf5",
    "full_hc_d0196_retention_cb975cfdd4105f97efa0b3afffe075cc.hdf5",
]



