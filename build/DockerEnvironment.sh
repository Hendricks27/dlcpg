#!/bin/bash


apt-get update
apt-get install -y software-properties-common

#add-apt-repository ppa:deadsnakes/ppa
#apt-get update

# apt-get install -y python3.7 python3-pip python-dev
apt-get install -y python3-pip python-dev

python3 -m pip install --upgrade pip
python3 -m pip install \
  deeplc==0.2.0 \
  joblib==1.1.0 \
  scikit-learn==1.0.2 \
  xgboost==1.5.2


chmod -R 0777 /code
chmod -R 0777 ./

mkdir /data
chmod -R 0777 /data

mkdir /data/config
chmod -R 0777 /data/config

mkdir /data/tmp
chmod -R 0777 /data/tmp



