#!/bin/bash

#BSUB -g /d.goldfarb/compute
#BSUB -G compute-d.goldfarb
#BSUB -a 'docker(wenjin27/deeplc:0.2.0)'
#BSUB -J wenjintest
#BSUB -R 'select[ gpuhost, mem>16000, tmp>20 ] rusage[mem=16000, tmp=20] span[hosts=1]'
#BSUB -q general
#BSUB -n 10
#BSUB -M 16GB
#BSUB -gpu 'num=1:gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes'
#BSUB -o /storage1/fs1/d.goldfarb/Active/Projects/Wenjin/code/dlcpg/train.%J.log
#BSUB -e /storage1/fs1/d.goldfarb/Active/Projects/Wenjin/code/dlcpg/train.%J.err



python3.7 run_full_mod.py

