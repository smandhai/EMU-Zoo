# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:56:02 2024

This code will cross-match cutouts with sources to ensure
every source is accounted for


@author: Soheb Mandhai
"""

import numpy as np
import os
import glob
import settings
import functools as ft

def extract_id(string,start_cut = None, end_cut = None):
    "Extracts source ID from a string"
    return "J"+string.split(start_cut)[-1].split(end_cut)[0]

#cutout_dir = [settings.prefix + "phase1/phase1_radio_cutouts_6x6/",
#              settings.prefix+"phase1/excluded_phase1_radio_cutouts_6x6/"]

cutout_dir = [settings.prefix+"beta_radio_cutouts_6x6"]

cat_dir = settings.dataloc+ ''.join(settings.cat_sub.split('/')[:-1])
files= glob.glob(cat_dir+"/AS101_Continuum_Component_Catalogue_*.csv")

total_reps = 0
rep_src_list = []
for i in range(len(files)):
    print(f"Loading file: {files[i].split('/')[-1]}")
    temp_f  = np.loadtxt(files[i],delimiter=',',dtype=str,skiprows=1,usecols=8)
    temp_src,temp_ind,uni_cnt = np.unique(temp_f,return_inverse=True,return_counts=True)
    rep_srcs = temp_src[np.where(uni_cnt>1)]
    total_reps += len(rep_srcs)
    print(f"Sources in file: {len(temp_f)}")
    print(f"Repeated Sources: {len(rep_srcs)}")
    if i ==0:
        total_src_list = temp_f
    if i !=0:
        total_src_list = np.append(total_src_list,temp_f)
    rep_src_list = np.append(rep_src_list,rep_srcs)
print(f"Total number of sources: {len(total_src_list)}")
#%%
cutout_files= [] 
for i in range(len(cutout_dir)):
    cutout_files.append(glob.glob(cutout_dir[i]+"/*.fits"))
    if i==0:
        cutout_IDs = list(map(ft.partial(extract_id,start_cut = 'J',end_cut='.fits'),cutout_files[-1]))
    else:
        temp_cutout_IDs = list(map(ft.partial(extract_id,start_cut = 'J',end_cut='.fits'),cutout_files[-1]))
        cutout_IDs = np.append(cutout_IDs,temp_cutout_IDs)
cutout_files = np.concatenate(cutout_files).astype(str)
cutout_argsort=np.argsort(cutout_IDs)

cutout_files = np.asarray(cutout_files)[cutout_argsort]
cutout_IDs = np.asarray(cutout_IDs)[cutout_argsort]
print(f"Total number of cutouts: {len(cutout_IDs)}")
print(f"Number of repeated sources in total: {total_reps}")
print(f"\nExpected number of missing matches = {len(total_src_list)-len(cutout_IDs)}")
np.savetxt("cutout_sources.dat",np.array([cutout_IDs, cutout_files]).T,fmt="%s",delimiter=',')
np.savetxt("duplicate_sources.dat",np.array([rep_src_list]).T,fmt="%s",delimiter=',')



# =============================================================================
# Actual matching of sources
# =============================================================================
miss_srcs = total_src_list[np.in1d(total_src_list, cutout_IDs, invert=True,assume_unique=False)]
print(f"Actual number of missing sources: {len(miss_srcs)}")
np.savetxt("missing_sources.dat",miss_srcs,fmt="%s",delimiter=',')
