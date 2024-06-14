# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:38:23 2024

This script will verify the sources align with the complexity analysis
+ allow subselections to be made

@author: Soheb Mandhai
"""

import numpy as np
import settings
import glob
import functools as ft
import shutil as sh
import utils as ut
import os
from astropy.io import ascii

#%%
"Completed Fields"
comp_fields= ["9287","9351","9434","9437","9442"]
#%%
root_dir = settings.prefix  
chosen_sources = "/complexity_sources/subject_list_PhaseI_RGZ-EMU_cpx_ab4000_majax_ab20.csv"
complexity_cat = "/complexity_sources/merged_10selavy_Gary_complexity.csv"
image_dir = settings.prefix+"/phase1_rel/phase1_3150"

"Extract component ID"
complexity_data = np.loadtxt(root_dir+complexity_cat,dtype=str,delimiter=',',skiprows=0)
complexity_data_headers = complexity_data[0] #Headers of the complexity data table
complexity_data = complexity_data[1:]#remove headers
"Component suffix - extract column with suffix values"
suffix_colname = 'suffix'
suffix_colind = np.where(complexity_data_headers==suffix_colname)[0][0]
"Component name - for matching values together"
name_colname= "component_name"
name_colind = np.where(complexity_data_headers==name_colname)[0][0]
"Component ID"
id_colname= "component_id"
id_colind = np.where(complexity_data_headers==id_colname)[0][0]
"Add strings together into a single list"
complexity_name = list(map(ut.add_strings,complexity_data.T[name_colind],complexity_data.T[suffix_colind]))
complex_id = list(map(ft.partial(ut.extract_srcs,split_by='_',add_prefix='',ind=0),complexity_data.T[id_colind]))
complexity_name = list(map(ut.add_strings,complex_id,complexity_name))


"Data"
chosen_data= np.loadtxt(root_dir+chosen_sources,dtype=str,delimiter=',',skiprows=0)

image_dir = glob.glob(image_dir+'/*')
image_dir_all = image_dir

"ID's of all the imagse"
image_srcs= list(map(ft.partial(ut.extract_srcs,split_by='SB',add_prefix='SB'),image_dir))
image_srcs_clean= list(map(ft.partial(ut.extract_srcs,split_by='_cross_',add_prefix='',ind=0),image_srcs))

"Cross-match image dir sources with complexity sources"
comp_images_match,comp_match_inds,img_match_inds = np.intersect1d(complexity_name,image_srcs_clean,return_indices=True)
"Reorder variables to match table ordering"
sorting = np.argsort(comp_match_inds)
comp_images_match = comp_images_match[sorting]
comp_match_inds = comp_match_inds[sorting]
img_match_inds = img_match_inds[sorting]

"Missing images - For debugging"
missing_img_ind = np.in1d(np.arange(len(image_srcs_clean)),img_match_inds,invert=True)
missing_imgs = np.asarray(image_srcs_clean)[missing_img_ind]

#%%
"Find complexity values"

comp_colname= "Complexity (bytes)"
comp_colind = np.where(complexity_data_headers==comp_colname)[0][0]
complexities = complexity_data.T[comp_colind]
complexities_match = np.asarray(complexities[comp_match_inds]).astype(float)

"Sort by complexities"
sorting = np.argsort(complexities_match)[::-1]
comp_images_match = comp_images_match[sorting]
comp_match_inds = comp_match_inds[sorting]
img_match_inds = img_match_inds[sorting]
complexities_match = complexities_match[sorting]

"Reorder images"
image_srcs_clean = np.asarray(image_srcs_clean)[img_match_inds]

"Image dirs reordered"
image_dir = np.asarray(image_dir)[img_match_inds]
  
#%%
"Apply image cut, remove sources that are unneeded"
select_sources = 3150
complexity_cut=4000
image_dir_cut = image_dir[complexities_match>=complexity_cut]
#%%
"Only keep sources with complete fields"
image_keep_cond = np.tile(False,len(image_dir_cut))

for field in comp_fields:
	temp_cond = list(map(ft.partial(ut.in_string,value='SB{}'.format(field)),image_dir_cut))
	image_keep_cond[np.where(np.asarray(temp_cond)==True)]=True
image_dir_cut = image_dir_cut[image_keep_cond]
#%%
"Delete unneeded sources"
count = 0
for image in image_dir_all:
	if image not in image_dir_cut:
		#print("Deleting Image: {}".format(image))
		os.remove(image)
		count+=1
		
print(f"Number of files deleted: {count}")

#%%
"Store remaining outputs to a folder"
meta_data = complexity_data[comp_match_inds]
ascii.write(meta_data,"metadata.csv",overwrite=True,delimiter=',',names=complexity_data_headers)
