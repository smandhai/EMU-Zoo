# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:32:45 2024

Check field completion against input list

@author: Soheb Mandhai
"""
import numpy as np
from glob import glob as gb
import settings
import functools as ft
import utils as ut


img_store = {"exclude":settings.prefix+"phase_1_hr_rel_v2/phase1_rel_hr_excluded_radio_cutouts/*",
			 "include":settings.prefix+"phase_1_hr_rel_v2/phase1_rel_hr_radio_cutouts/*" }

input_list = "..\sample_data\des\complexity_sources\subject_list_PhaseI_RGZ-EMU_cpx_ab4000_majax_ab20.csv"

total_src_list = []
"Iterate over included and excluded directories to note all files"
for key in img_store.keys():
	total_src_list += [gb(img_store[key])]

"List of sources including those that are excluded"
total_src_list = np.concatenate(total_src_list)

"Identify sources with _cross_ in their name"
src_keep_cond = list(map(ft.partial(ut.in_string,value='_cross_'),total_src_list))
total_src_list = total_src_list[src_keep_cond] #Only keep file names of interest
#%%
"Exctract field ID"
src_fields= list(map(ft.partial(ut.extract_srcs,split_by="SB"),total_src_list)) #Remove everything before SB
src_fields= list(map(ft.partial(ut.extract_srcs,split_by="J",ind=0),src_fields)) #Remove everything after SB
src_fields = np.asarray(src_fields).astype(int)

"Repeat the same for sources in the input list"
input_data= np.loadtxt(input_list,dtype=str,delimiter=',',skiprows=0)
input_src_list = input_data.T[0][1:]
input_fields =list(map(ft.partial(ut.extract_srcs,split_by="SB"),input_src_list)) #Remove everything before SB
input_fields= list(map(ft.partial(ut.extract_srcs,split_by="J",ind=0),input_fields)) #Remove everything after SB
input_fields = np.asarray(input_fields).astype(int)

"Get count of field entries in the input complexity"
input_unique = np.unique(input_fields,return_counts=True)
input_dic = {}

for i,field in enumerate(input_unique[0]):
	input_dic[field] = input_unique[1][i] 
"Get count of field entries in the src list"

src_unique = np.unique(src_fields,return_counts=True)
src_dic = {}
for i,field in enumerate(src_unique[0]):
	src_dic[field] = src_unique[1][i] 
#%%
"Calculate fractions of field completion"
for key in input_dic.keys():
	if key in src_dic.keys():
		frac = src_dic[key]/input_dic[key]
	else:
		frac = 0
	print("Fraction for field [{}]: {:.2f}".format(key,frac))