# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 17:22:46 2024

Utilities to conduct basic operations

@author: Soheb Mandhai
"""

import os
import shutil as sh

def to_abc(num):
	"Only works for two letters"
	reps = int(num/27)
	string = ''
	pre_num = 0
	for r in range(reps+1):
# 		print("String length: ",len(string))
# 		print(r)
		if (num>26):
			num-=26
			pre_num += 1
# 		print(num)
# 		print(pre_num)
		string += chr(ord('`')+num)
		if (pre_num > 0) & (len(string)>1):
			string = string.replace(string[-2],chr(ord('`')+1),pre_num)#Sets previous letter as a
		if pre_num > 26:
			pre_num = 0
				
	return string
#%%
def make_dir(dirname):
	if os.path.isdir(dirname)==False:
		os.mkdir(dirname)
		
def extract_sources(src_fname,field_split=' ',add_prefix='',ind=-1,return_suffix=True):
	temp_src = add_prefix + src_fname.split(field_split)[ind].split("_")[0].split(".")[0]
	if return_suffix:
		return temp_src,src_fname.split(temp_src)[-1].split(".")[0]
	return temp_src

def move_file(src_fname,new_dir='',split="\\"):
	if len(src_fname.split(split))>0:
		src_fname_source = src_fname.split(split)[-1] #If the source is nested in other directories
	else:
		src_fname_source = src_fname #If the file is already in the root directory
	sh.move(src_fname,new_dir+src_fname_source)
	
