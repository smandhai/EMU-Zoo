# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:53:48 2024

A simple post-processing script to move cross images

@author: Soheb Mandhai
"""

import glob
import settings
import shutil as sh
import functools as ft
import numpy as np
import os

def reduce_list(string,snippet=""):
	if len(string.split(snippet))>1:
		return string
	else:
		return
		
def copy_item(old_dir,new_dir):
	itm = old_dir.split("\\")[-1]
	new_dir += "\\"+itm
	sh.copy(old_dir,new_dir)

file_dir = settings.prefix+"/phase_1b_hr_rel/phase1b_for_release/"
new_dir = settings.prefix+"/phase_1b_hr_rel/final_release"

file_list = glob.glob(file_dir+"*")

reduced_list = list(map(ft.partial(reduce_list,snippet='_cross_'),file_list))
reduced_list = np.asarray(reduced_list)[np.where(np.asarray(reduced_list)!=None)[0]]

if os.path.isdir(new_dir)==False:
	os.mkdir(new_dir)
for i in range(len(reduced_list)):
	copy_item(reduced_list[i],new_dir)