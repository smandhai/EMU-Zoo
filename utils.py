# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 17:22:46 2024

Utilities to conduct basic operations

@author: Soheb Mandhai
"""

import os
import shutil as sh

def make_dir(dirname):
    if os.path.isdir(dirname)==False:
        os.mkdir(dirname)
        
def extract_sources(src_fname,field_split=''):
    temp_src = src_fname.split(field_split)[-1].split("_")[0] 
    return temp_src

def move_file(src_fname,new_dir='',split="\\"):
    if len(src_fname.split(split))>0:
        src_fname_source = src_fname.split(split)[-1] #If the source is nested in other directories
    else:
        src_fname_source = src_fname #If the file is already in the root directory
    sh.move(src_fname,new_dir+src_fname_source)
    
