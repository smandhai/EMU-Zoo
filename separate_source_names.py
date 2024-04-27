# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 16:36:54 2024

@author: Soheb
"""

import pandas as pd
import functools as ft
import settings
import utils as ut
import numpy as np

separate_by_field= True

data = pd.read_csv(settings.override_src)
srcs = list(map(ft.partial(ut.extract_sources,field_split='J',add_prefix="J",return_suffix=True)
				,data["File"].values))
SB_fields=  list(map(ft.partial(ut.extract_sources,field_split='J',ind=0,return_suffix=False),data["File"].values))

unique_SBfields = np.unique(SB_fields)

for field in unique_SBfields:
	cond = np.where(np.asarray(SB_fields)==field)
	temp_srcs= np.array(srcs).astype(str)[cond]
	np.savetxt("{}_".format(field)+settings.override_src.split("/")[-1],temp_srcs,fmt='%s' ,delimiter=',')
