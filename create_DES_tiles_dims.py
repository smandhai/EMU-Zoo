# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:19:51 2024

Collates the bounding boxes of all downloaded DES tiles into a single file. Used for referencing during EMU-ZOO image cutout creation.

@author: Soheb Mandhai
"""

import glob
from astropy.io import fits
import numpy as np

def find_id(tile_id):
	return tile_id.split("/")[-1].split("_")[0]
def return_bbox(tile_fname):
	tile_id = find_id(tile_fname)
	with fits.open(tile_fname) as tile:
		hdr = tile[1].header
		ra_min  = hdr["RAC2"] 
		ra_max = hdr["RAC1"]
		dec_max = hdr["DECC4"] 
		dec_min = hdr["DECC2"] 
		return tile_id, ra_max,ra_min,dec_max,dec_min
	
des_dir = "/mnt/shared/des/DES_Tiles_new"
tiles = glob.glob(des_dir + "/*_i.fits.fz")
#print(tiles)
#tile_ids = list(map(find_id,tiles))

coords= list(map(return_bbox,tiles))
#tab = np.c_[tile_ids,coords]
tab = np.array(coords)
np.savetxt("new_DES_tiles_dim.txt",tab,fmt="%s",delimiter=' ')
