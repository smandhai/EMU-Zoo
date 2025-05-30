# EMUZoo cutout generation code
# Author: Emma Alexander (emma.alexander@gmail.com)
# Modified by: Soheb Mandhai (soheb.mandhai@manchester.ac.uk)
# https://github.com/EmmaAlexander/EMU-Zoo

# NOTE: currently limited to EMU Pilot I survey fields area. Will needs tweaks for main survey data. 
# Also note: copied over from a jupyter notebook. Should hopefully still work.
# Original code (and data) can be found on AusSRC:
# scripts at /mnt/shared/home/ealexander/
# data at /mnt/shared/des

#-----------------------------------------------------------
#imports
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
from astroquery.skyview import SkyView
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob, os, sys
from astropy.nddata import Cutout2D,PartialOverlapError
import math
import cmasher as cmr
from astropy.visualization import lupton_rgb
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from matplotlib.patches import Ellipse
import settings
#import pandas as pd
import functools as ft
import utils as ut
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from importlib import reload
import warnings
from astropy.wcs import FITSFixedWarning
import h5py
reload(settings)

warnings.simplefilter("ignore",category=FITSFixedWarning)
#%%
#-----------------------------------------------------------
def ashinh_scale(array,zeropoint=0,scale=1):
	scaled=(array-zeropoint)*np.arcsinh(array*scale)/(array*np.arcsinh(scale))
	return scaled

def percentile(array,percent):
	val=np.percentile(array[np.isfinite(array)],percent)
	return val

def save_contours(contour,fname='source.h5',ppa=30,coords=(0,0),x_size=360,y_size=360): 
	"IN PROGRESS, TO BE COMPLETED"
	"Takes contour values and saves all contour levels to a h5 file"
	cs = contour
	with h5py.File(fname,'w') as f: #Open File
		f.attrs["pixels_per_arcmin"] = ppa #pixel per arcmin translation
		f.attrs["ra_cen"] = "{:.5f}".format(float(coords[0])) #central coordinate, (0,0) means unspecified coordinates.
		f.attrs["dec_cen"] = "{:.5f}".format(float(coords[1]))
		f.attrs["x_px_size"] = x_size
		f.attrs["y_px_size"] = y_size
		for i in range(len(cs.allsegs)):
# 			if len(cs.allsegs[i])<3:
# 				print(i)
			cs_cond = np.where(np.asarray(list(map(ft.partial(point_check,coord=(x_size/2,y_size/2)),cs.allsegs[i]))).astype(bool)==True)
			if len(cs_cond[0])>0:
				#print(i)
				cs_cond = cs_cond[0][0]
				f["Level_{}".format(i)] = cs.allsegs[i][cs_cond]
	return

def point_check(contour,coord=(0,0),search_pix = 5):
	"Checks if a point exists within the contour"
	point = Point(coord)
	poly=  Polygon(contour)
	#circle_p = Point(coord)
	exists = poly.contains(point)
	if exists ==False:	
		"Searching around to check if pixels have been missed"
		circle = point.buffer(search_pix)
		exists = circle.intersects(poly)
	
# 	exists = poly.contains(point)
# 	if exists ==False:
# 		"Search around the centre"
# 		point = Point((coord[0]+search_pix,coord[1]))
# 		poly=  Polygon(contour)
# 		exists = poly.contains(point)
# 		if exists:
# 			return exists		
# 		point = Point((coord[0]-search_pix,coord[1]))
# 		poly=  Polygon(contour)
# 		exists = poly.contains(point)	
# 		if exists:
# 			return exists		
# 		point = Point((coord[0],coord[1]+search_pix))
# 		poly=  Polygon(contour)
# 		exists = poly.contains(point)
# 		if exists:
# 			return exists		
# 		point = Point((coord[0],coord[1]-search_pix))
# 		poly=  Polygon(contour)
# 		exists = poly.contains(point)		
# 		if exists:
# 			return exists
	return exists
def masking(data,contours,mask_value=0,ppa=30,exclude=True,pixel_thresh=1):
	#global cs
	excluded_source=False
	x_size,y_size= data.shape
	pixels_per_arcmin = ppa
	# selection_window_x = (int(x_size/2 - pixels_per_arcmin),int(x_size/2 + pixels_per_arcmin))
	# selection_window_y = (int(y_size/2 - pixels_per_arcmin),int(y_size/2 + pixels_per_arcmin))
	# data = data[selection_window_x[0]:selection_window_x[1],
	#					   selection_window_y[0]:selection_window_y[1]]
	masked_data = np.array(data)
	masked_data[(masked_data<contours[0])]= 0 #Remove all background
	"Remove edges to ensure the contours can be complete connected"
	edge1 = masked_data[0:x_size+1,0]
	edge2 = masked_data[0:x_size+1,y_size-1]
	edge3 = masked_data[0,0:y_size+1]
	edge4 = masked_data[x_size-1,0:y_size+1]

	masked_data[0:x_size+1,0] = 0
	masked_data[0:x_size+1,y_size-1] = 0
	masked_data[0,0:y_size+1] = 0
	masked_data[x_size-1,0:y_size+1] = 0
	"Debug contour plotting"
	#plt.imshow(masked_data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
	cs = plt.contour(masked_data,levels=contours,colors='grey')
	
	"Find the average distance from the centre of the frame"
	r = np.asarray([np.quantile(np.linalg.norm(cs.allsegs[0][i]-(x_size/2,y_size/2),axis=1),0.6) for i in range(len(cs.allsegs[0]))])
	# r_ind_sorted=  np.argsort(r)
	# r = r[r_ind_sorted]
	#print(r)
	d_left = np.where(masked_data[int(x_size/2)-1:0:-1,int(y_size/2)] ==0)[0]
	d_right = np.where(masked_data[int(x_size/2)-1:int(x_size):,int(y_size/2)] ==0)[0]
	d_down =  np.where(masked_data[int(x_size/2),int(y_size/2)-1:0:-1] ==0)[0]
	d_up =  np.where(masked_data[int(x_size/2),int(y_size/2)-1:int(y_size)+1] ==0)[0]
	if d_left.size == 0:
		d_left = [0]
	if d_right.size == 0:
		d_right = [0]
	if d_up.size == 0:
		d_up = [0]
	if d_down.size==0:
		d_down = [0]
		
	
	extended_source =False
	try:
		d_mean = np.mean([d_left[0],d_right[0],d_down[0],d_up[0]])
# 		d_mean = np.amax([d_left[0]+d_right[0],d_down[0]+d_up[0]])
		#print(d_mean)
		if (d_mean > x_size/2)|(d_mean > y_size/2):
 			print("Potential extended source")
 			extended_source =True
		# print([d_left[0],d_right[0],d_down[0],d_up[0]])
	except:
		print("Potential extended source")
		extended_source =True #This is triggered if the source is likely large
		d_mean = np.mean([x_size,y_size])
	   # d_mean = np.mean([d_left[0],d_right[0],d_down[0],d_up[0]])
# 		d_mean = np.amax([d_left[0]+d_right[0],d_down[0]+d_up[0]])

 	# print("D_mean: ", d_mean)
 	# print(r)
	if len(r)>0:
# 		
# 		if len(r)>1:
#  			if len(r[r<ppa])>0:
# 				 r_cond = (r<np.quantile(r[r<ppa],0.5))&(r<ppa)
# 				 print(r[r_cond])
# 				 if (x_size>ppa)|(y_size>ppa)&(extended_source ==True):
#  					r_cond = r<np.quantile(r,0.5)
#  					print("we rolling")
#  			else:
# 				 r_cond = []
#  			if len(r[r_cond])==0:
# 				 r_cond = r<np.quantile(r,0.50)
# 				 r_cond = r==r.min()
# # 				print("we rolling again")
# # 			print(r[r_cond])
#  			cs_cond = np.where(r_cond)
# # 			if (len(cs_cond)>1) &(extended_source==False):
# # 				 print(12121)
# # 				 cs_cond = cs_cond[np.where(np.abs(r[cs_cond]-d_mean) ==np.abs(r[cs_cond]-d_mean).min())]
#  			cs_cond = np.where(r_cond)[-1]
#  			print(cs_cond)
#  			if len(cs_cond)>1:#If there's more than one radius, find the closest to the mean distance
# 			 cs_cond_selector =  np.where(np.abs(r[cs_cond]-d_mean) == np.abs(r[cs_cond]-d_mean).min() )[0]
# # 				print(cs_cond_selector)
# 				cs_cond = cs_cond[cs_cond_selector][0]
#  			else:
# 				cs_cond = cs_cond[0]
#  			print(r[r_cond])
#  			print(r[cs_cond])
# 		else:
#  			cs_cond = 0
		
		# print(r_ind_sorted)
		# print(cs.allsegs[0])
		cs_cond = np.where(np.asarray(list(map(ft.partial(point_check,coord=(x_size/2,y_size/2)),cs.allsegs[0]))).astype(bool)==True)
		if len(cs_cond[0]) >0:
			#print(cs_cond)
			cs_cond = cs_cond[0][0]
		else:
			#print(cs_cond)
			warnings.warn("Insuffient contours: Potentially due to unclosed contours. Data not masked.")
			excluded_source=True
			return masked_data, excluded_source
		x = cs.allsegs[0][cs_cond][:,0]
		y = cs.allsegs[0][cs_cond][:,1]
		#print(cs_cond)
		plt.plot(x,y)
		"Radius to boundary"
		#r= np.sqrt(x**2+y**2)
		#theta = np.arccos(x/r)
		#x_split = np.split(tiny,np.round(x).astype(int))
		#y_split = np.split(tiny,np.round(y).astype(int))
		#tiny.T[x.astype(int),y.astype(int)]=0 #Show boundary
		#plt.imshow(masked_data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
		mask_sorter = np.argsort(np.round(x).astype(int)) #Arrange values that need to be masked
		x_sorted = np.round(x).astype(int)[mask_sorter] #Order
		y_sorted = np.round(y).astype(int)[mask_sorter]
		unique_x_sorted = np.unique(x_sorted) #Find rows that need to be masked over
		"Mask out the main source row by row"
		pixel_bright = 0 #Initialise number of bright pixels
		excluded_source=True #Source is excluded unless it has a second contour
		exclude = True
		for x_ind in unique_x_sorted:
			cond = np.where(x_sorted==x_ind)
			#print(tiny.T[x_ind,y_sorted[cond].min():y_sorted[cond].max()+1])
			"Mask out the source"
			row = masked_data.T[x_ind,y_sorted[cond].min()-pixel_thresh:y_sorted[cond].max()+1+pixel_thresh] 
			cond_row =  np.where(row>=contours[1])[0]
			# print(row[cond_row])
			if len(cond_row)>=1:
				#print("Multi-contour source, do not exclude")
				excluded_source=False
				exclude=False
				"If there's a single bright pixel"
				pixel_bright += len(np.where(row>=contours[1])[0])
				
				#print(len(np.where(row>contours[1])[0]))
			masked_data.T[x_ind,y_sorted[cond].min():y_sorted[cond].max()+1] = mask_value
		#plt.imshow(tiny,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
		#tiny[tiny<radio_contours[1]]= 0
		#plt.imshow(tiny,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))

		"If there's a bright source found after the source has been masked out"
		#print("pixels ",pixel_bright)
		if pixel_bright <=pixel_thresh:
			"If there's a single bright pixel"
			exclude=True
			excluded_source=True
		else:
			exclude = False
			excluded_source =False
		if (len(np.where(masked_data>contours[1])[0]) !=0)&(exclude==True):
			print("Bright source found. Source is to be excluded")
			excluded_source=True
			
	else:
		print("No contours found")
		excluded_source=True
	plt.close() #Ensures the figure has been closed
	"Restore edges"
	masked_data[0:x_size+1,0] = edge1
	masked_data[0:x_size+1,y_size-1] = edge2
	masked_data[0,0:y_size+1] = edge3
	masked_data[x_size-1,0:y_size+1] = edge4
	return masked_data,excluded_source
#-----------------------------------------------------------

# =============================================================================
# OVERRIDE SETTINGS
# =============================================================================
override_src = settings.override_src #'J202254-540537' #"J202505-540405" 'J202254-540537'#None

#-----------------------------------------------------------


# Hard coded variables (data locations etc)

dataloc=settings.dataloc #askap data location
WISEfiles=glob.glob(settings.WISEfiles_dir) #WISE raw data tiles location
duplicate_sources = {} #A dictionary containing duplicate sources
SB= settings.SB#'9351'
add_suff= "" #Initialise extra suffix in case it needs to be added
dss = False #DSS not used by default
#hard coded image names due to inconsistent naming conventions but could be streamlined
image='image.i.SB'+SB+'.cont.taylor.0.restored.fits'

#note that there are two catalogue types: island and component
#cat=dataloc+'catalogues/AS101_Continuum_Island_Catalogue_'+SB+'.csv'
cat=dataloc+settings.cat_sub.format(SB,settings.component_number)
cat_type=settings.cat_type #change to true if using island catalogue

# folder for cutout data output if you want to save radio fits files
# not currently implemented due to data storage issues
radiooutloc=settings.radio_output+SB
if os.path.isdir(radiooutloc) ==False:
	os.system('mkdir '+radiooutloc+SB)

overlayloc=settings.overlayloc+SB
if os.path.isdir(overlayloc) ==False:
	os.system('mkdir '+overlayloc+SB)
overlay_suffix='.png'
#you can add to the above if you need different versions 
#e.g. change to '_v2.png'

arcmins=12. #maximum size of cutout in armins
# disclaimer: things may break if you change this...
# number of pixel dimensions
npix_edge=int(15*arcmins)
deg_edge=arcmins*0.025/3.

#-----------------------------------------------------------
# set up some plot parameters
# note that many of these are no longer used but kept for posterity

dpi=300
dpi_save = 300
my_dpi = 100#300
sub_font_size = 10
plt.rc('font', size=0.5)  
plt.rcParams.update({'lines.linewidth':0.8})
plt.rc('font', size=10)		  # controls default text sizes
infernocmap=plt.cm.inferno
infernocmap.set_bad('black',1)
magmacmap=plt.cm.magma
magmacmap.set_bad('black',1)
redcmap=plt.cm.Reds_r
redcmap.set_bad('black',1)
orangecmap=plt.cm.Oranges_r
orangecmap.set_bad('black',1)
greycmap=plt.cm.Greys_r
greycmap.set_bad('white',1)
viridis=plt.cm.viridis
viridis.set_bad('black',1)
gist_heat=plt.cm.gist_heat
gist_heat.set_bad('black',1)
greencmap=plt.cm.Greens_r
greencmap.set_bad('black',1)

twl_blue = cmr.get_sub_cmap('twilight_shifted', 0, 0.5)
twl_red = cmr.get_sub_cmap('twilight', 0.5, 1)
twl_red.set_bad(twl_red(0),1)
twl_red_r = cmr.get_sub_cmap('twilight_shifted', 0.5, 1)
twl_blue_r = cmr.get_sub_cmap('twilight',0,0.5)

#-----------------------------------------------------------

#read in full ASKAP image 
hdu=fits.open(dataloc+image)
image=hdu[0].data.squeeze()
header=hdu[0].header
wcs= WCS(hdu[0].header).celestial
#hdu.close()
#some useful parameters that don't get used here
pixscale=header['CDELT2']
bmaj_pix=header['BMAJ']/header['CDELT2']
bmin_pix=header['BMIN']/header['CDELT2']
bpa=header['BPA']

# read in catalogue
catalogue=np.genfromtxt(cat.format(SB),dtype='str',delimiter=',')
headers=catalogue[0,:]
data=catalogue[1:,:]
#below lines for if you'd like to check catalogue headers
#for i in range(0,45):
	#print(i, headers[i],data[0,i])

#sort data by major axis (optional, but data_sorted is array used going forward
#note that 19 needs changing to 31 for an island catalogue because they have different columns
#TODO: use header name rather than column index (pandas...?)
data_sorted = data[data[:,19].argsort()[::-1]] #componenet cat, major axis

#-----------------------------------------------------------

"DEPRECATED"
# =============================================================================
# #get the coordintes of the WISE images to know which to load and use
# WISEcoords=[]
# for file in WISEfiles:
#	 filename=file.split('/')[-1] #get the actual file name
#	 ra=float(filename[0:4])/10
#	 dec=-1*float(filename[5:8])/10
#	 WISEcoords.append(SkyCoord(ra,dec,unit=u.degree,frame='fk5'))
# WISEcoords=SkyCoord(np.asarray(WISEcoords))
# =============================================================================

WISEtiles=np.genfromtxt(settings.WISEtiles,dtype='str') #this file should be in the same location as the code
DEStiles=np.genfromtxt(settings.DEStiles,dtype='str') #this file should be in the same location as the code

#-----------------------------------------------------------
#%%

# =============================================================================
# Prepare output directories
# =============================================================================
ut.make_dir(settings.radio_output.split(settings.field_ref)[0]) #Create output directory
ut.make_dir(settings.exclusion_dir) # Create exclusion directory

#MAIN BIT
src_ind = 0
if cat_type=="island":
	src_ind = 6
	#_,find_src,_=np.intersect1d(data_sorted[:,6],override_src,assume_unique=True,return_indices=True)
	#find_src = np.where(np.asarray(data_sorted[:,6]).astype(str)==override_src)
elif cat_type=="catwise":
	src_ind = 2
	#_,find_src,_=np.intersect1d(data_sorted[:,2],override_src,assume_unique=True,return_indices=True)
elif cat_type=="component":
	src_ind = 7
	#_,find_src,_=np.intersect1d(data_sorted[:,7],override_src,assume_unique=True,return_indices=True)
	#find_src = np.where(np.asarray(data_sorted[:,7]).astype(str)==override_src)
elif cat_type=="component_new":
	src_ind = 8

add_suff= "" #Intialise suffix
#Override Source
if override_src != None:
	unique_srcs = False #Assume Unique srcs, by default every source is assumed to be unique

	if settings.use_file:
		#source_list = pd.read_csv(override_src,header=None)
		source_list = np.genfromtxt(override_src,dtype=str,delimiter=',')
		if len(source_list.shape)>1:
			override_src = list(source_list[0:,0])
			src_mode = list(source_list[0:,-1])
		else:
			override_src = list(source_list)
			src_mode = list(np.tile("",len(override_src)))
		if len(settings.override_src.split("duplicate"))>1:
			unique_srcs= True #Consider non-unique sources
			
	else:
		if type(override_src) == str: #Convert single source name to a list
			override_src = [override_src]
	_,find_src,_=np.intersect1d(data_sorted[:,src_ind],override_src,assume_unique=unique_srcs,return_indices=True)
	find_src = np.arange(len(data_sorted))[np.in1d(data_sorted[:,src_ind],override_src)]
	if len(find_src)>0:
		data_sorted = data_sorted[find_src]
	else:
		raise ValueError("Override value(s) not found")
	
	"Special exception for sources in the duplicate file"
	"Ensures that sources that are duplicated in the SAME field are accounted for only"


	"Find the count of sources that exist multiple times in this field"
	temp_srcs,temp_ind,temp_counts = np.unique(data_sorted[:,src_ind],return_index=True,return_counts=True)
	"Find where counts are >1 (i.e. duplicates). Only keep sources with duplicates"
	keep_srcs = data_sorted[:,src_ind][temp_ind[np.where(temp_counts>1)] ]
	if len(keep_srcs)>0:
		"Find the indices to with the catalogue to pair duplicates with"
		keep_inds = [np.where(data_sorted[:,src_ind]==keep_srcs[i]) for i in range(len(keep_srcs))]
		if len(keep_inds) ==0:
			raise ValueError("No duplicates found, nothing to do...")
	"Overwrite data_sorted to remove non-duplicate sources"
	if override_src ==  None:	
		override_src = ""
	if type(override_src) !=list:
		if len(override_src.split("duplicate"))>1:
			print("Special overwride file found... filtering out non-duplicates")
			data_sorted = data_sorted[np.concatenate(keep_inds).flatten()]
		
else:
	"If no overrides are found"
	src_mode = list(np.tile("",len(data_sorted)))
exclusion_list = [] #List of sources to be excluded

#%%
#loop over sources in the list

for i in range(0,len(data_sorted)):
	if i % 100 == 0:
		print("Current Iteration: ",i)
		#or use a better tracking method idk

	if settings.cat_type=="island":
		#get values from catalogue file
		# not all of these are used but I have left them in
		src = data_sorted[i,6]	
		n_components = data_sorted[i,7]
		ra_hms_cont = data_sorted[i,8]
		dec_dms_cont = data_sorted[i,9]
		ra_deg_cont = data_sorted[i,10]
		dec_deg_cont = data_sorted[i,11]
		maj_axis = float(data_sorted[i,13])
		min_axis = float(data_sorted[i,14])
		pos_ang = float(data_sorted[i,15])
		flux_int = float(data_sorted[i,16])
		flux_int_err = float(data_sorted[i,16])
		flux_peak = float(data_sorted[i,18])
		background_noise = float(data_sorted[i,20])
		n_pix  = float(data_sorted[i,30])
		solid_angle = float(data_sorted[i,31])
		x_cen = float(data_sorted[i,35])
		y_cen = float(data_sorted[i,36])
		rms_median = data_sorted[i,24]
	elif settings.cat_type=="catwise":
		src=data_sorted[i,2]
		ra_hms_cont	=data_sorted[i,3]
		dec_dms_cont=data_sorted[i,4]
		ra_deg_cont	=data_sorted[i,5]
		dec_deg_cont=data_sorted[i,6]	  
		background_noise =float(data_sorted[i,32]) 
		rms_median = data_sorted[i,32]	

	elif settings.cat_type=='component':
		#component catalogue being used
	
		src=data_sorted[i,7]
		ra_hms_cont	=data_sorted[i,8]
		dec_dms_cont=data_sorted[i,9]
		ra_deg_cont	=data_sorted[i,10]
		dec_deg_cont=data_sorted[i,11]	  
		background_noise =float(data_sorted[i,37]) 
		rms_median = data_sorted[i,32]

		#re-iterating todo: change to header names rather than index (likely to break)	
	elif settings.cat_type=='component_new':
		#component catalogue being used
	
		src=data_sorted[i,8]
		ra_hms_cont	=data_sorted[i,9]
		dec_dms_cont=data_sorted[i,10]
		ra_deg_cont	=data_sorted[i,11]
		dec_deg_cont=data_sorted[i,12]	  
		background_noise =float(data_sorted[i,38]) 
		rms_median = data_sorted[i,33]

		#re-iterating todo: change to header names rather than index (likely to break)	

	else:
		raise ValueError("Catalogue type is not valid")
	overwrite = settings.overwrite #Follow global overwrite settings
	filename=overlayloc+src+overlay_suffix
	filename_cross=overlayloc+src+'_cross_'+overlay_suffix
	src_coords =SkyCoord(ra_deg_cont,dec_deg_cont,frame='fk5',unit='deg')
	"Special exception for sources in the duplicate file"
	if (type(settings.override_src)!=list)&(settings.override_src!=None):	
		if (len(settings.override_src.split("duplicate"))>1):
			print("Handling duplicate found in the same field")
			overwrite = False
			if src in duplicate_sources.keys():
				duplicate_sources[src] +=1 
			else:
				duplicate_sources[src] = 0
			
			if (settings.cat_type=='island')|(settings.selavy_convention ==False):
				add_suff = "_{}".format(duplicate_sources[src])
			if settings.selavy_convention:
				add_suff = data_sorted[i,src_ind-1][-1]
			print("Source: {} [{} --> {}]".format(src,duplicate_sources[src],data_sorted[i,src_ind-1][-1]))
		else:
			if src in keep_srcs:
				print("Found repeated component ID, first value = _0")
				if src in duplicate_sources.keys():
					duplicate_sources[src] +=1 
				else:
					duplicate_sources[src] = 0
				if duplicate_sources[src] == 0:
					add_suff = '_0'
				print("Source: {} [{} --> {}]".format(src,duplicate_sources[src],data_sorted[i,src_ind-1][-1]))
	#print(data_sorted[i,src_ind-1][-1])
	if (settings.selavy_convention)&(settings.cat_type!="island"):
		add_suff = data_sorted[i,src_ind-1][-1]
	filename = filename.split(".png")[0]+add_suff+filename.split(".png")[-1]+".png"
	filename_DSS = filename.split(settings.field_ref)[0]+"dss_"+filename.split("/")[-1]
	filename_cross = filename.split(".png")[0]+'_cross_'+filename.split(".png")[-1]+".png"
	dir_to_save= settings.radio_output.split(settings.field_ref)[0].split(settings.prefix)[-1].split("/")
	dir_to_save =list(filter(None,dir_to_save))[0]
	
	# check to see if file already exists for this source
	# aka you shouldn't be able to overwrite existing files unless you tweak this!
	excluded_source =False #Check if the source should be excluded

	if (os.path.isfile(filename) == False)&(os.path.isfile(filename_DSS) == False)|overwrite:#|(os.path.isfile(filename) == True): #Second condition is for debugging
		print(i,SB,src)
		coords=SkyCoord(ra_deg_cont,dec_deg_cont,frame='fk5',unit=u.degree)   

		ra_max=float((coords.ra/u.degree)+deg_edge) 
		ra_min=float((coords.ra/u.degree)-deg_edge)
		dec_max=float((coords.dec/u.degree)+deg_edge) 
		dec_min=float((coords.dec/u.degree)-deg_edge)
		
		x_cen,y_cen= wcs.world_to_pixel(coords)
		
		xmin=int(x_cen-npix_edge)
		xmax=int(x_cen+npix_edge)
		ymin=int(y_cen-npix_edge)
		ymax=int(y_cen+npix_edge)
		
		#print(ra_max,ra_min,dec_max,dec_min)
		"Cutouts are created here"
		# find which DES and WISE tiles to use
		indices=np.where( (DEStiles[:,2].astype(float)>=ra_min) & (DEStiles[:,1].astype(float)<=ra_max) & (DEStiles[:,3].astype(float)>=dec_min) & (DEStiles[:,4].astype(float)<=dec_max))
		DES_tiles_to_use=DEStiles[indices,0]
		DES_tiles_to_use_coords=DEStiles[indices,:]

		wiseindices=np.where( (WISEtiles[:,6].astype(float)>=ra_min) & (WISEtiles[:,8].astype(float)<=ra_max) & (WISEtiles[:,13].astype(float)>=dec_min) & (WISEtiles[:,9].astype(float)<=dec_max))
		wise_tiles_to_use=WISEtiles[wiseindices,16]

		
		# get radio cutout
		radio_cutout = Cutout2D(image, position=(x_cen,y_cen), size=(2*npix_edge), wcs=wcs, mode='trim')
		# calculate contours
		# TAKE NOTE THIS IS WHERE TO TWEAK RADIO CONTOUR LEVELS
		contourexps=np.arange(start=0,stop=32,step=0.5) 
		#use step=1 for contours doubling each time, 0.5 for a factor of root 2 etc...
		contourmults=np.power(2,contourexps)
		#basecont=3.*background_noise/1000.
		#OR
		norm_background=np.quantile(np.random.normal(scale=background_noise/1000,size=int(1e7)),0.997)
		# basecont=max(min(norm_background,float(rms_median)/1e3*background_noise,0.00012),background_noise/1000)#,float(data_sorted[i,31])/1e6)#0.00012 #Median Value
		#basecont=max(min(norm_background,0.00012,float(rms_median)/100),background_noise/1000)#,float(data_sorted[i,31])/1e6)#0.00012 #Median Value
		"This snippet will ensure the image isn't overshadowed if it's noisy"
		min_thresh = np.quantile(radio_cutout.data[~np.isnan(radio_cutout.data)],0.99)
		basecont = min(max(norm_background,background_noise/1000),min_thresh)
		#basecont = 0.00012
		# else:
			# basecont=max(min(norm_background,float(data_sorted[i,32])/1e3*background_noise,0.00012),background_noise/1000)
		radio_contours = [basecont * i for i in contourmults]
	   
		radio_max=np.nanmax(radio_cutout.data)
		cont_cond = np.flatnonzero(radio_contours < radio_max)
		if len(cont_cond)>0:
			nconts=np.nanmax(cont_cond)
		else:
			nconts = settings.cont_limit
			excluded_source = True #The data is too noisy
		#Only choose a selected number of contours
		"Restrict number of shown contours"
		"Trim the radio contours"
		radio_contours = np.asarray(radio_contours)[radio_contours <= radio_max]
		# if len(radio_contours) > settings.cont_limit:
		#	 radio_hist,radio_bins = np.histogram(radio_contours,bins=settings.cont_limit) 
		#	 radio_contours = radio_bins+ (radio_bins[1]-radio_bins[0])/2		
		if len(radio_contours)>0:
			radio_contours = np.linspace(radio_contours.min(),radio_contours.max(),settings.cont_limit)
		else:
			#radio_contours = np.linspace(radio_contours.min(),radio_contours.max(),settings.cont_limit)
			radio_contours = np.asarray([radio_cutout.data[~np.isnan(radio_cutout.data)].max(),basecont])
		#radio_contours = np.quantile(radio_cutout.data,[0.97,0.98,0.99])
		"Log 2 increments - works well but solo contours still exist"
		nconts =settings.cont_limit
		radio_contours= np.logspace(np.log2(radio_contours.min()*0.6),np.log2(radio_contours.max()),nconts,base=2)+background_noise/1000
		
		#print(radio_contours, filename)
		# =============================================================================
		#		 MASKING TEST - Comment on during actual runs
		# =============================================================================
		radio_cutout_contours=  np.array(radio_cutout.data)
		radio_cutout_contours[radio_cutout_contours<radio_contours[1]]= 0
		x_size,y_size= radio_cutout_contours.shape
		pixels_per_arcmin = x_size/arcmins
		selection_window_x = (int(x_size/2 - pixels_per_arcmin),int(x_size/2 + pixels_per_arcmin))
		selection_window_y = (int(y_size/2 - pixels_per_arcmin),int(y_size/2 + pixels_per_arcmin))
		radio_cutout_window = radio_cutout_contours[selection_window_x[0]:selection_window_x[1],
							selection_window_y[0]:selection_window_y[1]]
		# =============================================================================
		# 	Contour extraction	
		# =============================================================================
		if settings.extract_contours:
			"Create contour"
			radio_cutout_ext = np.array(radio_cutout.data) #Create duplicate array for extraction
			radio_cutout_ext[0:,0] = 0 #Set the left boundary to 0
			radio_cutout_ext[0:,-1] = 0#Set the right boundary to 0
			radio_cutout_ext[0,0:] = 0#Set the upper boundary to 0
			radio_cutout_ext[-1,0:] = 0#Set the bottom boundary to 0
			ext_conts = plt.contour(radio_cutout_ext,levels=radio_contours,colors='grey')
			if os.path.isdir(settings.extract_contours_dir)==False:
				os.mkdir(settings.extract_contours_dir)
			save_contours(ext_conts, fname = settings.extract_contours_dir+"SB"+filename.split("SB")[-1].split(".png")[0]+".h5",coords=(ra_deg_cont,dec_deg_cont)
				 ,x_size= x_size,y_size=y_size)
#%%	 # =============================================================================
		#		Create smaller masking window to filter bright sources 
		# =============================================================================
		radio_cutout_tiny = Cutout2D(image, position=(x_cen,y_cen), size=(npix_edge/6*2), wcs=wcs, mode='trim') #Creats a 2x2' cutout
		if np.isnan(np.asarray(radio_cutout_tiny.data).min()):
			excluded_source=True
			warnings.warn("Invalid values found... source is likely on the edge of the detector image")
		else:
			
			masked_tiny,excluded_source = masking(radio_cutout_tiny.data,radio_contours,mask_value=settings.mask_value)
			#print(excluded_source)
			plt.imshow(masked_tiny,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
			"Remove single contours"
			if settings.remove_single_contours:
				if excluded_source:
					masked,excluded_source= masking(radio_cutout.data,radio_contours,mask_value=settings.mask_value,exclude=False)
					#print(excluded_source)
				else:
					masked,_= masking(radio_cutout.data,radio_contours,mask_value=settings.mask_value,exclude=False)
				plt.imshow(masked,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
			

#%%
		"Check if there is a source within the masked region"
		if (len(np.where(radio_cutout_window>0)[1]) ==0)|(excluded_source==True):
			"Reset filename to exclude source and place it in another directory"
			filename=  filename.split(dir_to_save)[0]+settings.exclusion_dir.split(settings.prefix)[-1]+filename.split(dir_to_save)[1]
			filename_cross =filename.split(".")[0] +"_cross_"+".png"#filename.split(".png")[-1] 
			exclusion_list.append(src)
			excluded_source = True
			if (os.path.isfile(filename) == True)&(overwrite==False):
				"If the file exists, proceed"
				continue
		# radio_cutout.data  = radio_cutout_contours
		
		contcolors=[]
		for c in range(0,nconts+1):
			"Gradually shift the colourmap incrementally based on contour level"
			contcolors.append(greycmap(0.5+(0.5*c/(nconts+1))))  

		R_list=[]
		G_list=[]
		B_list=[]
		if settings.skip_plotting:
			pass
		else:
			for j in range(0,len(DES_tiles_to_use[0])):
				#for j in ind:
				try:
					print("Loading Rhdu")
					Rhdu=fits.open(glob.glob(settings.DESfiles_dir.split("*")[0]+DES_tiles_to_use[0][j]+'*_i.fits*')[0])
				except IndexError:
					#raise IOError("DES File not found...{}".format(DES_tiles_to_use[0][j]+'*_i.fits*'))
					print("DES File not found...{}".format(DES_tiles_to_use[0][j]+'*_i.fits*'))
					dss=True
					continue
					#break
				except FileNotFoundError:
					#raise IOError("DES File not found...{}".format(DES_tiles_to_use[0][j]+'*_i.fits*'))
					print("DES File not found...{}".format(DES_tiles_to_use[0][j]+'*_i.fits*'))
					dss=True
					continue
					#break
				except:
					print("Something has gone wrong with the DES tiles...")
					dss=True
					continue
				finally:
					print("Continuing to next step...\n")
# 				print(1)
				R=Rhdu[1].data
				des_wcs=WCS(Rhdu[1].header)
				Rhdu.close()
# 				quit()
# 				print(2)
				Ghdu=fits.open(glob.glob(settings.DESfiles_dir.split("*")[0]+DES_tiles_to_use[0][j]+'*_r.fits*')[0])
				G=Ghdu[1].data
				Ghdu.close()
# 				print(3)
				Bhdu=fits.open(glob.glob(settings.DESfiles_dir.split("*")[0]+DES_tiles_to_use[0][j]+'*_g.fits*')[0])
				B=Bhdu[1].data
				Bhdu.close()
				R_cutout=Cutout2D(R,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
				G_cutout=Cutout2D(G,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
				B_cutout=Cutout2D(B,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
# 				print(4)
				"""Note -SM [30/04/2024]: The commented out exception does work but it also filters out viable sources"""
# 				try:
# 					R_cutout=Cutout2D(R,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='strict')
# 					G_cutout=Cutout2D(G,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='strict')
# 					B_cutout=Cutout2D(B,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='strict')
# 				except PartialOverlapError:
# 					print("Partial coverage")
# 					if excluded_source ==False:
# 						filename=  filename.split(dir_to_save)[0]+settings.exclusion_dir.split(settings.prefix)[-1]+filename.split(dir_to_save)[1]
# 						filename_cross =filename.split(".")[0] +"_cross_"+filename.split(".")[-1] 
# 						exclusion_list.append(src)
# 					excluded_source = True
# 					R_cutout=Cutout2D(R,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
# 					G_cutout=Cutout2D(G,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
# 					B_cutout=Cutout2D(B,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
# 				print(5)
				R_hdu=fits.PrimaryHDU(data=R_cutout.data, header=R_cutout.wcs.to_header())
				G_hdu=fits.PrimaryHDU(data=G_cutout.data, header=G_cutout.wcs.to_header())
				B_hdu=fits.PrimaryHDU(data=B_cutout.data, header=B_cutout.wcs.to_header())
# 				print(6)
				R_list.append(R_hdu)
				G_list.append(G_hdu)
				B_list.append(B_hdu)  
				print("Image created")
			if len(R_list)==0:
				print("something wrong")
				if len(DES_tiles_to_use)==1:
					dss= True
					print("No DES Images found, using DSS2 instead...")

					#continue
			#elif len(R_list)==0:
				else:
					print("Optical image obtained")
					#only one image so no need to mosaic
					R=R_list[0].data
					G=G_list[0].data
					B=B_list[0].data
					des_wcs=WCS(R_list[0].header)
			else:
				#need to combine them
				print("Combining DES Tiles")
				des_wcs, shape_out = find_optimal_celestial_wcs(R_list)
				R, footprint_R = reproject_and_coadd(R_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
				G, footprint_G = reproject_and_coadd(G_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
				B, footprint_B = reproject_and_coadd(B_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
				#If there are more than X% nan values, filter out the DES tiles and use DSS
				if len(R[np.isnan(R)])>0.2*len(R):
					dss=True
					print("Optical image is corrupted or incomplete")
			
			if dss:
				print("Downloading DSS Alternative")
				countdown = 5
				while countdown > 0 : #Keep trying if download fails
					try: #Download all the tiles
						R_hdu= SkyView.get_images(src_coords,survey=["DSS2 Red"],coordinates='J2000',radius=12*u.arcmin)[0]
						G_hdu= SkyView.get_images(src_coords,survey=["DSS"],coordinates='J2000',radius=12*u.arcmin)[0]
						B_hdu= SkyView.get_images(src_coords,survey=["DSS2 Blue"],coordinates='J2000',radius=12*u.arcmin)[0]
						G=np.zeros_like(R)#No G band present in DSS#G_hdu[0].data
						B=B_hdu[0].data- np.quantile(B_hdu[0].data,0.5)
						R=R_hdu[0].data - np.quantile(R_hdu[0].data,0.5)
						des_wcs=WCS(R_hdu[0].header)
						break #Leave while loop
					except: #If non-red tiles are unavailable
						R_hdu= SkyView.get_images(src_coords,survey=["DSS2 Red"],coordinates='J2000',radius=12*u.arcmin)[0]
						R=R_hdu[0].data - np.quantile(R_hdu[0].data,0.5)
						des_wcs=WCS(R_hdu[0].header)
						break #Leave while loop
					countdown -= 1
				
				filename = filename.split(settings.field_ref)[0]+"dss_"+filename.split("/")[-1]
				filename_cross = filename_cross.split(settings.field_ref)[0]+"dss_"+filename_cross.split("/")[-1]
			
			if dss:
				img = R
			else:
				img=lupton_rgb.make_lupton_rgb(R,G,B,Q=10,stretch=50,minimum=1)
			"Note - SM [27/04/2024]: Some sources have partial images, these need to be flagged up."
			
			fig = plt.figure(constrained_layout=False,figsize=(1024/my_dpi, 1024/my_dpi),dpi=my_dpi)#*0.55)# <- Uncomment this if you want to downscale the images
			# Set figure background as white
			fig.patch.set_facecolor('w')		
	
			ax1=plt.subplot(331,projection=radio_cutout.wcs,fc='grey')			
			ax2=plt.subplot(332,projection=radio_cutout.wcs,fc='grey') 
			ax3=plt.subplot(333,projection=radio_cutout.wcs,fc='grey')	
			ax4=plt.subplot(334,projection=radio_cutout.wcs,fc='grey')	
			ax5=plt.subplot(335,projection=radio_cutout.wcs,fc='grey')	
			ax6=plt.subplot(336,projection=radio_cutout.wcs,fc='grey')	
			ax7=plt.subplot(337,projection=radio_cutout.wcs,fc='grey')	
			ax8=plt.subplot(338,projection=radio_cutout.wcs,fc='grey')	
			ax9=plt.subplot(339,projection=radio_cutout.wcs,fc='grey')	
	
			ax1.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
			ax1.contour(radio_cutout.data,levels=radio_contours,colors='grey')  
			ax2.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
			ax2.contour(radio_cutout.data,levels=radio_contours,colors='grey')
			ax3.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
			ax3.contour(radio_cutout.data,levels=radio_contours,colors='grey')
			if dss ==False: #If DSS is not used, plot the DES image
				ax4.imshow(img,transform=ax4.get_transform(des_wcs),origin='lower') 
				ax5.imshow(img,transform=ax5.get_transform(des_wcs),origin='lower') 
				ax6.imshow(img,transform=ax6.get_transform(des_wcs),origin='lower')
			else:
				ax4.imshow(img,transform=ax4.get_transform(des_wcs),origin='lower',vmin=np.quantile(img,0.1),vmax=np.quantile(img,0.99))#,cmap=plt.cm.binary) 
				ax5.imshow(img,transform=ax5.get_transform(des_wcs),origin='lower',vmin=np.quantile(img,0.1),vmax=np.quantile(img,0.99))#,cmap=plt.cm.binary) 
				ax6.imshow(img,transform=ax6.get_transform(des_wcs),origin='lower',vmin=np.quantile(img,0.1),vmax=np.quantile(img,0.99))#,cmap=plt.cm.binary)
			ax4.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
			ax5.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
			ax6.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
			
			try:
				wise_list=[]
				for k in range(0,len(wise_tiles_to_use[0])):
					#get wise cutout 
					download_wise= False#Should wise images be downloaded
					wise_im=settings.WISEfiles_dir.split("*")[0]+str(wise_tiles_to_use[0][k])+'-w1-int-3.fits'
					if os.path.isfile(wise_im)==False:
						print("WISE Tile not found... moving on")
						try:
							print("Attempting to download WISE image")
							#wise_im = SkyView.get_images(src_coords,survey=["WISE 3.4"],coordinates='J2000',radius=12*u.arcmin)
							download_wise =  True
						except:
							print("Issue encountered... moving on")
							continue
					if download_wise==False:
						wise_hdu=fits.open(wise_im)
					else:
						#radio_coords =SkyCoord(radio_cutout.wcs.wcs.crval[0],radio_cutout.wcs.wcs.crval[1],frame='fk5',unit='deg')
						#wise_hdu= SkyView.get_images(radio_coords,survey=["WISE 3.4"],coordinates='J2000',radius=12*u.arcmin)[0]
						wise_hdu= SkyView.get_images(src_coords,survey=["WISE 3.4"],coordinates='J2000',radius=12*u.arcmin)[0]
										
					wise_data= wise_hdu[0].data
					wise_wcs= WCS(wise_hdu[0].header)
					wise_hdu.close()
					wise_cutout=Cutout2D(wise_data,position=coords,size=1.05*arcmins*u.arcmin,wcs=wise_wcs,mode='trim')
					wise_hdu=fits.PrimaryHDU(data=wise_cutout.data, header=wise_cutout.wcs.to_header())			   
					wise_list.append(wise_hdu)
				if len(wise_list)==0:
					print("something wrong")
				elif len(wise_list)==1:
					print("only one wise")
					#only one image so no need to mosaic
					wise_data=wise_list[0].data
					wise_wcs=WCS(wise_list[0].header)
					#print("Skipping source: {}".format(src))
					#continue
	
				else:
					print("combining wise")
					#need to combine them
					wise_wcs, wise_shape_out = find_optimal_celestial_wcs(wise_list)
					wise_data, footprint_wise = reproject_and_coadd(wise_list,wise_wcs,shape_out=wise_shape_out,reproject_function=reproject_interp)
	
				#wise_cutout_hdu=fits.PrimaryHDU(data=wise_data, header=wise_wcs.to_header()) 
				#wise_cutout_hdu.writeto(filename)
				if download_wise==False:
					ax7.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax7.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
					ax8.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax8.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
					ax9.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax9.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
				else:
					ax7.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax7.get_transform(wise_wcs),origin='lower',cmap=gist_heat,vmin=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.5),vmax=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.99))
					ax8.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax8.get_transform(wise_wcs),origin='lower',cmap=gist_heat,vmin=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.5),vmax=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.99))
					ax9.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax9.get_transform(wise_wcs),origin='lower',cmap=gist_heat,vmin=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.5),vmax=np.quantile(ashinh_scale(wise_data,zeropoint=2,scale=100),0.99))					
					
			except:
				print("wise failed")
				for k in range(0,len(wise_tiles_to_use[0])):
					#get wise cutout 
					wise_im=settings.WISEfiles_dir.split("*")[0]+str(wise_tiles_to_use[0][k])+'-w1-int-3.fits'
					if os.path.isfile(wise_im)==False:
						print("WISE Tile not found... moving on")
						continue
					wise_hdu=fits.open(wise_im)
					wise_data= wise_hdu[0].data
					wise_wcs= WCS(wise_hdu[0].header)
					wise_hdu.close()
					ax7.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax7.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
					ax8.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax8.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
					ax9.imshow(ashinh_scale(wise_data,zeropoint=2,scale=100),transform=ax9.get_transform(wise_wcs),origin='lower',cmap=gist_heat)
# 			except FileNotFoundError:
# 				print("WISE Tiles not found... moving on")
# 				continue
	
			ax7.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
			ax8.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
			ax9.contour(radio_cutout.data,levels=radio_contours,colors=contcolors)
	
	
			ax1.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax1.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax2.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax2.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax3.set_xlim(0,(2*npix_edge)-1)
			ax3.set_ylim(0,(2*npix_edge)-1)
			
			ax4.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax4.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax5.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax5.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax6.set_xlim(0,(2*npix_edge)-1)
			ax6.set_ylim(0,(2*npix_edge)-1)
			
			ax7.set_xlim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax7.set_ylim(0.75*npix_edge,(1.25*npix_edge)-1)
			ax8.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax8.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
			ax9.set_xlim(0,(2*npix_edge)-1)
			ax9.set_ylim(0,(2*npix_edge)-1)
	
	
			ax1.axis('off')
			ax2.axis('off')
			ax3.axis('off')
			ax4.axis('off')
			ax5.axis('off')
			ax6.axis('off')
			ax7.axis('off')
			ax8.axis('off')
			ax9.axis('off')
	
			plt.rc('font', size=10)	
			plt.annotate('Radio',xy=(0.08,0.75),xycoords='figure fraction',ha='center',va='center',rotation=90)
			plt.annotate('Optical',xy=(0.08,0.45),xycoords='figure fraction',ha='center',va='center',rotation=90)
			plt.annotate('Infrared',xy=(0.08,0.15),xycoords='figure fraction',ha='center',va='center',rotation=90)
			plt.annotate('Zoomed in',xy=(0.25,0.91),xycoords='figure fraction',ha='center')
			plt.annotate('Default',xy=(0.55,0.91),xycoords='figure fraction',ha='center')
			plt.annotate('Zoomed out',xy=(0.85,0.91),xycoords='figure fraction',ha='center')
			plt.subplots_adjust(left=0.1, bottom=0.01, right=0.99, top=0.9, hspace=0,wspace=0.02) # control space between figure and whitespace
			try:
				plt.savefig(filename,dpi=my_dpi)
				#plt.savefig(filename.split("png")[0]+"svg",dpi=my_dpi)
			except:
				print("Something went wrong when saving source figure: {}".format(src))
				continue
			
			xp,yp=utils.skycoord_to_pixel(coords,radio_cutout.wcs)
	
			#plt.rcParams.update({'lines.linewidth':1.2})
	
			ax1.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax1.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax1.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax1.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			ax2.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax2.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax2.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax2.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			ax3.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax3.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax3.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax3.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			ax4.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax4.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax4.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax4.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			ax5.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax5.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax5.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax5.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			ax6.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax6.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax6.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax6.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
			
			ax7.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax7.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax7.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax7.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
			
			ax8.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax8.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax8.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax8.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
			
			ax9.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
			ax9.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
			ax9.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
			ax9.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
	
			plt.savefig(filename_cross,dpi=my_dpi)
			#plt.savefig(filename_cross.split("png")[0]+"svg",dpi=my_dpi)
			plt.show()
			plt.close()
#%%
		# =============================================================================
		# Create standlone cutouts
		# =============================================================================
		if settings.create_cutout:
			print("Creating 6x6 cutout for source: {}".format(src))
			#Create file
			cutout_suffix = ""
			cutout_filename =settings.cutout_dir+"SB"+settings.SB+src+cutout_suffix+add_suff
			cutout_dir = settings.cutout_dir
			if excluded_source:
				"Dissect and append to the filename"
				cutout_filename = cutout_filename.split(settings.cutout_dir.split(settings.prefix)[-1])[0] +"excluded_" + settings.cutout_dir.split(settings.prefix)[-1] + cutout_filename.split(settings.cutout_dir.split(settings.prefix)[-1])[-1]
				cutout_dir=  cutout_filename.split(settings.field_ref)[0]
			ut.make_dir(cutout_dir) #Checks if directory already exists
			if settings.export_fits_cutout:
				scale_override= 1 #Should be 1 unless testing
				radio_cutout_small = Cutout2D(image, position=(x_cen,y_cen), size=(scale_override*npix_edge), wcs=wcs, mode='trim')
				hdu[0].data = radio_cutout_small.data
				hdu[0].header.update(radio_cutout_small.wcs.to_header())
				hdu.writeto(cutout_filename+".fits", overwrite=True)
			"Create cutout figure"
			if settings.skip_cutout_plotting ==False:
				fig = plt.figure(constrained_layout=False,figsize=(1024/my_dpi,1024/my_dpi))
				ax = plt.subplot(111,projection=radio_cutout.wcs,fc='grey')
				ax.imshow(radio_cutout.data,origin='lower',cmap=magmacmap,norm=colors.LogNorm(vmin=basecont/5, vmax=radio_max))
				ax.contour(radio_cutout.data,levels=radio_contours,colors='grey')
				ax.axis('off')
				ax.set_axis_off()
				ax.set_xlim(0.5*npix_edge,(1.5*npix_edge)-1)
				ax.set_ylim(0.5*npix_edge,(1.5*npix_edge)-1)
				if settings.use_cross:
					xp,yp=utils.skycoord_to_pixel(coords,radio_cutout.wcs)
					cutout_suffix = "_cross"
					ax.axhline(y=yp,xmin=0,xmax=0.45,c='w',linestyle=':')
					ax.axhline(y=yp,xmin=0.55,xmax=1,c='w',linestyle=':')
					ax.axvline(x=xp,ymin=0,ymax=0.45,c='w',linestyle=':')
					ax.axvline(x=xp,ymin=0.55,ymax=1,c='w',linestyle=':')
				plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
				plt.savefig(cutout_filename+".png",pad_inches=0)
				#plt.savefig(cutout_filename.split("png")[0]+"svg",dpi=my_dpi)
				plt.close(fig)

hdu.close()   #Close radio fits file
#%%
"Use numpy to save sources in the exclusion list - Removed temporarily"
#if len(exclusion_list)>0:
#	np.savetxt(settings.exclusion_list,exclusion_list)
#%%

# =============================================================================
# THE FOLLOWING IS INOPERABLE WITHOUT A PANDAS INSTALLATION
# =============================================================================
# print("Running clean-up and exclusion routine...\n")
# # =============================================================================
# # Creating exclusion list
# # =============================================================================
# if len(exclusion_list)>0: #If there are sources in the exclusion list
#	 table = pd.DataFrame(exclusion_list,columns=["Source_ID"]) #Create a new table
#	 table["Field"] = np.tile(settings.SB,len(exclusion_list)) #Assign field ID
#	 if os.path.isfile(settings.exclusion_list): #If file exists
#		 temp_table = pd.read_csv(settings.exclusion_list) #Load in existing file
#		 table= temp_table.concatenate(table) #Append current table to existing source list
#	 table.drop_duplicates()
#	 table.to_csv(settings.exclusion_list,index=False) #Save file

# # =============================================================================
# # Move excluded sources to their own designated folder
# # =============================================================================
#	 source_filenames= glob.glob(settings.radio_output+"*") #List of cutout files
	
#	 "Create a list of source files that have been created"
#	 source_list = list(map(ft.partial(ut.extract_sources,field_split=settings.SB),source_filenames)) #Extract source names from files
	
#	 "Find indices corresponding to the matched file location"
#	 _,match_index,_=np.intersect1d(source_list,table["Source_ID"],return_indices=True)
#	 excluded_source_list = source_list[match_index] #List of excluded source filenames
#	 excluded_source_filenames= source_filenames[match_index] #Source filenames
	
#	 "Move files"
	
#	 map(ft.partial(ut.extract_sources,new_dir=settings.exclusion_dir),excluded_source_filenames) #Extract source names from files
# print("Operation completed...\n")