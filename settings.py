# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:30:20 2024

Initialise directories for testing among other variable components.

@author: Soheb Mandhai
"""


# =============================================================================
# Field Settings
# =============================================================================
SB='9442'
component_number = 86
cat_type="component_new"#"island"#True#True
# =============================================================================
# Directories
# =============================================================================

#INPUTS
prefix= "../sample_data/des/"#"/mnt/shared/"
dataloc=prefix+"askap_data/"#'/mnt/shared/des/askap_data/' #askap data location
WISEfiles_dir=prefix+'WISEtiles/*.fits' #WISE raw data tiles location
DESfiles_dir=prefix + 'DEStiles/' #WISE raw data tiles location

if cat_type=='island':
    cat_sub = "catalogue/AS101_Continuum_Island_Catalogue_{}.csv"#'catalogues/AS101_Continuum_Component_Catalogue_{}_78.csv'
elif cat_type=='catwise':
    cat_sub= "catalogue/emu_with_catwise2020_desdr1_scos_desi_nocutoff.csv"
elif cat_type== "component":
    cat_sub = "catalogue/AS101_Continuum_Component_Catalogue_{}_{}.csv"#'catalogues/AS101_Continuum_Component_Catalogue_{}_78.csv'
elif cat_type== "component_new":
    cat_sub = "catalogues/AS101_Continuum_Component_Catalogue_{}_{}.csv"#'catalogues/AS101_Continuum_Component_Catalogue_{}_78.csv'
else:
    print("Valid catalogue filetype not specified in settings.")
WISEtiles=  'WISE_tiles.txt'
DEStiles= 'DES_tiles_dims.txt'

#OUTPUT
radio_output=prefix+'radio_cutouts/SB'
overlayloc = prefix+'radio_cutouts/SB'


"For excluded sources"
exclusion_list = prefix+"excluded_sources.csv"
exclusion_dir = prefix+"excluded_radio_cutouts/"

# =============================================================================
# Source Settings
# =============================================================================
use_file= False# If True, set override_src to the file name of the source list
override_src = ["J202715-553043","J202505-540405","J210240-552627","J202710-530638","J202254-540537"]#"J202505-540405"#'J202254-540537' #Default =None
override_src = prefix+'beta_test_subject_ids.txt'
override_src = None
override_src= "J210102-582918"
#Note: If you are using a file, ensure the source names are in the first column


# =============================================================================
# Contour Settings
# =============================================================================
cont_limit = 8

# =============================================================================
# Additional Options
# =============================================================================
create_cutout = True #Create cutout of default view?
export_fits_cutout=True #Create fits cutout?
cutout_dir = prefix+'radio_cutouts_6x6/'
use_cross = False #Put a cross_hair on 6x6 cutout sources
skip_plotting = True #Default: False
skip_cutout_plotting = False #Default: False
# =============================================================================
# ADVANCED OPTIONS - DO NOT CHANGE UNLESS YOU KNOW WHAT YOU'RE DOING
# =============================================================================
field_ref = 'SB' #Reference to field