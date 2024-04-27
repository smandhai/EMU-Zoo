# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:30:20 2024

Initialise directories for testing among other variable components.

@author: Soheb Mandhai
"""


# =============================================================================
# Field Settings
# =============================================================================
SB="9351"#"9351"#"9434"#'9442'#"9351"#
component_number="78"#78# 90#86#78# 86
# SB="9351"#"9434"#"9351"#"9434"#'9442'#"9351"#""
# component_number=78
cat_type="component_new"#"component_new"#"component_new"#"island"#True#True
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
radio_output=prefix+'phase1_radio_cutouts/SB'
overlayloc = prefix+'phase1_radio_cutouts/SB'


"For excluded sources"
exclusion_list = prefix+"excluded_sources.csv"
exclusion_dir = prefix+"excluded_radio_cutouts/"

# =============================================================================
# Source Settings
# =============================================================================
use_file= True# If True, set override_src to the file name of the source list
#override_src = ["J202715-553043","J202505-540405","J210240-552627","J202710-530638","J202254-540537"]#"J202505-540405"#'J202254-540537' #Default =None
override_src = prefix+'beta_test_subject_ids.txt'
#override_src = None
#override_src= "J203605-532105"#"J202710-530638"# "J203638-532628"#"J210102-582918"
#Note: If you are using a file, ensure the source names are in the first column
# override_src= ["J211525-531726","J211547-532405",]#"J211507-501139"##9434
# override_src= ["J204554-510432","J204221-494947","J210524-531430","J204539-505452","J205505-502604","J211447-510910"]#"J211507-501139"##9287
# override_src= ["J210337-490023","J204310-510020","J210951-531001","J210957-531003","J211536-485429"]#"J211507-501139"##9287
# override_src = "J211229-533041"#"J205308-535913"#"J205913-481721"
override_src = "duplicate_sources.dat" #<- Validation for duplicates
override_src = prefix + "complexity_sources/SB{}_subject_list_PhaseI_RGZ-EMU_cpx_ab4000_majax_ab20.csv".format(SB)
#override_src = prefix + "complexity_sources/subject_list_PhaseI_RGZ-EMU_cpx_ab4000_majax_ab20.csv".format(SB)
# =============================================================================
# Contour Settings
# =============================================================================
cont_limit = 8

# =============================================================================
# Additional Options
# =============================================================================
remove_single_contours = True #Remove single contours outside of target?
create_cutout = True #Create cutout of default view?
export_fits_cutout=True #Create fits cutout?
cutout_dir = prefix+'phase1_radio_cutouts_6x6/'
use_cross = True #Put a cross_hair on 6x6 cutout sources
skip_plotting = True #Default: False
skip_cutout_plotting = False #Default: False
mask_value = 100 #Masked value
duplicate_sources = "duplicate_sources.dat"
selavy_convention = True#Uses the selavy lettering convention BUT needs to have this in the catalogue file, column before component names
# =============================================================================
# ADVANCED OPTIONS - DO NOT CHANGE UNLESS YOU KNOW WHAT YOU'RE DOING
# =============================================================================
field_ref = 'SB' #Reference to field
overwrite = True #Overwrite existing files, default=False