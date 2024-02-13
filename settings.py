# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:30:20 2024

Initialise directories for testing among other variable components.

@author: Soheb Mandhai
"""

# =============================================================================
# Source Settings
# =============================================================================

override_src = ["J210240-552627","J202710-530638","J202505-540405","J202254-540537","J202715-553043"]#"J202505-540405"#'J202254-540537' #Default =None

# =============================================================================
# Field Settings
# =============================================================================
SB='9351'
island=True
# =============================================================================
# Directories
# =============================================================================

prefix= "../sample_data/"#"/mnt/shared/"
dataloc=prefix+"des/askap_data/"#'/mnt/shared/des/askap_data/' #askap data location
WISEfiles_dir=prefix+'des/WISEtiles/*.fits' #WISE raw data tiles location
DESfiles_dir=prefix + 'des/DEStiles/' #WISE raw data tiles location

if island:
    cat_sub = "catalogue/AS101_Continuum_Island_Catalogue_{}.csv"#'catalogues/AS101_Continuum_Component_Catalogue_{}_78.csv'
else:
    cat_sub = "catalogue/AS101_Continuum_Component_Catalogue_{}_78.csv"#'catalogues/AS101_Continuum_Component_Catalogue_{}_78.csv'

radio_output=prefix+'des/radio_cutouts/SB'
overlayloc = prefix+'des/radio_cutouts/SB'
WISEtiles=  'WISE_tiles.txt'
DEStiles= 'DES_tiles_dims.txt'

# =============================================================================
# Contour Settings
# =============================================================================
cont_limit = 6