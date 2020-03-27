# fuv_age.github.io

FUV as age indicator

full_analysis.py:
Combines all of my Jupyter notebooks such that I output all plots needed for the paper. 


sierchio_fuv.py:
Takes in Sierchio et al. 2014 catalog and uses astroquery to pull FUV magnitudes and errors from GALEX catlog. 

isaacson_fuv.py:
Takes in Isaacson an Fischer 2010 catalog and uses astroquery to pull FUV magnitudes and errors form GALEX. Note that it only pulls stars with HIP numbers. The rest I pulled seperately. 

ballering_fuv.py:
Takes in Ballering et al 2013 catalog and uses astroquery to pull FUV magnitudes and errors from GALEX.

fit_stars_only.py:
Fits age vs Q for just the field stars.

import_data.py:
Imports stars and the properties from Isaacson, Ballering, Sierchio, and Solar Twin catalogs

age_vs_rhk.py:
Plots age vs. RHK to demonstate age-activity relations

fuv_vs_b.py:
Plots two-color diagram of FUV-B vs B-V using the stars form for catalogs but does not use duplicate stars

mag_cuts.py:
Reduces the sample to stars within 4.3<Mv<5.3 (ie those with solar-like magnitudes.

some_mags.py:
Finds B-V, FUV-B and Q for all catalogs

individual_fits.py:
Fits the four stellar samples individually to determine outliers and then collects all non-outlier stars into arrays.
