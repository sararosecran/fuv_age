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
