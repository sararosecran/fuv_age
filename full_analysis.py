#This code contains everything needed to get age vs FUV fits
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import import_data
import age_vs_rhk
import fuv_vs_b
import mag_cuts
import some_mags
import individual_fits
import import_metallicities
import fit_stars_only
import fit_mg
import fit_stars_mg_clusters
import fuv_error
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#FUV-B vs B-V plot
os.chdir("txt")
#Import  Isaacson
i_hd, i_bmag, i_fuv, i_vmag, i_age, i_parallax = import_data.import_isaacson()
#Import Ballering
b_hd, b_vmag, b_bmag, b_fuv, b_age, b_parallax, b_rhk = import_data.import_ballering()
#plot Ballering age vs RHK
age_vs_rhk.rhk_plot(b_age, b_rhk)
#Import Sierchio
s_hd, s_vmag, s_bmag, s_fuv, s_age, s_parallax, s_fuv = import_data.import_sierchio()
#Import Solar Twin
twin_hd, twin_vmag, twin_bmag, twin_fuv, twin_age, twin_parallax = import_data.import_solar_twin()
#plot FUV-B vs B-V
fuv_vs_b.two_color_plot(i_hd, i_vmag, i_bmag, i_fuv,s_hd, s_vmag, s_bmag, s_fuv, b_hd, b_vmag, b_bmag, b_fuv, twin_hd, twin_vmag, twin_bmag, twin_fuv)
#absolute magnitude cuts
i_parallax, i_vmag,i_age, i_bmag, i_fuv, i_hd, b_parallax, b_vmag, b_age, b_bmag, b_fuv, b_hd, s_parallax, s_vmag, s_age, s_bmag, s_fuv, s_hd, twin_parallax, twin_vmag, twin_age, twin_bmag, twin_fuv, twin_hd = mag_cuts.mag_cut(i_parallax, i_vmag,i_age, i_bmag, i_fuv, i_hd, b_parallax, b_vmag, b_age, b_bmag, b_fuv, b_hd, s_parallax, s_vmag, s_age, s_bmag, s_fuv, s_hd, twin_parallax, twin_vmag, twin_age, twin_bmag, twin_fuv, twin_hd)
#Find B-V
i_bv, b_bv, s_bv, twin_bv = some_mags.b_v(i_bmag, i_vmag, b_bmag, b_vmag, s_bmag, s_vmag, twin_bmag, twin_vmag)
#find FUV-B
i_fuvb, b_fuvb, s_fuvb, twin_fuvb = some_mags.fuv_b(i_bmag, i_fuv, b_bmag, b_fuv, s_bmag, s_fuv, twin_bmag, twin_fuv)
#Find Q
i_Q, b_Q, s_Q, twin_Q = some_mags.Q(i_bmag, i_fuvb, i_bv, b_bmag, b_fuvb, b_bv, s_bmag, s_fuvb, s_bv, twin_bmag, twin_fuvb, twin_bv)
#Individual Field Star Sample Fits
all_field_x, all_field_y, all_field_bv, all_field_hd = individual_fits.combined_data(i_bv, i_Q, i_age, i_hd, i_fuvb, i_vmag, b_bv, b_Q, b_age, b_hd, b_fuvb, b_vmag, s_bv, s_Q, s_age, s_hd, s_fuvb, s_vmag, twin_bv, twin_Q, twin_age, twin_hd, twin_fuvb, twin_vmag)
#Get metallicities from Casagrande 2011
metal = import_metallicities.feH()
# Fit age vs Q for combine field star sample excluding outliers
fit_stars_only.field_stars(all_field_bv, all_field_x, all_field_y, metal)
# Fit age vs Q for moving groups only
fit_mg.moving_groups()
# Fit star, moving groups, and clusters
fit_stars_mg_clusters.fit_all(all_field_x, all_field_y, all_field_bv, metal)
# test how FUV contributes to error
fuv_error.fuv_test()

