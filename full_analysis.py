#This code contains everything needed to get age vs FUV fits
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
import math
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#FUV-B vs B-V plot
os.chdir("txt")
#Import  Isaacson
f = open('isaacson_new_hd.txt')
line = f.readlines()[1:]
f.close()
i_hd = np.array([])
i_bmag = np.array([])
i_vmag = np.array([])
i_fuv = np.array([])
i_age = np.array([])
i_parallax = np.array([])
for i in range(len(line)):
	i_hd = np.append(i_hd, str(line[i].split()[0]))
	i_bmag = np.append(i_bmag, str(line[i].split()[1]))
	i_vmag = np.append(i_vmag, str(line[i].split()[2]))
	i_age = np.append(i_age, str(line[i].split()[4]))
	i_fuv = np.append(i_fuv, str(line[i].split()[5]))
	i_parallax = np.append(i_parallax, float(line[i].split()[7]))
ignore = np.argwhere(i_fuv == 'n/a') #don't use data without galex fuv
i_hd = np.delete(i_hd, ignore)
i_bmag = np.delete(i_bmag, ignore)
i_fuv= np.delete(i_fuv, ignore)
i_vmag = np.delete(i_vmag, ignore)
i_age = np.delete(i_age, ignore)
i_parallax = np.delete(i_parallax,ignore)
ignore = np.argwhere(i_age == 'n/a') #don't use data without ages
i_hd = np.delete(i_hd, ignore)
i_bmag = np.delete(i_bmag, ignore)
i_fuv= np.delete(i_fuv, ignore)
i_vmag = np.delete(i_vmag, ignore)
i_age = np.delete(i_age, ignore)
i_parallax = np.delete(i_parallax,ignore)
ignore = np.argwhere(i_bmag == 'n/a') #don't use data without bmags
i_hd = np.delete(i_hd, ignore)
i_bmag = np.delete(i_bmag, ignore)
i_fuv= np.delete(i_fuv, ignore)
i_vmag = np.delete(i_vmag, ignore)
i_age = np.delete(i_age, ignore)
i_parallax = np.delete(i_parallax,ignore)
i_age = i_age.astype(np.float)
i_fuv = i_fuv.astype(np.float)
i_vmag = i_vmag.astype(np.float)
i_bmag = i_bmag.astype(np.float)
#Import Ballering
f = open('ballering_new_simbad.txt')
line = f.readlines()[1:]
f.close()
b_hd = np.array([])
b_bmag = np.array([])
b_vmag = np.array([])
b_fuv = np.array([])
b_age = np.array([])
b_parallax = np.array([])
b_rhk = np.array([])
for i in range(len(line)):
	b_hd = np.append(b_hd, str(line[i].split()[1]))
	b_age = np.append(b_age, str(line[i].split()[2]))
	b_bmag = np.append(b_bmag, str(line[i].split()[3]))
	b_vmag = np.append(b_vmag, str(line[i].split()[4]))
	b_fuv = np.append(b_fuv, str(line[i].split()[5]))
	b_rhk = np.append(b_rhk, str(line[i].split()[7]))
	b_parallax = np.append(b_parallax, float(line[i].split()[8]))
ignore = np.argwhere(b_fuv == 'n/a') #don't use data without galex fuv
b_hd = np.delete(b_hd, ignore)
b_vmag = np.delete(b_vmag, ignore)
b_bmag = np.delete(b_bmag, ignore)
b_fuv = np.delete(b_fuv, ignore)
b_age = np.delete(b_age, ignore)
b_parallax = np.delete(b_parallax, ignore)
b_rhk = np.delete(b_rhk, ignore)
ignore = np.argwhere(b_bmag == 'n/a') #don't use data without bmag
b_hd = np.delete(b_hd, ignore)
b_vmag = np.delete(b_vmag, ignore)
b_bmag = np.delete(b_bmag, ignore)
b_fuv = np.delete(b_fuv, ignore)
b_age = np.delete(b_age, ignore)
b_parallax = np.delete(b_parallax, ignore)
b_rhk = np.delete(b_rhk, ignore)
ignore = np.argwhere(b_age == 'n/a') #don't use data without ages
b_hd = np.delete(b_hd, ignore)
b_vmag = np.delete(b_vmag, ignore)
b_bmag = np.delete(b_bmag, ignore)
b_fuv = np.delete(b_fuv, ignore)
b_age = np.delete(b_age, ignore)
b_parallax = np.delete(b_parallax, ignore)
b_rhk = np.delete(b_rhk, ignore)
b_vmag = b_vmag.astype(np.float)
b_fuv = b_fuv.astype(np.float)
b_bmag = b_bmag.astype(np.float)
b_age = b_age.astype(np.float)
#plot Ballering age vs RHK
b_age_copy = np.copy(b_age)
ignore = np.argwhere(b_rhk == 'n/a')
b_age_copy = np.delete(b_age_copy, ignore)
b_rhk = np.delete(b_rhk, ignore)
b_rhk = b_rhk.astype(np.float)
fig = plt.figure(figsize=[10,7])
plt.scatter(b_rhk, b_age_copy)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('log$R^{\prime}_{HK}$', fontsize = 26)
plt.ylabel('Age (Gyr)', fontsize = 26)
plt.annotate('Age Versus Ca II Emission', xy=(0.6, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Ballering et al. (2013)', xy=(0.6, 0.90), xycoords='axes fraction', fontsize = 20)
#plt.show()
fig.savefig("f1a.png", bbox_inches='tight', dpi=300)
#Import Sierchio
f = open('sierchio_new_simbad.txt')
line = f.readlines()[1:]
f.close()
s_hd = np.array([])
s_bmag = np.array([])
s_vmag = np.array([])
s_fuv = np.array([])
s_age = np.array([])
s_parallax = np.array([])
for i in range(len(line)):
	s_hd = np.append(s_hd, str(line[i].split()[0]))
	s_bmag = np.append(s_bmag, str(line[i].split()[2]))
	s_vmag = np.append(s_vmag, str(line[i].split()[3]))
	s_age = np.append(s_age, str(line[i].split()[4]))
	s_fuv = np.append(s_fuv, str(line[i].split()[5]))
	s_parallax = np.append(s_parallax, float(line[i].split()[8]))
ignore = np.argwhere(s_fuv == 'n/a') #don't use data without galex fuv
s_hd = np.delete(s_hd, ignore)
s_vmag = np.delete(s_vmag, ignore)
s_bmag = np.delete(s_bmag, ignore)
s_fuv = np.delete(s_fuv, ignore)
s_age = np.delete(s_age, ignore)
s_parallax = np.delete(s_parallax, ignore)
ignore = np.argwhere(s_bmag== 'n/a') #don't use data without bmag
s_hd = np.delete(s_hd, ignore)
s_vmag = np.delete(s_vmag, ignore)
s_bmag = np.delete(s_bmag, ignore)
s_fuv = np.delete(s_fuv, ignore)
s_age = np.delete(s_age, ignore)
s_parallax = np.delete(s_parallax, ignore)
ignore = np.argwhere(s_age== 'n/a') #don't use data without ages
s_hd = np.delete(s_hd, ignore)
s_vmag = np.delete(s_vmag, ignore)
s_bmag = np.delete(s_bmag, ignore)
s_fuv = np.delete(s_fuv, ignore)
s_age = np.delete(s_age, ignore)
s_parallax = np.delete(s_parallax, ignore)
s_fuv = s_fuv.astype(np.float)
s_vmag = s_vmag.astype(np.float)
s_bmag = s_bmag.astype(np.float)
s_age = s_age.astype(np.float)
s_age = s_age / 1000.
#Import Solar Twin
f = open('solar_twin.txt')
line = f.readlines()[1:]
f.close()
twin_vmag = np.array([])
twin_bmag = np.array([])
twin_fuv = np.array([])
twin_hd = np.array([])
twin_age = np.array([])
twin_parallax = np.array([])
for i in range(len(line)):
	twin_hd = np.append(twin_hd, str(line[i].split()[1]))
	twin_bmag = np.append(twin_bmag, str(line[i].split()[2]))
	twin_vmag = np.append(twin_vmag, str(line[i].split()[3]))
	twin_age = np.append(twin_age, float(line[i].split()[5]))
	twin_fuv = np.append(twin_fuv, str(line[i].split()[8]))
	twin_parallax = np.append(twin_parallax, float(line[i].split()[10]))
ignore = np.argwhere(twin_fuv == 'n/a') #don't use data without galex fuv
twin_bmag = np.delete(twin_bmag, ignore)
twin_vmag = np.delete(twin_vmag, ignore)
twin_fuv = np.delete(twin_fuv, ignore)
twin_hd = np.delete(twin_hd, ignore)
twin_age = np.delete(twin_age, ignore)
twin_parallax = np.delete(twin_parallax, ignore)
ignore = np.argwhere(twin_bmag == 'n/a') #don't use data without bmags
twin_bmag = np.delete(twin_bmag, ignore)
twin_vmag = np.delete(twin_vmag, ignore)
twin_fuv = np.delete(twin_fuv, ignore)
twin_hd = np.delete(twin_hd, ignore)
twin_age = np.delete(twin_age, ignore)
twin_parallax = np.delete(twin_parallax, ignore)
twin_bmag = twin_bmag.astype(np.float)
twin_vmag = twin_vmag.astype(np.float)
twin_fuv = twin_fuv.astype(np.float)
#get rid of duplicates
full_cat_hd = np.copy(i_hd)
full_cat_vmag = np.copy(i_vmag)
full_cat_bmag = np.copy(i_bmag)
full_cat_fuv = np.copy(i_fuv)
hd_copy = np.copy(s_hd)
vmag_copy = np.copy(s_vmag)
bmag_copy = np.copy(s_bmag)
fuv_copy = np.copy(s_fuv)
ignore = np.array([])
for i in range(len(s_hd)):
    for j in range(len(full_cat_hd)):
        if (full_cat_hd[j] == s_hd[i]):
            ignore = np.append(ignore,i)
hd_copy = np.delete(hd_copy, ignore) #get rid of copies in Sierchio
vmag_copy = np.delete(vmag_copy, ignore)
bmag_copy = np.delete(bmag_copy, ignore)
fuv_copy = np.delete(fuv_copy,ignore)
full_cat_hd = np.append(full_cat_hd, hd_copy)
full_cat_vmag = np.append(full_cat_vmag, vmag_copy)
full_cat_bmag = np.append(full_cat_bmag, bmag_copy)
full_cat_fuv = np.append(full_cat_fuv, fuv_copy)
hd_copy = np.copy(b_hd)
vmag_copy = np.copy(b_vmag)
bmag_copy = np.copy(b_bmag)
fuv_copy = np.copy(b_fuv)
ignore = np.array([])
for i in range(len(b_hd)):
    for j in range(len(full_cat_hd)):
        if (full_cat_hd[j] == b_hd[i]):
            ignore = np.append(ignore,i)
hd_copy = np.delete(hd_copy, ignore) #get rid of copies in Ballering
vmag_copy = np.delete(vmag_copy, ignore)
bmag_copy = np.delete(bmag_copy, ignore)
fuv_copy = np.delete(fuv_copy,ignore)
full_cat_hd = np.append(full_cat_hd, hd_copy)
full_cat_vmag = np.append(full_cat_vmag, vmag_copy)
full_cat_bmag = np.append(full_cat_bmag, bmag_copy)
full_cat_fuv = np.append(full_cat_fuv, fuv_copy)
hd_copy = np.copy(twin_hd)
vmag_copy = np.copy(twin_vmag)
bmag_copy = np.copy(twin_bmag)
fuv_copy = np.copy(twin_fuv)
ignore = np.array([])
for i in range(len(twin_hd)):
    for j in range(len(full_cat_hd)):
        if (full_cat_hd[j] == twin_hd[i]):
            ignore = np.append(ignore,i)
hd_copy = np.delete(hd_copy, ignore) #get rid of copies in Solar Twin
vmag_copy = np.delete(vmag_copy, ignore)
bmag_copy = np.delete(bmag_copy, ignore)
fuv_copy = np.delete(fuv_copy,ignore)
full_cat_hd = np.append(full_cat_hd, hd_copy)
full_cat_vmag = np.append(full_cat_vmag, vmag_copy)
full_cat_bmag = np.append(full_cat_bmag, bmag_copy)
full_cat_fuv = np.append(full_cat_fuv, fuv_copy)
#fit
x = np.arange(0.5,0.9,0.001)
y = -29.701*x**2 + 48.159*x - 5.831
# Plot FUV-B vs B-V
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.scatter(full_cat_bmag - full_cat_vmag, full_cat_fuv - full_cat_bmag, color = 'black')
plt.plot(x, y, color='red')
plt.xlim(0.4,1.0)
plt.xlabel('B-V', fontsize = 26)
plt.ylabel('FUV-B', fontsize = 26)
plt.annotate('FUV-B vs. B-V', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.show()
fig.savefig("f4a.png", bbox_inches='tight', dpi=300)




#Individual Field Star Fits
#absolute magnitude cuts
i_parallax = i_parallax / 1000. # convert to arcsec
mv = np.array([])
for i in range(len(i_parallax)):
    mv = np.append(mv, i_vmag[i] + 5.0*(1.0 + np.log10(i_parallax[i])))
ignore = np.argwhere(mv <= 4.3) #cuts for Isaacson
i_age = np.delete(i_age, ignore)
i_bmag = np.delete(i_bmag, ignore)
i_vmag = np.delete(i_vmag, ignore)
i_fuv = np.delete(i_fuv, ignore)
i_parallax = np.delete(i_parallax, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
i_hd = np.delete(i_hd, ignore)
ignore = np.argwhere(mv >= 5.3)
i_age = np.delete(i_age, ignore)
i_bmag = np.delete(i_bmag, ignore)
i_vmag = np.delete(i_vmag, ignore)
i_fuv = np.delete(i_fuv, ignore)
i_parallax = np.delete(i_parallax, ignore)
i_hd = np.delete(i_hd, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
b_parallax = b_parallax / 1000. # convert to arcsec
mv = np.array([])
for i in range(len(b_parallax)):
    mv = np.append(mv, b_vmag[i] + 5.0*(1.0 + np.log10(b_parallax[i])))
ignore = np.argwhere(mv <= 4.3) #cuts for Ballering
b_age = np.delete(b_age, ignore)
b_bmag = np.delete(b_bmag, ignore)
b_vmag = np.delete(b_vmag, ignore)
b_fuv = np.delete(b_fuv, ignore)
b_parallax = np.delete(b_parallax, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
b_hd = np.delete(b_hd, ignore)
ignore = np.argwhere(mv >= 5.3)
b_age = np.delete(b_age, ignore)
b_bmag = np.delete(b_bmag, ignore)
b_vmag = np.delete(b_vmag, ignore)
b_fuv = np.delete(b_fuv, ignore)
b_hd = np.delete(b_hd, ignore)
b_parallax = np.delete(b_parallax, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
s_parallax = s_parallax / 1000. # convert to arcsec
mv = np.array([])
for i in range(len(s_parallax)):
    mv = np.append(mv, s_vmag[i] + 5.0*(1.0 + np.log10(s_parallax[i])))
ignore = np.argwhere(mv <= 4.3) #cuts for Sierchio
s_age = np.delete(s_age, ignore)
s_bmag = np.delete(s_bmag, ignore)
s_vmag = np.delete(s_vmag, ignore)
s_fuv = np.delete(s_fuv, ignore)
s_parallax = np.delete(s_parallax, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
s_hd = np.delete(s_hd, ignore)
ignore = np.argwhere(mv >= 5.3)
s_age = np.delete(s_age, ignore)
s_bmag = np.delete(s_bmag, ignore)
s_vmag = np.delete(s_vmag, ignore)
s_fuv = np.delete(s_fuv, ignore)
s_parallax = np.delete(s_parallax, ignore)
s_hd = np.delete(s_hd, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
twin_parallax = twin_parallax / 1000. # convert to arcsec
mv = np.array([])
for i in range(len(twin_parallax)):
    mv = np.append(mv, twin_vmag[i] + 5.0*(1.0 + np.log10(twin_parallax[i])))
ignore = np.argwhere(mv <= 4.3) #cuts for Solar_twin
twin_age = np.delete(twin_age, ignore)
twin_bmag = np.delete(twin_bmag, ignore)
twin_vmag = np.delete(twin_vmag, ignore)
twin_fuv = np.delete(twin_fuv, ignore)
twin_parallax = np.delete(twin_parallax, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
twin_hd = np.delete(twin_hd, ignore)
ignore = np.argwhere(mv >= 5.3)
twin_age = np.delete(twin_age, ignore)
twin_bmag = np.delete(twin_bmag, ignore)
twin_vmag = np.delete(twin_vmag, ignore)
twin_fuv = np.delete(twin_fuv, ignore)
twin_parallax = np.delete(twin_parallax, ignore)
twin_hd = np.delete(twin_hd, ignore)
#rhk = np.delete(rhk, ignore)
mv = np.delete(mv, ignore)
#find B-V
i_bv = np.array([]) #Isaacson
for i in range(len(i_bmag)):
    i_bv = np.append(i_bv, i_bmag[i] - i_vmag[i])
b_bv = np.array([]) #Ballering
for i in range(len(b_bmag)):
    b_bv = np.append(b_bv, b_bmag[i] - b_vmag[i])
s_bv = np.array([]) #Sierchio
for i in range(len(s_bmag)):
    s_bv = np.append(s_bv, s_bmag[i] - s_vmag[i])
twin_bv = np.array([]) #Solar Twin
for i in range(len(twin_bmag)):
    twin_bv = np.append(twin_bv, twin_bmag[i] - twin_vmag[i])
#find FUV-B
i_fuvb = np.array([]) #Isaacson
for i in range(len(i_bmag)):
    i_fuvb = np.append(i_fuvb, i_fuv[i] - i_bmag[i])
b_fuvb = np.array([]) #Ballering
for i in range(len(b_bmag)):
    b_fuvb = np.append(b_fuvb, b_fuv[i] - b_bmag[i])
s_fuvb = np.array([]) #Sierchio
for i in range(len(s_bmag)):
    s_fuvb = np.append(s_fuvb, s_fuv[i] - s_bmag[i])
twin_fuvb = np.array([]) #Solar Twin
for i in range(len(twin_bmag)):
    twin_fuvb = np.append(twin_fuvb, twin_fuv[i] - twin_bmag[i])
#Find Q
i_Q = np.array([]) #Isaacson
for i in range(len(i_bmag)):
    i_Q = np.append(i_Q, i_fuvb[i] + 29.701*(i_bv[i])**2 - 48.159*(i_bv[i]) + 5.831)
b_Q = np.array([]) #Ballering
for i in range(len(b_bmag)):
    b_Q = np.append(b_Q, b_fuvb[i] + 29.701*(b_bv[i])**2 - 48.159*(b_bv[i]) + 5.831)
s_Q = np.array([]) #Sierchio
for i in range(len(s_bmag)):
    s_Q = np.append(s_Q, s_fuvb[i] + 29.701*(s_bv[i])**2 - 48.159*(s_bv[i]) + 5.831)
twin_Q = np.array([]) #Solar Twin
for i in range(len(twin_bmag)):
    twin_Q = np.append(twin_Q, twin_fuvb[i] + 29.701*(twin_bv[i])**2 - 48.159*(twin_bv[i]) + 5.831)

#print 'number of field stars before outliers (there are duplicates)', len(i_Q) + len(b_Q) + len(s_Q) + len(twin_Q)



#0.553-0.603 B-V
print "0.553-0.603"
rang = np.argwhere((i_bv >= 0.553) & (i_bv <= 0.603))
total = np.array([]) #total is the number of stars within 0.553-0.753 range including outliers
total = np.append(total, len(rang))
i_xdata = np.array([])
i_ydata = np.array([])
i_hd_used = np.array([])
i_bv_used = np.array([])
i_fuvb_used = np.array([])
i_vmag_used = np.array([])
for i in range(len(rang)):
	i_xdata = np.append(i_xdata, i_Q[rang[i]])
	i_ydata = np.append(i_ydata, i_age[rang[i]])
	i_hd_used = np.append(i_hd_used, i_hd[rang[i]])
	i_bv_used = np.append(i_bv_used, i_bv[rang[i]])
	i_fuvb_used = np.append(i_fuvb_used, i_fuvb[rang[i]])
	i_vmag_used = np.append(i_vmag_used, i_vmag[rang[i]])
i_p = np.polyfit(i_xdata,np.log(i_ydata),1)
dist = np.abs((i_p[0] * i_xdata) - np.log(i_ydata) + i_p[1]) / np.abs(i_p[0])
i_outlier = np.argwhere(dist > 0.62)
i_fit_x = np.delete(i_xdata, i_outlier)
i_fit_y = np.delete(i_ydata, i_outlier)
i_fit_hd = np.delete(i_hd_used, i_outlier)
i_fit_bv = np.delete(i_bv_used, i_outlier)
i_fit_fuvb = np.delete(i_fuvb_used, i_outlier)
i_fit_vmag = np.delete(i_vmag_used, i_outlier)
print "test!!!!!"
for i in range(len(i_hd_used)):
	print i_hd_used[i], i_vmag_used[i]
print "Number in Isaacson sample", len(i_fit_x)
avg_bv = np.array([])
for i in range(len(i_fit_x)):
    avg_bv = np.append(avg_bv, i_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(i_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "i" + "\t" +"&"  + i_fit_hd[i] + "\t" +"&"  + str(round(i_fit_vmag[i],2)) + "\t" +"&" + str(round(i_fit_bv[i],2)) + "\t" +"&"  + str(round(i_fit_fuvb[i],2)) + "\t" +"&" + str(round(i_fit_x[i],2)) + "\t" +"&" + str(i_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.array([])
all_field_y = np.array([])
all_field_bv = np.array([])
all_field_hd = np.array([])
all_field_x = np.append(all_field_x, i_fit_x)
all_field_y = np.append(all_field_y, i_fit_y)
all_field_bv = np.append(all_field_bv, i_fit_bv)
all_field_hd = np.append(all_field_hd, i_fit_hd)
i_p = np.polyfit(i_fit_x,np.log(i_fit_y),1)
print "Isaacson", "a=", np.exp(i_p[1]), "b=", i_p[0]
slope, intercept, r_value, p_value, std_err = linregress(i_fit_x, i_fit_y)
print "Isaacson", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(i_xdata, i_ydata)
print "Isaacson", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(i_fit_x, i_fit_y)
print "Isaacson", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(i_xdata, i_ydata)
print "Isaacson", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((b_bv >= 0.553) & (b_bv <= 0.603))
total = np.append(total, len(rang))
b_xdata = np.array([])
b_ydata = np.array([])
b_hd_used = np.array([])
b_bv_used = np.array([])
b_fuvb_used = np.array([])
b_vmag_used = np.array([])
for i in range(len(rang)):
	b_xdata = np.append(b_xdata, b_Q[rang[i]])
	b_ydata = np.append(b_ydata, b_age[rang[i]])
	b_hd_used = np.append(b_hd_used, b_hd[rang[i]])
	b_bv_used = np.append(b_bv_used, b_bv[rang[i]])
	b_fuvb_used = np.append(b_fuvb_used, b_fuvb[rang[i]])
	b_vmag_used = np.append(b_vmag_used, b_vmag[rang[i]])
all_field_x = np.append(all_field_x, b_xdata)
all_field_y = np.append(all_field_y, b_ydata)
all_field_bv = np.append(all_field_bv, b_bv_used)
all_field_hd = np.append(all_field_hd, b_hd_used)
print "Number in Ballering sample", len(b_xdata)
avg_bv = np.array([])
for i in range(len(b_hd_used)):
    avg_bv = np.append(avg_bv, b_bv_used[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(b_hd_used)):
	filename = open('field_stars_fit.txt', "a")
	text = "b" + "\t" +"&"  + b_hd_used[i] + "\t" +"&"  + str(round(b_vmag_used[i],2)) + "\t" +"&" + str(round(b_bv_used[i],2)) + "\t" +"&"  + str(round(b_fuvb_used[i],2)) + "\t" +"&" + str(round(b_xdata[i],2)) + "\t" +"&" + str(b_ydata[i]) +  "\\"+ "\\" +  "\n"
	filename.write(text)
	filename.close()
b_p = np.polyfit(b_xdata,np.log(b_ydata),1)
print "Ballering", "a=", np.exp(b_p[1]), "b=", b_p[0]
slope, intercept, r_value, p_value, std_err = linregress(b_xdata, b_ydata)
print "Ballering", "r^2", r_value**2
rho, p = spearmanr(b_xdata, b_ydata)
print "Ballering", 'Spearman coefficient: ', rho
num = len(b_xdata)
stderr = 1.0 / math.sqrt(num - 3)
delta = 1.96 * stderr
lower = math.tanh(math.atanh(rho) - delta)
upper = math.tanh(math.atanh(rho) + delta)
print "Ballering", '95% confidence interval', lower, upper
rang = np.argwhere((s_bv >= 0.553) & (s_bv <= 0.603))
total = np.append(total, len(rang))
s_xdata = np.array([])
s_ydata = np.array([])
s_hd_used = np.array([])
s_bv_used = np.array([])
s_fuvb_used = np.array([])
s_vmag_used = np.array([])
for i in range(len(rang)):
	s_xdata = np.append(s_xdata, s_Q[rang[i]])
	s_ydata = np.append(s_ydata, s_age[rang[i]])
	s_hd_used = np.append(s_hd_used, s_hd[rang[i]])
	s_bv_used = np.append(s_bv_used, s_bv[rang[i]])
	s_fuvb_used = np.append(s_fuvb_used, s_fuvb[rang[i]])
	s_vmag_used = np.append(s_vmag_used, s_vmag[rang[i]])
s_p = np.polyfit(s_xdata,np.log(s_ydata),1)
dist = np.abs((s_p[0] * s_xdata) - np.log(s_ydata) + s_p[1]) / np.abs(s_p[0])
s_outlier = np.argwhere(dist > 0.4)
s_fit_x = np.delete(s_xdata, s_outlier)
s_fit_y = np.delete(s_ydata, s_outlier)
s_fit_hd = np.delete(s_hd_used, s_outlier)
s_fit_bv = np.delete(s_bv_used, s_outlier)
s_fit_fuvb = np.delete(s_fuvb_used, s_outlier)
s_fit_vmag = np.delete(s_vmag_used, s_outlier)
print "Number in Sierchio sample", len(s_fit_x)
avg_bv = np.array([])
for i in range(len(s_fit_x)):
    avg_bv = np.append(avg_bv, s_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(s_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "s" + "\t" +"&"  + s_fit_hd[i] + "\t" +"&"  + str(round(s_fit_vmag[i],2)) + "\t" +"&" + str(round(s_fit_bv[i],2)) + "\t" +"&"  + str(round(s_fit_fuvb[i],2)) + "\t" +"&" + str(round(s_fit_x[i],2)) + "\t" +"&" + str(s_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
s_p = np.polyfit(s_fit_x,np.log(s_fit_y),1)
all_field_x = np.append(all_field_x, s_fit_x)
all_field_y = np.append(all_field_y, s_fit_y)
all_field_bv = np.append(all_field_bv, s_fit_bv)
all_field_hd = np.append(all_field_hd, s_fit_hd)
print "Sierchio", "a=", np.exp(s_p[1]), "b=", s_p[0]
slope, intercept, r_value, p_value, std_err = linregress(s_fit_x, s_fit_y)
print "Sierchio", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(s_xdata, s_ydata)
print "Sierchio", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(s_fit_x, s_fit_y)
print "Sierchio", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(s_xdata, s_ydata)
print "Sierchio", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((twin_bv >= 0.553) & (twin_bv <= 0.603))
print "Number in Twin sample", len(rang)


#0.603-0.653
print "0.603-0.653"
rang = np.argwhere((i_bv >= 0.603) & (i_bv <= 0.653))
total = np.append(total, len(rang))
i_xdata = np.array([])
i_ydata = np.array([])
i_hd_used = np.array([])
i_bv_used = np.array([])
i_fuvb_used = np.array([])
i_vmag_used = np.array([])
for i in range(len(rang)):
	i_xdata = np.append(i_xdata, i_Q[rang[i]])
	i_ydata = np.append(i_ydata, i_age[rang[i]])
	i_hd_used = np.append(i_hd_used, i_hd[rang[i]])
	i_bv_used = np.append(i_bv_used, i_bv[rang[i]])
	i_fuvb_used = np.append(i_fuvb_used, i_fuvb[rang[i]])
	i_vmag_used = np.append(i_vmag_used, i_vmag[rang[i]])
i_p = np.polyfit(i_xdata,np.log(i_ydata),1)
dist = np.abs((i_p[0] * i_xdata) - np.log(i_ydata) + i_p[1]) / np.abs(i_p[0])
i_outlier = np.argwhere(dist > 0.515)
i_fit_x = np.delete(i_xdata, i_outlier)
i_fit_y = np.delete(i_ydata, i_outlier)
i_fit_hd = np.delete(i_hd_used, i_outlier)
i_fit_bv = np.delete(i_bv_used, i_outlier)
i_fit_fuvb = np.delete(i_fuvb_used, i_outlier)
i_fit_vmag = np.delete(i_vmag_used, i_outlier)
all_field_x = np.append(all_field_x, i_fit_x)
all_field_y = np.append(all_field_y, i_fit_y)
all_field_bv = np.append(all_field_bv, i_fit_bv)
all_field_hd = np.append(all_field_hd, i_fit_hd)
print "Number in Isaacson sample", len(i_fit_x)
avg_bv = np.array([])
for i in range(len(i_fit_x)):
    avg_bv = np.append(avg_bv, i_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(i_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "i" + "\t" +"&"  + i_fit_hd[i] + "\t" +"&"  + str(round(i_fit_vmag[i],2)) + "\t" +"&" + str(round(i_fit_bv[i],2)) + "\t" +"&"  + str(round(i_fit_fuvb[i],2)) + "\t" +"&" + str(round(i_fit_x[i],2)) + "\t" +"&" + str(i_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
i_p = np.polyfit(i_fit_x,np.log(i_fit_y),1)
print "Isaacson", "a=", np.exp(i_p[1]), "b=", i_p[0]
slope, intercept, r_value, p_value, std_err = linregress(i_fit_x, i_fit_y)
print "Isaacson", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(i_xdata, i_ydata)
print "Isaacson", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(i_fit_x, i_fit_y)
print "Isaacson", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(i_xdata, i_ydata)
print "Isaacson", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((b_bv >= 0.603) & (b_bv <= 0.653))
total = np.append(total, len(rang))
b_xdata = np.array([])
b_ydata = np.array([])
b_hd_used = np.array([])
b_bv_used = np.array([])
b_fuvb_used = np.array([])
b_vmag_used = np.array([])
for i in range(len(rang)):
	b_xdata = np.append(b_xdata, b_Q[rang[i]])
	b_ydata = np.append(b_ydata, b_age[rang[i]])
	b_hd_used = np.append(b_hd_used, b_hd[rang[i]])
	b_bv_used = np.append(b_bv_used, b_bv[rang[i]])
	b_fuvb_used = np.append(b_fuvb_used, b_fuvb[rang[i]])
	b_vmag_used = np.append(b_vmag_used, b_vmag[rang[i]])
b_p = np.polyfit(b_xdata,np.log(b_ydata),1)
dist = np.abs((b_p[0] * b_xdata) - np.log(b_ydata) + b_p[1]) / np.abs(b_p[0])
b_outlier = np.argwhere(dist > 0.45)
b_fit_x = np.delete(b_xdata, b_outlier)
b_fit_y = np.delete(b_ydata, b_outlier)
b_fit_hd = np.delete(b_hd_used, b_outlier)
b_fit_bv = np.delete(b_bv_used, b_outlier)
b_fit_fuvb = np.delete(b_fuvb_used, b_outlier)
b_fit_vmag = np.delete(b_vmag_used, b_outlier)
all_field_x = np.append(all_field_x, b_fit_x)
all_field_y = np.append(all_field_y, b_fit_y)
all_field_bv = np.append(all_field_bv, b_fit_bv)
all_field_hd = np.append(all_field_hd, b_fit_hd)
print "Number in Ballering sample", len(b_fit_x)
avg_bv = np.array([])
for i in range(len(b_fit_x)):
    avg_bv = np.append(avg_bv, b_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(b_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "b" + "\t" +"&"  + b_fit_hd[i] + "\t" +"&"  + str(round(b_fit_vmag[i],2)) + "\t" +"&" + str(round(b_fit_bv[i],2)) + "\t" +"&"  + str(round(b_fit_fuvb[i],2)) + "\t" +"&" + str(round(b_fit_x[i],2)) + "\t" +"&" + str(b_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
b_p = np.polyfit(b_fit_x,np.log(b_fit_y),1)
print "Ballering", "a=", np.exp(b_p[1]), "b=", b_p[0]
slope, intercept, r_value, p_value, std_err = linregress(b_fit_x, b_fit_y)
print "Ballering", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(b_xdata, b_ydata)
print "Ballering", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(b_fit_x, b_fit_y)
print "Ballering", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(b_xdata, b_ydata)
print "Ballering", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((s_bv >= 0.603) & (s_bv <= 0.653))
total = np.append(total, len(rang))
s_xdata = np.array([])
s_ydata = np.array([])
s_hd_used = np.array([])
s_bv_used = np.array([])
s_fuvb_used = np.array([])
s_vmag_used = np.array([])
for i in range(len(rang)):
	s_xdata = np.append(s_xdata, s_Q[rang[i]])
	s_ydata = np.append(s_ydata, s_age[rang[i]])
	s_hd_used = np.append(s_hd_used, s_hd[rang[i]])
	s_bv_used = np.append(s_bv_used, s_bv[rang[i]])
	s_fuvb_used = np.append(s_fuvb_used, s_fuvb[rang[i]])
	s_vmag_used = np.append(s_vmag_used, s_vmag[rang[i]])
s_p = np.polyfit(s_xdata,np.log(s_ydata),1)
dist = np.abs((s_p[0] * s_xdata) - np.log(s_ydata) + s_p[1]) / np.abs(s_p[0])
s_outlier = np.argwhere(dist > 0.34)
s_fit_x = np.delete(s_xdata, s_outlier)
s_fit_y = np.delete(s_ydata, s_outlier)
s_fit_hd = np.delete(s_hd_used, s_outlier)
s_fit_bv = np.delete(s_bv_used, s_outlier)
s_fit_fuvb = np.delete(s_fuvb_used, s_outlier)
s_fit_vmag = np.delete(s_vmag_used, s_outlier)
print "Number in Sierchio sample", len(s_fit_x)
avg_bv = np.array([])
for i in range(len(s_fit_x)):
    avg_bv = np.append(avg_bv, s_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(s_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "s" + "\t" +"&"  + s_fit_hd[i] + "\t" +"&"  + str(round(s_fit_vmag[i],2)) + "\t" +"&" + str(round(s_fit_bv[i],2)) + "\t" +"&"  + str(round(s_fit_fuvb[i],2)) + "\t" +"&" + str(round(s_fit_x[i],2)) + "\t" +"&" + str(s_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, s_fit_x)
all_field_y = np.append(all_field_y, s_fit_y)
all_field_bv = np.append(all_field_bv, s_fit_bv)
all_field_hd = np.append(all_field_hd, s_fit_hd)
s_p = np.polyfit(s_fit_x,np.log(s_fit_y),1)
print "Sierchio", "a=", np.exp(s_p[1]), "b=", s_p[0]
slope, intercept, r_value, p_value, std_err = linregress(s_fit_x, s_fit_y)
print "Sierchio", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(s_xdata, s_ydata)
print "Sierchio", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(s_fit_x, s_fit_y)
print "Sierchio", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(s_xdata, s_ydata)
print "Sierchio", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((twin_bv >= 0.603) & (twin_bv <= 0.653))
total = np.append(total, len(rang))
twin_xdata = np.array([])
twin_ydata = np.array([])
twin_hd_used = np.array([])
twin_bv_used = np.array([])
twin_fuvb_used = np.array([])
twin_vmag_used = np.array([])
for i in range(len(rang)):
	twin_xdata = np.append(twin_xdata, twin_Q[rang[i]])
	twin_ydata = np.append(twin_ydata, twin_age[rang[i]])
	twin_hd_used = np.append(twin_hd_used, twin_hd[rang[i]])
	twin_bv_used = np.append(twin_bv_used, twin_bv[rang[i]])
	twin_fuvb_used = np.append(twin_fuvb_used, twin_fuvb[rang[i]])
	twin_vmag_used = np.append(twin_vmag_used, twin_vmag[rang[i]])
twin_p = np.polyfit(twin_xdata,np.log(twin_ydata),1)
dist = np.abs((twin_p[0] * twin_xdata) - np.log(twin_ydata) + twin_p[1]) / np.abs(twin_p[0])
twin_outlier = np.argwhere(dist > 0.54)
twin_fit_x = np.delete(twin_xdata, twin_outlier)
twin_fit_y = np.delete(twin_ydata, twin_outlier)
twin_fit_hd = np.delete(twin_hd_used, twin_outlier)
twin_fit_bv = np.delete(twin_bv_used, twin_outlier)
twin_fit_fuvb = np.delete(twin_fuvb_used, twin_outlier)
twin_fit_vmag = np.delete(twin_vmag_used, twin_outlier)
twin_fit_bv = np.delete(twin_bv_used, twin_outlier)
print "Number in Twin sample", len(twin_fit_x)
avg_bv = np.array([])
for i in range(len(twin_fit_x)):
    avg_bv = np.append(avg_bv, twin_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(twin_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "l" + "\t" +"&"  + twin_fit_hd[i] + "\t" +"&"  + str(round(twin_fit_vmag[i],2)) + "\t" +"&" + str(round(twin_fit_bv[i],2)) + "\t" +"&"  + str(round(twin_fit_fuvb[i],2)) + "\t" +"&" + str(round(twin_fit_x[i],2)) + "\t" +"&" + str(twin_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, twin_fit_x)
all_field_y = np.append(all_field_y, twin_fit_y)
all_field_bv = np.append(all_field_bv, twin_fit_bv)
all_field_hd = np.append(all_field_hd, twin_fit_hd)
twin_p = np.polyfit(twin_fit_x,np.log(twin_fit_y),1)
print "Twin", "a=", np.exp(twin_p[1]), "b=", twin_p[0]
slope, intercept, r_value, p_value, std_err = linregress(twin_fit_x, twin_fit_y)
print "Twin", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(twin_xdata, twin_ydata)
print "Twin", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(twin_fit_x, twin_fit_y)
print "Twin", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(twin_xdata, twin_ydata)
print "Twin", 'Spearman coefficient (w/ outliers): ', rho_2
for i in range(len(twin_fit_x)):
    #filename = open('avg_exclude.txt', "a")
    filename = open('test_twin.txt', "a")
    text = twin_hd_used[i] + "\t" + str(twin_bv_used[i]) + "\t" + str(twin_fit_x[i]) + "\t" + str(twin_fit_y[i]) + "\n"
    filename.write(text)
    filename.close()

#0.653-0.703
print "0.653-0.703"
rang = np.argwhere((i_bv >= 0.653) & (i_bv <= 0.703))
total = np.append(total, len(rang))
i_xdata = np.array([])
i_ydata = np.array([])
i_hd_used = np.array([])
i_bv_used = np.array([])
i_fuvb_used = np.array([])
i_vmag_used = np.array([])
for i in range(len(rang)):
	i_xdata = np.append(i_xdata, i_Q[rang[i]])
	i_ydata = np.append(i_ydata, i_age[rang[i]])
	i_hd_used = np.append(i_hd_used, i_hd[rang[i]])
	i_bv_used = np.append(i_bv_used, i_bv[rang[i]])
	i_fuvb_used = np.append(i_fuvb_used, i_fuvb[rang[i]])
	i_vmag_used = np.append(i_vmag_used, i_vmag[rang[i]])
i_p = np.polyfit(i_xdata,np.log(i_ydata),1)
dist = np.abs((i_p[0] * i_xdata) - np.log(i_ydata) + i_p[1]) / np.abs(i_p[0])
i_outlier = np.argwhere(dist > 0.493)
i_fit_x = np.delete(i_xdata, i_outlier)
i_fit_y = np.delete(i_ydata, i_outlier)
i_fit_hd = np.delete(i_hd_used, i_outlier)
i_fit_bv = np.delete(i_bv_used, i_outlier)
i_fit_fuvb = np.delete(i_fuvb_used, i_outlier)
i_fit_vmag = np.delete(i_vmag_used, i_outlier)
all_field_x = np.append(all_field_x, i_fit_x)
all_field_y = np.append(all_field_y, i_fit_y)
all_field_bv = np.append(all_field_bv, i_fit_bv)
all_field_hd = np.append(all_field_hd, i_fit_hd)
print "Number in Isaacson sample", len(i_fit_x)
avg_bv = np.array([])
for i in range(len(i_fit_x)):
    avg_bv = np.append(avg_bv, i_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(i_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "i" + "\t" +"&"  + i_fit_hd[i] + "\t" +"&"  + str(round(i_fit_vmag[i],2)) + "\t" +"&" + str(round(i_fit_bv[i],2)) + "\t" +"&"  + str(round(i_fit_fuvb[i],2)) + "\t" +"&" + str(round(i_fit_x[i],2)) + "\t" +"&" + str(i_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
i_p = np.polyfit(i_fit_x,np.log(i_fit_y),1)
print "Isaacson", "a=", np.exp(i_p[1]), "b=", i_p[0]
slope, intercept, r_value, p_value, std_err = linregress(i_fit_x, i_fit_y)
print "Isaacson", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(i_xdata, i_ydata)
print "Isaacson", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(i_fit_x, i_fit_y)
print "Isaacson", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(i_xdata, i_ydata)
print "Isaacson", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((b_bv >= 0.653) & (b_bv <= 0.703))
total = np.append(total, len(rang))
b_xdata = np.array([])
b_ydata = np.array([])
b_hd_used = np.array([])
b_bv_used = np.array([])
b_fuvb_used = np.array([])
b_vmag_used = np.array([])
for i in range(len(rang)):
	b_xdata = np.append(b_xdata, b_Q[rang[i]])
	b_ydata = np.append(b_ydata, b_age[rang[i]])
	b_hd_used = np.append(b_hd_used, b_hd[rang[i]])
	b_bv_used = np.append(b_bv_used, b_bv[rang[i]])
	b_fuvb_used = np.append(b_fuvb_used, b_fuvb[rang[i]])
	b_vmag_used = np.append(b_vmag_used, b_vmag[rang[i]])
b_p = np.polyfit(b_xdata,np.log(b_ydata),1)
dist = np.abs((b_p[0] * b_xdata) - np.log(b_ydata) + b_p[1]) / np.abs(b_p[0])
b_outlier = np.argwhere(dist > 0.29)
b_fit_x = np.delete(b_xdata, b_outlier)
b_fit_y = np.delete(b_ydata, b_outlier)
b_fit_hd = np.delete(b_hd_used, b_outlier)
b_fit_bv = np.delete(b_bv_used, b_outlier)
b_fit_fuvb = np.delete(b_fuvb_used, b_outlier)
b_fit_vmag = np.delete(b_vmag_used, b_outlier)
print "Number in Ballering sample", len(b_fit_x)
avg_bv = np.array([])
for i in range(len(b_fit_x)):
    avg_bv = np.append(avg_bv, b_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(b_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "b" + "\t" +"&"  + b_fit_hd[i] + "\t" +"&"  + str(round(b_fit_vmag[i],2)) + "\t" +"&" + str(round(b_fit_bv[i],2)) + "\t" +"&"  + str(round(b_fit_fuvb[i],2)) + "\t" +"&" + str(round(b_fit_x[i],2)) + "\t" +"&" + str(b_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, b_fit_x)
all_field_y = np.append(all_field_y, b_fit_y)
all_field_bv = np.append(all_field_bv, b_fit_bv)
all_field_hd = np.append(all_field_hd, b_fit_hd)
b_p = np.polyfit(b_fit_x,np.log(b_fit_y),1)
print "Ballering", "a=", np.exp(b_p[1]), "b=", b_p[0]
slope, intercept, r_value, p_value, std_err = linregress(b_fit_x, b_fit_y)
print "Ballering", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(b_xdata, b_ydata)
print "Ballering", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(b_fit_x, b_fit_y)
print "Ballering", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(b_xdata, b_ydata)
print "Ballering", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((s_bv >= 0.653) & (s_bv <= 0.703))
total = np.append(total, len(rang))
s_xdata = np.array([])
s_ydata = np.array([])
s_hd_used = np.array([])
s_bv_used = np.array([])
s_fuvb_used = np.array([])
s_vmag_used = np.array([])
for i in range(len(rang)):
	s_xdata = np.append(s_xdata, s_Q[rang[i]])
	s_ydata = np.append(s_ydata, s_age[rang[i]])
	s_hd_used = np.append(s_hd_used, s_hd[rang[i]])
	s_bv_used = np.append(s_bv_used, s_bv[rang[i]])
	s_fuvb_used = np.append(s_fuvb_used, s_fuvb[rang[i]])
	s_vmag_used = np.append(s_vmag_used, s_vmag[rang[i]])
s_p = np.polyfit(s_xdata,np.log(s_ydata),1)
dist = np.abs((s_p[0] * s_xdata) - np.log(s_ydata) + s_p[1]) / np.abs(s_p[0])
s_outlier = np.argwhere(dist > 0.4)
s_fit_x = np.delete(s_xdata, s_outlier)
s_fit_y = np.delete(s_ydata, s_outlier)
s_fit_hd = np.delete(s_hd_used, s_outlier)
s_fit_bv = np.delete(s_bv_used, s_outlier)
s_fit_fuvb = np.delete(s_fuvb_used, s_outlier)
s_fit_vmag = np.delete(s_vmag_used, s_outlier)
print "Number in Sierchio sample", len(s_fit_x)
avg_bv = np.array([])
for i in range(len(s_fit_x)):
    avg_bv = np.append(avg_bv, s_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(s_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "s" + "\t" +"&"  + s_fit_hd[i] + "\t" +"&"  + str(round(s_fit_vmag[i],2)) + "\t" +"&" + str(round(s_fit_bv[i],2)) + "\t" +"&"  + str(round(s_fit_fuvb[i],2)) + "\t" +"&" + str(round(s_fit_x[i],2)) + "\t" +"&" + str(s_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, s_fit_x)
all_field_y = np.append(all_field_y, s_fit_y)
all_field_bv = np.append(all_field_bv, s_fit_bv)
all_field_hd = np.append(all_field_hd, s_fit_hd)
s_p = np.polyfit(s_fit_x,np.log(s_fit_y),1)
print "Sierchio", "a=", np.exp(s_p[1]), "b=", s_p[0]
slope, intercept, r_value, p_value, std_err = linregress(s_fit_x, s_fit_y)
print "Sierchio", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(s_xdata, s_ydata)
print "Sierchio", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(s_fit_x, s_fit_y)
print "Sierchio", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(s_xdata, s_ydata)
print "Sierchio", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((twin_bv >= 0.653) & (twin_bv <= 0.703))
total = np.append(total, len(rang))
twin_xdata = np.array([])
twin_ydata = np.array([])
twin_hd_used = np.array([])
twin_bv_used = np.array([])
twin_fuvb_used = np.array([])
twin_vmag_used = np.array([])
for i in range(len(rang)):
	twin_xdata = np.append(twin_xdata, twin_Q[rang[i]])
	twin_ydata = np.append(twin_ydata, twin_age[rang[i]])
	twin_hd_used = np.append(twin_hd_used, twin_hd[rang[i]])
	twin_bv_used = np.append(twin_bv_used, twin_bv[rang[i]])
	twin_fuvb_used = np.append(twin_fuvb_used, twin_fuvb[rang[i]])
	twin_vmag_used = np.append(twin_vmag_used, twin_vmag[rang[i]])
twin_p = np.polyfit(twin_xdata,np.log(twin_ydata),1)
dist = np.abs((twin_p[0] * twin_xdata) - np.log(twin_ydata) + twin_p[1]) / np.abs(twin_p[0])
twin_outlier = np.argwhere(dist > 0.42)
twin_fit_x = np.delete(twin_xdata, twin_outlier)
twin_fit_y = np.delete(twin_ydata, twin_outlier)
twin_fit_hd = np.delete(twin_hd_used, twin_outlier)
twin_fit_bv = np.delete(twin_bv_used, twin_outlier)
twin_fit_fuvb = np.delete(twin_fuvb_used, twin_outlier)
twin_fit_vmag = np.delete(twin_vmag_used, twin_outlier)
print "Number in Twin sample", len(twin_fit_x)
avg_bv = np.array([])
for i in range(len(twin_fit_x)):
    avg_bv = np.append(avg_bv, twin_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(twin_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "l" + "\t" +"&"  + twin_fit_hd[i] + "\t" +"&"  + str(round(twin_fit_vmag[i],2)) + "\t" +"&" + str(round(twin_fit_bv[i],2)) + "\t" +"&"  + str(round(twin_fit_fuvb[i],2)) + "\t" +"&" + str(round(twin_fit_x[i],2)) + "\t" +"&" + str(twin_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, twin_fit_x)
all_field_y = np.append(all_field_y, twin_fit_y)
all_field_bv = np.append(all_field_bv, twin_fit_bv)
all_field_hd = np.append(all_field_hd, twin_fit_hd)
twin_p = np.polyfit(twin_fit_x,np.log(twin_fit_y),1)
print "Twin", "a=", np.exp(twin_p[1]), "b=", twin_p[0]
slope, intercept, r_value, p_value, std_err = linregress(twin_fit_x, twin_fit_y)
print "Twin", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(twin_xdata, twin_ydata)
print "Twin", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(twin_fit_x, twin_fit_y)
print "Twin", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(twin_xdata, twin_ydata)
print "Twin", 'Spearman coefficient (w/ outliers): ', rho_2
for i in range(len(twin_fit_x)):
    filename = open('test_twin.txt', "a")
    text = twin_hd_used[i] + "\t" + str(twin_bv_used[i]) + "\t" + str(twin_fit_x[i]) + "\t" + str(twin_fit_y[i]) + "\n"
    filename.write(text)
    filename.close()
#0.703-0.753
print "0.703-0.753"
rang = np.argwhere((i_bv >= 0.703) & (i_bv <= 0.753))
total = np.append(total, len(rang))
i_xdata = np.array([])
i_ydata = np.array([])
i_hd_used = np.array([])
i_bv_used = np.array([])
i_fuvb_used = np.array([])
i_vmag_used = np.array([])
for i in range(len(rang)):
	i_xdata = np.append(i_xdata, i_Q[rang[i]])
	i_ydata = np.append(i_ydata, i_age[rang[i]])
	i_hd_used = np.append(i_hd_used, i_hd[rang[i]])
	i_bv_used = np.append(i_bv_used, i_bv[rang[i]])
	i_fuvb_used = np.append(i_fuvb_used, i_fuvb[rang[i]])
	i_vmag_used = np.append(i_vmag_used, i_vmag[rang[i]])
i_p = np.polyfit(i_xdata,np.log(i_ydata),1)
dist = np.abs((i_p[0] * i_xdata) - np.log(i_ydata) + i_p[1]) / np.abs(i_p[0])
i_outlier = np.argwhere(dist > 0.64)
i_fit_x = np.delete(i_xdata, i_outlier)
i_fit_y = np.delete(i_ydata, i_outlier)
i_fit_hd = np.delete(i_hd_used, i_outlier)
i_fit_bv = np.delete(i_bv_used, i_outlier)
i_fit_fuvb = np.delete(i_fuvb_used, i_outlier)
i_fit_vmag = np.delete(i_vmag_used, i_outlier)
print "Number in Isaacson sample", len(i_fit_x)
avg_bv = np.array([])
for i in range(len(i_fit_x)):
    avg_bv = np.append(avg_bv, i_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(i_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "i" + "\t" +"&"  + i_fit_hd[i] + "\t" +"&"  + str(round(i_fit_vmag[i],2)) + "\t" +"&" + str(round(i_fit_bv[i],2)) + "\t" +"&"  + str(round(i_fit_fuvb[i],2)) + "\t" +"&" + str(round(i_fit_x[i],2)) + "\t" +"&" + str(i_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, i_fit_x)
all_field_y = np.append(all_field_y, i_fit_y)
all_field_bv = np.append(all_field_bv, i_fit_bv)
all_field_hd = np.append(all_field_hd, i_fit_hd)
i_p = np.polyfit(i_fit_x,np.log(i_fit_y),1)
print "Isaacson", "a=", np.exp(i_p[1]), "b=", i_p[0]
slope, intercept, r_value, p_value, std_err = linregress(i_fit_x, i_fit_y)
print "Isaacson", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(i_xdata, i_ydata)
print "Isaacson", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(i_fit_x, i_fit_y)
print "Isaacson", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(i_xdata, i_ydata)
print "Isaacson", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((b_bv >= 0.703) & (b_bv <= 0.753))
total = np.append(total, len(rang))
b_xdata = np.array([])
b_ydata = np.array([])
b_hd_used = np.array([])
b_bv_used = np.array([])
b_fuvb_used = np.array([])
b_vmag_used = np.array([])
for i in range(len(rang)):
	b_xdata = np.append(b_xdata, b_Q[rang[i]])
	b_ydata = np.append(b_ydata, b_age[rang[i]])
	b_hd_used = np.append(b_hd_used, b_hd[rang[i]])
	b_bv_used = np.append(b_bv_used, b_bv[rang[i]])
	b_fuvb_used = np.append(b_fuvb_used, b_fuvb[rang[i]])
	b_vmag_used = np.append(b_vmag_used, b_vmag[rang[i]])
b_p = np.polyfit(b_xdata,np.log(b_ydata),1)
dist = np.abs((b_p[0] * b_xdata) - np.log(b_ydata) + b_p[1]) / np.abs(b_p[0])
b_outlier = np.argwhere(dist > 0.3)
b_fit_x = np.delete(b_xdata, b_outlier)
b_fit_y = np.delete(b_ydata, b_outlier)
b_fit_hd = np.delete(b_hd_used, b_outlier)
b_fit_bv = np.delete(b_bv_used, b_outlier)
b_fit_fuvb = np.delete(b_fuvb_used, b_outlier)
b_fit_vmag = np.delete(b_vmag_used, b_outlier)
print "Number in Ballering sample", len(b_fit_x)
avg_bv = np.array([])
for i in range(len(b_fit_x)):
    avg_bv = np.append(avg_bv, b_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(b_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "b" + "\t" +"&"  + b_fit_hd[i] + "\t" +"&"  + str(round(b_fit_vmag[i],2)) + "\t" +"&" + str(round(b_fit_bv[i],2)) + "\t" +"&"  + str(round(b_fit_fuvb[i],2)) + "\t" +"&" + str(round(b_fit_x[i],2)) + "\t" +"&" + str(b_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, b_fit_x)
all_field_y = np.append(all_field_y, b_fit_y)
all_field_bv = np.append(all_field_bv, b_fit_bv)
all_field_hd = np.append(all_field_hd, b_fit_hd)
b_p = np.polyfit(b_fit_x,np.log(b_fit_y),1)
print "Ballering", "a=", np.exp(b_p[1]), "b=", b_p[0]
slope, intercept, r_value, p_value, std_err = linregress(b_fit_x, b_fit_y)
print "Ballering", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(b_xdata, b_ydata)
print "Ballering", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(b_fit_x, b_fit_y)
print "Ballering", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(b_xdata, b_ydata)
print "Ballering", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((s_bv >= 0.703) & (s_bv <= 0.753))
total = np.append(total, len(rang))
s_xdata = np.array([])
s_ydata = np.array([])
s_hd_used = np.array([])
s_bv_used = np.array([])
s_fuvb_used = np.array([])
s_vmag_used = np.array([])
for i in range(len(rang)):
	s_xdata = np.append(s_xdata, s_Q[rang[i]])
	s_ydata = np.append(s_ydata, s_age[rang[i]])
	s_hd_used = np.append(s_hd_used, s_hd[rang[i]])
	s_bv_used = np.append(s_bv_used, s_bv[rang[i]])
	s_fuvb_used = np.append(s_fuvb_used, s_fuvb[rang[i]])
	s_vmag_used = np.append(s_vmag_used, s_vmag[rang[i]])
s_p = np.polyfit(s_xdata,np.log(s_ydata),1)
dist = np.abs((s_p[0] * s_xdata) - np.log(s_ydata) + s_p[1]) / np.abs(s_p[0])
s_outlier = np.argwhere(dist > 0.3)
s_fit_x = np.delete(s_xdata, s_outlier)
s_fit_y = np.delete(s_ydata, s_outlier)
s_fit_hd = np.delete(s_hd_used, s_outlier)
s_fit_bv = np.delete(s_bv_used, s_outlier)
s_fit_fuvb = np.delete(s_fuvb_used, s_outlier)
s_fit_vmag = np.delete(s_vmag_used, s_outlier)
print "Number in Sierchio sample", len(s_fit_x)
avg_bv = np.array([])
for i in range(len(s_fit_x)):
    avg_bv = np.append(avg_bv, s_fit_bv[i])
print "avg. B-V", np.mean(avg_bv)
for i in range(len(s_fit_hd)):
	filename = open('field_stars_fit.txt', "a")
	text = "s" + "\t" +"&"  + s_fit_hd[i] + "\t" +"&"  + str(round(s_fit_vmag[i],2)) + "\t" +"&" + str(round(s_fit_bv[i],2)) + "\t" +"&"  + str(round(s_fit_fuvb[i],2)) + "\t" +"&" + str(round(s_fit_x[i],2)) + "\t" +"&" + str(s_fit_y[i]) +  "\\" + "\\" + "\n"
	filename.write(text)
	filename.close()
all_field_x = np.append(all_field_x, s_fit_x)
all_field_y = np.append(all_field_y, s_fit_y)
all_field_bv = np.append(all_field_bv, s_fit_bv)
all_field_hd = np.append(all_field_hd, s_fit_hd)
s_p = np.polyfit(s_fit_x,np.log(s_fit_y),1)
print "Sierchio", "a=", np.exp(s_p[1]), "b=", s_p[0]
slope, intercept, r_value, p_value, std_err = linregress(s_fit_x, s_fit_y)
print "Sierchio", "r^2 (w/o outliers)", r_value**2
slope, intercept, r_value_2, p_value, std_err = linregress(s_xdata, s_ydata)
print "Sierchio", "r^2 (w/ outliers)", r_value_2**2
rho, p = spearmanr(s_fit_x, s_fit_y)
print "Sierchio", 'Spearman coefficient (w/o outliers): ', rho
rho_2, p = spearmanr(s_xdata, s_ydata)
print "Sierchio", 'Spearman coefficient (w/ outliers): ', rho_2
rang = np.argwhere((twin_bv >= 0.703) & (twin_bv <= 0.753))
print "Number in Twin sample", len(rang)
#print 'number of field stars including outliers and copies: ', total


#Get [Fe/H] from Casagrande 2011
#Import  metallicities
f = open('combined_metallicities.txt')
line = f.readlines()[1:]
f.close()
hd_casa = np.array([])
metal = np.array([])
for i in range(len(line)):
    hd_casa = np.append(hd_casa, str(line[i].split()[0]))
    metal = np.append(metal, float(line[i].split()[1]))
#print metallicity histogram.
nans = np.argwhere(np.isnan(metal) == True)
fe_h_plot = np.delete(metal, nans)
print "Range of metallicities: "
print np.min(fe_h_plot), np.max(fe_h_plot)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel( "[Fe/H]", fontsize=26)
plt.annotate('Distribution of metallicities for', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('stars used in FUV-age calibration', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.hist(fe_h_plot, bins='auto')
plt.savefig("f3a.png", bbox_inches='tight',dpi=300)
print '[Fe/H]', np.shape(np.where(np.logical_and(fe_h_plot>=-0.2, fe_h_plot<=0.2)))
print 'nans', len(nans)

#Fit only Field Stars
#0.553-0.603
print "All Fields 0.553-0.603"
rang = np.argwhere((all_field_bv>= 0.553) & (all_field_bv<= 0.603))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
fe_copy = metal[rang]
x_copy = x_copy.flatten()
y_copy = y_copy.flatten()
fe_copy = fe_copy.flatten()
low = np.array([])
med = np.array([])
high = np.array([])
no_good = np.array([])
for i in range(len(rang)):
    if ((fe_copy[i] >=-0.8) and (fe_copy[i]<-0.2)):
        low = np.append(low, i)
    if ((fe_copy[i] >=-0.2) and (fe_copy[i]<0.1)):
        med = np.append(med, i)
    if ((fe_copy[i] >=0.1) and (fe_copy[i]<0.5)):
        high = np.append(high, i)
    if (np.isnan(fe_copy[i]) == True):
        no_good = np.append(no_good, i)
low = low.astype(int)
med = med.astype(int)
high = high.astype(int)
no_good = no_good.astype(int)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
print 'TEST!'
print p
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='.', s=120, label = r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='+', s=120, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=120, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=120, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
#plt.scatter(all_field_x[rang],np.log(all_field_y[rang]),color='black',marker='.', s=70)
plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
plt.plot(x, math.log(22.685) + x*2.190, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(39.671) + x*2.581, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(28.375) + x*2.336, color='black', linestyle='-.', label='Sierchio et al.')
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gigayears)", fontsize=26)
plt.annotate('0.553 $\leq$ B-V $\leq$ 0.603', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.6,0.0)
plt.ylim(-3.5,3)
plt.savefig("f8a.png", bbox_inches='tight',dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
print "All Fields 0.603-0.653"
rang = np.argwhere((all_field_bv>= 0.603) & (all_field_bv<= 0.653))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
fe_copy = metal[rang]
x_copy = x_copy.flatten()
y_copy = y_copy.flatten()
fe_copy = fe_copy.flatten()
low = np.array([])
med = np.array([])
high = np.array([])
no_good = np.array([])
for i in range(len(rang)):
    if ((fe_copy[i] >=-0.8) and (fe_copy[i]<-0.2)):
        low = np.append(low, i)
    if ((fe_copy[i] >=-0.2) and (fe_copy[i]<0.1)):
        med = np.append(med, i)
    if ((fe_copy[i] >=0.1) and (fe_copy[i]<0.5)):
        high = np.append(high, i)
    if (np.isnan(fe_copy[i]) == True):
        no_good = np.append(no_good, i)
low = low.astype(int)
med = med.astype(int)
high = high.astype(int)
no_good = no_good.astype(int)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
x = np.arange(-3.1,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='.', s=120, label = r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='+', s=120, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=120, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=120, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
plt.plot(x, math.log(30.609) + x*1.914, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(39.571) + x*2.130, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(40.097) + x*2.140, color='black', linestyle='-.', label='Sierchio et al.')
plt.plot(x, math.log(44.328) + x*2.247, color='black', linestyle='-', label='Lorenzo et al.')
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.603 $\leq$ B-V $\leq$ 0.653', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.5,0.0)
plt.ylim(-3.5,3)
plt.savefig("f8b.png", bbox_inches='tight',dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
print "All Fields 0.653-0.703"
rang = np.argwhere((all_field_bv>= 0.653) & (all_field_bv<= 0.703))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
fe_copy = metal[rang]
x_copy = x_copy.flatten()
y_copy = y_copy.flatten()
fe_copy = fe_copy.flatten()
low = np.array([])
med = np.array([])
high = np.array([])
no_good = np.array([])
for i in range(len(rang)):
    if ((fe_copy[i] >=-0.8) and (fe_copy[i]<-0.2)):
        low = np.append(low, i)
    if ((fe_copy[i] >=-0.2) and (fe_copy[i]<0.1)):
        med = np.append(med, i)
    if ((fe_copy[i] >=0.1) and (fe_copy[i]<0.5)):
        high = np.append(high, i)
    if (np.isnan(fe_copy[i]) == True):
        no_good = np.append(no_good, i)
low = low.astype(int)
med = med.astype(int)
high = high.astype(int)
no_good = no_good.astype(int)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
x = np.arange(-3.1,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='.', s=120, label = r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='+', s=120, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=120, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=120, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
plt.plot(x, math.log(31.016) + x*1.699, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(31.181) + x*1.665, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(39.556) + x*1.909, color='black', linestyle='-.', label='Sierchio et al.')
plt.plot(x, math.log(19.510) + x*1.110, color='black', linestyle='-', label='Lorenzo et al.')
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.653 $\leq$ B-V $\leq$ 0.703', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-3.1,-0.25)
plt.ylim(-3.5,3)
plt.savefig("f8c.png", bbox_inches='tight',dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
print "All Fields 0.703-0.753"
rang = np.argwhere((all_field_bv>= 0.703) & (all_field_bv<= 0.753))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
fe_copy = metal[rang]
x_copy = x_copy.flatten()
y_copy = y_copy.flatten()
fe_copy = fe_copy.flatten()
low = np.array([])
med = np.array([])
high = np.array([])
no_good = np.array([])
for i in range(len(rang)):
    if ((fe_copy[i] >=-0.8) and (fe_copy[i]<-0.2)):
        low = np.append(low, i)
    if ((fe_copy[i] >=-0.2) and (fe_copy[i]<0.1)):
        med = np.append(med, i)
    if ((fe_copy[i] >=0.1) and (fe_copy[i]<0.5)):
        high = np.append(high, i)
    if (np.isnan(fe_copy[i]) == True):
        no_good = np.append(no_good, i)
low = low.astype(int)
med = med.astype(int)
high = high.astype(int)
no_good = no_good.astype(int)
p = np.polyfit(x_copy,np.log(y_copy),1, full = True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='.', s=120, label = r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='+', s=120, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=120, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=120, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
plt.plot(x, math.log(19.550) + x*1.227, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(35.289) + x*1.648, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(31.079) + x*1.539, color='black', linestyle='-.', label='Sierchio et al.')
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.703 $\leq$ B-V $\leq$ 0.753', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.7,0.0)
plt.ylim(-3.5,3)
plt.savefig("f8d.png", bbox_inches='tight', dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)



#Fit just moving groups
f = open('kevin_mg_fits.txt')
line = f.readlines()[1:]
f.close()
c0 = np.array([])
c1 = np.array([])
c2 = np.array([])
c3 = np.array([])
c4 = np.array([])
c5 = np.array([])
bv_low = np.array([])
bv_high = np.array([])
mg_age = np.array([])
for i in range(len(line)):
    c0 = np.append(c0, float(line[i].split()[1]))
    c1 = np.append(c1, float(line[i].split()[2]))
    c2 = np.append(c2, float(line[i].split()[3]))
    c3 = np.append(c3, float(line[i].split()[4]))
    c4 = np.append(c4, float(line[i].split()[5]))
    c5 = np.append(c5, float(line[i].split()[6]))
    bv_low = np.append(bv_low, float(line[i].split()[7]))
    bv_high = np.append(bv_high, float(line[i].split()[8]))
    mg_age = np.append(mg_age, float(line[i].split()[9]))
def fuv_b(c0,c1,c2,c3,c4,c5,bv_low,bv_high,bv):
    if ((bv > bv_low) & (bv < bv_high)):
        fuvb = c0 + c1*(bv) + c2*(bv**2) + c3*(bv**3) + c4*(bv**4) + c5*(bv**5)
    else:
        fuvb = 0.0001
    return fuvb
print "Moving Groups 0.553-0.603"
bv_mg = 0.578 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
p = np.polyfit(Q_mg,np.log(age_dupe),1, full = True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'rms', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.553 $\leq$ B-V $\leq$ 0.603', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlim(-2.5,0.0)
plt.ylim(-3.5,3)
plt.savefig("f9a.png", bbox_inches='tight', dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(Q_mg, age_dupe)
print "r^2", r_value**2
rho, p = spearmanr(Q_mg, age_dupe)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(Q_mg)
print "Moving Groups 0.603-0.653"
bv_mg = 0.628 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
p = np.polyfit(Q_mg,np.log(age_dupe),1, full = True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'rms', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.603 $\leq$ B-V $\leq$ 0.653', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlim(-2.5,0.0)
plt.ylim(-3.5,3)
plt.savefig("f9b.png", bbox_inches='tight', dpi=300)
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(Q_mg, age_dupe)
print "r^2", r_value**2
rho, p = spearmanr(Q_mg, age_dupe)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(Q_mg)
print "Moving Groups 0.653-0.703"
bv_mg = 0.678 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
p = np.polyfit(Q_mg,np.log(age_dupe),1, full = True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'rms', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.653 $\leq$ B-V $\leq$ 0.703', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlim(-2.5,0.0)
plt.ylim(-3.5,3)
plt.savefig("f9c.png", bbox_inches='tight')
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(Q_mg, age_dupe)
print "r^2", r_value**2
rho, p = spearmanr(Q_mg, age_dupe)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(Q_mg)
print "Moving Groups 0.703-0.753"
bv_mg = 0.728 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
p = np.polyfit(Q_mg,np.log(age_dupe),1, full=True)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'rms', rms
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.703 $\leq$ B-V $\leq$ 0.753', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlim(-3.0,0.0)
#plt.ylim(-3.5,3)
plt.savefig("f9d.png", bbox_inches='tight')
plt.show()
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(Q_mg, age_dupe)
print "r^2", r_value**2
rho, p = spearmanr(Q_mg, age_dupe)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(Q_mg)



#Fit field + moving group + cluster stars
#Import Moving Groups
f = open('kevin_mg_fits.txt')
line = f.readlines()[1:]
f.close()
c0 = np.array([])
c1 = np.array([])
c2 = np.array([])
c3 = np.array([])
c4 = np.array([])
c5 = np.array([])
bv_low = np.array([])
bv_high = np.array([])
mg_age = np.array([])
for i in range(len(line)):
    c0 = np.append(c0, float(line[i].split()[1]))
    c1 = np.append(c1, float(line[i].split()[2]))
    c2 = np.append(c2, float(line[i].split()[3]))
    c3 = np.append(c3, float(line[i].split()[4]))
    c4 = np.append(c4, float(line[i].split()[5]))
    c5 = np.append(c5, float(line[i].split()[6]))
    bv_low = np.append(bv_low, float(line[i].split()[7]))
    bv_high = np.append(bv_high, float(line[i].split()[8]))
    mg_age = np.append(mg_age, float(line[i].split()[9]))
def fuv_b(c0,c1,c2,c3,c4,c5,bv_low,bv_high,bv):
    if ((bv > bv_low) & (bv < bv_high)):
        fuvb = c0 + c1*(bv) + c2*(bv**2) + c3*(bv**3) + c4*(bv**4) + c5*(bv**5)
    else:
        fuvb = 0.0001
    return fuvb
#0.553-0.603
print "All stars 0.553-0.603"
rang = np.argwhere((all_field_bv>= 0.553) & (all_field_bv<= 0.603))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
hyades_bv = coma_bv = 0.578 # get cluster Q
hyades_fuvb = 4.47 + 10.08*hyades_bv
hyades_Q = hyades_fuvb + 29.701*(hyades_bv)**2 - 48.159*(hyades_bv)+ 5.831
hyades_age = 0.625
coma_fuvb = 3.84 + 11.32*coma_bv
coma_Q = coma_fuvb + 29.701*(coma_bv)**2 - 48.159*(coma_bv)+ 5.831
coma_age = 0.500
bv_mg = 0.578 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
	fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
	Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
x_copy = np.append(x_copy, hyades_Q)
x_copy = np.append(x_copy, coma_Q)
x_copy = np.append(x_copy, Q_mg)
y_copy = np.append(y_copy, hyades_age)
y_copy = np.append(y_copy, coma_age)
y_copy = np.append(y_copy, age_dupe)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.scatter(all_field_x[rang],np.log(all_field_y[rang]),color='black',marker='.', s=120, label = 'Field Stars')
print 'test', p
plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
plt.plot(x, math.log(22.685) + x*2.190, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(39.671) + x*2.581, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(28.375) + x*2.336, color='black', linestyle='-.', label='Sierchio et al.')
plt.scatter(coma_Q,np.log(coma_age),color='blue', label = 'Coma Ber Cluster', marker = "+", s = 150)
plt.scatter(hyades_Q,np.log(hyades_age),color='green', label = 'Hyades Cluster', marker = "+", s = 150)
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.553 $\leq$ B-V $\leq$ 0.603', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.5,0.0)
plt.ylim(-3.5,3)
plt.savefig("f6a.png", bbox_inches='tight', dpi=300)
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
#0.603-0.653
print "All stars 0.603-0.653"
rang = np.argwhere((all_field_bv>= 0.603) & (all_field_bv<= 0.653))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
bv_mg = 0.628 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
x_copy = np.append(x_copy, Q_mg)
y_copy = np.append(y_copy, age_dupe)
p = np.polyfit(x_copy,np.log(y_copy),1, full = True)
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.scatter(all_field_x[rang],np.log(all_field_y[rang]),color='black',marker='.', s=120, label = 'Field Stars')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red',label = 'Final Fit')
plt.plot(x, math.log(30.609) + x*1.914, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(39.571) + x*2.130, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(40.097) + x*2.140, color='black', linestyle='-.', label='Sierchio et al.')
plt.plot(x, math.log(44.328) + x*2.247, color='black', linestyle='-', label='Lorenzo et al.')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.603 $\leq$ B-V $\leq$ 0.653', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
#plt.annotate('RMS = 0.42', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 14)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.5,-0.25)
plt.ylim(-3.5,3)
plt.savefig("f6b.png", bbox_inches='tight', dpi=300)
#plt.show()
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
#0.653-0.703
print "All stars 0.653-0.703"
rang = np.argwhere((all_field_bv>= 0.653) & (all_field_bv<= 0.703))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
bv_mg = 0.678 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
x_copy = np.append(x_copy, Q_mg)
y_copy = np.append(y_copy, age_dupe)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.scatter(all_field_x[rang],np.log(all_field_y[rang]),color='black',marker='.', s=120, label = 'Field Stars')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red',label = 'Final Fit')
plt.plot(x, math.log(31.016) + x*1.699, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(31.181) + x*1.665, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(39.556) + x*1.909, color='black', linestyle='-.', label='Sierchio et al.')
plt.plot(x, math.log(19.510) + x*1.110, color='black', linestyle='-', label='Lorenzo et al.')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.653 $\leq$ B-V $\leq$ 0.703', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-2.5,-0.25)
plt.ylim(-3.5,3)
plt.savefig("f6c.png", bbox_inches='tight', dpi=300)
#plt.show()
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)
#0.703-0.753
print "All stars 0.703-0.753"
rang = np.argwhere((all_field_bv>= 0.703) & (all_field_bv<= 0.753))
avg_bv = np.array([])
for i in range(len(rang)):
    avg_bv = np.append(avg_bv, all_field_bv[rang[i]])
print "avg. B-V", np.mean(avg_bv)
bv_mg = 0.728 # get moving group Q
fuvb_mg = np.array([])
Q_mg = np.array([])
for i in range(len(c0)):
    fuvb_mg = np.append(fuvb_mg, fuv_b(c0[i], c1[i], c2[i], c3[i], c4[i], c5[i], bv_low[i], bv_high[i], bv_mg))
    Q_mg = np.append(Q_mg, float(fuvb_mg[i]) + 29.701*(bv_mg)**2 - 48.159*(bv_mg)+ 5.831)
test = np.argwhere(fuvb_mg == 0.0001)
age_dupe = np.copy(mg_age)
fuvb_mg = np.delete(fuvb_mg, test)
Q_mg = np.delete(Q_mg, test)
age_dupe = np.delete(age_dupe, test)
x_copy = all_field_x[rang]
y_copy = all_field_y[rang]
x_copy = np.append(x_copy, Q_mg)
y_copy = np.append(y_copy, age_dupe)
p = np.polyfit(x_copy,np.log(y_copy),1, full=True)
x = np.arange(-3.0,0.0,0.001)
fig = plt.figure(figsize=[10,7])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.scatter(all_field_x[rang],np.log(all_field_y[rang]),color='black',marker='.', s=120, label = 'Field Stars')
plt.plot(x, p[0][1] + (p[0][0] * x),color='red',label = 'Final Fit')
plt.plot(x, math.log(19.550) + x*1.227, color='black', linestyle='dashed', label='Isaacson $\&$ Fischer')
plt.plot(x, math.log(35.289) + x*1.648, color='black', linestyle='dotted', label='Ballering et al.')
plt.plot(x, math.log(31.079) + x*1.539, color='black', linestyle='-.', label='Sierchio et al.')
plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
plt.xlabel("Q", fontsize=26)
plt.ylabel("Log Age (Log Gyr)", fontsize=26)
plt.annotate('0.703 $\leq$ B-V $\leq$ 0.753', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.legend(loc=4, fontsize=15)
plt.xlim(-3.0,-0.25)
plt.ylim(-3.5,3)
plt.savefig("f6d.png", bbox_inches='tight', dpi=300)
plt.show()
rms =  np.sqrt(p[1][0] / len(y_copy))
print 'residual = ', p[1][0]
print 'rms of fit = ', rms
print "a=", np.exp(p[0][1]), "b=", p[0][0]
slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
print "r^2", r_value**2
rho, p = spearmanr(x_copy, y_copy)
print 'Spearman coefficient: ', rho
print 'number in range: ', len(x_copy)

print 'number of field stars without outliers and copies: ', len(all_field_bv)





#fits using moving groups, field stars, and cluster stars
def fit_1(q):
    return np.log(28.449634253325204) + (2.3571960506583163 * q)
def fit_2(q):
    return np.log(34.823542221518814) + (2.0483178580525685 * q)
def fit_3(q):
    return np.log(33.49950280094755) + (1.7530621740118566 * q)
def fit_4(q):
    return np.log(25.629847730865514) + (1.5073484444345047 * q)


est_age = np.array([])
hd_age = np.array([])
bv_age = np.array([])
Q_age = np.array([])
age_age = np.array([])
low_age = np.array([])
high_age = np.array([])
for i in range(len(all_field_bv)):
	if ((all_field_bv[i] < 0.603) and (all_field_bv[i] > 0.553)):
		est_age = np.append(est_age, fit_1(all_field_x[i]))
		bv_age = np.append(bv_age, all_field_bv[i])
		Q_age = np.append(Q_age, all_field_x[i])
		age_age = np.append(age_age, all_field_y[i])
		low_age = np.append(low_age, fit_1(all_field_x[i] - 0.4))
		high_age = np.append(high_age, fit_1(all_field_x[i] + 0.4))
	if ((all_field_bv[i] < 0.653) and (all_field_bv[i] > 0.603)):
		est_age = np.append(est_age, fit_2(all_field_x[i]))
		bv_age = np.append(bv_age, all_field_bv[i])
		Q_age = np.append(Q_age, all_field_x[i])
		age_age = np.append(age_age, all_field_y[i])
		low_age = np.append(low_age, fit_2(all_field_x[i] - 0.4))
		high_age = np.append(high_age, fit_2(all_field_x[i] + 0.4))
	if ((all_field_bv[i] < 0.703) and (all_field_bv[i] > 0.653)):
		est_age = np.append(est_age, fit_3(all_field_x[i]))
		bv_age = np.append(bv_age, all_field_bv[i])
		Q_age = np.append(Q_age, all_field_x[i])
		age_age = np.append(age_age, all_field_y[i])
		low_age = np.append(low_age, fit_3(all_field_x[i] - 0.4))
		high_age = np.append(high_age, fit_3(all_field_x[i] + 0.4))
	if ((all_field_bv[i] < 0.753) and (all_field_bv[i] > 0.703)):
		est_age = np.append(est_age, fit_4(all_field_x[i]))
		bv_age = np.append(bv_age, all_field_bv[i])
		Q_age = np.append(Q_age, all_field_x[i])
		age_age = np.append(age_age, all_field_y[i])
		low_age = np.append(low_age, fit_4(all_field_x[i] - 0.4))
		high_age = np.append(high_age, fit_4(all_field_x[i] + 0.4))

x = np.arange(0,18, 0.0001)
y = x
#x1 = np.arange(0,4, 0.03)
#x2 = np.arange(0,35,0.5)
#y_low = x1*(2.7/0.2)
#y_high = x2*(0.9/4.8)
#fig = plt.figure(figsize=[10,7])
#plt.scatter(np.exp(est_age),all_field_y,color='black')
#plt.scatter(np.exp(low_age), all_field_y , color = 'green', marker = '+', label = 'Q - 0.4')
#plt.scatter(np.exp(high_age), all_field_y, color = 'blue', marker = '+', label = 'Q + 0.4')
#plt.plot(x1, y_low, '.', color = 'green', label = 'Q - 0.4')
#plt.plot(x2, y_high, '.', color = 'blue', label = 'Q + 0.4')
#plt.plot(x, y, color = 'red', label = '1:1')
#plt.ylabel('Quoted Age', fontsize = 20)
#plt.ylim(0,15)
#plt.xlim(0,33)
#plt.xlabel('Estimated Age', fontsize = 20)
#plt.legend()
#plt.savefig("f7a.png", bbox_inches='tight')
#plt.show()

#fig = plt.figure(figsize=[10,7])
#plt.scatter(np.exp(est_age),all_field_y,color='black')
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#plt.ylabel('Literature-Quoted Age', fontsize = 26)
#plt.plot(x, y, color = 'red', label = '1:1')
#plt.ylim(0,15)
#plt.xlim(0,33)
#plt.xlabel('FUV Estimated Age', fontsize = 26)
#plt.legend(fontsize=15)
#plt.ylim(0,2)
#plt.xlim(0,2)
#plt.savefig("f8a.png", bbox_inches='tight', dpi=300)
#plt.show()







#moving group test: Input values from Cochrane 2019 MG FUV-B and age. Compare to my estimates from just MG fits and overal fits
f = open('mg_fuv_age.txt')
line = f.readlines()[1:]
f.close()
mg_name = np.array([])
mg_18_bv = np.array([])
mg_18_fuvb = np.array([])
mg_18_age = np.array([])
for i in range(len(line)):
    mg_name = np.append(mg_name, str(line[i].split()[0]))
    mg_18_bv = np.append(mg_18_bv, float(line[i].split()[1]))
    mg_18_fuvb = np.append(mg_18_fuvb, float(line[i].split()[2]))
    mg_18_age = np.append(mg_18_age, float(line[i].split()[3]))

ignore = np.argwhere(mg_18_bv >0.651) #don't use data without galex fuv
mg_name = np.delete(mg_name, ignore)
mg_18_bv = np.delete(mg_18_bv, ignore)
mg_18_age = np.delete(mg_18_age, ignore)
mg_18_fuvb = np.delete(mg_18_fuvb, ignore)



#estimated age using only moving group fits
def fit_1_mg(q):
    return np.log(80.355) + (3.182 * q)
def fit_2_mg(q):
    return np.log(31.110) + (2.293 * q)
def fit_3_mg(q):
    return np.log(32.982) + (2.169 * q)
def fit_4_mg(q):
    return np.log(22.678) + (1.934 * q)

Q_18_mg = np.array([])
mg_18_est_age = np.array([])
star_18_est_age = np.array([])
for i in range(len(mg_name)): #find q
	Q_18_mg = np.append(Q_18_mg, float(mg_18_fuvb[i]) + 29.701*(mg_18_bv[i])**2 - 48.159*(mg_18_bv[i])+ 5.831)

for i in range(len(mg_name)): #find estimated age given only moving group fit
	if ((0.553 <= mg_18_bv[i]) and (0.603 >= mg_18_bv[i])):
		age_test = fit_1_mg(Q_18_mg[i])
	if ((0.603 <= mg_18_bv[i]) and (0.653 >= mg_18_bv[i])):
		age_test = fit_2_mg(Q_18_mg[i])
		mg_18_est_age = np.append(mg_18_est_age, age_test)
	if ((0.653 <= mg_18_bv[i]) and (0.703 >= mg_18_bv[i])):
		age_test = fit_3_mg(Q_18_mg[i])
	if ((0.703 <= mg_18_bv[i]) and (0.753 >= mg_18_bv[i])):
		age_test = fit_4_mg(Q_18_mg[i])


#fits using moving groups, field stars, and cluster stars
def fit_1_star(q):
	return np.log(27.497) + (2.337 * q)
def fit_2_star(q):
	return np.log(34.904) + (2.030 * q)
def fit_3_star(q):
	return np.log(32.731) + (1.705 * q)
def fit_4_star(q):
	return np.log(21.147) + (1.282 * q)
	
for i in range(len(mg_name)): #find estimated age given only field star fit
	if ((0.553 <= mg_18_bv[i]) and (0.603 >= mg_18_bv[i])):
		age_test = fit_1_star(Q_18_mg[i])
	if ((0.603 <= mg_18_bv[i]) and (0.653 >= mg_18_bv[i])):
		age_test = fit_2_star(Q_18_mg[i])
		star_18_est_age = np.append(star_18_est_age, age_test)
	if ((0.653 <= mg_18_bv[i]) and (0.703 >= mg_18_bv[i])):
		age_test = fit_3_star(Q_18_mg[i])
	if ((0.703 <= mg_18_bv[i]) and (0.753 >= mg_18_bv[i])):
		age_test = fit_4_star(Q_18_mg[i])


all_18_est_age = np.array([])
for i in range(len(mg_name)): #find estimated age given form mg + field + cluster stars
	if ((0.553 <= mg_18_bv[i]) and (0.603 >= mg_18_bv[i])):
		age_test = fit_1(Q_18_mg[i])
	if ((0.603 <= mg_18_bv[i]) and (0.653 >= mg_18_bv[i])):
		age_test = fit_2(Q_18_mg[i])
		all_18_est_age = np.append(all_18_est_age, age_test)
	if ((0.653 <= mg_18_bv[i]) and (0.703 >= mg_18_bv[i])):
		age_test = fit_3(Q_18_mg[i])
	if ((0.703 <= mg_18_bv[i]) and (0.753 >= mg_18_bv[i])):
		age_test = fit_4(Q_18_mg[i])





#fig, ax = plt.subplots(figsize=[10,7])
#fig = plt.figure(figsize=[10,7])
#cmap = plt.cm.get_cmap('viridis')
#plt.scatter(np.exp(mg_18_est_age),mg_18_age,marker = '+', c='black', label = 'Moving Group Fit')
#guess = plt.scatter(np.exp(star_18_est_age),mg_18_age, marker = 'v', c= 'black', label = 'Field Star Fit')
#plt.scatter(np.exp(all_18_est_age),mg_18_age, c = 'black', label = 'Field, Moving Group, and Cluster Star Fit')
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#plt.ylabel('Literature-Quoted Age', fontsize = 20)
#plt.plot(x, y, color = 'red', label = '1:1')
#plt.ylim(0,10)
#plt.xlim(0,9)
#plt.xlabel('FUV Estimated Age', fontsize = 20)
#plt.legend()
#plt.savefig("f12a.png", bbox_inches='tight')
#plt.show()




#Test FUV systematic error
low_2_est_age_1 = np.array([]) # -0.02 FUV
high_2_est_age_1 = np.array([]) # +0.02 FUV
low_5_est_age_1 = np.array([]) # -0.05 FUV
high_5_est_age_1 = np.array([]) # -0.05 FUV
low_2_est_age_2 = np.array([]) # -0.02 FUV
high_2_est_age_2 = np.array([]) # +0.02 FUV
low_5_est_age_2 = np.array([]) # -0.05 FUV
high_5_est_age_2 = np.array([]) # -0.05 FUV
low_2_est_age_3 = np.array([]) # -0.02 FUV
high_2_est_age_3 = np.array([]) # +0.02 FUV
low_5_est_age_3 = np.array([]) # -0.05 FUV
high_5_est_age_3 = np.array([]) # -0.05 FUV
low_2_est_age_4 = np.array([]) # -0.02 FUV
high_2_est_age_4 = np.array([]) # +0.02 FUV
low_5_est_age_4 = np.array([]) # -0.05 FUV
high_5_est_age_4 = np.array([]) # -0.05 FUV
test_age = np.arange(0.25, 17.75, 0.25)


for i in range(len(test_age)):
	Q_test_1 = (np.log(test_age[i]) - np.log(28.449634253325204)) / 2.3571960506583163
	Q_test_2 = (np.log(test_age[i]) - np.log(34.823542221518814)) / 2.0483178580525685
	Q_test_3 = (np.log(test_age[i]) - np.log(33.49950280094755)) / 1.7530621740118566
	Q_test_4 = (np.log(test_age[i]) - np.log(25.629847730865514)) / 1.5073484444345047
	low_2_est_age_1 = np.append(low_2_est_age_1, fit_1(Q_test_1 - 0.02))
	high_2_est_age_1 = np.append(high_2_est_age_1, fit_1(Q_test_1 + 0.02))
	low_5_est_age_1 = np.append(low_5_est_age_1, fit_1(Q_test_1 - 0.05))
	high_5_est_age_1 = np.append(high_5_est_age_1, fit_1(Q_test_1 + 0.05))
	low_2_est_age_2 = np.append(low_2_est_age_2, fit_2(Q_test_2 - 0.02))
	high_2_est_age_2 = np.append(high_2_est_age_2, fit_2(Q_test_2 + 0.02))
	low_5_est_age_2 = np.append(low_5_est_age_2, fit_2(Q_test_2 - 0.05))
	high_5_est_age_2 = np.append(high_5_est_age_2, fit_2(Q_test_2 + 0.05))
	low_2_est_age_3 = np.append(low_2_est_age_3, fit_3(Q_test_3 - 0.02))
	high_2_est_age_3 = np.append(high_2_est_age_3, fit_3(Q_test_3 + 0.02))
	low_5_est_age_3 = np.append(low_5_est_age_3, fit_3(Q_test_3 - 0.05))
	high_5_est_age_3 = np.append(high_5_est_age_3, fit_3(Q_test_3 + 0.05))
	low_2_est_age_4 = np.append(low_2_est_age_4, fit_4(Q_test_4 - 0.02))
	high_2_est_age_4 = np.append(high_2_est_age_4, fit_4(Q_test_4 + 0.02))
	low_5_est_age_4 = np.append(low_5_est_age_4, fit_4(Q_test_4 - 0.05))
	high_5_est_age_4 = np.append(high_5_est_age_4, fit_4(Q_test_4 + 0.05))



fig = plt.figure(figsize=[10,7])
plt.plot(test_age, np.exp(low_2_est_age_2), linestyle='solid', color = '#377eb8', label = 'FUV $\pm$ 0.02')

plt.plot(test_age, np.exp(high_2_est_age_2), linestyle='solid', color = '#377eb8')
plt.plot(test_age, np.exp(low_5_est_age_2), linestyle='solid', color = '#ff7f00', label = 'FUV $\pm$ 0.05')
plt.plot(test_age, np.exp(high_5_est_age_2), linestyle='solid', color = '#ff7f00')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Theoretical Age', fontsize = 26)
plt.annotate('Error in FUV', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
plt.annotate('$B-V = 0.65$', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
plt.plot(x, y, color = '#4daf4a', linestyle= 'dashed',  label = '1:1')
plt.ylim(0,14.0)
plt.xlim(0,14.0)
plt.xlabel('FUV Estimated Age', fontsize = 26)
plt.legend(loc = 4, fontsize=15)
#plt.xticks(np.arange(0, 1, 0.1))
plt.savefig("f10.png", bbox_inches='tight', dpi=300)
plt.show()


#plot fit parameters a vs average b for field+moving group+cluster fits
a = [26.402, 36.053, 32.181, 23.509]
b = [2.314, 2.085, 1.711, 1.408]
avg_bv = [0.585, 0.633, 0.680, 0.724]
quad_fit_a = np.polyfit(avg_bv,a,2)
x = np.arange(0.55,0.75,0.001)
fig = plt.figure(figsize=[10,7])
plt.scatter(avg_bv, a, color = '#377eb8', label = 'FUV $\pm$ 0.02')
plt.plot(x, quad_fit_a[2] + (quad_fit_a[1] * x) + (quad_fit_a[0]*x**2),color='red')
print quad_fit_a[2], quad_fit_a[1], quad_fit_a[0]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('$a$', fontsize = 26)
plt.annotate('Fit Parameter $a$ vs. avg. B-V', xy=(0.05, 0.05), xycoords='axes fraction', fontsize = 20)
#plt.annotate('$B-V = 0.65$', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 14)
#plt.ylim(0,17.5)
plt.xlim(0.55,0.75)
plt.xlabel('Average B-V', fontsize = 26)
plt.savefig("f7a.png", bbox_inches='tight', dpi=300)
quad_fit_b = np.polyfit(avg_bv,b,1)
print quad_fit_b[1], quad_fit_b[0]
fig = plt.figure(figsize=[10,7])
plt.scatter(avg_bv, b, color = '#377eb8', label = 'FUV $\pm$ 0.02')
plt.plot(x, quad_fit_b[1] + (quad_fit_b[0]*x),color='red')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('$b$', fontsize = 26)
plt.annotate('Fit Parameter $b$ vs. avg. B-V', xy=(0.05, 0.05), xycoords='axes fraction', fontsize = 20)
#plt.annotate('$B-V = 0.65$', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 14)
#plt.ylim(0,17.5)
plt.xlim(0.55,0.75)
plt.xlabel('Average B-V', fontsize = 26)
plt.savefig("f7b.png", bbox_inches='tight', dpi=300)
