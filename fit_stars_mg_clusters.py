# fits age vs Q for field stars, moving groups, and clusters.
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
import math
from scipy.optimize import curve_fit
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def fuv_b(c0,c1,c2,c3,c4,c5,bv_low,bv_high,bv):
    if ((bv > bv_low) & (bv < bv_high)):
        fuvb = c0 + c1*(bv) + c2*(bv**2) + c3*(bv**3) + c4*(bv**4) + c5*(bv**5)
    else:
        fuvb = 0.0001
    return fuvb

def two_lines_1(x, a, b, c, d):
    out = np.empty_like(x)
    mask = x < -0.83
    out[mask] = a*x[mask] + b
    out[~mask] = c*x[~mask] + d
    return out

def line_one(x, a, b):
    return a*x + b


def fit_all(all_field_x, all_field_y, all_field_bv, metal):
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
	fe_copy = metal[rang]
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
	pw0 = (2.0,-4.0,0.2,1.0)
	pw, cov = curve_fit(two_lines_1, x_copy, np.log(y_copy), pw0)
	print 'rms: ', np.sqrt(np.mean((np.log(y_copy)-two_lines_1(x_copy, *pw))**2))
	b1 = pw[1]
	b2 = pw[3]
	m1 = pw[0]
	m2 = pw[2]
	#find intersect
	x_intersect = (b1-b2) / (m2-m1)
	y_intersect = (m1*x_intersect) + b1
	x = np.arange(-3.5,0.0,0.0001)
	out = np.empty_like(x)
	mask = x < x_intersect
	out[mask] = m1*x[mask] + b1
	out[~mask] = m2*x[~mask] + b2
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.plot(x, out, '-', color = 'red', label = 'Final Fit')
	plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='o', s=150, label= r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
	plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='v', s=150, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
	plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=150, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
	plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=150, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
	plt.scatter(coma_Q,np.log(coma_age),color='blue', label = 'Coma Ber Cluster', marker = "+", s = 150)
	plt.scatter(hyades_Q,np.log(hyades_age),color='green', label = 'Hyades Cluster', marker = "+", s = 150)
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.55 $\leq$ B-V $<$ 0.60', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0.0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f7a.png", bbox_inches='tight', dpi=300)
	print "m1=", m1, "b1=", np.exp(b1), "m2=", m2, "b2=", np.exp(b2)
	print 'intersect ', (x_intersect, y_intersect)
	#Find R^2
	pw_rms_1 = np.delete(pw,[2,3])
	pw_rms_2 = np.delete(pw,[0,1])
	x_fit = np.array([])
	y_fit = np.array([])
	for i in range(len(x_copy)):
		if x_copy[i] < x_intersect:
			x_fit = np.append(x_fit, x_copy[i])
			y_fit = np.append(y_fit, y_copy[i])
	residuals_1 = np.log(y_fit)- line_one(x_fit, pw_rms_1[0], pw_rms_1[1])
	ss_tot_1 = np.sum((np.log(y_fit)-np.mean(np.log(y_fit)))**2)
	x_fit = np.array([])
	y_fit = np.array([])
	for i in range(len(x_copy)):
		if x_copy[i] > x_intersect:
			x_fit = np.append(x_fit, x_copy[i])
			y_fit = np.append(y_fit, y_copy[i])
	residuals_2 = np.log(y_fit)- line_one(x_fit, pw_rms_2[0], pw_rms_2[1])
	ss_res = np.sum(residuals_1**2) + np.sum(residuals_2**2)
	ss_tot_2 = np.sum((np.log(y_fit)-np.mean(np.log(y_fit)))**2)
	ss_tot = ss_tot_1 + ss_tot_2
	r_squared = 1 - (ss_res / ss_tot)
	print "R^2: ", r_squared
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
	fe_copy = metal[rang]
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
	pw0 = (2.0,-4.0,0.2,1.0)
	pw, cov = curve_fit(two_lines_1, x_copy, np.log(y_copy), pw0)
	print 'rms: ', np.sqrt(np.mean((np.log(y_copy)-two_lines_1(x_copy, *pw))**2))
	b1 = pw[1]
	b2 = pw[3]
	m1 = pw[0]
	m2 = pw[2]
	#find intersect
	x_intersect = (b1-b2) / (m2-m1)
	y_intersect = (m1*x_intersect) + b1
	x = np.arange(-3.5,0.0,0.001)
	out = np.empty_like(x)
	mask = x < x_intersect
	out[mask] = m1*x[mask] + b1
	out[~mask] = m2*x[~mask] + b2
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.plot(x, out, '-', color = 'red', label = 'Final Fit')
	plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='o', s=150, label= r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
	plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='v', s=150, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
	plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=150, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
	plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=150, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.60 $\leq$ B-V $<$ 0.65', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	#plt.annotate('RMS = 0.42', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 14)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0.0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f7b.png", bbox_inches='tight', dpi=300)
	print "m1=", m1, "b1=", np.exp(b1), "m2=", m2, "b2=", np.exp(b2)
	print 'intersect ', (x_intersect, y_intersect)
	#Find R^2
	pw_rms_1 = np.delete(pw,[2,3])
	pw_rms_2 = np.delete(pw,[0,1])
	x_fit = np.array([])
	y_fit = np.array([])
	for i in range(len(x_copy)):
		if x_copy[i] < x_intersect:
			x_fit = np.append(x_fit, x_copy[i])
			y_fit = np.append(y_fit, y_copy[i])
	residuals_1 = np.log(y_fit)- line_one(x_fit, pw_rms_1[0], pw_rms_1[1])
	ss_tot_1 = np.sum((np.log(y_fit)-np.mean(np.log(y_fit)))**2)
	x_fit = np.array([])
	y_fit = np.array([])
	for i in range(len(x_copy)):
		if x_copy[i] > x_intersect:
			x_fit = np.append(x_fit, x_copy[i])
			y_fit = np.append(y_fit, y_copy[i])
	residuals_2 = np.log(y_fit)- line_one(x_fit, pw_rms_2[0], pw_rms_2[1])
	ss_res = np.sum(residuals_1**2) + np.sum(residuals_2**2)
	ss_tot_2 = np.sum((np.log(y_fit)-np.mean(np.log(y_fit)))**2)
	ss_tot = ss_tot_1 + ss_tot_2
	r_squared = 1 - (ss_res / ss_tot)
	print "R^2: ", r_squared
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
	x = np.arange(-3.5,0.0,0.001)
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='o', s=150, label= r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
	plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='v', s=150, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
	plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=150, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
	plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=150, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
	plt.plot(x, p[0][1] + (p[0][0] * x),color='red',label = 'Final Fit')
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', label = 'Moving Groups', marker = "+", s = 150)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.65 $\leq$ B-V $\leq$ 0.70', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0.0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f7c.png", bbox_inches='tight', dpi=300)
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
	return
