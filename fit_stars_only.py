# fits age vs Q for just the field stars.
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
import math
from scipy.optimize import curve_fit
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def two_lines_1(x, a, b, c, d):
    out = np.empty_like(x)
    mask = x < -0.83
    out[mask] = a*x[mask] + b
    out[~mask] = c*x[~mask] + d
    return out

def line_one(x, a, b):
    return a*x + b


def field_stars(all_field_bv, all_field_x, all_field_y, metal):
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
	#fit two lines
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
	x = np.arange(-3.1,0.0,0.0001)
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
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.55 $\leq$ B-V $<$ 0.60', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	#plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f5a.png", bbox_inches='tight',dpi=300)
	plt.show()
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
	#fit two lines
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
	x = np.arange(-3.1,0.0,0.001)
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
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.60 $\leq$ B-V $<$ 0.65', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f5b.png", bbox_inches='tight',dpi=300)
	plt.show()
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
	plt.scatter(all_field_x[rang[low]],np.log(all_field_y[rang[low]]),color='#377eb8',marker='o', s=150, label = r'$-0.8\leq [\mathrm{Fe}/\mathrm{H}] <-0.2$')
	plt.scatter(all_field_x[rang[med]],np.log(all_field_y[rang[med]]),color='#ff7f00',marker='v', s=150, label = r'$-0.2\leq [\mathrm{Fe}/\mathrm{H}] <0.1$')
	plt.scatter(all_field_x[rang[high]],np.log(all_field_y[rang[high]]),color='#4daf4a',marker='d', s=150, label = r'$0.1\leq [\mathrm{Fe}/\mathrm{H}] <0.5$')
	plt.scatter(all_field_x[rang[no_good]],np.log(all_field_y[rang[no_good]]),color='#f781bf',marker='*', s=150, label = r'no $[\mathrm{Fe}/\mathrm{H}]$')
	plt.plot(x, p[0][1] + (p[0][0] * x),color='red', label = 'Final Fit')
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.65 $\leq$ B-V $\leq$ 0.70', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Field Stars Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.legend(loc=4, fontsize=15)
	plt.xlim(-3.1,0)
	plt.ylim(-3.5,2.5)
	plt.savefig("f5c.png", bbox_inches='tight',dpi=300)
	plt.show()
	print "a=", np.exp(p[0][1]), "b=", p[0][0]
	slope, intercept, r_value, p_value, std_err = linregress(x_copy, y_copy)
	print "r^2", r_value**2
	rho, p = spearmanr(x_copy, y_copy)
	print 'Spearman coefficient: ', rho
	print 'number in range: ', len(x_copy)
