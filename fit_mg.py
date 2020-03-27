# fits age vs Q for moving groups.#Fit just moving groups
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

def moving_groups():

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
	rms =  np.sqrt(p[1][0] / len(age_dupe))
	print 'rms', rms
	x = np.arange(-3.0,0.0,0.001)
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 400)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.55 $\leq$ B-V $<$ 0.60', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xlim(-2.5,0.0)
	plt.ylim(-2.0,2.5)
	plt.savefig("f6a.png", bbox_inches='tight', dpi=300)
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
	rms =  np.sqrt(p[1][0] / len(age_dupe))
	print 'rms', rms
	x = np.arange(-3.0,0.0,0.001)
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 400)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.60 $\leq$ B-V $<$ 0.65', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xlim(-2.5,0.0)
	plt.ylim(-2.0,2.5)
	plt.savefig("f6b.png", bbox_inches='tight', dpi=300)
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
	rms =  np.sqrt(p[1][0] / len(age_dupe))
	print 'rms', rms
	x = np.arange(-3.0,0.0,0.001)
	fig = plt.figure(figsize=[10,7])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.plot(x, p[0][1] + (p[0][0] * x),color='red')
	plt.scatter(Q_mg, np.log(age_dupe) ,color='m', marker = "+", s = 400)
	plt.xlabel("Q", fontsize=26)
	plt.ylabel("Log Age (Log Gyr)", fontsize=26)
	plt.annotate('0.65 $\leq$ B-V $\leq$ 0.70', xy=(0.05, 0.95), xycoords='axes fraction', fontsize = 20)
	plt.annotate('Moving Groups Only', xy=(0.05, 0.90), xycoords='axes fraction', fontsize = 20)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xlim(-2.5,0.0)
	plt.ylim(-2.0,2.5)
	plt.savefig("f6c.png", bbox_inches='tight')
	plt.show()
	print "a=", np.exp(p[0][1]), "b=", p[0][0]
	slope, intercept, r_value, p_value, std_err = linregress(Q_mg, age_dupe)
	print "r^2", r_value**2
	rho, p = spearmanr(Q_mg, age_dupe)
	print 'Spearman coefficient: ', rho
	print 'number in range: ', len(Q_mg)
	return
