# fits the four stellar samples individually to determine outliers and then collects all non-outlier stars into arrays.
import numpy as np
from scipy.stats import spearmanr, linregress
import math


def combined_data(i_bv, i_Q, i_age, i_hd, i_fuvb, i_vmag, b_bv, b_Q, b_age, b_hd, b_fuvb, b_vmag, s_bv, s_Q, s_age, s_hd, s_fuvb, s_vmag, twin_bv, twin_Q, twin_age, twin_hd, twin_fuvb, twin_vmag):

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
	i_outlier = np.argwhere(dist > 0.62) #these outlier cutoffs were determine in jupyter notebooks
	i_fit_x = np.delete(i_xdata, i_outlier)
	i_fit_y = np.delete(i_ydata, i_outlier)
	i_fit_hd = np.delete(i_hd_used, i_outlier)
	i_fit_bv = np.delete(i_bv_used, i_outlier)
	i_fit_fuvb = np.delete(i_fuvb_used, i_outlier)
	i_fit_vmag = np.delete(i_vmag_used, i_outlier)
	print "Number in Isaacson sample ignoring outliers", len(i_fit_x)
	avg_bv = np.array([])
	for i in range(len(i_fit_x)):
		avg_bv = np.append(avg_bv, i_fit_bv[i])
	print "avg. B-V ignoring outliers", np.mean(avg_bv)
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
	print "Number in Ballering sample excluding outliers ", len(b_xdata)
	avg_bv = np.array([])
	for i in range(len(b_hd_used)):
		avg_bv = np.append(avg_bv, b_bv_used[i])
	print "avg. B-V", np.mean(avg_bv)
	for i in range(len(b_hd_used)):
		filename = open('field_stars_fit.txt', "a")
		text = "b" + "\t" +"&"  + b_hd_used[i] + "\t" +"&"  + str(round(b_vmag_used[i],2)) + "\t" +"&" + str(round(b_bv_used[i],2)) + "\t" +"&"  + str(round(b_fuvb_used[i],2)) + "\t" +"&" + str(round(b_xdata[i],2)) + "\t" +"&" + str(b_ydata[i]) +  "\\"+ "\\" +  "\n"
		filename.write(text)
		filename.close()
	b_p	= np.polyfit(b_xdata,np.log(b_ydata),1)
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
	print "Number in Sierchio sample exclusing outliers", len(s_fit_x)
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
	i_p	= np.polyfit(i_xdata,np.log(i_ydata),1)
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
	print "Number in Isaacson sample excluding outliers", len(i_fit_x)
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
	print "Number in Ballering sample excluding outliers", len(b_fit_x)
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
	print "Number in Sierchio sample excluding outliers", len(s_fit_x)
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
	print "Number in Twin sample excluding outliers ", len(twin_fit_x)
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
	print "Number in Isaacson sample excluding outliers", len(i_fit_x)
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
	print "Number in Ballering sample excluding outliers", len(b_fit_x)
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
	print "Number in Sierchio sample excluding outliers", len(s_fit_x)
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
	print "Number in Twin sample excluding outliers ", len(twin_fit_x)
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

	return all_field_x, all_field_y, all_field_bv, all_field_hd
