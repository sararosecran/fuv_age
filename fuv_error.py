# tests the systematic error contributed by FUV measurements.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


#fits using moving groups, field stars, and cluster stars
def fit_1(q):
    return np.log(28.449634253325204) + (2.3571960506583163 * q)
def fit_2(q):
    return np.log(34.823542221518814) + (2.0483178580525685 * q)
def fit_3(q):
    return np.log(33.49950280094755) + (1.7530621740118566 * q)
def fit_4(q):
    return np.log(25.629847730865514) + (1.5073484444345047 * q)


def fuv_test():
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
	
	x = np.arange(0,18, 0.0001)
	y = x

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
	plt.savefig("f9a.png", bbox_inches='tight', dpi=300)
	plt.show()
	return
