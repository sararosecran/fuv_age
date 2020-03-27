#gets metallicities from Casagrande 2011 to be used in the histogram showing the [Fe/H] spread and in fitting plots.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def feH():
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

	return metal 


