import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def rhk_plot(b_age, b_rhk):
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
	return
