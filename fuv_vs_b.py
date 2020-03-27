# code plots two-color diagram of FUV-B vs B-V using the stars form for catalogs but does not use duplicate stars

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def two_color_plot(i_hd, i_vmag, i_bmag, i_fuv,s_hd, s_vmag, s_bmag, s_fuv, b_hd, b_vmag, b_bmag, b_fuv, twin_hd, twin_vmag, twin_bmag, twin_fuv):

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
	#plt.show()
	fig.savefig("f4a.png", bbox_inches='tight', dpi=300)
	return
