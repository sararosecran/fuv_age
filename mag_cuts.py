# reduces the sample to stars within 4.3<Mv<5.3 (ie those with solar-like magnitudes.
import numpy as np

def mag_cut(i_parallax, i_vmag,i_age, i_bmag, i_fuv, i_hd, b_parallax, b_vmag, b_age, b_bmag, b_fuv, b_hd, s_parallax, s_vmag, s_age, s_bmag, s_fuv, s_hd, twin_parallax, twin_vmag, twin_age, twin_bmag, twin_fuv, twin_hd):
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
	mv = np.delete(mv, ignore)
	i_hd = np.delete(i_hd, ignore)
	ignore = np.argwhere(mv >= 5.3)
	i_age = np.delete(i_age, ignore)
	i_bmag = np.delete(i_bmag, ignore)
	i_vmag = np.delete(i_vmag, ignore)
	i_fuv = np.delete(i_fuv, ignore)
	i_parallax = np.delete(i_parallax, ignore)
	i_hd = np.delete(i_hd, ignore)

	mv = np.delete(mv, ignore)
	b_parallax = b_parallax / 1000. # convert to arcsec
	mv = np.array([])
	for	i in range(len(b_parallax)):
		mv = np.append(mv, b_vmag[i] + 5.0*(1.0 + np.log10(b_parallax[i])))
	ignore = np.argwhere(mv <= 4.3) #cuts for Ballering
	b_age = np.delete(b_age, ignore)
	b_bmag = np.delete(b_bmag, ignore)
	b_vmag = np.delete(b_vmag, ignore)
	b_fuv = np.delete(b_fuv, ignore)
	b_parallax = np.delete(b_parallax, ignore)
	mv = np.delete(mv, ignore)
	b_hd = np.delete(b_hd, ignore)
	ignore = np.argwhere(mv >= 5.3)
	b_age = np.delete(b_age, ignore)
	b_bmag = np.delete(b_bmag, ignore)
	b_vmag = np.delete(b_vmag, ignore)
	b_fuv = np.delete(b_fuv, ignore)
	b_hd = np.delete(b_hd, ignore)
	b_parallax = np.delete(b_parallax, ignore)
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
	mv = np.delete(mv, ignore)
	s_hd = np.delete(s_hd, ignore)
	ignore = np.argwhere(mv >= 5.3)
	s_age = np.delete(s_age, ignore)
	s_bmag = np.delete(s_bmag, ignore)
	s_vmag = np.delete(s_vmag, ignore)
	s_fuv = np.delete(s_fuv, ignore)
	s_parallax = np.delete(s_parallax, ignore)
	s_hd = np.delete(s_hd, ignore)
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
	mv = np.delete(mv, ignore)
	twin_hd = np.delete(twin_hd, ignore)
	ignore = np.argwhere(mv >= 5.3)
	twin_age = np.delete(twin_age, ignore)
	twin_bmag = np.delete(twin_bmag, ignore)
	twin_vmag = np.delete(twin_vmag, ignore)
	twin_fuv = np.delete(twin_fuv, ignore)
	twin_parallax = np.delete(twin_parallax, ignore)
	twin_hd = np.delete(twin_hd, ignore)
	mv = np.delete(mv, ignore)
	return i_parallax, i_vmag,i_age, i_bmag, i_fuv, i_hd, b_parallax, b_vmag, b_age, b_bmag, b_fuv, b_hd, s_parallax, s_vmag, s_age, s_bmag, s_fuv, s_hd, twin_parallax, twin_vmag, twin_age, twin_bmag, twin_fuv, twin_hd
