# code imports stars and the properties from Isaacson, Ballering, Sierchio, and Solar Twin catalogs

import numpy as np

def import_isaacson():
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
	return i_hd, i_bmag, i_fuv, i_vmag, i_age, i_parallax

def import_ballering():
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
	return b_hd, b_vmag, b_bmag, b_fuv, b_age, b_parallax, b_rhk

def import_sierchio():
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
	return s_hd, s_vmag, s_bmag, s_fuv, s_age, s_parallax, s_fuv

def import_solar_twin():
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
	return twin_hd, twin_vmag, twin_bmag, twin_fuv, twin_age, twin_parallax


