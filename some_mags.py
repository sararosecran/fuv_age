# finds B-V, FUV-B and Q for all catalogs
import numpy as np

#find B-V
def b_v(i_bmag, i_vmag, b_bmag, b_vmag, s_bmag, s_vmag, twin_bmag, twin_vmag):
	i_bv = np.array([]) #Isaacson
	for i in range(len(i_bmag)):
		i_bv = np.append(i_bv, i_bmag[i] - i_vmag[i])
	b_bv = np.array([]) #Ballering
	for i in range(len(b_bmag)):
		b_bv = np.append(b_bv, b_bmag[i] - b_vmag[i])
	s_bv = np.array([]) #Sierchio
	for i in range(len(s_bmag)):
		s_bv = np.append(s_bv, s_bmag[i] - s_vmag[i])
	twin_bv = np.array([]) #Solar Twin
	for i in range(len(twin_bmag)):
		twin_bv = np.append(twin_bv, twin_bmag[i] - twin_vmag[i])
	return (i_bv, b_bv, s_bv, twin_bv)

#find FUV-B
def fuv_b(i_bmag, i_fuv, b_bmag, b_fuv, s_bmag, s_fuv, twin_bmag, twin_fuv):
	i_fuvb = np.array([]) #Isaacson
	for i in range(len(i_bmag)):
		i_fuvb = np.append(i_fuvb, i_fuv[i] - i_bmag[i])
	b_fuvb = np.array([]) #Ballering
	for i in range(len(b_bmag)):
		b_fuvb = np.append(b_fuvb, b_fuv[i] - b_bmag[i])
	s_fuvb = np.array([]) #Sierchio
	for i in range(len(s_bmag)):
		s_fuvb = np.append(s_fuvb, s_fuv[i] - s_bmag[i])
	twin_fuvb = np.array([]) #Solar Twin
	for i in range(len(twin_bmag)):
		twin_fuvb = np.append(twin_fuvb, twin_fuv[i] - twin_bmag[i])
	return (i_fuvb, b_fuvb, s_fuvb, twin_fuvb)

#find Q
def Q(i_bmag, i_fuvb, i_bv, b_bmag, b_fuvb, b_bv, s_bmag, s_fuvb, s_bv, twin_bmag, twin_fuvb, twin_bv):
	i_Q = np.array([]) #Isaacson
	for i in range(len(i_bmag)):
		i_Q = np.append(i_Q, i_fuvb[i] + 29.701*(i_bv[i])**2 - 48.159*(i_bv[i]) + 5.831)
	b_Q = np.array([]) #Ballering
	for i in range(len(b_bmag)):
		b_Q = np.append(b_Q, b_fuvb[i] + 29.701*(b_bv[i])**2 - 48.159*(b_bv[i]) + 5.831)
	s_Q = np.array([]) #Sierchio
	for i in range(len(s_bmag)):
		s_Q = np.append(s_Q, s_fuvb[i] + 29.701*(s_bv[i])**2 - 48.159*(s_bv[i]) + 5.831)
	twin_Q = np.array([]) #Solar Twin
	for i in range(len(twin_bmag)):
		twin_Q = np.append(twin_Q, twin_fuvb[i] + 29.701*(twin_bv[i])**2 - 48.159*(twin_bv[i]) + 5.831)
	return (i_Q, b_Q, s_Q, twin_Q)
