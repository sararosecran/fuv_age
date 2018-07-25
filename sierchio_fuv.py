from astroquery.mast import Observations
from astroquery.mast import Catalogs
import numpy as np

f = open('sierchio_new.txt', "r")
line = f.readlines()[1:]
f.close()
name = np.array([])
for i in range(len(line)):
	name = np.append(name, str(line[i].split()[0]))

for i in range(len(name)-1):
	temp_name = "HD" + name[i]
	catalogData = Catalogs.query_object(temp_name, radius=0.004 , catalog="Galex") #searches within 0.2 arcmin
	length = len(catalogData)
	if length > 1:
		fuv_array = np.array([])
		err_array = np.array([])
		for j in range(length):
			if catalogData[j][5] != 1:
				fuv_array = np.append(fuv_array, catalogData[j][10])
				err_array = np.append(err_array, catalogData[j][11])
		fuv_mag = np.mean(fuv_array)
		fuv_mag_err = np.mean(err_array)
	if length < 1:
		fuv_mag = 'n/a'
		fuv_mag_err = 'n/a'
	if length == 1:
		if catalogData[0][5] == 1:
			fuv_mag = 'n/a'
			fuv_mag_err = 'n/a'
		else:
			fuv_mag = catalogData[0][10]
			fuv_mag_err = catalogData[0][11]

	test_file = open('test_file.txt', "a")
	text = str(temp_name) + "\t" +  str(fuv_mag) + "\t" + str(fuv_mag_err) + "\n"
	test_file.write(text)
	test_file.close()


