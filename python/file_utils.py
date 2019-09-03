import os
import re
import matplotlib.pyplot as plt

def find_prefixed_file(directory, prefix):
	for i in os.listdir(directory):
		if re.search('^'+prefix, i) is not None:
			return os.path.join(directory, i)
	return None

def find_prefixed_files(directory, prefix):
	files = []
	for i in os.listdir(directory):
		if re.search('^'+prefix, i) is not None:
			files.append(os.path.join(directory, i))
	return files

def find_dicom_directory(directory):
	dcm_directory = None
	ct_prefix = 'CT'
	dose_prefix = 'RTDOSE'

	for i in os.listdir(directory):
		if os.path.isdir(os.path.join(directory, i)):
			dcm_directory = os.path.join(directory, i)

	return dcm_directory

def implay(cube):
	plt.ion()
	plt.figure()
	plt.show()

	for i in range(cube.shape[2]):
		plt.clf()
		plt.title('Frame %d' % (i+1))
		plt.imshow(cube[:, :, i])
		plt.colorbar()
		plt.pause(0.01)
	
	plt.ioff()
