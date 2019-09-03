import pydicom
import numpy as np
import pickle

from file_utils import find_prefixed_files, find_dicom_directory

def load_rtdose_files(files):
	dosegrids = []
	voxelsize = None
	for i, fname in enumerate(files):
		data = pydicom.dcmread(fname)
		dosegrids.append(data.pixel_array.astype(float))
		dosegrids[i] *= data.DoseGridScaling
	
		voxelsize = data.PixelSpacing[0] * data.PixelSpacing[1] * data.SliceThickness
	
	voxelsize /= 1000. #Rescale from mm^3 to cm^3
	return dosegrids, voxelsize

def check_organ_dose(mask, dosegrid, voxelsize, dosechecks=[5, 10, 15, 20]):
	print(np.max(dosegrid))
	dosemask = dosegrid[mask['Mask']]
	maxdose = np.max(dosemask)
	meandose = np.sum(dosemask) / np.sum(mask['Mask'])
	volume = np.sum(mask['Mask']) * voxelsize
	intdose = meandose * volume

	dosevols = [np.sum(dosemask >= j) * voxelsize for j in dosechecks]

	return maxdose, meandose, volume, intdose, dosevols

def check_all_organs(masks, dosegrid, voxelsize, dosechecks=[5, 10, 15, 20]):
	for mask in masks:
		mx, mn, vl, itd, dv = check_organ_dose(mask, dosegrid, voxelsize, dosechecks)
		outstr = '%20s,%g,%g,%g,%g' % (mask['Name'], mx, mn, vl, itd)
		for i in dv:
			outstr += ',%g' % (i / vl)
		print(outstr)

if __name__=='__main__':
	directory = '../data/AA'
	dcm_directory = find_dicom_directory(directory)

	rtdose_files = find_prefixed_files(dcm_directory, 'RTDOSE')
	dosegrids, voxelsize = load_rtdose_files(rtdose_files)

	with open('../data/AA/masks.pickle', 'rb') as infile:
		masks = pickle.load(infile)
	
	check_all_organs(masks, np.sum(np.array(dosegrids), axis=0).transpose(2, 1, 0), voxelsize)
