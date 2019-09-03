import numpy as np
import pydicom
import os
import pickle

from file_utils import find_prefixed_file, find_dicom_directory, implay, find_prefixed_files
from sys import exit

basic_mask_dicts = [
	{'NameStrings': ['other', 'organs'], 'GV': False, 'Stationary': False, 'CardiacOutput': -1},
	{'NameStrings': ['lung', 'total'], 'GV': False, 'Stationary': False, 'CardiacOutput': 0.025},
	{'NameStrings': ['aorta'], 'GV': True, 'Stationary': False, 'CardiacOutput': 1.},
	{'NameStrings': ['pa'], 'GV': True, 'Stationary': False, 'CardiacOutput': 1.},
	{'NameStrings': ['vc'], 'GV': True, 'Stationary': False, 'CardiacOutput': 1.},
	{'NameStrings': ['thoracic', 'spine'], 'GV': False, 'Stationary': True, 'CardiacOutput': 0.},
]

def find_matching_contour_idx(contours, name_strings):
	for i, nm in enumerate(contours['ROIName']):
		lnm = nm.lower()
		found = True
		for j in name_strings:
			if not j.lower() in lnm:
				found = False
		if found:
			return i
	
	return -1

def get_conversion_grids(ct_info, dose_info, dim_vol, dim_dos):
	#Get voxel dimensions
	dim_voxd = np.array([dose_info.PixelSpacing[0], dose_info.PixelSpacing[1], dose_info.SliceThickness])
	dim_voxv = np.array([ct_info.PixelSpacing[0], ct_info.PixelSpacing[1], ct_info.SliceThickness])

	#Get image corners
	corner_d = np.array([dose_info.ImagePositionPatient])
	corner_v = np.array([ct_info.ImagePositionPatient])

	posd = np.transpose(np.array(np.meshgrid(
		np.arange(dim_dos[0]),
		np.arange(dim_dos[1]),
		np.arange(dim_dos[2]))), (1, 2, 3, 0))
	
	posv = (corner_d - corner_v + posd * dim_voxd) / dim_voxv
	posv = posv.astype(int)

	valid_voxel = np.logical_and(posv[:,:,:,0] >= 0, posv[:,:,:,1] >= 0)
	valid_voxel = np.logical_and(valid_voxel, posv[:,:,:,2] >= 0)
	valid_voxel = np.logical_and(valid_voxel, posv[:,:,:,0] < dim_vol[1])
	valid_voxel = np.logical_and(valid_voxel, posv[:,:,:,1] < dim_vol[0])
	valid_voxel = np.logical_and(valid_voxel, posv[:,:,:,2] < dim_vol[2])

	return posv[valid_voxel], posd[valid_voxel]

def layer_size(mask):
	num = np.sum(mask, axis=2)
	num = num[num > 0]
	return np.sum(num) / len(num)

def get_first_CT_frame(ct_infos):
	first = ct_infos[0]
	for i in ct_infos:
		if float(i.ImagePositionPatient[2]) < float(first.ImagePositionPatient[2]):
			first = i
	return first

'''
@param contours - The result from structure_loading.load_structures
@param ct_info - Header information from a CT file
@param dose_info - Header information from a DOSE file
@param mask_dicts - Information about which masks to include
	Entries in mask_dicts have the format:
		NameStrings - search for ROIs in contours which have these name strings (use lowercase)
		GV - True if the organ is a great vessel
		Stationary - True if the organ is stationary (thoracic spine)
		Cardiac Output - Percentage of total cardiac output (0-1)
'''
def mask_generation(
	contours, 
	ct_infos,
	dose_info,
	mask_dicts=basic_mask_dicts):

	ct_info = get_first_CT_frame(ct_infos)

	dim_vol = np.array([ct_info.Rows, ct_info.Columns, len(ct_infos)])
	dim_dos = np.array([dose_info.Columns, dose_info.Rows, dose_info.NumberOfFrames])

	print(dim_dos, dim_vol)

	posv, posd = get_conversion_grids(ct_info, dose_info, dim_vol, dim_dos)

	posv = np.ravel_multi_index(posv.transpose(), dim_vol)
	posd = np.ravel_multi_index(posd.transpose(), dim_dos)

	masks = []
	other_organs_ind = -1
	used_voxels = None
	#mask_dicts.reverse()
	for i, dct in enumerate(mask_dicts):
		contour_idx = find_matching_contour_idx(contours, dct['NameStrings'])
		mdict = {}
		mdict['Name'] = contours['ROIName'][contour_idx]
		mdict['GV'] = dct['GV']
		mdict['Stationary'] = dct['Stationary']
		mdict['CardiacOutput'] = dct['CardiacOutput']

		print('Creating mask for organ %s' % mdict['Name'])

		seg = contours['Segmentation'][contour_idx]
		seg = seg.flatten()
		mdict['Mask'] = np.zeros(dim_dos.prod(), dtype=bool)
		mdict['Mask'][posd] = seg[posv]
		mdict['Mask'] = mdict['Mask'].reshape(dim_dos)
		print(np.sum(mdict['Mask']))

		mdict['LayerSize'] = layer_size(mdict['Mask'])
		masks.append(mdict)

		if masks[i]['CardiacOutput'] == -1:
			other_organs_ind = i
		else:
			if used_voxels is None:
				used_voxels = np.copy(masks[i]['Mask'])
			else:
				used_voxels = np.logical_or(used_voxels, masks[i]['Mask'])
	
	#Now remove duplicated voxels in other organs
	if other_organs_ind != -1:
		masks[other_organs_ind]['Mask'] = np.logical_and(
			masks[other_organs_ind]['Mask'], np.logical_not(used_voxels))

	return masks

if __name__=='__main__':
	with open('../data/AA/contours.pickle', 'rb') as infile:
		contours = pickle.load(infile)
	
	directory = '../data/AA'
	try:
		dcm_directory = find_dicom_directory(directory)
		ct_prefix = 'CT'
		dose_prefix = 'RTDOSE'
		
		ct_infos = [pydicom.dcmread(f) for f in find_prefixed_files(dcm_directory, ct_prefix)]
		dose_info = pydicom.dcmread(find_prefixed_file(dcm_directory, dose_prefix))
	except Exception as ex:
		print(type(ex), ex)
		print('Could not load in ct/dose info')
		exit(0)
	
	masks = mask_generation(contours, ct_infos, dose_info)
	with open('../data/AA/masks.pickle', 'wb') as outfile:
		pickle.dump(masks, outfile)
	implay(masks[0]['Mask'].astype(int))
