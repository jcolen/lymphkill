import numpy as np
import pydicom
import os
import pickle

from structure_loading import find_prefixed_file
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

'''
Rescale a ct image into the size of a dosemap image
'''
def voxelvolume_to_dosemap(volume, ct_info, dose_info):
	#Get grid dimensions
	dim_dos = np.array([dose_info.Rows, dose_info.Columns, dose_info.NumberOfFrames])
	dim_vol = volume.shape
	
	#Get voxel dimensions
	dim_voxd = np.array([dose_info.PixelSpacing[0], dose_info.PixelSpacing[1], dose_info.SliceThickness])
	dim_voxv = np.array([ct_info.PixelSpacing[0], ct_info.PixelSpacing[1], ct_info.SliceThickness])

	#Get image corners
	corner_d = np.array([dose_info.ImagePositionPatient])
	corner_v = np.array([ct_info.ImagePositionPatient])

	dosemap = np.zeros(dim_dos)

	x, y, z = np.meshgrid(
		np.arange(dim_dos[0]),
		np.arange(dim_dos[1]),
		np.arange(dim_dos[2])
	)

	a = np.ravel_multi_index((x.flatten(), y.flatten(), z.flatten()), dim_dos)

	xyz = np.transpose(np.array([y, x, z]), (1, 2, 3, 0))

	pos = np.empty(xyz.shape)
	pos[:, :, :] = (corner_d - corner_v + xyz * dim_voxd) / dim_voxv

	pos = pos.reshape([-1, 3])
	xyz = xyz.reshape([-1, 3])

	valid_voxel = np.logical_and(pos[:, 0] >= 0, pos[:, 1] >= 0)
	valid_voxel = np.logical_and(valid_voxel, pos[:, 2] >= 0)
	valid_voxel = np.logical_and(valid_voxel, pos[:, 0] < dim_vol[1])
	valid_voxel = np.logical_and(valid_voxel, pos[:, 1] < dim_vol[0])
	valid_voxel = np.logical_and(valid_voxel, pos[:, 2] < dim_vol[2])

	pos = pos.astype(int)
	pos = pos[valid_voxel]
	xyz = xyz[valid_voxel]

	xyz = np.ravel_multi_index(xyz.transpose(), [dim_dos[1], dim_dos[0], dim_dos[2]])
	pos = np.ravel_multi_index(pos.transpose(), dim_vol)

	dosemap = dosemap.flatten()
	dosemap[xyz] = volume.flatten()[pos]
	dosemap = dosemap.reshape(dim_dos)

	return dosemap

def layer_size(mask):
	num = np.sum(mask, axis=2)
	num = num[num > 0]
	return np.sum(num) / len(num)

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
	ct_info,
	dose_info,
	mask_dicts=basic_mask_dicts):

	masks = []

	other_organs_ind = -1
	used_voxels = None
	
	for i, dct in enumerate(mask_dicts):
		contour_idx = find_matching_contour_idx(contours, dct['NameStrings'])
		mdict = {}
		mdict['Name'] = contours['ROIName'][contour_idx]
		mdict['GV'] = dct['GV']
		mdict['Stationary'] = dct['Stationary']
		mdict['CardiacOutput'] = dct['CardiacOutput']

		print('Creating mask for organ %s' % mdict['Name'])

		mdict['Mask'] = voxelvolume_to_dosemap(
			contours['Segmentation'][contour_idx],
			ct_info=ct_info,
			dose_info=dose_info).astype(bool)
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
		dcm_directory = None
		ct_prefix = 'CT'
		dose_prefix = 'RTDOSE'

		for i in os.listdir(directory):
			if os.path.isdir(os.path.join(directory, i)):
				dcm_directory = os.path.join(directory, i)
		
		ct_info = pydicom.dcmread(find_prefixed_file(dcm_directory, ct_prefix))
		dose_info = pydicom.dcmread(find_prefixed_file(dcm_directory, dose_prefix))
	except Exception as ex:
		print(type(ex), ex)
		print('Could not load in ct/dose info')
		exit(0)
	
	masks = mask_generation(contours, ct_info, dose_info)
	with open('../data/AA/masks.pickle', 'wb') as outfile:
		pickle.dump(masks, outfile)
