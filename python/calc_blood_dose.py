import pydicom
import pickle
import os
import numpy as np

from file_utils import find_dicom_directory, find_prefixed_files, load_rtdose_files
from plan_info import get_beam_info

def calc_beam_doses(beam_dose, tvox, layer_size, time_per_beam, gated):
	beam_time = time_per_beam if not gated else time_per_beam * 3.
	frac_dose = np.zeros(beam_dose.shape)
	t = 0
	while t < time_per_beam:
		if not gated or (gated and t % 4. < 1.33):
			frac_dose += beam_dose
		t += tvox
		frac_dose = np.roll(frac_dose, layer_size, axis=1)

'''
@param masks - The masks structure returned by mask_generation
@param time_per_beam - The beam on time for each beam
@param dosegrids - The grids loaded by each RTDOSE file
@param fracs - The number of fractions to compute
@param gated - True if the plan is gated
'''
def calc_blood_dose(masks, time_per_beam, dosegrids, fracs=1, gated=False):
	print('Mask size: %s\tDose size: %s' % (masks[0]['Mask'].shape, dosegrids[0].shape))

	#Find which voxels are part of the body and which are part of contoured organs
	dosed_voxels = np.sum(np.array(dosegrids), axis=0).astype(bool)
	masked_voxels = np.sum(np.array([mask['Mask'] for mask in masks]), axis=0).astype(bool)
	unmasked_voxels = np.logical_and(dosed_voxels, ~masked_voxels)
	
	layer_sizes = np.sum(dosed_voxels, axis=(0, 1))
	nvoxels = np.sum(layer_sizes)
	layer_size = np.floor(nvoxels / np.sum(layer_sizes > 0))
	blood_voxels = 3 * nvoxels	#Assume abdomen is ~1/3 of the body

	'''
	Basic calculation of blood flow rate
	blood density = 5000 cm^3 / blood_voxels
	blood flow rate = 5000 cm^3 / 30 s
	voxel flow rate = blood flow rate * cardiac output / (layer size * blood density)
	This simplifies to VFL = CO * blood_voxels / (30 * layer size)
	In GVs, density is ~n times higher, so VFL is divided by n
	Thus, time_voxel = n * 30 * layer size / (co * blood_voxels)
	'''
	gv_density = 8	#GV density factor
	h2h = 30.0 		#Heart to heart time in seconds
	min_tvox = 0.01 	#Maximum blood velocity of 25 cm/s

	remainingCO = 1.0
	other_ind = -1

	for j, mask in enumerate(masks):
		if mask['CardiacOutput'] == -1:
			other_ind = j
		else:
			mask['TimeVoxel'] = h2h * mask['LayerSize'] / (mask['CardiacOutput'] * blood_voxels)
			if mask['GV']:
				mask['TimeVoxel'] *= gv_density
			else:
				remainingCO -= mask['CardiacOutput']
			mask['TimeVoxel'] = max(mask['TimeVoxel'], min_tvox)
	
	nunmasked = np.sum(unmasked_voxels)
	remaining_layers = np.sum(np.sum(unmasked_voxels, axis=(0, 1)) > 0)
	remaining_layer_size = nunmasked / remaining_layers

	if other_ind >= 0:
		nother = np.sum(masks[other_ind]['Mask'])
		masks[other_ind]['CardiacOutput'] = remainingCO * nother / (nother + nunmasked)
		masks[other_ind]['TimeVoxel'] = h2h * masks[other_ind]['LayerSize'] /\
			(masks[other_ind]['CardiacOutput'] * blood_voxels)
		remainingCO -= masks[other_ind]['CardiacOutput']
		masks[other_ind]['TimeVoxel'] = max(min_tvox, masks[other_ind]['TimeVoxel'])
	
	remaining_tvox = h2h * remaining_layer_size / (remainingCO * blood_voxels)
	remaining_tvox = max(min_tvox, remaining_tvox)

	dosegrids = np.array(dosegrids)

	
	padding_factor = 2 * gv_density + 2
	for i, mask in enumerate(masks):
		mask_nvoxels = np.sum(mask['Mask'])
		mask['BeamDose'] = np.zeros([dosegrids.shape[0], blood_voxels])
		if mask['GV']:
			sub_size = int(np.floor(blood_voxels / padding_factor))
			sub_mask = np.zeros([dosegrids.shape[0], sub_size])
			sub_mask[:, :mask_nvoxels] = dosegrids[:, mask['Mask']].reshape([dosegrids.shape[0], -1])
			mask['BeamDose'][:, :gv_density*sub_size] = np.repeat(sub_mask, 8, axis=1)
		else:
			mask['BeamDose'][:, :mask_nvoxels] = dosegrids[:, mask['Mask']].reshape(\
				[dosegrids.shape[0], -1])

		if not mask['Stationary']:
			mask['BeamDose'] *= mask['TimeVoxel']
		print(mask['BeamDose'].shape)

	remaining_dose = dosegrids[:, unmasked_voxels] * remaining_tvox
	remaining_beam_dose = np.zeros([dosegrids.shape[0], blood_voxels])
	remaining_beam_dose[:, :nunmasked] = remaining_dose.reshape([dosegrids.shape[0], -1])

	print('Precalculating total beam doses')

	for mask in masks:
		print('Calculating mask %s' % mask['Name'])
		if mask['Stationary']:
			mask['FracDose'] = mask['BeamDose'] * time_per_beam
		else:
			mask['FracDose'] = calc_beam_doses(
			mask['BeamDose'], mask['TimeVoxel'], int(mask['LayerSize']), time_per_beam, gated)
	remaining_frac_dose = calc_beam_doses(
		remaining_dose, remaining_tvox, int(remaining_layer_size), time_per_beam, gated)	

	blood = np.zeros(blood_voxels)
	for day in range(len(fracs)):
		print('Beginning dose for Day %d' % day)
		for i in range(frac_dose.shape[0]):
			print('\tBeam %d of %d' % (i, len(dosegrids)))
			blood += remaining_frac_dose[i, :]
			blood = np.random.permutation(blood)
			for mask in masks:
				blood += mask['FracDose'][i, :]
				blood = np.random.permuation(blood)
	
	return blood


if __name__=='__main__':
	directory = '../data/AA'
	dcm_directory = find_dicom_directory(directory)
	rtdose_files = find_prefixed_files(dcm_directory, 'RTDOSE')
	dosegrids = load_rtdose_files(rtdose_files)
	print(len(dosegrids))

	with open(os.path.join(directory, 'masks.pickle'), 'rb') as infile:
		masks = pickle.load(infile)

	total_mu, active_beams, time_per_beam = get_beam_info(directory)

	blood_voxels = calc_blood_dose(masks, time_per_beam, dosegrids)
	print('Done calculating blood')

	import matplotlib.pyplot as plt
	plt.hist(blood_voxels)
	plt.show()
