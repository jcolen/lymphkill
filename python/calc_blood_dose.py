import pydicom
import pickle
import os
import numpy as np
import functools
import multiprocessing as mp
from time import time

from file_utils import find_dicom_directory, find_prefixed_files, load_rtdose_files
from plan_info import get_beam_info

def calc_beam_doses(mask, time_per_beam, gated):
	print('Calculating mask %s:\tTV: %g' % (mask['Name'], mask['TimeVoxel']))
	tm = time()
	if mask['Stationary']:
		return mask['BeamDose'] * time_per_beam

	beam_time = time_per_beam if not gated else time_per_beam * 3.
	beam_dose = mask['BeamDose']
	tvox = mask['TimeVoxel']
	layer_size = int(mask['LayerSize'])
	frac_dose = np.zeros(beam_dose.shape)
	t = 0
	while t < time_per_beam:
		if not gated or (gated and t % 4. < 1.33):
			frac_dose += beam_dose
		t += tvox
		frac_dose = np.roll(frac_dose, layer_size, axis=1)

	mask['FracDose'] = frac_dose
	print('\tTook %g s' % (time() - t))

'''
@param masks - The masks structure returned by mask_generation
@param time_per_beam - The beam on time for each beam
@param dosegrids - The grids loaded by each RTDOSE file
@param fracs - The number of fractions to compute
@param gated - True if the plan is gated
'''
def calc_blood_dose(masks, time_per_beam, dosegrids, 
					fracs=1, 
					gated=False, 
					use_multiprocessing=False):
	t = time()
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
		elif mask['CardiacOutput'] == 0:
			mask['TimeVoxel'] = 1.
		else:
			mask['TimeVoxel'] = h2h * mask['LayerSize'] / (mask['CardiacOutput'] * blood_voxels)
			if mask['GV']:
				mask['TimeVoxel'] *= gv_density
			else:
				remainingCO -= mask['CardiacOutput']
			mask['TimeVoxel'] = max(mask['TimeVoxel'], min_tvox)
	
	nunmasked = np.sum(unmasked_voxels)
	remaining_layers = np.sum(np.sum(unmasked_voxels, axis=(0, 1)) > 0)
	masks.append({'Name': 'Remaining'})
	masks[-1]['LayerSize'] = nunmasked / remaining_layers

	if other_ind >= 0:
		nother = np.sum(masks[other_ind]['Mask'])
		masks[other_ind]['CardiacOutput'] = remainingCO * nother / (nother + nunmasked)
		masks[other_ind]['TimeVoxel'] = h2h * masks[other_ind]['LayerSize'] /\
			(masks[other_ind]['CardiacOutput'] * blood_voxels)
		remainingCO -= masks[other_ind]['CardiacOutput']
		masks[other_ind]['TimeVoxel'] = max(min_tvox, masks[other_ind]['TimeVoxel'])
	
	masks[-1]['CardiacOutput'] = remainingCO
	remaining_tvox = h2h * masks[-1]['LayerSize'] / (remainingCO * blood_voxels)
	masks[-1]['TimeVoxel'] = max(min_tvox, remaining_tvox)
	masks[-1]['Mask'] = unmasked_voxels
	masks[-1]['Stationary'] = False
	masks[-1]['GV'] = False

	dosegrids = np.array(dosegrids)

	padding_factor = 2 * gv_density + 2
	for i, mask in enumerate(masks):
		print('Mask %s:\n\tCO: %g\n\tTV: %g' % (mask['Name'], mask['CardiacOutput'], mask['TimeVoxel']))
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

	print('Time to set up doses: %g' % (time() - t))
	print('Precalculating total beam doses')

	if use_multiprocessing:
		print('Using %d workers' % os.cpu_count())
		pool = mp.Pool(os.cpu_count())
		pool.map(functools.partial(calc_beam_doses, 
								   time_per_beam=time_per_beam, 
								   gated=gated), 
				 masks)
		pool.close()
	else:
		for mask in masks:
			calc_beam_doses(mask, time_per_beam, gated)

	blood = np.zeros(blood_voxels)
	for day in range(len(fracs)):
		print('Beginning dose for Day %d' % day)
		for i in range(frac_dose.shape[0]):
			print('\tBeam %d of %d' % (i, len(dosegrids)))
			for mask in masks:
				blood += mask['FracDose'][i, :]
				blood = np.random.permuation(blood)
	
	return blood


if __name__=='__main__':
	directory = '../data/AA'
	dcm_directory = find_dicom_directory(directory)
	rtdose_files = find_prefixed_files(dcm_directory, 'RTDOSE')
	dosegrids = load_rtdose_files(rtdose_files)

	with open(os.path.join(directory, 'masks.pickle'), 'rb') as infile:
		masks = pickle.load(infile)

	total_mu, active_beams, time_per_beam = get_beam_info(directory)

	blood_voxels = calc_blood_dose(masks, time_per_beam, dosegrids, use_multiprocessing=True)
	print('Done calculating blood')

	import matplotlib.pyplot as plt
	plt.hist(blood_voxels)
	plt.show()

	with open(os.path.join(directory, 'blood.pickle'), 'wb') as outfile:
		pickle.dump(outfile, blood_voxels)
