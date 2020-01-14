import numpy as np
import pandas as pd
import pickle
import os
import argparse

from lymphkill.calc_blood_kill import kill_contributions, calc_alpha_beta_frac

from scipy.io import loadmat

def process_dosestruct(mat):
	dose = mat['doseData']
	dose1frac = mat['doseData_1frac']

	print(dose.dtype.names)
	names = [n[0][0] for n in dose[['Name']][0]]
	pretx = [p[0][0][0] for p in dose[['PreTxLYA']][0]]
	day = [d[0][0][0] for d in dose[['Day']][0]]
	measured = [m[0][0][0] for m in dose[['Measured']][0]]

	blbs = [np.array(b[0][0]) for b in dose[['BinLowerBounds']][0]]
	bcs = [np.array(b[0][0]) for b in dose[['BinCounts']][0]]
	tvs = [np.sum(bc) for bc in bcs]
	
	blbs1 = [np.array(b[0][0]) for b in dose1frac[['BinLowerBounds']][0]]
	bcs1 = [np.array(b[0][0]) for b in dose1frac[['BinCounts']][0]]

	df = pd.DataFrame({
		'Name': names,
		'PreTxLYA': pretx,
		'Day': day,
		'Measured': measured,
		'BinLowerBounds': blbs,
		'BinCounts': bcs,
		'BinLowerBounds_1frac': blbs1,
		'BinCounts_1frac': bcs1,
		'TotalVoxels': tvs,
	})
	return df
	
def process_replenish(mat):
	return {'Variable': mat['RegenRate_var'], 'Constant': mat['RegenRate_const']}

def get_cell_kill_histograms(dose):
	alpha, beta = calc_alpha_beta_frac(0.2188, 1)	#Using fractionated dose with n=1 for LQ
	kills = []
	bnds = []
	for i in range(len(dose)):
		counts = dose['BinCounts'][i]
		edges = dose['BinLowerBounds'][i]
		kcs = kill_contributions(counts, edges, alpha, beta, 1)
		kcs /= dose['TotalVoxels'][i] #Percent kill at each bin
		kcs *= dose['PreTxLYA'][i]
		kills.append(kcs)
		bnds.append(edges)

	bnds = np.array(bnds[0])	#Assume all bounds are the same, which they are
	kills = np.array(kills)

	df =  pd.DataFrame({'Name': dose['Name'], 'PreTxLYA': dose['PreTxLYA']})
	for i, lb in enumerate(bnds):
		df.insert(len(df.columns), 'KC: %g Gy' % lb, kills[:, i])
	df.insert(len(df.columns), 'Total Kill', np.sum(kills, axis=1))
	return df

def get_dose_histograms(dose):
	counts = []
	bnds = []
	for i in range(len(dose)):
		counts.append(dose['BinCounts'][i] / dose['TotalVoxels'][i])
	bnds = np.array(dose['BinLowerBounds'][0])
	counts = np.array(counts)
	df = pd.DataFrame({'Name': dose['Name']})
	for i, lb in enumerate(bnds):
		df.insert(len(df.columns), 'Percent %g Gy' % lb, counts[:, i])
	return df

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('matfile', type=str)
	parser.add_argument('-r', '--replenish', type=str, default='')
	args = parser.parse_args()

	dosemat = loadmat(args.matfile)
	dose = process_dosestruct(dosemat)
	print(dose)
	with open('%s.pickle' % args.matfile[:-4], 'wb') as out:
		pickle.dump(dose, out)

	kill_hists = get_cell_kill_histograms(dose)
	print(kill_hists)
	kill_hists.to_csv('../data/cell_kill_histograms.csv')

	get_dose_histograms(dose).to_csv('../data/dose_histograms.csv')

	if len(args.replenish) > 0:
		replenish = loadmat(args.replenish)
		rep_py = process_replenish(replenish)
		print(rep_py)
		with open('%s.pickle' % args.replenish[:-4], 'wb') as out:
			pickle.dump(rep_py, out)

		

