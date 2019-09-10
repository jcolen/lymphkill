import pydicom
import pickle
import os
import numpy as np

def calc_alpha_beta_frac(k05, fracs):
	x1 = 5.
	x2 = 0.5
	k1 = 0.992
	beta = fracs / (x2 - x1) * (np.log(1 - k1) / x1 - np.log(1 - k05) / x2)
	alpha = -np.log(1 - k1) / x1 - beta * x1 / fracs
	return alpha, beta

def calc_kill_frac(counts, edges, k05=0.2044, fracs=5):
	#Raw kill
	alpha, beta = calc_alpha_beta_frac(k05, fracs)
	killed = np.sum(kill_contributions(counts, edges[:-1], alpha, beta, fracs))
	percent = killed / np.sum(counts)

	return percent

def kill_contributions(counts, edges, alpha, beta, fracs):
	return  counts * (1. - np.exp(-fracs * edges * (alpha + beta * edges)))

def regeneration(percent, regen_rate, day=25):
	day = min(day, 105)
	percent -= (day - 25) * regen_rate if day > 25 else 0
	return percent

if __name__=='__main__':
	directory = '../data/AA'
	with open(os.path.join(directory, 'blood_hist.pickle'), 'rb') as infile:
		counts, edges = pickle.load(infile)
	percent = calc_kill_frac(counts, edges)
	percent = regeneration(percent, 0.001, day=25)
	print('Total Percent Kill:\t%g' % percent)
