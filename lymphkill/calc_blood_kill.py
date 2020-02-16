import pydicom
import pickle
import os
import numpy as np
import argparse

replenish_variable = np.array([
	[0.5, 0.0025],
	[1.0, 0.0014],
	[1.5, 0.0004],
	[2.0, -0.0005]])

replenish_const = 0.0006

'''
Calculate alpha and beta for the fractionated LQ model
Parameters:
	k05 - The kill percentage at 0.5 Gy
	fracs - The number of fractions
Returns:
	alpha, beta
'''
def calc_alpha_beta_frac(k05, fracs):
	x1 = 5.
	x2 = 0.5
	k1 = 0.992
	beta = 1 / (fracs * (x1 - x2)) * (np.log(1 - k05) / x2 - np.log( 1 - k1) / x1)
	alpha = -np.log(1 - k1) / (fracs * x1) - beta * x1
	return alpha, beta

'''
Calculate the raw blood cell kill using the fractionated LQ model
Parameters:
	counts, edges - A histogram of the blood matrix
	k05 - The kill percentage at 0.5 Gy
	fracs - The number of fractions
Returns:
	percent - The kill percentage
'''
def calc_kill_frac(counts, edges, k05=0.2044, fracs=5):
	#Raw kill
	alpha, beta = calc_alpha_beta_frac(k05, fracs)
	killed = np.sum(kill_contributions(counts, edges[:-1], alpha, beta, fracs))
	percent = killed / np.sum(counts)

	return percent

'''
Calculate the kill contributions for each dose range
Parameters:
	counts, edges - A histogram of the blood matrix
	alpha, beta - alpha and beta for the fractionated LQ model
	fracs - The number of fractions
Returns:
	The contribution to the total kill for each bin in the blood dose histogram
'''
def kill_contributions(counts, edges, alpha, beta, fracs):
	return  counts * (1. - np.exp(-fracs * edges * (alpha + beta * edges)))

'''
Apply regeneration to the percent kill figure
Parameters:
	percent - The raw kill percentage
	regen_rate - The regeneration rate in percentage points per day
	day - The measurement day (found that regeneration levels off at day 105)
Returns:
	The percent kill including regeneration
'''
def regeneration(percent, regen_rate, day=25):
	day = min(day, 105)
	percent -= (day - 25) * regen_rate if day > 25 else 0
	return percent

'''
Determine the natural cell death at a given day
Parameters:
	pretx - The Pre-Tx LYA in cells/L * 1e9
	day - The measurement day
Returns:
	The LYA drop due to natural cell death at day
	TODO change this so it is nonzero
'''
def natural_cell_death(pretx, day):
	return 0

'''
Determine the regenerated LYA at a given day
The current regeneration function is:
	regenerated = standard + max(0, 1 - (a * pretx + b) * day^2 + (c * pretx + d) * day^2
Parameters:
	pretx - The Pre-Tx LYA in cells/L * 1e9
	day - The measurement day
	standard - The standard cell regeneration rate
	a, b, c, d - The coefficients given above
Returns:
	The LYA increase due to regeneration, due to release from lymphoid organs and
	standard cell regeneration
'''
def regeneration_curve(pretx, day, standard=0, a=1, b=1, c=1, d=1):
	return standard + max(0, 1 - (a * pretx + b) * day * day + day + (c * pretx + d))

'''
Determine the killed LYA at a given day
The current function is:
	killed = pretx * percent * max(0, 1 - (a * pretx + b) * day^2)
	This should be refined
Parameters:
	percent - The kill percent at day = infinity
	pretx - The Pre-Tx LYA in cells/L * 1e9
	day - The measurement day
	a, b - The coefficients given above
'''
def killed_LYA_day(percent, pretx, day, a=1, b=1):
	return pretx * percent * max(0, 1 - (a * pretx + b) * day * day)

'''
Determine total post-treatment LYA level, including regeneration
posttx = pretx * (1 - kill_percent(day)) 
	+ regeneration(pretx, day) 
	- natural_cell_death(pretx, day)

Parameters:
	counts, edges - A histogram of the blood matrix
	pretx - The Pre-Tx LYA in cells/L * 1e9
	day - The measurement day
Return:
	A Post-Tx LYA estimate
'''
def post_treatment_LYA(counts, edges, pretx, day):
	percent = calc_kill_frac(counts, edges)
	return pretx - killed_LYA_day(percent, pretx, day) + \
		regeneration_curve(pretx, day) - \
		natural_cell_death(pretx, day)
	
if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('directory', type=str, help='The patient directory to look in')
	parser.add_argument('-d', '--day', type=int, nargs='+', default=[5, 30, 180], 
		help='Measurement day(s)')
	parser.add_argument('-p', '--pretx', type=float, default=None, help='Pre Treatment LYA')
	args = parser.parse_args()
	
	with open(os.path.join(args.directory, 'blood_hist.pickle'), 'rb') as infile:
		counts, edges = pickle.load(infile)
	percent = calc_kill_frac(counts, edges)

	if args.pretx is None:
		regenRate = replenish_const 
	else:
		ridx = np.argwhere(replenish_variable[:, 0] < args.pretx)[-1]
		regenRate = replenish_variable[ridx, 1]
	
	for day in args.day:
		#regenRate = -0.0090
		percent = regeneration(percent, regenRate, day=day)
		print('Percent Kill after %d Days:\t%g' % (day, percent))
