import os
import pickle
import argparse

import numpy as np
import pandas as pd

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--dosestruct', type=str, default='../data/dose_struct_fix.pickle')
	args = parser.parse_args()

	with open(args.dosestruct, 'rb') as infile:
		ds = pickle.load(infile)

	print(ds.dtypes)
