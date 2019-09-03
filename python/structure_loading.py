import pydicom
import os
import re
import pickle

import numpy as np
from matplotlib.path import Path

def find_prefixed_file(directory, prefix):
	for i in os.listdir(directory):
		if prefix in i:
			return os.path.join(directory, i)
	return None

def load_dicom_imageheaders(directory, siUID):
	#Loop over dicom files in directory
	imgheaders = []
	for i in os.listdir(directory):
		try:
			info = pydicom.dcmread(os.path.join(directory, i))
			if info.StudyInstanceUID == siUID and re.search('^RT', info.Modality) is None:
				imgheaders.append(info)
		except:
			continue
	imgheaders.sort(key=lambda x: int(x.InstanceNumber))
	return imgheaders

def build_affine_transformation(headers):
	N = len(headers)
	dr, dc = float(headers[0].PixelSpacing[0]), float(headers[0].PixelSpacing[1])
	F = np.zeros([3, 2])
	F[:, 0] = headers[0].ImageOrientationPatient[:3]
	F[:, 1] = headers[0].ImageOrientationPatient[3:6]

	T1 = np.array(headers[0].ImagePositionPatient, dtype='float')
	TN = np.array(headers[-1].ImagePositionPatient, dtype='float')
	k = (T1 - TN) / (1 - N)

	A = np.array([
		[F[0, 0] * dr, F[0, 1] * dc, k[0], T1[0]],
		[F[1, 0] * dr, F[1, 1] * dc, k[1], T1[1]],
		[F[2, 0] * dr, F[2, 1] * dc, k[2], T1[2]],
		[0, 0, 0, 1]])
	
	return A


def read_structures(rtssheader, imageheaders):
	xfm = build_affine_transformation(imageheaders)
	dim_min = np.array([0, 0, 0, 1])
	dim_max = np.array([int(imageheaders[0].Columns)-1, int(imageheaders[1].Rows)-1, len(imageheaders)-1, 1])

	template = np.zeros([imageheaders[0].Columns, imageheaders[1].Rows, len(imageheaders)])

	roi_contour_sequence = rtssheader.ROIContourSequence
	roi_structureset_sequence = rtssheader.StructureSetROISequence

	nrois = len(roi_contour_sequence)
	contours = {
		'ROIName': [None] * nrois, 
		'Points': [None] * nrois, 
		'VoxPoints': [None] * nrois, 
		'Segmentation': [None] * nrois
	}

	for i in range(nrois):
		roi_contour = roi_contour_sequence[i]
		roi_number = roi_contour.ReferencedROINumber
		roi_structureset = next((x for x in roi_structureset_sequence if x.ROINumber == roi_number), None)
		contours['ROIName'][i] = roi_structureset.ROIName
		contours['Segmentation'][i] = template

		print(contours['ROIName'][i])

		try:
			contour_sequence = roi_contour.ContourSequence
			segments = [None] * len(contour_sequence)
			#Loop through segments
			for j in range(len(contour_sequence)):
				contour = contour_sequence[j]
				if contour.ContourGeometricType == 'CLOSED_PLANAR':
					#Read segment points
					segments[j] = np.reshape(contour.ContourData,
						[contour.NumberOfContourPoints, 3])

					#Make lattice
					d4seg = np.ones([segments[j].shape[0], 4])
					d4seg[:, :3] = segments[j]
					d4seg = d4seg.transpose()
					d4start = np.append(segments[j][0], 1)
					points = np.linalg.solve(xfm, d4seg)
					start = np.linalg.solve(xfm, d4start)

					minvox = np.maximum(np.floor(np.min(points, axis=1)), dim_min)
					maxvox = np.minimum(np.ceil(np.max(points, axis=1)), dim_max)

					minvox[2] = round(start[2])
					maxvox[2] = round(start[2])
				
					minvox = minvox.astype(int)
					maxvox = maxvox.astype(int)

					x, y, z = np.meshgrid(
						np.arange(minvox[0], maxvox[0]+1),
						np.arange(minvox[1], maxvox[1]+1),
						np.arange(minvox[2], maxvox[2]+1))

					points = np.matmul(
						xfm,
						np.array([x.flatten(), y.flatten(), z.flatten(), np.ones(x.size)]))
						
					#Make binary image
					segpath = Path(segments[j][:, :2])
					inpoly = segpath.contains_points(points[:2].transpose())
					inpoly = inpoly.reshape(x.shape)
					inpoly = np.transpose(inpoly, (1, 0, 2))
					contours['Segmentation'][i][minvox[0]:maxvox[0]+1, minvox[1]:maxvox[1]+1, minvox[2]:maxvox[2]+1] = inpoly 

			contours['Points'][i] = np.vstack(segments)
			d4p = np.ones([contours['Points'][i].shape[0], 4])
			d4p[:, :3] = contours['Points'][i]
			d4p = d4p.transpose()
			contours['VoxPoints'][i] = np.linalg.solve(xfm, d4p)
			contours['VoxPoints'][i] = contours['VoxPoints'][i][:3]
		except Exception as ex:
			print(type(ex), ex)
			print('Loading contour %d: %s failed' % (i, contours['ROIName'][i]))

	return contours

def load_structures(directory, struct_prefix='RTSTRUCT', dose_prefix='RTDOSE', ct_prefix='CT'):
	rtstruct_file = find_prefixed_file(directory, struct_prefix)
	ct_file = find_prefixed_file(directory, ct_prefix)
	rtdose_file = find_prefixed_file(directory, dose_prefix)

	rtssheader = pydicom.dcmread(rtstruct_file)
	ct_info = pydicom.dcmread(ct_file)
	rtd_info = pydicom.dcmread(rtdose_file)

	print(ct_file)
	print(rtdose_file)

	imageheaders = load_dicom_imageheaders(directory, rtssheader.StudyInstanceUID)

	print('Found Organs:')
	for struct in rtssheader.StructureSetROISequence:
		print(struct.ROIName)

	print('Loading in contours')
	contours = read_structures(rtssheader, imageheaders)

	return contours 


if __name__=='__main__':
	contours = load_structures('../data/AA/clinical')
	for i in contours['ROIName']:
		print(i)
	outfile = open('../data/AA/contours.pickle', 'wb')
	pickle.dump(contours, outfile)
	outfile.close()
