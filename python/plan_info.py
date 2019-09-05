import pydicom
from file_utils import find_dicom_directory, find_prefixed_file

def get_beam_info(directory):
	dcm_directory = find_dicom_directory(directory)
	rtplan_file = find_prefixed_file(dcm_directory, 'RTPLAN')
	rtplan = pydicom.dcmread(rtplan_file)

	total_mu = 0
	active_beams = 0
	total_time = 0.
	for i, beam in enumerate(rtplan.FractionGroupSequence[0].ReferencedBeamSequence):
		total_mu += beam.BeamMeterset
		doserate = rtplan.BeamSequence[i].ControlPointSequence[0].DoseRateSet
		total_time += beam.BeamMeterset / doserate * 60
		if beam.BeamMeterset > 0:
			active_beams += 1	

	time_per_beam = total_time / active_beams
	
	return total_mu, active_beams, time_per_beam

if __name__=='__main__':
	total_mu, active_beams, time_per_beam = get_beam_info('../data/AA')	
	print('Total MU: %d\nActive Beams: %d\nTime Per Beam: %g' % \
		(total_mu, active_beams, time_per_beam))
