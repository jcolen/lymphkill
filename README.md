# Blood Dose and Lymphocyte Kill Simulation

This is a repository containing the code referenced in \<CITE PAPER UPON SUBMISSION\>. The folder `matlabcode` the matlab scripts which were used for the original run of code. The folder `lymphkill` contains those same scripts but rewritten in Python with the goal of readability and better code structure. 

The Python code requires the following libraries to be installed:
- numpy
- matplotlib
- pickle
- functools
- pandas
- pydicom

They will be installed automatically with running:

`python setup.py install`

Because the goal is to replace all Matlab scripts with Python, we will not outline the structure of the Matlab code. The python source files are loacted in the `lymphkill` subdirectory. The scripts, which can be run standalone from the command line, are located in the `scripts` subdirectory, and are installed with `setup.py`. They are:
- `calc_blood_dose` - Calculation of the dose absorbed by blood during a treatment
  - Output stored in `blood_dose.pickle`, `blood_hist.pickle`
- `calc_blood_kill` - Calculation of the blood cell kill, given the absorbed dose
- `static_organ_info` - Calculation of static organ doses
- `mask_generation` - Generation of mask structures for use in blood dose calculation
  - Output stored in `masks.pickle`
- `plan_info` - Parse RTPLAN files for information such as number of beams and beam on time
- `run_pipeline` - Given a zipped file containing patient CT, dose, plan, and structure set info (in dicom format), run the entire pipeline of scripts
- `spreadsheet_processing` - Data processing for summary spreadsheets
- `structure_loading` - Loading of physician contoured organs into boolean masks
  - Output stored in `contours.pickle`

Each script takes a single command line argument. `spreadsheet_processing.py` accepts a path to a csv file which has a summary of patient and prediction information. 

The remaining scripts accept a path to a patient directory. All generated pickle files will be stored in this directory. It is assumed that this directory has a subdirectory which contains all of the dicom files needed (CT, RTDOSE, RTPLAN, RTSTRUCT). To run the full pipeline, simply place the zip file containing all of this information in that directory and execute:

`python run_pipeline.py <path_to_directory>`

This will unzip the patient information and perform all of the processing. Structures generated in intermediate steps will be stored as pickle files. Any other scripts can be run as, for example:

`python structure_loading.py <path_to_directory>`
