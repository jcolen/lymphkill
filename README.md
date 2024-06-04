# Blood Dose and Lymphocyte Kill Simulation

This is a repository containing code originally written to compute the absorbed dose to circulating lymphocytes during ratdiation therapy treatments and predict radiation-induced immune suppression (RIIS). ***This is a preliminary version of the code!*** The final version, used in the manuscript [Predicting radiation-induced immune suppression in lung cancer patients treated with stereotactic body radiation therapy](https://aapm.onlinelibrary.wiley.com/doi/full/10.1002/mp.17181) is maintained by Cam Nguyen.

The repo contains three folders. The folder `lymphkill` contains Python code for processing and analyzing patient treatment plan information. The folder `matlabcode` contains MATLAB scripts which were used for the original implementation of this project and are preserved here for the sake of completeness. 

The Python code requires the following libraries to be installed:
- numpy
- scipy
- matplotlib
- pickle
- functools
- pandas
- pydicom

They will be installed automatically with running:

`python setup.py install`

The Python source files are located in the `lymphkill` subdirectory. The scripts, which can be run standalone from the command line, are located in the `scripts` subdirectory, and are installed with `setup.py`. They are:
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
