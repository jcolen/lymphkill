from distutils.core import setup

setup(
	name='LymphKill',
	version='0.1dev',
	packages=['lymphkill'],
	description='Lymphocyte kill estimation',
	author='Jonathan Colen',
	author_email='jcolen19@gmail.com',
	url='https://github.com/jcolen/lymphkill',
	install_requires=['numpy',
			  'scipy',
			  'matplotlib',
			  'pandas',
			  'pydicom',],
	scripts=['scripts/calc_blood_dose',
			 'scripts/calc_blood_kill',
			 'scripts/mask_generation',
			 'scripts/plan_info',
			 'scripts/run_pipeline',
			 'scripts/spreadsheet_processing',
			 'scripts/static_organ_info',
			 'scripts/structure_loading']
)
