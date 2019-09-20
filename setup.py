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
					  'matplotlib',
					  'pickle',
					  'functools',
					  'pandas',
					  'pydicom']
)
