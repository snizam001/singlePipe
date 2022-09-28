from setuptools import find_packages,setup

requirements = open("requirements.txt").read().strip().split("\n")

exec(open('singlePipe/__init__.py').read())

setup(
	name='singlePipe',
	packages=find_packages(),
	entry_points={
		"console_scripts": ['singlePipe = singlePipe.run:main']
		},
	version=__version__,
	description="pipeline for the analysis of single cell datasets (Granja et al. (Nat. Biotech 2019)",
	author="Sheikh Nizamuddin and ..",
	license="MIT",
	include_package_data=True,
	package_data={"greenPipe": ["data/*"]},
#	scripts=['rscripts/*'],
#	package_rscript={"greenCUTRUN": ["rscripts/*"]},
#	package_dir={"":'greenCUTRUN'},
	install_requires=requirements,
#	url='https://github.com/snizam001/singlePipe'
) 
