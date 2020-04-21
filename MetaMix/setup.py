#!/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
# ----------------------------------------------------------------------------------------------------------------------
# 888b     d888          888             888b     d888 d8b
# 8888b   d8888          888             8888b   d8888 Y8P
# 88888b.d88888          888             88888b.d88888
# 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.0'


from setuptools import setup, find_packages

setup(
	name='MetaMix',
	version='1.0',
	author='JungWon Park(KOSST)',
	author_email='kosstommani@macrogen.com',
	packages=find_packages(include=['bin']),
	include_package_data=True,
	zip_safe=False,
	install_requires=[
		'Click',
		'humanize',
		'bs4',
		'jinja2',
	],
	entry_points={
		'console_scripts': [
			'MetaMix = bin.MetaMix:metamix',
			'PreMA = bin.PreMA:main',
			'CoFI = bin.CoFI:main',
			'theCups = bin.theCups:main',
		],
	},
)
