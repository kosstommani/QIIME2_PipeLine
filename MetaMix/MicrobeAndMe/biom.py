# ----------------------------------------------------------------------------------------------------------------------
#                                      888b     d888          888             888b     d888 d8b
#                                      8888b   d8888          888             8888b   d8888 Y8P
#                                      88888b.d88888          888             88888b.d88888
# 88888b.d88b.   8888b.  88888b.d88b.  888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b     "88b 888 "888 "88b 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 .d888888 888  888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 888  888 888  888  888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888 "Y888888 888  888  888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.0'

from click import echo, secho
import os
from SpoON.util import check_file_type, run_cmd, check_run_cmd
from CoFI.biom import get_biom_summarize_cmd


class Biom:
	def __init__(self, kargs: dict, p_verbose: bool = True, p_exit: bool = True):
		"""

		:param kargs:
				r_dada2_path:
				biom_path:
				taxonomy_path:
				db_name:
				# metadata:
		"""
		self.verbose = p_verbose
		self.exit = p_exit
		self.r_dada2_path = kargs['r_dada2_path']
		self.biom_path = kargs['biom_path']
		self.taxonomy_path = kargs['taxonomy_path']
		self.db_name = kargs['db_name']
		self.metadata = None

	@property
	def asvs_table_file(self):
		return os.path.join(self.r_dada2_path, self.asvs_table_name)

	@property
	def tax_assignment_file(self):
		return os.path.join(self.taxonomy_path, f'blast_{self.db_name}', self.taxonomy_name)

	@property
	def base_biom_file(self):
		return os.path.join(self.biom_path, self.biom_name)

	@property
	def taxonomy_biom_file(self):
		splitext = os.path.splitext(self.biom_name)
		return os.path.join(self.biom_path, f'{splitext[0]}.blast_{self.db_name}{splitext[1]}')

	@property
	def biom_table_file(self):
		return os.path.join(self.biom.path, self.biom_table_name)

	@property
	def taxonomy_biom_table_file(self):
		splitext = os.path.splitext(self.biom_table_name)
		return os.path.join(self.biom_path, f'{splitext[0]}.bast_{self.db_name}{splitext[1]}')

	@property
	def recipe_file(self):
		if self.db_name is None:
			recipe_name = 'biom_recipe'
		else:
			recipe_name = f'biom.{self.db_name}.recipe'
		return os.path.join(self.biom_path, recipe_name)

	def check_biom_dir(self, p_verbose=False):
		if check_file_type(self.biom_path, 'exists'):
			pass
		else:
			os.mkdir(self.biom_path)
			if p_verbose is True:
				secho('>>> BIOM 디렉터리 생성', fg='cyan')
				echo(self.biom_path)

	def run_cmd(self, cmd, true_meg, false_meg):
		run = run_cmd(cmd)
		check_run_cmd(
			{
				'run': run,
				'true_meg': true_meg,
				'false_meg': false_meg,
			}, self.exit
		)
		return run

	def write_recipe(self, p_cmd: str):
		with open(self.recipe_file, 'a') as o_recipe:
			o_recipe.write(p_cmd)
			o_recipe.write('\n')
			o_recipe.write('\n')

	def convert_biom(self, p_no_run: bool = False) -> object:
		"""
		
		:param p_no_run: recipe 파일만 생성시 
		:return: 
		"""
		if self.verbose is True:
			true_meg = 'BIOM 생성'
			false_meg = 'biom convert'
		else:
			true_meg = None
			false_meg = None

		biom_convert_cmd = \
			'biom convert ' \
			f'-i {self.asvs_table_file} ' \
			f'-o {self.base_biom_file} ' \
			'--to-json'
		self.write_recipe(biom_convert_cmd)
		if p_no_run is False:
			run = self.run_cmd(biom_convert_cmd, true_meg, false_meg)
			return run

	def add_metadata(self, p_no_run: bool = False, p_add_sample_metadata: bool = False) -> object:
		"""
		
		:param p_no_run: recipe 파일만 생성시
		:param p_add_sample_metadata: sample metadata 추가시. only NGS
		:return: 
		"""
		if self.verbose is True:
			true_meg = 'Add Taxonomy 완료'
			false_meg = 'biom add-metadata'
		else:
			true_meg = None
			false_meg = None

		biom_add_metadata_cmd = \
			'biom add-metadata ' \
			f'-i {self.base_biom_file} ' \
			f'-o {self.taxonomy_biom_file} ' \
			f'--observation-metadata-fp {self.tax_assignment_file} ' \
			f'{f"--sample-metadata-fp {self.metadata} " if p_add_sample_metadata else ""}' \
			'--observation-header ASVsID,taxonomy,confidence,count ' \
			'--sc-separated taxonomy ' \
			'--output-as-json'
		self.write_recipe(biom_add_metadata_cmd)
		if p_no_run is False:
			run = self.run_cmd(biom_add_metadata_cmd, true_meg, false_meg)
			return run

	def convert_base_table(self, p_no_run: bool = False) -> object:
		"""
		
		:param p_no_run: recipe 파일만 생성시 
		:return: 
		"""
		if self.verbose is True:
			true_meg = 'Convert Table 완료'
			false_meg = 'biom convert'
		else:
			true_meg = None
			false_meg = None
		biom_table_cmd = \
			'biom convert ' \
			f'-i {self.taxonomy_biom_file} ' \
			f'-o {self.biom_table_file} ' \
			'--to-tsv'
		self.write_recipe(biom_table_cmd)
		if p_no_run is False:
			run = self.run_cmd(biom_table_cmd, true_meg, false_meg)
			return run

	def convert_taxonomy_table(self, p_no_run: bool = False) -> object:
		"""
		
		:param p_no_run: recipe 파일만 생성시
		:return:
		"""
		if self.verbose is True:
			true_meg = 'Convert Table 완료'
			false_meg = 'biom convert'
		else:
			true_meg = None
			false_meg = None
		biom_table_cmd = \
			'biom convert ' \
			f'-i  {self.taxonomy_biom_file} ' \
			f'-o {self.taxonomy_biom_table_file} ' \
			'--to-tsv ' \
			'--header-key taxonomy'
		self.write_recipe(biom_table_cmd)
		if p_no_run is False:
			run = self.run_cmd(biom_table_cmd, true_meg, false_meg)
			return run


class BiomForMicrobeAndMe(Biom):
	def __init__(self, kargs: dict, p_verbose: bool = True, p_exit: bool = True):
		super().__init__(kargs, p_verbose, p_exit)
		self.asvs_table_name = 'ASVs.tsv'
		self.biom_name = 'ASVs.biom'
		self.biom_table_name = 'ASVs_table.tsv'
		self.taxonomy_name = 'otus_rep_tax_assignments.txt'


class BiomForNGS(Biom):
	def __init__(self, kargs: dict, p_verbose: bool = True, p_exit: bool = True):
		super().__init__(kargs, p_verbose, p_exit)
		self.metadata = kargs['metadata']
		self.asvs_table_name = 'all_ASVs.tsv'
		self.biom_name = 'all_ASVs.biom'
		self.biom_table_name = 'all_ASVs_table.tsv'
		self.taxonomy_name = 'otus_rep_tax_assignments.txt'

	def make_summary_file(self, p_no_run: bool = False) -> object:
		if self.verbose is True:
			true_meg = 'BIOM Summary 생성'
			false_meg = 'biom summarize-table'
		else:
			true_meg = None
			false_meg = None
		summary_file = f'{os.path.splitext(self.base_biom_file)[0]}_summary.txt'
		summary_cmd = get_biom_summarize_cmd(self.base_biom_file, summary_file)
		self.write_recipe(summary_cmd)
		if p_no_run is False:
			run = self.run_cmd(summary_cmd, true_meg, false_meg)
			return run
