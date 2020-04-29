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

import sys
import os
from click import echo, secho
from SpoON.util import parse_config
from MicrobeAndMe.data_structure import Microbe, ASVsMerger, FilesCollector

# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']


def all_pipeline(kargs):
    microbe = Microbe(kargs, mode='all', analysis_type='MicrobeAndMe')
    microbe.run_prema()
    if (kargs['copy_rawdata'] is False) and (kargs['no_prema'] is False):
        exit()
    microbe.make_analysis_number_dir()
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_dada2_using_r()
    microbe.check_jobs()
    microbe.assign_taxonomy()
    microbe.make_biom_by_db()
    microbe.run_alpha_diversity()
    microbe.summrize_taxa()


def dada2_using_r(kargs):
    # microbe = Microbe(kargs, mode='r_dada2', analysis_type='MicrobeAndMe')
    microbe = Microbe(kargs, mode='r_dada2', analysis_type='NGS')
    microbe.make_analysis_number_dir()
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_dada2_using_r()
    microbe.check_jobs()


def merge_asvs_r(kargs):
    asv_merge = ASVsMerger(kargs)
    asv_merge.read_order_number_file()
    asv_merge.glob_samples()
    asv_merge.check_rds_files()
    all_rds_dir_path = asv_merge.copy_rds_files('rds_files')
    asv_merge.run_merge_asvs(all_rds_dir_path)


def collect_files(kargs):
    collection = FilesCollector(kargs)
    collection.read_order_number_file()
    collection.glob_samples()
    collection.check_files()
    collection.copy_files('collected_files')


def biom(kargs):
    microbe = Microbe(kargs, mode='biom', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.make_biom_by_db()


def taxonomy(kargs):
    microbe = Microbe(kargs, mode='taxonomy', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.assign_taxonomy()


def alpha_diversity(kargs):
    microbe = Microbe(kargs, mode='alpha_diversity', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_alpha_diversity()


def summarize_taxa(kargs):
    microbe = Microbe(kargs, mode='summarize_taxa', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.summrize_taxa()


def score():
    from MicrobeAndMe.score import AndMe
    path = '/garnet/Analysis/BI/AmpliconMetaGenome_MAM/HN00122513/Analysis_1/35.8/Summarized_Taxa'
    and_me = AndMe(
        {
            'bact_file': '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/items_and_bacteria.yaml',
            'dist_db_file': '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/MicrobeAndMe_Ref_V2.1.yaml',
            'phylum_file': os.path.join(path, 'ASVs.NCBI_16S_L2.txt'),
            'family_file': os.path.join(path, 'ASVs.NCBI_16S_L5.txt'),
            'genus_file': os.path.join(path, 'ASVs.NCBI_16S_L6.txt'),
            'species_file': os.path.join(path, 'ASVs.NCBI_16S_L7.txt'),
            'sample_name': '35.8'
        })
    and_me.read_bacteria_file()
    and_me.read_dist_db_file()
    and_me.read_genus_file()
    and_me.read_species_file()
    and_me.compute_intestinal_health()
    and_me.compute_wellness()
    # DiversityIndex('/garnet/Analysis/BI/AmpliconMetaGenome_MAM/HN00122513/Analysis_1/35.8/Alpha_Diversity/adiv.txt')


def insert():
    pass

