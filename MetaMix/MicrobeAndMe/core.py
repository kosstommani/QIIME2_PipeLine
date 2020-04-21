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
    microbe = Microbe(kargs, mode='all')
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
    microbe = Microbe(kargs, mode='r_dada2')
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
    microbe = Microbe(kargs, mode='biom')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.make_biom_by_db()


def taxonomy(kargs):
    microbe = Microbe(kargs, mode='taxonomy')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.assign_taxonomy()


def alpha_diversity(kargs):
    microbe = Microbe(kargs, mode='alpha_diversity')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_alpha_diversity()


def summarize_taxa(kargs):
    microbe = Microbe(kargs, mode='summarize_taxa')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.summrize_taxa()


def score():
    pass


def insert():
    pass
