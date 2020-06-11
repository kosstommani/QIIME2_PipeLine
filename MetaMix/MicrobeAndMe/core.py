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

from click import secho
from SpoON.util import parse_config
from MicrobeAndMe.data_structure import Microbe, ASVsMerger, FilesCollector
from SpoON.run_time import print_mam_metamix_run_time
import time

# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']


def all_pipeline(kargs):
    start_time = time.time()
    microbe = Microbe(kargs, mode='all', analysis_type='MicrobeAndMe')
    if (kargs['copy_rawdata'] is False) and (kargs['no_prema'] is False):
        microbe.run_prema()
        exit()
    elif (kargs['copy_rawdata'] is True) and (kargs['no_prema'] is False):
        microbe.run_prema()
    elif (kargs['copy_rawdata'] is False) and (kargs['no_prema'] is True):
        pass
    elif (kargs['copy_rawdata'] is True) and (kargs['no_prema'] is True):
        secho('Error: --copy_rawdata 와 --no_prema 옵션을 같이 사용할 수 없습니다.', fg='red')
    else:
        raise ValueError(f'예기치 않은 오류입니다. copy_rawdata:{kargs["copy_rawdata"]}, no_prema:{kargs["no_prema"]}')
    microbe.make_analysis_number_dir()
    microbe.set_path()
    microbe.make_sample_list()
    sample_end_time = time.time()
    microbe.run_dada2_using_r()
    microbe.check_jobs()
    dada2_end_time = time.time()
    microbe.assign_taxonomy()
    taxonomy_end_time = time.time()
    microbe.make_biom_by_db()
    biom_end_time = time.time()
    microbe.run_alpha_diversity()
    diverstiy_end_time = time.time()
    microbe.summarize_taxa()
    summarize_end_time = time.time()
    microbe.run_score()
    score_end_time = time.time()
    microbe.make_info()
    info_end_time = time.time()
    microbe.insert_score()
    insert_end_time = time.time()
    microbe.find_score()
    end_time = time.time()
    pipeline_end_time = {
        'make_sample': sample_end_time,
        'dada2': dada2_end_time,
        'assign_taxonomy': taxonomy_end_time,
        'biom': biom_end_time,
        'alpha_diversity': diverstiy_end_time,
        'summarize_taxa': summarize_end_time,
        'score': score_end_time,
        'info': info_end_time,
        'insert_score': insert_end_time,
    }
    print_mam_metamix_run_time(start_time, pipeline_end_time, end_time)


def dada2_using_r(kargs):
    microbe = Microbe(kargs, mode='r_dada2', analysis_type='MicrobeAndMe')
    # microbe = Microbe(kargs, mode='r_dada2', analysis_type='NGS')
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
    microbe.summarize_taxa()


def score(kargs):
    microbe = Microbe(kargs, mode='score', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_score()


def db_insert(kargs):
    microbe = Microbe(kargs, mode='db_insert', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.insert_score()


def db_find(kargs):
    microbe = Microbe(kargs, mode='db_find', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.find_score()


def db_del(kargs):
    microbe = Microbe(kargs, mode='db_del', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.del_score()


def info(kargs):
    microbe = Microbe(kargs, mode='info', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.make_info()


def delete(kargs):
    microbe = Microbe(kargs, mode='delete', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_delete()


def qdel(kargs):
    microbe = Microbe(kargs, mode='qdel', analysis_type='MicrobeAndMe')
    microbe.set_path()
    microbe.make_sample_list()
    microbe.run_qdel()

