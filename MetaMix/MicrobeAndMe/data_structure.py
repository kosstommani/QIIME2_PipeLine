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

import os
import sys
import shutil
from click import echo, secho, style, echo_via_pager, pause
from CoFI.data_structure import RunData
from CoFI.data_handler import check_duplication
from SpoON.util import parse_config, run_cmd, check_run_cmd, glob_dir, \
    read_data, check_file_type, check_jobs_in_queue, launcher_cmd, parse_job_id, \
    read_csv, check_db_version, get_samples_info, read_csv_to_dict, write_sample_info, \
    read_yaml, run_cmd_qsub, parse_job_id, check_jobs_in_queue, print_qsub_status_message, \
    make_sample_list, make_analysis_dir
from theCups.data_structure import Path
from MicrobeAndMe.dada2 import run_dada2_r, run_merge_asvs
from MicrobeAndMe.score import AndMe
from MicrobeAndMe.biom import BiomForMicrobeAndMe, BiomForNGS
from DBcontrol.db_microbe_and_me import DBMicrobeAndMe, DBDownloadSite
from collections import Counter
from tqdm import tqdm
import time
import hashlib
import hmac

# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']
ANALYSIS_DIR_NAME = CONFIG['Analysis_Dir_Name']
ITEMS_AND_BACTERIA_FILE = CONFIG['mamMetaMix']['items_and_bacteria']
DISTRIBUTION_DB = CONFIG['mamMetaMix']['dist_db']
TAX_ID = CONFIG['mamMetaMix']['tax_id']
DB = CONFIG['mamMetaMix']['DB']
DB_DOWNLOAD_SITE = CONFIG['mamMetaMix']['DB_DownloadSite']
HASH_KEY = CONFIG['mamMetaMix']['HashKey']


class PathMAM(Path):
    def sample_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name)

    def r_dada2_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, self.analysis_dir_name['r_dada2'])

    def taxonomy_assignment_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, self.analysis_dir_name['taxonomy_assignment'])

    def biom_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, self.analysis_dir_name['biom'])

    def alpha_diversity_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, 'Alpha_Diversity')

    def summarized_taxa_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, 'Summarized_Taxa')

    def score_path(self, p_sample_name: str) -> str:
        return os.path.join(self.analysis_number_path, p_sample_name, 'Score')


class Microbe:
    def __init__(self, kargs: dict, mode, analysis_type):
        """
        Microbe&Me 서비스를 위한 분석을 진행합니다.

        :param kargs:
        :param mode: [all | r_dada2 | taxonomy | biom | alpha_diversity |
                      summarize_taxa | score | insert | delete | qdel] 분석 명령어
        :param analysis_type: [NGS | MicrobeAndMe]
        """
        self.analysis_type = analysis_type
        if mode == 'all':
            self._all(kargs)
        elif mode == 'r_dada2':
            self._r_dada2(kargs)
        elif mode == 'taxonomy':
            self._taxonomy(kargs)
        elif mode == 'biom':
            self._biom(kargs)
        elif mode == 'alpha_diversity':
            self._alpha_diversity(kargs)
        elif mode == 'summarize_taxa':
            self._summarize_taxa(kargs)
        elif mode == 'score':
            self._score(kargs)
        elif mode == 'db_insert':
            self._db_insert(kargs)
        elif mode == 'db_find':
            self._db_find(kargs)
        elif mode == 'db_del':
            self._db_del(kargs)
        elif mode == 'info':
            self._info(kargs)
        elif mode == 'delete':
            self._delete(kargs)
        elif mode == 'qdel':
            self._qdel(kargs)
        else:
            raise ValueError(f'Microbe mode 설정 오류: {mode}')

    def _all(self, kargs):
        """
        PreMA -> R_DADA2 -> Taxonomy Assignment -> BIOM -> (작성 중)
        전체 파이프라인을 실행하기 위해 필요한 값들을 설정한다.
        :param kargs:
                    rawdata_base_path:
                    prema_target_dir_suffix:
                    copy_rawdata:
                    index_kit:
                    sample_read:
                    no_prema:
                    analysis_base_path:
                    analysis_dir_base_name:
                    cofi_target_dir_suffix:
                    r1_suffix:
                    r2_suffix:
                    queue:
                    slots:
                    no_queue:
                    order_number:
                    include:
                    exclude:
        :return:
        """
        self.rawdata_base_path = kargs['rawdata_base_path']
        self.prema_target_dir_suffix = kargs['prema_target_dir_suffix']
        self.copy_rawdata = kargs['copy_rawdata']
        self.index_kit = kargs['index_kit']
        self.sample_read = kargs['sample_read']
        self.no_prema = kargs['no_prema']
        self.analysis_base_path = kargs['analysis_base_path']
        self.analysis_dir_base_name = kargs['analysis_dir_base_name']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.queue = kargs['queue']
        self.slots = kargs['slots']
        self.no_queue = kargs['no_queue']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.database = kargs['database']
        self.table = kargs['table']
        self.print_record = False
        # 설정
        self._set_slots()
        self.i_path = None
        self.l_nt_run_list = None
        self.analysis_number = None

    def _r_dada2(self, kargs: dict):
        """
        R_DADA2 실행에 필요한 값들을 설정한다.
        :param kargs:
                    analysis_base_path:
                    analysis_dir_base_name:
                    cofi_target_dir_suffix:
                    r1_suffix:
                    r2_suffix:
                    queue:
                    slots:
                    no_queue:
                    order_number:
                    include:
                    exclude:
                NGS :
                    trim_left_f:
                    trim_left_r:
                    trunc_len_f:
                    trunc_len_r:
                    maxn:
                    maxee_f
                    maxee_r:
                    trunc_q:
                    min_overlap:
                    chimera_method:
        :return:
        """
        self.analysis_number = None
        self.analysis_base_path = kargs['analysis_base_path']
        self.analysis_dir_base_name = kargs['analysis_dir_base_name']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.queue = kargs['queue']
        self.slots = kargs['slots']
        self.no_queue = kargs['no_queue']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        if self.analysis_type == 'NGS':
            self.trim_left_f = kargs['trim_left_f']
            self.trim_left_r = kargs['trim_left_r']
            self.trunc_len_f = kargs['trunc_len_f']
            self.trunc_len_r = kargs['trunc_len_r']
            self.max_n = kargs['maxn']
            self.maxee_f = kargs['maxee_f']
            self.maxee_r = kargs['maxee_r']
            self.trunc_q = kargs['truncq']
            self.min_overlap = kargs['min_overlap']
            self.chimera_method = kargs['chimera_method']
        elif self.analysis_type == 'MicrobeAndMe':
            self.trim_left_f = None
            self.trim_left_r = None
            self.trunc_len_f = None
            self.trunc_len_r = None
            self.max_n = None
            self.maxee_f = None
            self.maxee_r = None
            self.trunc_q = None
            self.min_overlap = None
            self.chimera_method = None
        else:
            raise ValueError(f'self.analysis_type: {self.analysis_type}. 적절하지 않은 값입니다.')
        # 설정
        self._set_slots()
        self.i_path = None
        self.l_nt_run_list = None

    def _taxonomy(self, kargs: dict):
        """
        Taxonomy Assignment 실행에 필요한 값들을 설정한다.
        :param kargs:
                    analysis_number:
                    analysis_base_path:
                    analysis_dir_base_name:
                    cofi_target_dir_suffix:
                    r1_suffix:
                    r2_suffix:
                    queue:
                    no_queue:
                    order_number:
                    include:
                    exclude:
        :return:
        """
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.queue = kargs['queue']
        self.no_queue = kargs['no_queue']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.database = kargs['database']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _biom(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.order_number = kargs['order_number']
        if self.analysis_type == 'MicrobeAndMe':
            self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
            self.r1_suffix = kargs['r1_suffix']
            self.r2_suffix = kargs['r2_suffix']
            self.queue = kargs['queue']
            self.no_queue = kargs['no_queue']
            self.include = kargs['include']
            self.exclude = kargs['exclude']
        elif self.analysis_type == 'NGS':
            self.no_queue = True
            self.metadata = kargs['metadata']
        self.database = kargs['database']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _alpha_diversity(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.queue = kargs['queue']
        self.no_queue = kargs['no_queue']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _summarize_taxa(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.queue = kargs['queue']
        self.no_queue = kargs['no_queue']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.database = kargs['database']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _score(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _db_insert(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _db_find(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.table = kargs['table']
        self.print_record = kargs['print_record']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _db_del(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.table = kargs['table']
        self.print_record = kargs['print_record']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _info(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _delete(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.t_dir_name = kargs['dir']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    def _qdel(self, kargs):
        self.analysis_number = kargs['analysis_number']
        self.analysis_base_path = kargs['analysis_base_path']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.target_dir = kargs['target']
        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    @property
    def max_slots(self):
        if self.queue == 'bi3m.q':
            return 50
        elif self.queue == 'meta.q':
            return 144
        elif self.queue == 'meta.q@denovo06':
            return 48
        elif self.queue == 'meta.q@denovo11':
            return 144

    @property
    def l_db(self):
        """
        마이크로브앤미 전용
        :return: 
        """
        if self.database == 'all':
            return ['NCBI_16S', 'NCBI_Probiotics21']
        else:
            return [self.database]

    def _set_slots(self):
        """
        Queue 당 사용할 수 있는 slots 의 수를 지정한다.
        bi3m.q : core 56개 - cm456(Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz) -
        meta.q@denovo06 : core 48개 - denovo06(Intel(R) Xeon(R) CPU E7540  @ 2.00GHz)

                   Passmark CPU Mark    CPU Value(CPUmark/$Price)    Price
        Gold 6123:            27,522	                     6.34	 $4,339.90(Total for 2 CPUs (2020-03-18))
        E7540    :             7,946                       124.38    $63.88(2018-12-31)

        :return:
        """
        if self.queue == 'bi3m.q':
            if self.slots is None:
                self.slots = self.max_slots
        elif self.queue == 'meta.q':
            if self.slots is None:
                self.slots = self.max_slots
        elif self.queue == 'meta.q@denovo06':
            if self.slots is None:
                self.slots = self.max_slots
        elif self.queue == 'meta.q@denovo11':
            if self.slots is None:
                self.slots = self.max_slots
        else:
            if self.no_queue is False:
                secho('Error: Queue 지정이 잘못되었습니다.', fg='red', blink=True, err=True)
                echo(f'Queue: {self.queue}', err=True)
                exit(1)

    def make_analysis_number_dir(self):
        if self.analysis_number is None:
            analysis_path = make_analysis_dir(self.analysis_base_path, self.order_number, self.analysis_dir_base_name)
            self.analysis_number = os.path.basename(analysis_path)
        else:
            raise ValueError(f'self.analysis_number의 값이 있습니다. ({self.analysis_number})')

    def set_path(self):
        """
        분석시 생성되는 디렉터리 경로를 설정한다.
        make_analysis_number_dir() 메소드를 먼저 실행하거나 self.analysis_number에 값을 설정해야한다.

        :return:
        """
        if self.analysis_number is None:
            raise ValueError('self.analysis_number의 값이 None입니다. make_analysis_number_dir method를 실행하세요.')
        else:
            self.i_path = PathMAM(self.analysis_base_path, self.order_number, self.analysis_number)

    @property
    def order_dir_path(self):
        return os.path.join(self.analysis_base_path, self.order_number)

    @property
    def run_dir_path(self):
        """
        Analysis --> RawData 디렉터리에 Run 디렉터리에서 1개만 선택되도록 지정한다.
        MicrobeAndMe 분석에 적용된다.
        :return:
        """
        return glob_dir(self.order_dir_path, self.cofi_target_dir_suffix, 'only')

    @property
    def sample_count(self):
        return sum([len(nt_run.samples) for nt_run in self.l_nt_run_list])

    def make_sample_list(self):
        """
        분석을 진행할 시료들의 목록을 작성한다.

        :return: l_nt_run_list - named tuple인 RunData 를 원소로 가지는 리스트
                 RunData.path
                        .run_info
                        .samples: named tuple인 SampleList을 원소로 가지는 리스트
                            SampleList.path: 경로
                                      .name: 시료명
        """
        if self.analysis_type == 'NGS':
            pattern = self.cofi_target_dir_suffix
        elif self.analysis_type == 'MicrobeAndMe':
            pattern = f'RawData*/{os.path.basename(self.run_dir_path)}'
        else:
            secho('Error: self.analysis_type 지정이 잘못되었습니다.', fg='red', blink=True, err=True)
            echo(f'self.analysis_type: {self.analysis_type}')
            exit()

        self.l_nt_run_list = make_sample_list(
            self.analysis_base_path,
            self.order_number,
            pattern,
            None,
        )
        if self.include is not None:
            l_include = self.include.strip().split()
            new_l_nt_run_list = list()
            for nt_run in self.l_nt_run_list:
                l_sample = list()
                for sample in nt_run.samples:
                    if sample.name in l_include:
                        l_sample.append(sample)
                new_run_data = RunData(path=nt_run.path, run_info=nt_run.run_info, samples=l_sample)
                new_l_nt_run_list.append(new_run_data)
            self.l_nt_run_list = new_l_nt_run_list
            echo('\n')  # 빈 줄
            secho('+++ 시료 목록 수정: include +++', fg='yellow', blink=True)
            check_duplication(self.l_nt_run_list)
        elif self.exclude is not None:
            l_exclude = self.exclude.strip().split()
            new_l_nt_run_list = list()
            for nt_run in self.l_nt_run_list:
                l_index = list()
                for index, sample in enumerate(nt_run.samples):
                    if sample.name in l_exclude:
                        l_index.append(index)
                l_index.reverse()
                _ = [nt_run.samples.pop(index) for index in l_index]
                new_run_data = RunData(path=nt_run.path, run_info=nt_run.run_info, samples=nt_run.samples)
                new_l_nt_run_list.append(new_run_data)
            self.l_nt_run_list = new_l_nt_run_list
            echo('\n')  # 빈 줄
            secho('--- 시료 목록 수정: exclude ---', fg='yellow', blink=True)
            check_duplication(self.l_nt_run_list)

    def get_prema_cmd(self):
        global PREMA
        prema_cmd = \
            f'{PREMA} ' \
            'CORE ' \
            f'--rawdata_base_path {self.rawdata_base_path} ' \
            f'--analysis_base_path {self.analysis_base_path} ' \
            f'--target_dir_suffix "{self.prema_target_dir_suffix}" ' \
            f'--index_kit {self.index_kit} ' \
            '--microbe_and_me ' \
            f'{"-c " if self.copy_rawdata else ""}' \
            f'{self.order_number}'
        return prema_cmd

    def run_prema(self):
        prema_cmd = self.get_prema_cmd()
        run = run_cmd(prema_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'PreMA CORE 완료',
                'false_meg': 'PreMA CORE',
            }, p_stdout=False, p_stderr=False, p_exit=True,
        )

    def set_queue_mode(self, p_count: int):
        if self.queue_mode == 'auto':
            if p_count >= 50:
                self.slots = 1
            elif p_count >= 25:
                self.slots = 2
            elif p_count > 1:
                self.slots = int((self.max_slots - p_count) / p_count) + 1
            elif p_count == 1:
                self.slots = 10
            else:
                raise ValueError('p_count는 양수.')

    def run_dada2_using_r(self):
        """

        :return:
        """
        if self.analysis_type == 'MicrobeAndMe':
            if len(self.l_nt_run_list) > 1:
                secho('Error: Run Dir 의 개수가 1이 아닙니다.', fg='red', blink=True, err=True)
                echo([nt_run.path for nt_run in self.l_nt_run_list], err=True)
                exit(1)
        # self.set_queue_mode(self.sample_count)
        d_return_code = dict()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            run_dir = nt_run.path.split(self.order_number)[1]
            for nt_sample in tqdm(nt_run.samples, desc='R_DADA2', ascii=True, leave=False):
                if self.analysis_type == 'MicrobeAndMe':
                    return_code = run_dada2_r(
                        {
                            'order_number': self.order_number,
                            'analysis_number': self.analysis_number,
                            'run_dir': run_dir,
                            'analysis_base_path': self.analysis_base_path,
                            'sample_name': nt_sample.name,
                            'queue': self.queue,
                            'no_queue': self.no_queue,
                            'slots': self.slots,
                        }, self.analysis_type
                    )
                elif self.analysis_type == 'NGS':
                    return_code = run_dada2_r(
                        {
                            'order_number': self.order_number,
                            'analysis_number': self.analysis_number,
                            'run_dir': run_dir,
                            'analysis_base_path': self.analysis_base_path,
                            'sample_name': nt_sample.name,
                            'queue': self.queue,
                            'no_queue': self.no_queue,
                            'slots': self.slots,
                            'trim_left_f': self.trim_left_f,
                            'trim_left_r': self.trim_left_r,
                            'trunc_len_f': self.trunc_len_f,
                            'trunc_len_r': self.trunc_len_r,
                            'max_n': self.max_n,
                            'maxee_f': self.maxee_f,
                            'maxee_r': self.maxee_r,
                            'trunc_q': self.trunc_q,
                            'min_overlap': self.min_overlap,
                            'chimera_method': self.chimera_method,
                        }, self.analysis_type
                    )
                else:
                    raise ValueError(f'self.analysis_type: {self.analysis_type}. 적절하지 않은 값입니다.')
                d_return_code[nt_sample.name] = return_code

        if all(d_return_code.values()):
            if self.no_queue is True:
                secho(f'>>> DADA2 작업 - {len(d_return_code.values())} 완료', fg='cyan')
            else:
                secho(f'>>> Queue 작업 등록 - {len(d_return_code.values())} 완료', fg='cyan')
        else:
            failed_job = [key for key, value in d_return_code.items() if value != 0]
            success_job = [key for key, value in d_return_code.items() if value == 0]
            if len(failed_job) == 0:
                failed_text = f'{len(failed_job)} 실패'
            else:
                failed_text = style(f'{len(failed_job)} 실패', fg='red', blink=True)
            success_text = style(f'{len(success_job)} 완료', fg='cyan')
            if self.no_queue is True:
                echo(f'>>> DADA2 작업 - {success_text}, {failed_text}')
            else:
                echo(f'>>> Queue 작업 등록 - {success_text}, {failed_text}')
                echo()  # 빈 줄

    def check_jobs(self):
        """
        R_DADA2 전용        
        :return:
        """
        if self.no_queue is False:
            l_job_id = list()
            while self.sample_count != len(l_job_id):
                l_job_id = list()
                l_qsub_log = glob_dir(self.i_path.analysis_number_path, '*_qsub.log', 'many')
                if l_qsub_log is None:
                    secho('Error: qsub.log 파일이 없습니다.', fg='red', blink=True, err=True)
                    echo(self.i_path.analysis_number_path)
                for log in l_qsub_log:
                    with open(log, 'r') as o_log:
                        # log_text : Your job 5125107 ("Job_MBS191230NE001") has been submitted
                        job_id = o_log.read().strip().replace('Your job ', '').split()[0]
                        l_job_id.append(int(job_id))

            echo('!! Queue 에 등록된 작업은 mamMetaMix 프로그램과 연결이 끊기더라도 동작합니다.')
            echo('============ R_DADA2 : Queue 작업 상태 ============'.center(60))
            check_jobs_in_queue(l_job_id)
            echo()  # 빈 줄

    def assign_taxonomy(self):
        l_cmd = list()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Assign_Taxonomy', ascii=True, leave=False):
                asv_fasta_file = os.path.join(self.i_path.r_dada2_path(nt_sample.name), 'ASVs.fasta')
                taxonomy_dir_path = self.i_path.taxonomy_assignment_path(nt_sample.name)
                os.makedirs(taxonomy_dir_path)
                for db in self.l_db:
                    cofi_taxonomy_cmd = \
                        f'{COFI} TAXONOMY None None ' \
                        f'-i {asv_fasta_file} ' \
                        f'--out_dir {taxonomy_dir_path} ' \
                        f'--queue {self.queue} ' \
                        f'--blast_db {db}'
                    l_cmd.append(cofi_taxonomy_cmd)
        echo()  # 빈 줄
        launcher_cmd(l_cmd, 'Taxonomy Assignment', 50, False)

    def make_biom(self):
        """
        서버 사용
        
        :return:
        """
        if self.analysis_type == 'MicrobeAndMe':
            for db in self.l_db:
                self.make_biom_mam(db)
        elif self.analysis_type == 'NGS':
            for db in self.database:
                self.make_biom_ngs(db)

    @staticmethod
    def check_job_for_biom(l_run_process: list, base_biom_status: bool) -> None:
        failed_job = [run for run in l_run_process if run.returncode != 0]
        success_job = [run for run in l_run_process if run.returncode == 0]
        if base_biom_status is True:
            secho('>>> BASE BIOM 파일이 존재합니다.', fg='yellow')
            secho('\t---> BASE BIOM 파일 생성이 생략됩니다.', fg='magenta')
        if len(failed_job) == 0:
            failed_text = f'{len(failed_job)} 실패'
        else:
            failed_text = style(f'{len(failed_job)} 실패', fg='red', blink=True)
        success_text = style(f'{len(success_job)} 완료', fg='cyan')
        echo(f'>>> BIOM 생성(작업 3개) - {success_text}, {failed_text}')
        if len(failed_job) != 0:
            for run in l_run_process:
                check_run_cmd(
                    {
                        'run': run,
                        'true_meg': None,
                        'false_meg': 'BIOM 생성',
                    }, p_exit=False,
                )
            exit()
        else:
            echo()  # 빈 줄

    def make_biom_mam(self, p_db):
        """
        마이크로브앤미 처리용

        :param p_db:
        :return:
        """
        l_run_process = list()
        base_biom_status = False
        echo('>>> BIOM 파일 생성 시작')
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc=f'BIOM({p_db})', ascii=True, leave=False):
                i_biom = BiomForMicrobeAndMe(
                    {
                        'r_dada2_path': self.i_path.r_dada2_path(nt_sample.name),
                        'biom_path': self.i_path.biom_path(nt_sample.name),
                        'taxonomy_path': self.i_path.taxonomy_assignment_path(nt_sample.name),
                        'db_name': p_db,
                    }, p_verbose=False, p_exit=False,
                )
                i_biom.check_biom_dir()
                if check_file_type(i_biom.base_biom_file, 'exists') is True:
                    base_biom_status = True
                else:
                    l_run_process.append(i_biom.convert_biom())  # Base BIOM
                l_run_process.append(i_biom.add_metadata())  # Taxonomy BIOM
                l_run_process.append(i_biom.convert_taxonomy_table())
        self.check_job_for_biom(l_run_process, base_biom_status)

    def make_biom_ngs(self, p_db: str) -> None:
        """
        NGS 처리용

        :param p_db:
        :param p_metadata
        :return:
        """
        l_run_process = list()
        base_biom_status = False
        echo(f'>>> BIOM 파일 생성 시작: {p_db}')
        i_biom = BiomForNGS(
            {
                'r_dada2_path': self.i_path.r_dada2_path(''),
                'biom_path': self.i_path.biom_path(''),
                'taxonomy_path': self.i_path.taxonomy_assignment_path(''),
                'db_name': p_db,
                'metadata': self.metadata,
            }, p_verbose=False, p_exit=True,
        )
        i_biom.check_biom_dir()
        if check_file_type(i_biom.base_biom_file, 'exists') is True:
            base_biom_status = True
        else:
            l_run_process.append(i_biom.convert_biom())  # Base BIOM
            echo('    Summary File 생성 추가')
            l_run_process.append(i_biom.make_summary_file())  # Summary File for Base BIOM
        l_run_process.append(i_biom.add_metadata(p_add_sample_metadata=True))  # Taxonomy BIOM
        l_run_process.append(i_biom.convert_taxonomy_table())
        self.check_job_for_biom(l_run_process, base_biom_status)

    def make_biom_using_queue(self, p_db: str, p_no_base: bool = False):
        if self.analysis_type == 'MicrobeAndMe':
            Biom = BiomForMicrobeAndMe
        elif self.analysis_type == 'NGS':
            Biom = BiomForNGS
        else:
            raise ValueError(f'self.analysis_type: {self.analysis_type}')
        l_qsub_log = list()
        l_returncode = list()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc=f'BIOM({p_db})', ascii=True, leave=False):
                i_biom = Biom(
                    {
                        'r_dada2_path': self.i_path.r_dada2_path(nt_sample.name),
                        'biom_path': self.i_path.biom_path(nt_sample.name),
                        'taxonomy_path': self.i_path.taxonomy_assignment_path(nt_sample.name),
                        'db_name': p_db,
                    }, p_verbose=False, p_exit=False,
                )
                i_biom.check_biom_dir()
                if p_no_base is False:
                    i_biom.convert_biom()  # Base BIOM
                i_biom.add_metadata()  # Taxonomy BIOM
                i_biom.convert_taxonomy_table()
                l_returncode.append(run_cmd_qsub(
                    {
                        'stream': 'y',
                        'queue': self.queue,
                        'interpreter': '/bin/bash',
                        'cmd': i_biom.recipe_file,
                        'out_path': i_biom.biom_path,
                    }, p_verbose=False, p_exit=False, p_qsub_base_name=f'qsub.{p_db}',
                ))
                l_qsub_log.append((i_biom.biom_path, f'qsub.{p_db}.log'))
        print_qsub_status_message(l_returncode)
        l_job_id = list()
        for biom_path, log in tqdm(l_qsub_log, desc='Job ID', ascii=True):
            l_job_id.append(parse_job_id(biom_path, log, 'only', p_verbose=False))
        check_jobs_in_queue(l_job_id)

    def make_biom_by_db(self):
        """
        mamMetaMix - Microbe&Me 
        :return:
        """
        echo('>>> BIOM 생성 시작')
        if self.analysis_type == 'MicrobeAndMe':
            l_db = self.l_db
        elif self.analysis_type == 'NGS':
            l_db = self.database
        else:
            raise ValueError(f'self.analysis_type: {self.analysis_type}')
        for count, db in enumerate(l_db, 1):
            if self.no_queue is True:
                self.make_biom(db)
            else:
                if count == 1:
                    self.make_biom_using_queue(db)
                else:
                    self.make_biom_using_queue(db, p_no_base=True)
        secho('>>> BIOM 생성 완료', fg='cyan')
        echo()  # 빈 줄

    def run_alpha_diversity(self):
        l_run_process = list()
        l_returncode = list()
        l_alpha_diversity_path = list()
        echo('>>> Alpha Diversity 생성 시작')
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc=f'Alpha Diversity', ascii=True, leave=False):
                biom = os.path.join(self.i_path.biom_path(nt_sample.name), 'ASVs.biom')
                alpha_diversity_path = self.i_path.alpha_diversity_path(nt_sample.name)
                os.makedirs(alpha_diversity_path)
                adiv_file = os.path.join(alpha_diversity_path, 'adiv.txt')
                cmd = \
                    'alpha_diversity.py ' \
                    f'-i {biom} ' \
                    '-m shannon,simpson ' \
                    f'-o {adiv_file}'
                recipe_file = os.path.join(alpha_diversity_path, 'alpha_diversity.recipe')
                with open(recipe_file, 'w') as o_recipe:
                    o_recipe.write(cmd)
                    o_recipe.write('\n')
                if self.no_queue is True:
                    run = run_cmd(cmd)
                    l_run_process.append(run)
                else:
                    l_returncode.append(run_cmd_qsub(
                        {
                            'stream': 'y',
                            'queue': self.queue,
                            'interpreter': '/bin/bash',
                            'cmd': f'{recipe_file}',
                            'out_path': alpha_diversity_path,
                        }, p_verbose=False, p_exit=False,
                    ))
                    l_alpha_diversity_path.append(alpha_diversity_path)
        if self.no_queue is True:
            failed_job = [run for run in l_run_process if run.returncode != 0]
            success_job = [run for run in l_run_process if run.returncode == 0]
            if len(failed_job) == 0:
                failed_text = f'{len(failed_job)} 실패'
            else:
                failed_text = style(f'{len(failed_job)} 실패', fg='red', blink=True)
            success_text = style(f'{len(success_job)} 완료', fg='cyan')
            echo(f'>>> Alpha Diversity - {success_text}, {failed_text}')
            if len(failed_job) != 0:
                for run in l_run_process:
                    check_run_cmd(
                        {
                            'run': run,
                            'true_meg': None,
                            'false_meg': 'alpha_diversity',
                        }, p_exit=False,
                    )
                exit()
            else:
                echo()  # 빈 줄
        else:
            print_qsub_status_message(l_returncode)
            l_job_id = list()
            for adiv_path in tqdm(l_alpha_diversity_path, desc='Job ID', ascii=True):
                l_job_id.append(parse_job_id(adiv_path, 'qsub.log', 'only', p_verbose=False))
            check_jobs_in_queue(l_job_id)
        secho('>>> Alpha Diversity 완료', fg='cyan')
        echo()  # 빈 줄

    def summarize_taxa(self):
        l_run_process = list()
        l_returncode = list()
        l_qsub_log = list()
        echo('>>> Summarize Taxa 시작')
        if self.analysis_type == 'MicrobeAndMe':
            Biom = BiomForMicrobeAndMe
        elif self.analysis_type == 'NGS':
            Biom = BiomForNGS
        else:
            raise ValueError(f'self.analysis_type: {self.analysis_type}')
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Summarize Taxa', ascii=True, leave=False):
                summarized_taxa_path = self.i_path.summarized_taxa_path(nt_sample.name)
                os.makedirs(summarized_taxa_path)
                biom_dir_path = self.i_path.biom_path(nt_sample.name)
                for db in self.l_db:
                    i_biom = Biom(
                        {
                            'r_dada2_path': None,
                            'biom_path': biom_dir_path,
                            'taxonomy_path': None,
                            'db_name': db,
                        }
                    )
                    cmd = \
                        'summarize_taxa.py ' \
                        f'-i {i_biom.taxonomy_biom_file} ' \
                        f'-o {summarized_taxa_path} ' \
                        '-L 2,5,6,7'
                    recipe_file = os.path.join(summarized_taxa_path, f'summarize_taxa.{db}.recipe')
                    with open(recipe_file, 'w') as o_recipe:
                        o_recipe.write(cmd)
                        o_recipe.write('\n')
                    if self.no_queue is True:
                        run = run_cmd(cmd)
                        l_run_process.append(run)
                    else:
                        l_returncode.append(run_cmd_qsub(
                            {
                                'stream': 'y',
                                'queue': self.queue,
                                'interpreter': '/bin/bash',
                                'cmd': f'{recipe_file}',
                                'out_path': summarized_taxa_path,
                            }, p_verbose=False, p_exit=False, p_qsub_base_name=f'qsub.{db}'
                        ))
                        l_qsub_log.append((summarized_taxa_path, f'qsub.{db}.log'))
        if self.no_queue is True:
            failed_job = [run for run in l_run_process if run.returncode != 0]
            success_job = [run for run in l_run_process if run.returncode == 0]
            if len(failed_job) == 0:
                failed_text = f'{len(failed_job)} 실패'
            else:
                failed_text = style(f'{len(failed_job)} 실패', fg='red', blink=True)
            success_text = style(f'{len(success_job)} 완료', fg='cyan')
            echo(f'>>> Summarize Taxa - {success_text}, {failed_text}')
            if len(failed_job) != 0:
                for run in l_run_process:
                    check_run_cmd(
                        {
                            'run': run,
                            'true_meg': None,
                            'false_meg': 'summarize_taxa',
                        }, p_exit=False,
                    )
                exit()
            else:
                echo()  # 빈 줄
        else:
            print_qsub_status_message(l_returncode)
            l_job_id = list()
            for summarized_taxa_path, log in tqdm(l_qsub_log, desc='Job ID', ascii=True):
                l_job_id.append(parse_job_id(summarized_taxa_path, log, 'only', p_verbose=False))
            check_jobs_in_queue(l_job_id)
        secho('>>> Summarize Taxa 완료', fg='cyan')
        echo()  # 빈 줄

    def run_score(self):
        global ITEMS_AND_BACTERIA_FILE, DISTRIBUTION_DB
        for nt_run in self.l_nt_run_list:
            for nt_sample in nt_run.samples:
                echo()  # 빈 줄
                secho(f'=== Score 계산 : {nt_sample.name} ===', fg='yellow', bold=True)
                and_me = AndMe(
                    {
                        'bact_file': ITEMS_AND_BACTERIA_FILE,
                        'dist_db_file': DISTRIBUTION_DB,
                        'tax_id_file': TAX_ID,
                        'phylum_file': os.path.join(self.i_path.summarized_taxa_path(nt_sample.name),
                                                    'ASVs.blast_NCBI_16S_L2.txt'),
                        'family_file': os.path.join(self.i_path.summarized_taxa_path(nt_sample.name),
                                                    'ASVs.blast_NCBI_16S_L5.txt'),
                        'genus_file': os.path.join(self.i_path.summarized_taxa_path(nt_sample.name),
                                                   'ASVs.blast_NCBI_16S_L6.txt'),
                        'species_file': os.path.join(self.i_path.summarized_taxa_path(nt_sample.name),
                                                     'ASVs.blast_NCBI_16S_L7.txt'),
                        'alpha_diversity_file': os.path.join(self.i_path.alpha_diversity_path(nt_sample.name),
                                                             'adiv.txt'),
                        'probiotics19_species_file': os.path.join(self.i_path.summarized_taxa_path(nt_sample.name),
                                                                  'ASVs.blast_NCBI_Probiotics21_L7.txt'),
                        'sample_name': nt_sample.name,
                        'out_path': self.i_path.sample_path(nt_sample.name),
                    }
                )
                and_me.read_bacteria_file()
                and_me.read_dist_db_file()
                and_me.read_tax_id_file()
                and_me.read_phylum_file()
                and_me.read_family_file()
                and_me.read_genus_file()
                and_me.read_species_file()
                and_me.make_score_dir()
                and_me.read_probiotics19_species_file()
                and_me.compute_intestinal_health()
                and_me.compute_wellness()
                and_me.compute_more_bowel_environments()
                and_me.compute_supplement()
                and_me.compute_diversity_index()
                and_me.compute_supplement_probiotics21()
                # and_me.print_scores()
                and_me.save_scores()

    @staticmethod
    def connect_db():
        global DB
        db_microbe = DBMicrobeAndMe(
            {
                'ip': DB['ip'],
                'user': DB['user'],
                'password': DB['password'],
                'db': DB['db'],
                'port': 3306,
            }
        )
        return db_microbe

    @staticmethod
    def connect_db_download_site():
        global DB_DOWNLOAD_SITE
        db_site = DBDownloadSite(
            {
                'ip': DB_DOWNLOAD_SITE['ip'],
                'user': DB_DOWNLOAD_SITE['user'],
                'password': DB_DOWNLOAD_SITE['password'],
                'db': DB_DOWNLOAD_SITE['db'],
                'port': DB_DOWNLOAD_SITE['port'],
            }
        )
        return db_site

    def insert_score(self):
        date_fmt = '%Y-%m-%d %H:%M:%S'
        analysis_time = time.strftime(date_fmt, time.localtime())

        # Score 파일 및 Sample Info 파일 존재 확인
        true_status_points = 0
        false_status_points = 0
        false_name_points = list()
        true_status_probiotics = 0
        false_status_probiotics = 0
        false_name_probiotics = list()
        true_status_info = 0
        false_status_info = 0
        false_name_info = list()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                points_csv = os.path.join(self.i_path.sample_path(nt_sample.name), 'Score',
                                          f'{nt_sample.name}_points.csv')
                probiotics_csv = os.path.join(self.i_path.sample_path(nt_sample.name), 'Score',
                                              f'{nt_sample.name}_probiotics.csv')
                sample_info_csv = os.path.join(self.i_path.sample_path(nt_sample.name), f'{nt_sample.name}_Info.csv')
                if check_file_type(points_csv, 'exists') is True:
                    true_status_points += 1
                else:
                    false_status_points += 1
                    false_name_points.append(nt_sample.name)
                if check_file_type(probiotics_csv, 'exists') is True:
                    true_status_probiotics += 1
                else:
                    false_status_probiotics += 1
                    false_name_probiotics.append(nt_sample.name)
                if check_file_type(sample_info_csv, 'exists') is True:
                    true_status_info += 1
                else:
                    false_status_info += 1
                    false_name_info.append(nt_sample.name)
        # points.csv 파일 존재 유무 확인 메시지 출력
        if false_status_points > 0:
            false_points_text = style(f'실패: {false_status_points}', fg='red', blink=True)
        else:
            false_points_text = style(f'실패: {false_status_points}', fg='cyan')
        true_points_text = style(f'>>> points.csv 파일 확인 - 완료: {true_status_points}', fg='cyan')
        echo(f'{true_points_text}, {false_points_text}')
        if false_status_points > 0:
            echo(false_name_points)
        # probiotics.csv 파일 존재 유무 확인 메시지 출력
        if false_status_probiotics > 0:
            false_probiotics_text = style(f'실패: {false_status_probiotics}', fg='red', blink=True)
        else:
            false_probiotics_text = style(f'실패: {false_status_probiotics}', fg='cyan')
        true_probiotics_text = style(f'>>> probiotics.csv 파일 확인 - 완료: {true_status_probiotics}', fg='cyan')
        echo(f'{true_probiotics_text}, {false_probiotics_text}')
        if false_status_probiotics > 0:
            echo(false_name_probiotics)
        # info.csv 파일 존재 유무 확인 메시지 출력
        if false_status_info > 0:
            false_info_text = style(f'실패: {false_status_info}', fg='red', blink=True)
        else:
            false_info_text = style(f'실패: {false_status_info}', fg='cyan')
        true_info_text = style(f'>>> info.csv 파일 확인 - 완료: {true_status_info}', fg='cyan')
        echo(f'{true_info_text}, {false_info_text}')
        if false_status_info > 0:
            echo(false_name_info)
        # 프로그램 종료 여부 결정
        if false_status_points + false_status_probiotics + false_status_info > 0:
            secho('--> 미존재 파일들을 확인하세요.', fg='magenta', bold=True)
            exit()

        # DB - 삽입
        secho('>>> DB 삽입 진행', bold=True)
        db_microbe = self.connect_db()
        db_microbe.start_transaction()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                # DB version 확인
                try:
                    points_csv = os.path.join(self.i_path.sample_path(nt_sample.name), 'Score',
                                              f'{nt_sample.name}_points.csv')
                    probiotics_csv = os.path.join(self.i_path.sample_path(nt_sample.name), 'Score',
                                                  f'{nt_sample.name}_probiotics.csv')
                    l_points = read_csv(points_csv, p_verbose=False)
                    l_probiotics = read_csv(probiotics_csv, p_verbose=False)
                    db_version_points = check_db_version(l_points)
                    db_version_probiotics = check_db_version(l_probiotics)
                    if db_version_points is False:
                        secho('Error: points.csv 파일에 여러 개의 DB 버전이 존재합니다.', fg='red', err=True)
                        echo(points_csv)
                    if db_version_probiotics is False:
                        secho('Error: probiotics.csv 파일에 여러 개의 DB 버전이 존재합니다.', fg='red', err=True)
                        echo(probiotics_csv)
                    if all([db_version_points, db_version_probiotics]) is False:
                        exit(1)
                    else:
                        db_version = l_points[0][-1]
                except Exception as err:
                    secho('Error: DB Version 확인 중 Exception 발생', fg='red', blink=True, err=True)
                    echo(err)
                    db_microbe.rollback()
                # Sample Info 파일 읽기
                sample_info_csv = os.path.join(self.i_path.sample_path(nt_sample.name), f'{nt_sample.name}_Info.csv')
                d_sample_info = read_csv_to_dict(sample_info_csv)[0]
                # 데이터 삽입
                try:
                    db_microbe.insert_analysis_version((0, nt_sample.name, db_version, analysis_time, 0))
                    db_microbe.insert_sample_info(d_sample_info, self.order_number, time.time())
                    db_microbe.insert_point(points_csv)
                    db_microbe.insert_biome(probiotics_csv)
                except KeyError as err:
                    secho('Error: 데이터베이스에 데이터 입력 중 KeyError 발생', fg='red', blink=True, err=True)
                    echo(err)
                    db_microbe.rollback()
                except Exception as err:
                    secho('Error: 데이터베이스에 데이터 입력 중 Exception 발생', fg='red', blink=True, err=True)
                    echo(err)
                    db_microbe.rollback()
        db_microbe.commit()
        db_microbe.close()

    @staticmethod
    def make_signature(self):
        global HASH_KEY
        hash_key = bytes(HASH_KEY, 'utf-8')
        code = bytes('TEST', 'utf-8')
        o_hash = hmac.new(hash_key, code, hashlib.sha256)
        signature = o_hash.hexdigest()
        return signature

    def find_score(self):
        d_bact_data = read_yaml(ITEMS_AND_BACTERIA_FILE)
        report_version = d_bact_data['Items_and_Bacteria']['report_version']
        db_microbe = self.connect_db()
        l_checked_sample = list()
        l_no_sample = list()
        l_wrong_sample = list()
        l_checked_sample_analysis = list()
        l_no_sample_analysis = list()
        l_wrong_sample_analysis = list()
        l_checked_sample_point = list()
        l_no_sample_point = list()
        l_wrong_sample_point = list()
        l_checked_sample_microbiome = list()
        l_no_sample_microbiome = list()
        l_wrong_sample_micorbiome = list()
        l_checked_analysis_version = list()
        l_no_analysis_version = list()
        l_wrong_analysis_version = list()
        d_db_data = dict()
        secho('>>> DB 레코드 검색 진행', bold=True)
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                # Sample Info 파일 읽기
                sample_info_csv = os.path.join(self.i_path.sample_path(nt_sample.name), f'{nt_sample.name}_Info.csv')
                d_sample_info = read_csv_to_dict(sample_info_csv, p_verbose=False)[0]
                sample = db_microbe.select_sample(d_sample_info['client_name'],
                                                  d_sample_info['Birth'],
                                                  d_sample_info['Sex'])
                sample_analysis = db_microbe.select_sample_analysis(d_sample_info['KitId'])
                sample_point = db_microbe.select_sample_point(d_sample_info['KitId'], report_version)
                sample_microbiome = db_microbe.select_sample_microbiome(d_sample_info['KitId'], report_version)
                analysis_version = db_microbe.select_analysis_version(d_sample_info['KitId'], report_version)
                d_db_data[nt_sample.name] = {'sample': sample,
                                             'sample_analysis': sample_analysis,
                                             'sample_point': sample_point,
                                             'sample_microbiome': sample_microbiome,
                                             'analysis_version': analysis_version}
                if len(sample) == 1:
                    l_checked_sample.append(nt_sample.name)
                elif len(sample) == 0:
                    l_no_sample.append(nt_sample.name)
                else:
                    l_wrong_sample.append(nt_sample.name)
                if len(sample_analysis) == 1:
                    l_checked_sample_analysis.append(nt_sample.name)
                elif len(sample_analysis) == 0:
                    l_no_sample_analysis.append(nt_sample.name)
                else:
                    l_wrong_sample_analysis.append(nt_sample.name)
                if len(sample_point) == 30:
                    l_checked_sample_point.append(nt_sample.name)
                elif len(sample_point) == 0:
                    l_no_sample_point.append(nt_sample.name)
                else:
                    l_wrong_sample_point.append(nt_sample.name)
                if len(sample_microbiome) == 21:
                    l_checked_sample_microbiome.append(nt_sample.name)
                elif len(sample_microbiome) == 0:
                    l_no_sample_microbiome.append(nt_sample.name)
                else:
                    l_wrong_sample_micorbiome.append(nt_sample.name)
                if len(analysis_version) == 1:
                    l_checked_analysis_version.append(nt_sample.name)
                elif len(analysis_version) == 0:
                    l_no_analysis_version.append(nt_sample.name)
                else:
                    l_wrong_analysis_version.append(nt_sample.name)
        db_microbe.close()

        def print_count(p_table_text, p_l_checked, p_l_no, p_l_wrong):
            secho(f'{p_table_text}\n'
                  f'    - 확인: {len(p_l_checked)}', fg='cyan')
            secho(f'    - 데이터없음: {len(p_l_no)}',
                  fg='cyan' if len(p_l_no) == 0 else 'yellow',
                  blink=False if len(p_l_no) == 0 else True)
            if len(p_l_no) != 0:
                echo(f'    {p_l_no}')
            secho(f'    - 문제있음: {len(p_l_wrong)}',
                  fg='cyan' if len(p_l_wrong) == 0 else 'red',
                  blink=False if len(p_l_wrong) == 0 else True)
            if len(p_l_wrong) != 0:
                echo(f'    {p_l_wrong}')

        # 데이터베이스 검색 요약 결과 출력
        echo()  # 빈 줄
        secho('<<<< DB Table 확인 결과 >>>>', fg='cyan', bold=True)
        if (self.table == 'all') or (self.table == 'sample'):
            print_count('sample Table(레코드 개수 = 1)',
                        l_checked_sample, l_no_sample, l_wrong_sample)
        if (self.table == 'all') or (self.table == 'sample_analysis'):
            print_count('sample_analysis Table(레코드 개수 = 1)',
                        l_checked_sample_analysis, l_no_sample_analysis, l_wrong_sample_analysis)
        if (self.table == 'all') or (self.table == 'sample_point'):
            print_count('sample_point Table(레코드 개수 = 30)',
                        l_checked_sample_point, l_no_sample_point, l_wrong_sample_point)
        if (self.table == 'all') or (self.table == 'sample_microbiome'):
            print_count('sample_microbiome Table(레코드 개수 = 21)',
                        l_checked_sample_microbiome, l_no_sample_microbiome, l_wrong_sample_micorbiome)
        if (self.table == 'all') or (self.table == 'analysis_version'):
            print_count('analysis_version Table(레코드 개수 = 1)',
                        l_checked_analysis_version, l_no_analysis_version, l_wrong_analysis_version)

        # 검색된 레코드 출력
        if self.print_record is True:
            pause('레코드을 확인하려면 아무키나 누르세요.')
            text = ''
            if self.table == 'all':
                l_table = ['sample', 'sample_analysis', 'analysis_version', 'sample_point', 'sample_microbiome']
            else:
                l_table = [self.table]
            for sample in d_db_data.keys():
                text += style(sample, fg='green', bold=True) + '\n'
                for table in l_table:
                    text += style(f'DB Table: {table} - 레코드: {len(d_db_data[sample][table])}개\n', fg='yellow')
                    text += '\n'.join([str(x) for x in sorted(d_db_data[sample][table])])
                    text += '\n'
                text += '\n'
            echo_via_pager(text)

        # 문제 있는 시료명 반환
        if self.table == 'all':
            return list(set(l_no_sample
                            + l_wrong_sample
                            + l_no_sample_analysis
                            + l_wrong_sample_analysis
                            + l_no_sample_point
                            + l_wrong_sample_point
                            + l_no_sample_microbiome
                            + l_wrong_sample_micorbiome
                            + l_no_analysis_version
                            + l_wrong_analysis_version))
        elif self.table == 'sample':
            return list(set(l_no_sample + l_wrong_sample))
        elif self.table == 'sample_analysis':
            return list(set(l_wrong_sample_analysis + l_no_sample_point))
        elif self.table == 'analysis_version':
            return list(set(l_no_analysis_version + l_wrong_analysis_version))
        elif self.table == 'sample_point':
            return list(set(l_no_sample_point + l_wrong_sample_point))
        elif self.table == 'sample_microbiome':
            return list(set(l_no_sample_microbiome + l_wrong_sample_micorbiome))
        else:
            raise ValueError(f'허용되지 않는 값입니다. self.table: {self.table}')

    def del_score(self):
        d_bact_data = read_yaml(ITEMS_AND_BACTERIA_FILE)
        report_version = d_bact_data['Items_and_Bacteria']['report_version']
        l_excluded_sample = self.find_score()
        secho(f'>>> 삭제 제외 시료 - {len(l_excluded_sample)}', fg='yellow', bold=True)
        secho(str(l_excluded_sample), fg='yellow')
        db_microbe = self.connect_db()
        db_microbe.start_transaction()
        secho('>>> DB 레코드 삭제 진행', bold=True)
        deleted_count = 0
        skipped_count = 0
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                if nt_sample.name in l_excluded_sample:
                    skipped_count += 1
                    continue
                else:
                    # Sample Info 파일 읽기
                    sample_info_csv = os.path.join(self.i_path.sample_path(nt_sample.name),
                                                   f'{nt_sample.name}_Info.csv')
                    d_sample_info = read_csv_to_dict(sample_info_csv, p_verbose=False)[0]
                    if (self.table == 'all') or (self.table == 'sample'):
                        db_microbe.del_sample(d_sample_info['client_name'],
                                              d_sample_info['Birth'],
                                              d_sample_info['Sex'])
                    if (self.table == 'all') or (self.table == 'sample_analysis'):
                        db_microbe.del_sample_analysis(d_sample_info['KitId'])
                    if (self.table == 'all') or (self.table == 'sample_point'):
                        db_microbe.del_sample_point(d_sample_info['KitId'], report_version)
                    if (self.table == 'all') or (self.table == 'sample_microbiome'):
                        db_microbe.del_sample_microbiome(d_sample_info['KitId'], report_version)
                    if (self.table == 'all') or (self.table == 'analysis_version'):
                        db_microbe.del_analysis_version(d_sample_info['KitId'], report_version)
                    deleted_count += 1
        db_microbe.commit()
        db_microbe.close()
        deleted_text = style(f'>>> DB 레코드 삭제 - 완료: {deleted_count}', fg='cyan')
        skipped_text = style(f'제외: {skipped_count}', fg='cyan' if skipped_count == 0 else 'yellow')
        echo(f'{deleted_text}, {skipped_text}')

    def make_info(self):
        samples_info_file = os.path.join(self.i_path.order_number_path, 'SamplesInfo.csv')
        if check_file_type(samples_info_file, 'exists'):
            pass
        else:
            secho('Warning: SampleInfo.csv 파일이 없습니다.', fg='yellow')
            secho('\t--> Sample Info를 불러옵니다.', fg='magenta')
            get_samples_info(self.order_number, samples_info_file)
        l_samples_info = read_csv_to_dict(samples_info_file)
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Summarize Taxa', ascii=True, leave=False):
                for d_info in l_samples_info:
                    if nt_sample.name == d_info['KitId']:
                        # 성별 데이터 변환
                        if d_info['Sex'] == '남성':
                            d_info['Sex'] = 'male'
                        elif d_info['Sex'] == '여성':
                            d_info['Sex'] = 'female'
                        else:
                            secho('SampleInfo.csv 파일의 성별 항목에 사용할 수 없는 값이 존재합니다.', fg='red', err=True)
                            echo(f'Sex: {d_info["Sex"]}')
                        sample_info_file = os.path.join(self.i_path.sample_path(nt_sample.name),
                                                        f'{nt_sample.name}_Info.csv')
                        write_sample_info(d_info, sample_info_file)
        secho('>>> Sample Info 생성 완료', fg='cyan')

    def run_delete(self):
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                for dir_name in tqdm(self.t_dir_name, desc='Target Dir', ascii=True, leave=False):
                    if dir_name == 'Alpha_Diversity':
                        shutil.rmtree(self.i_path.alpha_diversity_path(nt_sample.name), True)
                    elif dir_name == 'BIOM':
                        shutil.rmtree(self.i_path.biom_path(nt_sample.name), True)
                    elif dir_name == 'R_DADA2':
                        shutil.rmtree(self.i_path.r_dada2_path(nt_sample.name), True)
                    elif dir_name == 'Score':
                        shutil.rmtree(self.i_path.score_path(nt_sample.name), True)
                    elif dir_name == 'Summarized_Taxa':
                        shutil.rmtree(self.i_path.summarized_taxa_path(nt_sample.name), True)
                    elif dir_name == 'Taxonomy_Assignment':
                        shutil.rmtree(self.i_path.taxonomy_assignment_path(nt_sample.name), True)
                    else:
                        raise ValueError(f'dir_name({dir_name})이 값이 알맞지 않습니다.')

    def run_qdel(self):
        l_job_id = list()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            for nt_sample in tqdm(nt_run.samples, desc='Sample', ascii=True, leave=False):
                if self.target_dir == 'R_DADA2':
                    target_path = self.i_path.analysis_number_path
                    pattern = f'{nt_sample.name}_qsub.log'
                    job_id = parse_job_id(target_path, pattern, 'only', p_mode='last', p_verbose=False)
                    l_job_id.append(job_id)
                elif self.target_dir == 'Taxonomy_Assignment':
                    target_path = self.i_path.taxonomy_assignment_path(nt_sample.name)
                    pattern = f'*/qsub.log'
                    job_id = parse_job_id(target_path, pattern, 'many', p_mode='all', p_verbose=False)
                    l_job_id.extend(job_id)
                elif self.target_dir == 'BIOM':
                    target_path = self.i_path.biom_path(nt_sample.name)
                    pattern = 'qsub.*.log'
                    job_id = parse_job_id(target_path, pattern, 'many', p_mode='last', p_verbose=False)
                    l_job_id.extend(job_id)
                elif self.target_dir == 'Alpha_Diversity':
                    target_path = self.i_path.alpha_diversity_path(nt_sample.name)
                    pattern = 'qsub.log'
                    job_id = parse_job_id(target_path, pattern, 'many', p_mode='last', p_verbose=False)
                    l_job_id.extend(job_id)
                elif self.target_dir == 'Summarized_Taxa':
                    target_path = self.i_path.summarized_taxa_path(nt_sample.name)
                    pattern = 'qsub.*.log'
                    job_id = parse_job_id(target_path, pattern, 'many', p_mode='last', p_verbose=False)
                    l_job_id.extend(job_id)
                else:
                    raise ValueError(f'self.target_dir의 값이 허용되지 않는 값이 입력되었습니다. '
                                     f'self.target_dir: {self.target_dir}')
        secho(f'>>> Job ID 확인: {len(l_job_id)}', fg='cyan')
        echo(l_job_id)
        for ele_id in tqdm(l_job_id, desc='qdel', ascii=True):
            cmd = f'qdel {ele_id}'
            run_cmd(cmd)
        secho('>>> qDEL 작업 완료', fg='cyan')


class ASVsMerger:
    def __init__(self, kargs):
        self.analysis_base_path = kargs['analysis_base_path']
        self.order_number_file = kargs['order_number_file']
        self.out_dir = kargs['out_dir']
        self.rds_name = kargs['rds_name']
        self.dada2_dir_name = kargs['dada2_dir_name']
        self.queue = kargs['queue']
        self.no_queue = kargs['no_queue']

        # 설정
        self.l_order_analysis_number = None
        self.d_i_path = None
        self.l_nt_analysis_list = None

    def read_order_number_file(self):
        """
        ASVs 파일들을 읽어 들이기 위해서 해당 파일이 존재하는 수주번호와 분석디렉터리번호를 읽어들인다.

        파일 구조
        수주번호 (탭) 분석디렉터리번호
        ex)
        HN0012345   Analysis_1
        HN0054321   Analysis_1
        HN0054321   Analysis_2

        :return:
        """
        self.l_order_analysis_number = read_data(self.order_number_file, '\t')  # [[수주번호, 분석번호]]
        secho('>>> Order Number List 파일 읽기 완료', fg='cyan')
        order_num_count = Counter([f'{order_number}/{analysis_number}'
                                   for order_number, analysis_number in self.l_order_analysis_number])
        if max(order_num_count.values()) == 1:
            secho('>>> 수주번호 중복 확인 완료', fg='cyan')
            secho(f'수주번호 개수: {len(Counter([order_number for order_number, _ in self.l_order_analysis_number]))}')
            secho(f'분석번호 개수: {len(self.l_order_analysis_number)}', fg='yellow')
        else:
            secho('Error: 수주번호&분석번호에 중복이 있습니다.', fg='read', blink=True, err=True)
            for target, count in order_num_count.items():
                if count > 1:
                    echo(f'{target}: {count}')
            exit()

    def make_sample_list(self, p_order_number, p_analysis_number):
        l_nt_run_list = make_sample_list(
            self.analysis_base_path,
            p_order_number,
            f'{p_analysis_number}',
            None,
        )
        return l_nt_run_list

    def glob_samples(self):
        """

        l_nt_analysis_list 예시
        RunData(
            path='/garnet/Tools/Amplicon_MetaGenome/MicrobeAndMe/HN00124490/Analysis_1',
            run_info='Analysis_1',
            samples=
                [
                    SampleList(
                        path='/garnet/Tools/Amplicon_MetaGenome/MicrobeAndMe/HN00124490/Analysis_1/MBS200224NP01',
                        name='MBS200224NP01'),
                    SampleList(
                        path='/garnet/Tools/Amplicon_MetaGenome/MicrobeAndMe/HN00124490/Analysis_1/MBS200218NP002',
                        name='MBS200218NP002'),
                ]
        )

        :return:
        """
        l_nt_analysis_list = list()
        for order_number, analysis_number in self.l_order_analysis_number:
            l_sample_list_from_rawdata = self.make_sample_list(order_number, 'RawData*/*')
            l_sample_names_from_rawdata = list()
            for nt_run in l_sample_list_from_rawdata:
                for sample in nt_run.samples:
                    l_sample_names_from_rawdata.append(sample.name)
            l_sample_list_from_analysis = self.make_sample_list(order_number, analysis_number)
            count_removed_dir = 0
            nt_analysis = l_sample_list_from_analysis[0]  # 원소개수 1개!
            for index in range(len(nt_analysis.samples)-1, -1, -1):
                if nt_analysis.samples[index].name in l_sample_names_from_rawdata:
                    pass
                else:
                    nt_analysis.samples.pop(index)
                    count_removed_dir += 1
            secho(f'>>> 시료명 아닌 디렉터리 제거: {count_removed_dir}', fg='cyan')
            secho(f'시료 개수: {len(nt_analysis.samples)}', fg='yellow')
            l_nt_analysis_list.extend(l_sample_list_from_analysis)
        echo()  # 빈 줄
        secho('------ 전체 시료 확인 ------', fg='magenta')
        check_duplication(l_nt_analysis_list)
        self.l_nt_analysis_list = l_nt_analysis_list

    def check_rds_files(self):
        global ANALYSIS_DIR_NAME
        l_not_exist_rds = list()
        rds_file_count = 0
        for nt_analysis in self.l_nt_analysis_list:
            for nt_sample in nt_analysis.samples:
                rds_file = os.path.join(nt_sample.path, ANALYSIS_DIR_NAME['r_dada2'], self.rds_name)
                if check_file_type(rds_file, 'exists'):
                    rds_file_count += 1
                else:
                    l_not_exist_rds.append(rds_file)
        if len(l_not_exist_rds) != 0:
            secho(f'Error: {self.rds_name} 파일이 없는 시료가 있습니다.', fg='red', blink=True, err=True)
            secho(f'{l_not_exist_rds}', err=True)
            exit(1)
        else:
            secho(f'>>> {self.rds_name} 파일 확인: {rds_file_count}', fg='cyan')

    def copy_rds_files(self, p_dir_name) -> str:
        if check_file_type(self.out_dir, 'exists'):
            if check_file_type(self.out_dir, 'isdir'):
                pass
            else:
                secho(f'Error: {self.out_dir}은 디렉터리가 아닙니다.', fg='red', blink=True, err=True)
                secho(f'\t--out_dir 옵션에 디렉터리명(경로)을 입력하세요', fg='magenta')
        else:
            os.mkdir(self.out_dir)
            secho('>>> Out Dir 생성 완료', fg='cyan')
            echo(f'{self.out_dir}')

        all_rds_dir_path = os.path.join(self.out_dir, p_dir_name)
        os.mkdir(all_rds_dir_path)
        l_err = list()
        done_rds_file_count = 0
        for nt_analysis in tqdm(self.l_nt_analysis_list, desc='Analysis', ascii=True):
            for nt_sample in tqdm(nt_analysis.samples, desc='Sample', ascii=True, leave=False):
                ori_rds_file = os.path.join(nt_sample.path, ANALYSIS_DIR_NAME['r_dada2'], self.rds_name)
                des_rds_file = os.path.join(all_rds_dir_path, f'{nt_sample.name}_{self.rds_name}')
                if check_file_type(des_rds_file, 'exists'):
                    l_err.append(FileExistsError(f'{ori_rds_file} is existed'))
                try:
                    shutil.copy2(ori_rds_file, des_rds_file)
                except (FileNotFoundError, PermissionError) as err:
                    l_err.append(err)
                except Exception as err:
                    l_err.append(err)
                else:
                    done_rds_file_count += 1

        done_text = style(f'>>> {self.rds_name} 복사 - 완료: {done_rds_file_count}', fg='cyan')
        failed_text = style(f'실패: {len(l_err)}',
                            fg='red' if len(l_err) != 0 else '',
                            blink=True if len(l_err) != 0 else False)
        echo(f'{done_text}, {failed_text}')
        if len(l_err) != 0:
            [echo(err) for err in l_err]
            exit()
        return all_rds_dir_path

    def run_merge_asvs(self, p_rds_dir):
        run_merge_asvs(p_rds_dir, self.rds_name, self.queue, self.no_queue)


class FilesCollector(ASVsMerger):
    def __init__(self, kargs):
        self.analysis_base_path = kargs['analysis_base_path']
        self.order_number_file = kargs['order_number_file']
        self.out_dir = kargs['out_dir']
        self.rds_name = kargs['file_name']
        self.dada2_dir_name = None

        # 설정
        self.l_order_analysis_number = None
        self.d_i_path = None
        self.l_nt_analysis_list = None

    def check_files(self):
        self.check_rds_files()

    def copy_files(self, p_dir_name) -> str:
        dir_path = self.copy_rds_files(p_dir_name)
        return dir_path


