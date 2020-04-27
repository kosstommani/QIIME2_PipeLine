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
from click import echo, secho, style, progressbar
from CoFI.core import make_sample_list, make_analysis_dir
from CoFI.data_structure import RunData
from CoFI.data_handler import check_duplication
from SpoON.util import parse_config, run_cmd, check_run_cmd, glob_dir, \
    read_data, check_file_type, check_jobs_in_queue, launcher_cmd
from SpoON.run_time import compute_run_time
from theCups.data_structure import Path
from MicrobeAndMe.dada2 import run_dada2_r, run_merge_asvs
from collections import Counter, defaultdict
from tqdm import tqdm


# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']
ANALYSIS_DIR_NAME = CONFIG['Analysis_Dir_Name']


class PathMAM(Path):
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


class Microbe:
    def __init__(self, kargs: dict, mode, analysis_type):
        """
        Microbe&Me 서비스를 위한 분석을 진행합니다.

        :param kargs:
        :param mode: [all | r_dada2 | taxonomy | biom] 분석 명령어
        :param analysis_type: [NGS | MicrobeAndMe]
        """
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
        else:
            raise ValueError(f'Microbe mode 설정 오류: {mode}')
        self.analysis_type = analysis_type

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
                    queue_mode:
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
        self.queue_mode = kargs['queue_mode']
        self.slots = kargs['slots']
        self.no_queue = kargs['no_queue']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.database = kargs['database']

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
                    queue_mode:
                    slots:
                    no_queue:
                    order_number:
                    include:
                    exclude:
        :return:
        """
        self.analysis_number = None
        self.analysis_base_path = kargs['analysis_base_path']
        self.analysis_dir_base_name = kargs['analysis_dir_base_name']
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.queue = kargs['queue']
        self.queue_mode = kargs['queue_mode']
        self.slots = kargs['slots']
        self.no_queue = kargs['no_queue']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        
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
        self.cofi_target_dir_suffix = kargs['cofi_target_dir_suffix']
        self.r1_suffix = kargs['r1_suffix']
        self.r2_suffix = kargs['r2_suffix']
        self.order_number = kargs['order_number']
        self.include = kargs['include']
        self.exclude = kargs['exclude']
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
        self.include = kargs['include']
        self.exclude = kargs['exclude']
        self.database = kargs['database']

        # 설정
        self.i_path = None
        self.l_nt_run_list = None

    @property
    def max_slots(self):
        if self.queue == 'bi3m.q':
            return 50
        elif self.queue == 'meta.q':
            return 40

    @property
    def l_db(self):
        if self.database == 'all':
            return ['NCBI_16S', 'NCBI_Probiotics']
        else:
            return [self.database]

    def _set_slots(self):
        """
        Queue 당 사용할 수 있는 slots 의 수를 지정한다.
        bi3m.q : core 56개 - cm456(Intel(R) Xeon(R) Gold 6132 CPU @ 2.60GHz) -
        meta.q : core 48개 - denovo06(Intel(R) Xeon(R) CPU E7540  @ 2.00GHz)

                   Passmark CPU Mark    CPU Value(CPUmark/$Price)    Price
        Gold 6123:            27,522	                     6.34	 $4,339.90(Total for 2 CPUs (2020-03-18))
        E7540    :             7,946                       124.38    $63.88(2018-12-31)

        :return: 
        """
        if self.queue_mode == 'auto':
            if self.slots is None:
                self.slots = 2
        elif self.queue_mode == 'fix':
            if self.queue == 'bi3m.q':
                if self.slots is None:
                    self.slots = self.max_slots
            elif self.queue == 'meta.q':
                if self.slots is None:
                    self.slots = self.max_slots
            else:
                if self.no_queue is False:
                    secho('Error: Queue 지정이 잘못되었습니다.', fg='red', blink=True, err=True)
                    echo(f'Queue: {self.queue}', err=True)
                    exit(1)
        else:
            secho('Error: --queue_mode 옵션이 잘못되었습니다.', fg='red', blink=True, err=True)
            echo(f'--queue_mode: {self.queue_mode}')
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
        elif self.analysis_type == 'MicrboeAndMe':
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
            f'--target_dir_suffix {self.prema_target_dir_suffix} ' \
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

        :param p_mode: ['MicrobeAndMe', 'NGS]
        :return:
        """
        if self.analysis_type == 'MicrobeAndMe':
            if len(self.l_nt_run_list) > 1:
                secho('Error: Run Dir 의 개수가 1이 아닙니다.', fg='red', blink=True, err=True)
                echo([nt_run.path for nt_run in self.l_nt_run_list], err=True)
                exit(1)
        self.set_queue_mode(self.sample_count)
        d_return_code = dict()
        for nt_run in tqdm(self.l_nt_run_list, desc='Run Dir', ascii=True):
            run_dir = nt_run.path.split(self.order_number)[1]
            for nt_sample in tqdm(nt_run.samples, desc='R_DADA2', ascii=True, leave=False):
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
                    }
                )
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

    def make_biom(self, p_db):
        l_run_process = list()
        base_biom_status = False
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
                if check_file_type(i_biom.base_biom_file, 'exists') is True:
                    base_biom_status = True
                else:
                    l_run_process.append(i_biom.convert_biom())  # Base BIOM
                l_run_process.append(i_biom.add_metadata())  # Taxonomy BIOM
                l_run_process.append(i_biom.convert_taxonomy_table())
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

    def make_biom_by_db(self):
        for db in self.l_db:
            self.make_biom(db)

    def run_alpha_diversity(self):
        l_run_process = list()
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
                run = run_cmd(cmd)
                l_run_process.append(run)
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

    def summrize_taxa(self):
        l_run_process = list()
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
                        f'-L 6,7'
                    run = run_cmd(cmd)
                    l_run_process.append(run)
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

    # def set_path(self):
    # self.d_i_path = defaultdict(list)
    # if self.l_order_analysis_number is None:
    #     raise ValueError('self.l_order_number 에 값이 없습니다. ASVsMerger.read_order_number_file() 를 실행하세요.')
    # for order_number, analysis_number in self.l_order_analysis_number:
    #     i_path = PathMAM(self.analysis_base_path, order_number, analysis_number)
    #     self.d_i_path[order_number].append(i_path)

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
            l_nt_analysis_list.extend(self.make_sample_list(order_number, analysis_number))
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

    def copy_rds_files(self, p_dir_name):
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

        done_text = style(f'{done_rds_file_count}완료', fg='cyan')
        failed_text = style(f'{len(l_err)}실패',
                            fg='red' if len(l_err) != 0 else '',
                            blink=True if len(l_err) != 0 else False)
        secho(f'>>> {self.rds_name} 복사: {done_text}, {failed_text}')
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

    def copy_files(self, p_dir_name):
        self.copy_rds_files(p_dir_name)


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
        # self.metadata = kargs['metadata']
        self.asvs_table_name = 'ASVs.tsv'
        self.biom_name = 'ASVs.biom'
        self.biom_table_name = 'ASVs_table.tsv'
        self.taxonomy_name = 'otus_rep_tax_assignments.txt'

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
        return os.path.join(self.biom_path, f'{splitext[0]}.{self.db_name}{splitext[1]}')

    @property
    def biom_table_file(self):
        return os.path.join(self.biom.path, self.biom_table_name)

    @property
    def taxonomy_biom_table_file(self):
        splitext = os.path.splitext(self.biom_table_name)
        return os.path.join(self.biom_path, f'{splitext[0]}.{self.db_name}{splitext[1]}')

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

    def convert_biom(self) -> object:
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
        run = self.run_cmd(biom_convert_cmd, true_meg, false_meg)
        return run

    def add_metadata(self) -> object:
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
            '--observation-header ASVsID,taxonomy,confidence,count ' \
            '--sc-separated taxonomy ' \
            '--output-as-json'
        self.write_recipe(biom_add_metadata_cmd)
        run = self.run_cmd(biom_add_metadata_cmd, true_meg, false_meg)
        return run

    def convert_base_table(self) -> object:
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
        run = self.run_cmd(biom_table_cmd, true_meg, false_meg)
        return run

    def convert_taxonomy_table(self) -> object:
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
        run = self.run_cmd(biom_table_cmd, true_meg, false_meg)
        return run
