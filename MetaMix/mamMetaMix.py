#!/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
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

import click
from pprint import pprint
from SpoON.util import check_order_number_system, parse_config


# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']
RAWDATA_BASE_PATH = CONFIG['base_path']['rawdata_base_path']
ANALYSIS_BASE_PATH = CONFIG['mamMetaMix']['analysis_base_path']
PREMA_TARGET_DIR_SUFFIX = CONFIG['PreMA_CORE']['target_dir_suffix']
R1_SUFFIX = CONFIG['CoFI']['R1_suffix']
R2_SUFFIX = CONFIG['CoFI']['R2_suffix']
SAMPLE_READ = 20000  # DADA2 이후 10,000개 Reads
Q30 = CONFIG['PreMA_CORE']['q30']
N_BASE = CONFIG['PreMA_CORE']['n_base']
N_READ = CONFIG['PreMA_CORE']['n_read']
ADAPTER_TRIM_TOOL = CONFIG['PreMA_CORE']['adapter_trim_tool']
INDEX_KIT = CONFIG['PreMA_CORE']['index_kit']
SAMPLE_NAME_MODE = CONFIG['PreMA_CORE']['sample_name_mode']
COFI_TARGET_DIR_SUFFIX = CONFIG['CoFI']['target_dir_suffix']
ANALYSIS_DIR_BASE_NAME = CONFIG['CoFI']['analysis_dir_base_name']
POOL_WORKER = 10
# 옵션 기본값
INDEX_KIT_LIST = ['auto', 'Nextera', 'TruSeq']
DATABASE_LIST = ['all', 'NCBI_16S', 'NCBI_Probiotics21']
QUEUE_LIST = CONFIG['Queue_List']
DEFAULT_QUEUE = CONFIG['default_queue']
DB_TABLE = CONFIG['mamMetaMix']['db_table']
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def main(**kargs):
    """
    \b
                              ___  ___     _       ___  ____
                              |  \\/  |    | |      |  \\/  (_)
     _ __ ___   __ _ _ __ ___ | .  . | ___| |_ __ _| .  . |___  __
    | '_ ` _ \\ / _` | '_ ` _ \\| |\\/| |/ _ \\ __/ _` | |\\/| | \\ \\/ /
    | | | | | | (_| | | | | | | |  | |  __/ || (_| | |  | | |>  <
    |_| |_| |_|\\__,_|_| |_| |_\\_|  |_/\\___|\\__\\__,_\\_|  |_/_/_/\\_\\

    작성자 : 박정원(JungWon Park, KOSST)
    
    \b
    Microbe&Me 서비스 파이프라인
    mamMetaMix - ALL
        ├─── PreMA CORE
        ├─── R_DADA2
        ├─── TAXONOMY
        ├─── BIOM
        ├  ─ ALPHA_DIVERSITY
        ├  ─ SUMMARIZE_TAXA
        ├─── SCORE
        ├─── INFO
        ├─── DB_INSERT
        └─── DB_FIND
    \b
    DADA2 Utility
        └─── COLLECT
        └─── MERGE_ASVs
    \b
    Utility
        └─── DEL
        └─── qDEL
        └─── DB_FIND
        └─── DB_DEL
    """
    pass


@main.command('ALL', short_help='파이프라인')
@click.argument('order_number')
@click.option('--rawdata_base_path', '-rp',
              default=RAWDATA_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='Raw Data 기본 경로')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--prema_target_dir_suffix', '-pt',
              default=PREMA_TARGET_DIR_SUFFIX,
              show_default=True,
              help='PreMA의 target_dir_suffix. Raw Data(Demultiplexing)')
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--sample_read', '-sr',
              default=SAMPLE_READ,
              show_default=True,
              help='시료별 생산될 Read의 개수.')
@click.option('--index_kit', '-ik',
              default=INDEX_KIT,
              show_default=True,
              type=click.Choice(INDEX_KIT_LIST))
@click.option('--no_prema',
              is_flag=True,
              help='PreMA 실행 안함.')
@click.option('--database', '-db',
              default='all',
              show_default=True,
              type=click.Choice(DATABASE_LIST))
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--slots', '-s',
              type=click.INT,
              help='Queue 작업시 사용할 슬롯(CPU)개수. '
                   'bi3m.q: 56, meta.q: 144, meta.q@denovo06:48, meta.q@denovo11: 144. '
                   '사용 가능한 CPU 개수 초과시 최대 개수 할당.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
@click.option('--table', '-t',
              default='all',
              show_default=True,
              type=click.Choice(DB_TABLE),
              help='검색할 데이터베이스 테이블 선택.')
@click.option('--print_record', '-p',
              is_flag=True,
              help='DB에 등록된 레코드 출력.')
def all(**kargs):
    from MicrobeAndMe.core import all_pipeline
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    pprint(kargs)
    all_pipeline(kargs)


@main.command('R_DADA2', short_help='R DADA2')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--slots', '-s',
              type=click.INT,
              help='Queue 작업시 사용할 슬롯(CPU)개수. '
                   'bi3m.q: 56, meta.q: 144, meta.q@denovo06:48, meta.q@denovo11: 144. '
                   '사용 가능한 CPU 개수 초과시 최대 개수 할당.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def r_dada2(**kargs):
    from MicrobeAndMe.core import dada2_using_r
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    dada2_using_r(kargs)


@main.command('MERGE_ASVs', short_help='ASVs 통합')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--order_number_file', '-of',
              required=True,
              help='ASVs 파일들을 통합할 수주들의 목록. '
                   '화이트스페이스 기준으로 구분. '
                   'ex) HN00124490	Analysis_1')
@click.option('--dada2_dir_name', '-dn',
              default='R_DADA2',
              show_default=True,
              help='DADA2 디렉터리 이름.')
@click.option('--out_dir', '-o',
              required=True,
              help='결과를 저장할 디렉터리명. 해당 디렉터리 경로에 --dada2_dir_name 옵션의 디렉터리 생성.'
                   '디렉터리 존재 유무 확인 후 없으면 생성.')
@click.option('--rds_name', '-rn',
              default='seqtab_noChimera.rds',
              show_default=True,
              type=click.Choice(['seqtab_noChimera.rds', 'seqtab.rds']),
              help='ASVs 정보를 저장한 rds 파일의 이름.')
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
def merge(**kargs):
    from MicrobeAndMe.core import merge_asvs_r
    merge_asvs_r(kargs)


@main.command('COLLECT', short_help='여러 수주들의 R_DADA2 디렉터리 파일 복사')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--order_number_file', '-of',
              required=True,
              help='ASVs 파일들을 통합할 수주들의 목록. 화이트스페이스 기준으로 구분.')
@click.option('--out_dir', '-o',
              required=True,
              help='결과를 저장할 디렉터리명. 해당 디렉터리 경로에 --dada2_dir_name 옵션의 디렉터리 생성.'
                   '디렉터리 존재 유무 확인 후 없으면 생성.')
@click.option('--file_name', '-fn',
              help=' 파일의 이름.')
def collect_files(**kargs):
    from MicrobeAndMe.core import collect_files
    collect_files(kargs)


@main.command('TAXONOMY', short_help='Taxnonomy: BLAST')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--database', '-db',
              default='all',
              show_default=True,
              type=click.Choice(DATABASE_LIST))
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def taxonomy(**kargs):
    from MicrobeAndMe.core import taxonomy
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    taxonomy(kargs)


@main.command('BIOM', short_help='BIOM 파일 생성')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--database', '-db',
              default='all',
              show_default=True,
              type=click.Choice(DATABASE_LIST))
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def biom(**kargs):
    from MicrobeAndMe.core import biom
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    biom(kargs)


@main.command('ALPHA_DIVERSITY', short_help='Alpha Diversity 계산(Shannon, Simpson)')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def alpha_diversity(**kargs):
    from MicrobeAndMe.core import alpha_diversity
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    alpha_diversity(kargs)


@main.command('SUMMARIZE_TAXA', short_help='Level별 Taxonomy Abundance 파일 생성')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--database', '-db',
              default='all',
              show_default=True,
              type=click.Choice(DATABASE_LIST))
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='DADA2(R)를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def summarize_taxa(**kargs):
    from MicrobeAndMe.core import summarize_taxa
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    summarize_taxa(kargs)


@main.command('SCORE', short_help='Score 생성')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def score(**kargs):
    from MicrobeAndMe.core import score
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    score(kargs)


@main.command('DB_INSERT', short_help='DB 삽입')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def db_insert(**kargs):
    from MicrobeAndMe.core import db_insert
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    db_insert(kargs)


@main.command('DB_FIND', short_help='DB 검색')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
@click.option('--table', '-t',
              default='all',
              show_default=True,
              type=click.Choice(DB_TABLE),
              help='검색할 데이터베이스 테이블 선택.')
@click.option('--print_record', '-p',
              is_flag=True,
              help='DB에 등록된 레코드 출력.')
def db_find(**kargs):
    from MicrobeAndMe.core import db_find
    # 수주번호 형식 확인                                                                                       `
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    db_find(kargs)
    

@main.command('DB_DEL', short_help='DB 삽입')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
@click.option('--table', '-t',
              default='all',
              show_default=True,
              type=click.Choice(DB_TABLE),
              help='삭제할 데이터베이스 테이블 선택.')
@click.option('--print_record', '-p',
              is_flag=True,
              help='DB에 등록된 레코드 출력.')
def db_del(**kargs):
    from MicrobeAndMe.core import db_del
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    db_del(kargs)


@main.command('INFO', short_help='고객정보 파일 생성')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def info(**kargs):
    from MicrobeAndMe.core import info
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    info(kargs)


@main.command('DEL', short_help='분석 디렉터리 삭제')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--dir', '-d',
              type=click.Choice(['Alpha_Diversity', 'BIOM', 'R_DADA2', 'Score',
                                 'Summarized_Taxa', 'Taxonomy_Assignment']),
              multiple=True,
              required=True,
              help='삭제할 디렉터리. 다중선택 가능. ex) -d BIOM -d Score')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def delete(**kargs):
    from MicrobeAndMe.core import delete
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    delete(kargs)


@main.command('qDEL', short_help='Queue Job 삭제')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--cofi_target_dir_suffix', '-pt',
              default=COFI_TARGET_DIR_SUFFIX,
              show_default=True,
              help='CoFI의 target_dir_suffix. Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--target', '-t',
              type=click.Choice(['R_DADA2', 'Taxonomy_Assignment', 'BIOM', 'Alpha_Diversity', 'Summarized_Taxa']),
              help='Queue에 등록된 작업번호 중 선택된 목표 명령어에 의해 생성된 작업번호를 삭제')
@click.option('--include',
              help='분석에 포함할 시료명.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명.'
                   'ex) --exclude "D E F"')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def qdel(**kargs):
    from MicrobeAndMe.core import qdel
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    qdel(kargs)


if __name__ == '__main__':
    main()
