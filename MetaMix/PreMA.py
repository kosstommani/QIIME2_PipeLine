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
#  _____         _____ _____
# |  _  |___ ___|     |  _  |
# |   __|  _| -_| | | |     |
# |__|  |_| |___|_|_|_|__|__|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------


__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.2'

import click
import os
from pprint import pprint
from time import time
from SpoON.util import parse_config, delete_dir, run_dir_tree, check_file_type, \
                       check_order_number_system
from SpoON.run_time import print_prema_run_time
from PreMA.core import run_core, check_custom_name


# 기본값 설정
CONFIG = parse_config()
RAWDATA_BASE_PATH = CONFIG['base_path']['rawdata_base_path']
ANALYSIS_BASE_PATH = CONFIG['base_path']['analysis_base_path']
TARGET_DIR_SUFFIX = CONFIG['PreMA_CORE']['target_dir_suffix']
R1_SUFFIX = CONFIG['PreMA_CORE']['R1_suffix']
R2_SUFFIX = CONFIG['PreMA_CORE']['R2_suffix']
SAMPLE_READ = CONFIG['PreMA_CORE']['sample_read']
Q30 = CONFIG['PreMA_CORE']['q30']
N_BASE = CONFIG['PreMA_CORE']['n_base']
N_READ = CONFIG['PreMA_CORE']['n_read']
ADAPTER_TRIM_TOOL = CONFIG['PreMA_CORE']['adapter_trim_tool']
INDEX_KIT = CONFIG['PreMA_CORE']['index_kit']
SAMPLE_NAME_MODE = CONFIG['PreMA_CORE']['sample_name_mode']
TREE_LEVEL = CONFIG['PreMA_DIR']['tree_level']
# 분석 디렉터리
ANALYSIS_TARGET_DIR_SUFFIX = CONFIG['CoFI']['target_dir_suffix']
ANALYSIS_R1_SUFFIX = CONFIG['CoFI']['R1_suffix']
ANALYSIS_R2_SUFFIX = CONFIG['CoFI']['R2_suffix']

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def main(**kargs):
    """
    \b
    888b     d888          888             888b     d888 d8b
    8888b   d8888          888             8888b   d8888 Y8P
    88888b.d88888          888             88888b.d88888
    888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
    888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
    888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
    888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
    888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
     _____         _____ _____
    |  _  |___ ___|     |  _  |
    |   __|  _| -_| | | |     |
    |__|  |_| |___|_|_|_|__|__|
    \b
    Amplicon MetaGenome Anlaysis
    분석 파이프라인
    MetaMix - PreMA
    """
    pass


@main.command('CORE', short_help='OTU분석을 위한 준비 작업')
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--rawdata_base_path', '-rp',
              default=RAWDATA_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='Raw Data 기본 경로')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='분석 기본 경로')
@click.option('--target_dir_suffix', '-t',
              default=TARGET_DIR_SUFFIX,
              show_default=True,
              help='Run 디렉터리 검색조건')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True,
              help='Read1 접미사')
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True,
              help='Read2 접미사')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_copy_fastqc', '-nqc',
              is_flag=True,
              help='FASTQC 데이터를 복사하지 않음.')
@click.option('--sample_read', '-sr',
              default=SAMPLE_READ,
              show_default=True,
              help='시료별 생산될 Read의 개수.')
@click.option('--q30', '-q30',
              default=Q30,
              show_default=True,
              nargs=2, type=click.IntRange(0, 100),
              help='Q30 기준값')
@click.option('--n_base', '-nb',
              default=N_BASE,
              show_default=True,
              nargs=2, type=click.FloatRange(0, 1),
              help='n base의 비율(%)')
@click.option('--n_read', '-nr',
              default=N_READ,
              show_default=True,
              help='')
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']),
              help='')
@click.option('--trim_tail',
              default=None,
              show_default=True,
              type=click.INT,
              help='fastp 전용. Read의 끝에서부터 제거할 base의 개수. '
                   'ex) 301bp를 251bp를 만들고자 할 경우 --trim_tail 50')
@click.option('--index_kit', '-ik',
              default=INDEX_KIT,
              show_default=True,
              type=click.Choice(['auto', 'Nextera', 'TruSeq']),
              help='')
@click.option('--sample_name_mode', '-nm',
              default=SAMPLE_NAME_MODE,
              show_default=True,
              type=click.Choice(['all', 'all_no', 'no_zero', 'no_mid']),
              help='시료명 변경에 관한 모드 설정. '
                   'all: zero, middle zero 모두 추가. '
                   'all_no: no_zero, no_mid 모두 적용. '
                   'no_zero - 문자형 숫자에 0을 추가하지 않음. '
                   'no_mid - "."를 구분자로 하는 문자열에서 구분자를 '
                   '기준으로 문자형 숫자가 있는 경우 0을 추가하지 않음.'.format(SAMPLE_NAME_MODE))
@click.option('--custom_sample_name', '-csn',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='변경할 시료명과 그룹 정보가 담긴 파일명. 각 행이 "#"으로 시작하면 주석 처리. '
                   '구조: #SampleID   New_name   Group1   Group2. 자세한 방법은 개발자에게 문의.')
@click.option('--no_multiqc',
              is_flag=True,
              help='MultiQC 실햄 안함. 시료의 개수가 많을 경우 수행 시간이 길어짐.')
@click.option('--my_story',
              is_flag=True,
              help='myBiomeStory - fastp QC 기준 다름.')
@click.option('--microbe_and_me',
              is_flag=True,
              help='MicrobeAndMe - fastp QC 기준 다름.')
def core(**kargs):
    """
    \b
    OTU분석을 진행하기 위한 준비 작업을 진행합니다.
    시료명 변경, STAT Check, QC Report 생성, FASTQ Copy, Adapter Trim
    """
    # 수주번호 형식 확인
    start_time = time()
    if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
        pass

    pprint(kargs)
    if all([kargs['my_story'], kargs['microbe_and_me']]):
        click.secho('Option Error: --my_stroy 와 --microbe_and_me 옵션 모두 적용 불가', fg='red')
    if check_file_type(os.path.join(kargs['rawdata_base_path'], kargs['order_number']), 'exists'):
        main_run_time, integ_run_time, multiqc_run_time = run_core(kargs)
        end_time = time()
        print_prema_run_time(start_time, main_run_time, integ_run_time, multiqc_run_time, end_time)
    else:
        click.secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')


@main.command('DIR', short_help='디렉터리 구조 확인 및 삭제')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='')
@click.option('--level', '-l',
              default=TREE_LEVEL,
              show_default=True,
              type=int,
              help='목표 디렉터리의 출력 깊이.'.format(TREE_LEVEL))
def dir_del(**kargs):
    """
    목표 위치의 디렉터리 구조 확인 및 삭제
    """
    # 수주번호 형식 확인
    # if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
    #     pass

    if check_file_type(os.path.join(*(kargs['analysis_base_path'], kargs['order_number'])), 'exists'):
        run_dir_tree(kargs['analysis_base_path'], kargs['order_number'], kargs['level'])
        previous_level = kargs['level']
        cmd_answer = None
        while cmd_answer != 'Q':
            cmd_answer = click.prompt('L: 출력 깊이 설정 / D: 삭제 / Q: 나가기\n명령어?', type=click.Choice(['L', 'D', 'Q']))
            if cmd_answer == 'D':
                answer = click.confirm('삭제하시겠습니까?')
                if answer:
                    answer2 = click.confirm(click.style('진짜? 정말?', fg='yellow', bold=True, underline=True))
                    if answer2:
                        delete_dir(kargs['analysis_base_path'], kargs['order_number'])
                        exit()
                else:
                    exit()
            elif cmd_answer == 'L':
                level = click.prompt('출력깊이? (이전값: {previous})'.format(previous=previous_level), type=click.INT)
                run_dir_tree(kargs['analysis_base_path'], kargs['order_number'], level)
                previous_level = level
    else:
        click.secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')


@main.command('CHECK_NAME', short_help='시료명 변경 오류 확인')
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--custom_sample_name', '-csn',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='변경할 시료명이 담긴 파일의 위치. 각 행이 "#"으로 시작하면 주석 처리. '
                   '구조: 현재시료명[공백 또는 탭]변경할 시료명. ex) TEST    Amp')
@click.option('--sample_name_mode', '-nm',
              default=SAMPLE_NAME_MODE,
              show_default=True,
              type=click.Choice(['all', 'all_no', 'no_zero', 'no_mid']),
              help='시료명 변경에 관한 모드 설정. '
                   'all: zero, middle zero 모두 추가. '
                   'all_no: no_zero, no_mid 모두 적용. '
                   'no_zero - 문자형 숫자에 0을 추가하지 않음. '
                   'no_mid - "."를 구분자로 하는 문자열에서 구분자를 '
                   '기준으로 문자형 숫자가 있는 경우 0을 추가하지 않음.'.format(SAMPLE_NAME_MODE))
@click.option('--rawdata-base-path', '-rp',
              default=RAWDATA_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='')
@click.option('--target_dir_suffix', '-t',
              default=TARGET_DIR_SUFFIX,
              show_default=True,
              help='')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True,
              help='')
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True,
              help='')
def check_name(**kargs):
    """
    \b
    시료명 변경 파일의 오류 여부를 확인합니다.
    시료명 중복 여부, 누락 여부, 최종 예상 시료명 출력.
    """
    # 수주번호 형식 확인
    if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
        pass

    check_custom_name(kargs)


@main.command('TRIM', short_help='특정길이로 Length Trim')
@click.argument('order_number')
@click.option('--trim_length', '-tl',
              required=True,
              type=click.INT,
              help='원하는 Read의 길이. ex) Read의 길이가 301bp 일 경우, --trim_length 250 이면 Read1 Read2 의 길이가 '
                   '각각 250bp임.')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='')
@click.option('--target_dir_suffix', '-t',
              default=ANALYSIS_TARGET_DIR_SUFFIX,
              show_default=True,
              help='')
@click.option('--sample_suffix', '-s',
              help='검출하고자 하는 시료의 접미사. ex) V3.4 ITS3.4')
@click.option('--R1_suffix', '-r1',
              default=ANALYSIS_R1_SUFFIX,
              show_default=True,
              help='')
@click.option('--R2_suffix', '-r2',
              default=ANALYSIS_R2_SUFFIX,
              show_default=True,
              help='')
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']),
              help='')
def length_trim(**kargs):
    """
    Adapter Trimming 이후에 Length Trimming을 진행할 경우 사용할 수 있다.

    \b
    MetaMix - Length Trimming 2가지 방법
      1. PreMA CORE --trim_tail 사용(fastp 옵션 이용).
         제거할 base의 길이 입력.(제거전 길이에 대해서 상대적으로 적용.)
      2. PreMA TRIM command 사용(in house code).
         최종 길이 입력.(절대적)
    """
    start_time = time()
    from CoFI.core import make_sample_list
    from PreMA.core import run_length_trim
    from SpoON.fastq_handler import exist_fastq
    from SpoON.run_time import print_legnth_trim_time
    # 수주번호 형식 확인
    if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
        pass
    pprint(kargs)
    if check_file_type(os.path.join(kargs['analysis_base_path'], kargs['order_number']), 'exists'):
        l_nt_run_list = make_sample_list(
            kargs['analysis_base_path'],
            kargs['order_number'],
            kargs['target_dir_suffix'],
            kargs['sample_suffix'])
        exist_fastq(l_nt_run_list, kargs)
    else:
        click.secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')
    r1_suffix = kargs['r1_suffix'].format(trim=kargs['adapter_trim_tool'])
    r2_suffix = kargs['r2_suffix'].format(trim=kargs['adapter_trim_tool'])
    run_length_trim(l_nt_run_list, kargs['trim_length'], r1_suffix, r2_suffix)
    end_time = time()
    print_legnth_trim_time(start_time, end_time)


@main.command('MAKE', short_help='분석 디렉터리 구조 생성(작성 예정)')
def make(**kargs):
    """
    고객으로부터 전달받은 RawData(fastq)를 분석에 적합한 구조 생성 및 QC(Adapter Trim 등)
    """
    pass


if __name__ == '__main__':
    main()
