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
#  _____     _____ _____
# |     |___|   __|     |
# |   --| . |   __|-   -|
# |_____|___|__|  |_____|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.1.8'

# ------------------------------------
# Ver. 1.1.6  2020.04.28
# Taxonomy - DB 다중 선택 기능
# ------------------------------------
# Ver. 1.1.7 2020.04.29
# Taxonomy - remove keep : Archaea 추가
# ------------------------------------
# Ver. 1.1.8 2020.05.19
# R_DADA2 추가


import click
import os
from pprint import pprint
from time import time
from humanize import intcomma
from SpoON.util import parse_config, check_file_type, check_order_number_system


# 기본값 설정
CONFIG = parse_config()
# pprint(CONFIG)
ANALYSIS_BASE_PATH = CONFIG['base_path']['analysis_base_path']
TARGET_DIR_SUFFIX = CONFIG['CoFI']['target_dir_suffix']
ANALYSIS_DIR_BASE_NAME = CONFIG['CoFI']['analysis_dir_base_name']
R1_SUFFIX = CONFIG['CoFI']['R1_suffix']
R2_SUFFIX = CONFIG['CoFI']['R2_suffix']
ADAPTER_TRIM_TOOL = CONFIG['CoFI']['adapter_trim_tool']
READ_LENGTH = CONFIG['CoFI']['read_length']
# FLASH
PERCENT_COMBINED_ALARM = CONFIG['CoFI_FLASH']['percent_combined_alarm']
MIN_OVERLAP = CONFIG['CoFI_FLASH']['min_overlap']
MAX_OVERLAP = CONFIG['CoFI_FLASH']['max_overlap']
LENGTH_FILTER = CONFIG['CoFI_FLASH']['length_filter']
MODIFY_HEADER = CONFIG['CoFI_FLASH']['modify_header']
POOLED_SAMPLE_FILE = CONFIG['CoFI_FLASH']['pooled_sample_file']
# OTU - CD-HIT-OTU
OTU_CUTOFF = CONFIG['CoFI_OTU']['otu_cutoff']
TRIM_CUTOFF = CONFIG['CoFI_OTU']['trim_cutoff']
CHECK_CHIMERA = CONFIG['CoFI_OTU']['check_chimera']
PCR_ERROR = CONFIG['CoFI_OTU']['pcr_error']
# OTU - Closed-Reference OTU Picking
CLOSED_JOBS_TO_START = CONFIG['CoFI_OTU']['closed_otu']['jobs_to_start']
ASSEMBLY_SUFFIX = CONFIG['CoFI_OTU']['closed_otu']['assembly_suffix']
DATABASE_LIST = list(CONFIG['CoFI_OTU']['closed_otu']['database'].keys())
# TAXONOMY
BLAST_DATABASE_LIST = list(CONFIG['CoFI_TAXONOMY']['blast_db'].keys())
UCLUST_DATABASE_LIST = list(CONFIG['CoFI_TAXONOMY']['uclust_db'].keys())
BLAST_READ_PER_JOB = CONFIG['CoFI_TAXONOMY']['queue']['read_per_job']
BLAST_NT_MAX_JOB = CONFIG['CoFI_TAXONOMY']['queue']['nt_max_job']
BLAST_QUEUE = CONFIG['CoFI_TAXONOMY']['queue']['name']

# DADA2
RANDOM_READ = CONFIG['CoFI_DADA2']['random_read']
P_TRUNC_LEN_F = CONFIG['CoFI_DADA2']['p_trunc_len_f']
P_TRUNC_LEN_R = CONFIG['CoFI_DADA2']['p_trunc_len_r']
P_TRIM_LEFT_F = CONFIG['CoFI_DADA2']['p_trim_left_f']
P_TRIM_LEFT_R = CONFIG['CoFI_DADA2']['p_trim_left_r']

QUEUE_LIST = CONFIG['Queue_List']
DEFAULT_QUEUE = CONFIG['default_queue']
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
     _____     _____ _____
    |     |___|   __|     |
    |   --| . |   __|-   -|
    |_____|___|__|  |_____|

    \b
    Amplicon MetaGenome Anlaysis
    분석 파이프라인
    MetaMix - CoFI

    \b
    \033[1;32;40m<OTU PipeLine>\033[m
    FLASH_HIT - FLASH(Pooled & FASTQ)
              - CD-HIT-OTU
    FLASH(Pooled & FASTA) ── CLOSED

    \b
    \033[1;32;40m<ASVs PipeLine>\033[m
    DADA2
     ├─── Q2_DADA2
     └─── R_DADA2

    \b
    \033[1;32;40m<분석 순서>\033[m
    FLASH_HIT(FLASH ── CD-HIT-OTU)
       ├  ─  QIIME1 : ALIGNMENT  ── PHYLOGENY
       ├  ─  QIIME2 : MAFFT_FastTree ── DIVERSITY
       ├─── TAXONOMY
       └─── BIOM
    \b
    Q2_DADA2
       ├─── MAFFT_FastTree ── DIVERSITY
       └─── TAXONOMY
    \b
    R_DADA2

    """
    pass


@main.command('INFO', short_help='info.txt 파일 생성')
@click.argument('order_number')
@click.option('--output_path', '-o',
              type=click.Path(exists=True, dir_okay=True))
@click.option('--sample_count', '-s',
              type=click.INT,
              help='분석에 사용된 시료의 개수')
@click.option('--target_region', '-r',
              type=click.Choice(['Bakt_341F-805R', 'ITS_3F-4R', '27F-Bact338']),
              help='')
def info(**kargs):
    from SpoON.util import make_info_file
    if kargs['target_region'] == 'Bakt_341F-805R':
        target_region = 'V3V4(Bakt_341F-805R)'
    elif kargs['target_region'] == 'ITS_3F-4R':
        target_region = 'ITS(ITS_3F-4R)'
    elif kargs['target_region'] == '27F-Bact338':
        target_region = 'V1V2(27F-Bact338)'
    else:
        target_region = 'None'
    make_info_file(kargs['order_number'], kargs['output_path'], kargs['sample_count'], target_region)


@main.command('PRIMER', short_help='Primer 정보 출력')
@click.option('--region', '-r',
              default='all',
              show_default=True,
              type=click.Choice(['all', 'Archaea_16S', 'Bacteria_16S', 'Fungi_ITS', '18S', 'COI']))
def primer(**kargs):
    """
    primer.yaml 파일에 등록된 Primer 정보들을 출력합니다.
    """
    def print_list(p_primer, p_region):
        if p_region == 'all':
            for region in p_primer:
                click.echo('-' * 60)
                click.secho(region, fg='cyan')
                for name in p_primer[region]:
                    if p_primer[region][name].get('default'):
                        click.secho(name, fg='yellow')
                        pprint(p_primer[region][name])
                        click.echo()
                    else:
                        click.echo(name)
                        pprint(p_primer[region][name])
                        click.echo()
        else:
            for name in p_primer[p_region]:
                if p_primer[p_region][name].get('default'):
                    click.secho(name, fg='yellow')
                    pprint(p_primer[p_region][name])
                    click.echo()
                else:
                    click.echo(name)
                    pprint(p_primer[p_region][name])
                    click.echo()

    dic_primer = parse_config('primer')
    if kargs['region'] == 'all':
        print_list(dic_primer, 'all')
    elif kargs['region'] == 'Archaea_16S':
        click.secho('Archaea_16S', fg='cyan')
        click.echo('-' * 60)
        print_list(dic_primer, 'Archaea_16S')
    elif kargs['region'] == 'Bacteria_16S':
        click.secho('Bacteria_16S', fg='cyan')
        click.echo('-' * 60)
        print_list(dic_primer, 'Bacteria_16S')
    elif kargs['region'] == 'Fungi_ITS':
        click.secho('Fungi_ITS', fg='cyan')
        click.echo('-' * 60)
        print_list(dic_primer, 'Fungi_ITS')
    elif kargs['region'] == '18S':
        click.secho('18S', fg='cyan')
        click.echo('-' * 60)
        print_list(dic_primer, '18S')
    elif kargs['region'] == 'COI':
        click.secho('COI', fg='cyan')
        click.echo('-' * 60)
        print_list(dic_primer, 'COI')
    else:
        raise ValueError(click.style('region 정보 오류 : {}'.format(kargs['region']), fg='red'))


@main.command('FLASH', short_help='Read Assembly')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--target_dir_suffix', '-t',
              default=TARGET_DIR_SUFFIX,
              show_default=True)
@click.option('--sample_suffix', '-s',
              help='검출하고자 하는 시료의 접미사. ex) V3.4 ITS3.4')
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']))
@click.option('--read_length', '-l',
              default=READ_LENGTH,
              show_default=True,
              type=click.IntRange(51, 301))
@click.option('--no_use_worker', '-noW',
              is_flag=True,
              help='FLASH 실행시 Worker를 사용하지 않음. 시료가 많아 Error 발생하는 경우 사용.')
@click.option('--target_size', '-ts',
              required=True,
              type=int,
              help='목표영역의 크기. --min_overlap 또는 --max_overlap의 값이 None일 경우, 해당 값을 자동 계산. '
                   '계산 방법은 --min_overlap, --max_overlap 옵션 참조.')
@click.option('--min_overlap', '-m',
              default=MIN_OVERLAP,
              show_default=True,
              help='None일 경우 --target_size를 이용하여 자동 계산. '
                   '수식: (read_length * 2) - Target_Size - 10'.format(MIN_OVERLAP))
@click.option('--max_overlap', '-M',
              default=MAX_OVERLAP,
              show_default=True,
              help='None일 경우 target_size를 이용하여 자동 계산. '
                   '수식: (read_length * 2) - Target_Size + 10.')
@click.option('--percent_combined_alarm', '-pca',
              default=PERCENT_COMBINED_ALARM,
              show_default=True,
              type=click.FloatRange(0, 100))
@click.option('--no_compress',
              is_flag=True,
              help='FLASH의 결과파일을 압축하지 않음.')
@click.option('--no_csv', '-noC',
              is_flag=True,
              help='CSV형식의 FLASH Summary 파일 생성 안함.')
@click.option('--no_excel', '-noE',
              is_flag=True,
              help='엑셀형식의 FLASH Summary 파일 생성 안함.')
@click.option('--no_Hcsv', '-noH',
              is_flag=True,
              help='Hcsv형식의 FLASH Summary 파일 생성 안함.')
@click.option('--length_filter', '-lf',
              default=LENGTH_FILTER,
              show_default=True,
              type=click.BOOL,
              required=True,
              help='Assembled Read에 대한 Length Fitler 적용 여부. '
                   '상한값, 하한값 중 하나만 적용하고자 할 경우 -lM, -lm 에서 no 입력.')
@click.option('--length_min', '-lm',
              # required=True,
              help='Length Filter의 최소길이. 최소길이미만 제거(<). 미적용시 no 입력.')
@click.option('--length_max', '-lM',
              # required=True,
              help='Length Filter의 최대길이. 최대길이초과 제거(>). 미적용시 no 입력.')
@click.option('--modify_header', '-mh',
              default=MODIFY_HEADER,
              show_default=True,
              type=click.BOOL,
              help='FASTQ의 Hader 정보를 ">시료명_SeqNum"으로 변경.')
@click.option('--pooled_sample', '-ps',
              default=POOLED_SAMPLE_FILE,
              show_default=True,
              type=click.BOOL,
              help='시료들의 Assembled Read 파일들을 하나의 파일로 통합. '
                   '미적용시 Clustering(OTU Picking):CD-HIT-OTU 분석시 별도 진행 필요.')
@click.option('--stat',
              default=True,
              show_default=True,
              type=click.BOOL,
              help='Filtered Assembled Read(fastq)에 대한 STAT 생성.'
                   '--out_seq_type 옵션이 FASTA일 경우 Assembled Read에 대한 STAT 생성(Filtered 아님).')
@click.option('--out_seq_type', '-st',
              default='FASTQ',
              show_default=True,
              type=click.Choice(['FASTQ', 'FASTA']),
              help='Filtered Read 파일과 Pooled Sample 파일의 형식을 지정.'
                   '--length_filter 과 --pooled_sample 이 모두 False 일 경우, 해당 옵션 적용 안됨.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def flash(**kargs):
    """
    \b
    FLASH 프로그램을 이용하여 Read Assembly를 진행합니다.
    Read Assembly, Summary 파일 생성, Length Filter, Pooled Sample File 생성.

    \b
    목표 영역에 대한 권장 설정 값
    Bakt(341F-805R) - Target Size: 465bp, Length Filter: 400, 500
                    : -ts 465 -lm 400 -lM 500
    ITS(3F-4R) - Target Size: 430bp, Length Filter: 300, 500
                    : -ts 430 -lm 300 -lM 500
                    : -ts 430 -m 100 -M 200 -lm 300 -lM 500
                    : -ts 430 -m 10 -M 200 -lm 300 -lM 500
    다른 Primer Set 정보에 대해서는 아래의 명령어를 이용하여 정보를 확인하세요.
        CoFI PRIMER --help

    \b
    Closed-Reference OTU Picking 진행시
    pooled_sample 방식
        : -ts 465 -lm 400 -lM 500 -st FASTA --no_compress
    by_sample 방식
        : -ts 465 -lm 400 -lM 500 -ps False -st FASTA --no_compress
    """
    start_time = time()
    from CoFI.core import flash_pipeline
    from SpoON.run_time import print_cofi_flash_run_time

    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass

    if kargs['length_filter'] is True:
        if kargs['length_min']:
            if kargs['length_min'].isdigit():
                pass
            elif kargs['length_min'] == 'no':
                pass
            else:
                click.echo('Error: --length_min 옵션의 입력값이 잘못되었습니다. (숫자 or no)')
                exit()
        elif kargs['length_min'] is None:
            click.echo('Error: --length_min 옵션을 입력하세요.')
            exit()

        if kargs['length_max']:
            if kargs['length_max'].isdigit():
                pass
            elif kargs['length_max'] == 'no':
                pass
            else:
                click.echo('Error: --length_max 옵션의 입력값이 잘못되었습니다. (숫자 or no)')
                exit()
        elif kargs['length_max'] is None:
            click.echo('Error: --length_max 옵션을 입력하세요.')
            exit()

    pprint(kargs)
    if check_file_type(os.path.join(kargs['analysis_base_path'], kargs['order_number']), 'exists'):
        _, assembly_run_time = flash_pipeline(kargs)
        end_time = time()
        print_cofi_flash_run_time(start_time, assembly_run_time, end_time)
    else:
        click.secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')


@main.command('LENGTH_FILTER', short_help='FLASH Excel & Length Filter')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']))
@click.option('--flash_log_suffix',
              default='*FLASH.log',
              show_default=True,
              help='FLASH Log을 저장한 파일의 접미사.')
@click.option('--mode', '-m',
              required=True,
              type=click.Choice(['excel', 'length_filter', 'all']),
              help='excel: FLASH*.log 파일들을 엑셀로 정리. '
                   'length_filter: Assembled Read들에 Length Filtering 적용. '
                   'all: excel & length fitler 모두 적용.')
@click.option('--assembled_dir', '-ad',
              required=False,
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='Read_Assembly의 경로.')
@click.option('--target_size', '-ts',
              type=int,
              help='csv, Hcsv, Excel 출력용.')
@click.option('--percent_combined_alarm', '-pca',
              default=PERCENT_COMBINED_ALARM,
              show_default=True,
              type=click.FloatRange(0, 100))
@click.option('--compression',
              default=True,
              show_default=True,
              type=click.BOOL,
              help='FLASH 결과 파일들의 압축 여부.')
@click.option('--length_min', '-lm',
              help='Length Filter의 최소길이. 최소길이미만 제거(<). 미적용시 no 입력.')
@click.option('--length_max', '-lM',
              help='Length Filter의 최대길이. 최대길이초과 제거(>). 미적용시 no 입력.')
@click.option('--modify_header', '-mh',
              default=MODIFY_HEADER,
              show_default=True,
              type=click.BOOL,
              help='FASTQ의 Hader 정보를 ">시료명_SeqNum"으로 변경.'.format(MODIFY_HEADER))
@click.option('--pooled_sample', '-ps',
              default=POOLED_SAMPLE_FILE,
              show_default=True,
              type=click.BOOL,
              help='시료들의 Assembled Read 파일들을 하나의 파일로 통합. '
                   '미적용시 Clustering(OTU Picking):CD-HIT-OTU 분석시 별도 진행 필요.'.format(POOLED_SAMPLE_FILE))
@click.option('--stat',
              default=True,
              show_default=True,
              type=click.BOOL,
              help='Filtered Assembled Read(fastq)에 대한 STAT 생성.'
                   '--out_seq_type 옵션이 FASTA일 경우 Assembled Read에 대한 STAT 생성(Filtered 아님).')
@click.option('--out_seq_type', '-st',
              default='FASTQ',
              show_default=True,
              type=click.Choice(['FASTQ', 'FASTA']),
              help='Filtered Read 파일과 Pooled Sample 파일의 형식을 지정.'
                   '--length_filter 과 --pooled_sample 이 모두 False 일 경우, 해당 옵션 적용 안됨.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def length_filter(**kargs):
    """
    \b
    FLASH 실행 이후, FLASH Log 파일들의 정보를 정리하여 Excel 파일을 생성한다.
    Length Filtering 을 수행한다.

    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI LENGTH_FILTER HN001XXXXX Analysis_1
        --mode [] --target_size [] --length_min [] --lenght_max []
    \b
    \033[1;32;40m<실행방법2>
    --assembled_dir 옵션 사용\033[m
    CoIF LENGTH_FILTER HN001XXXXX None --assembled_dir 디렉터리경로/Read_Assembly
        --mode [] --target_size [] --length_min [] --length_max []

    \b
    \033[1;35;40m<기존 데이터 삭제> (자동삭제)\033[m
    rm F_*.fastq.gz *.fastq.gz.log STAT.bash soft_link.bash STAT.txt
    rm HN00*
    find -type l | xargs rm

    \b
    \033[1;35;40m<주의>
    시료의 개수가 255개를 넘으면, Excel chart 생성 제약으로 인해 chart에 데이터가 모두 표현이 안 된다.\033[m

    \b
    Closed-Reference OTU Picking 진행시
    pooled_sample 방식
     -ad {PATH}/Read_Assembly -on {수주번호} -m all -ts 465 -lm 400 -lM 500 -st FASTA --compression False

    """
    start_time = time()
    from CoFI.core import length_filter_pipeline
    from SpoON.run_time import print_cofi_legnth_filter_run_time
    pprint(kargs)

    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if kargs['order_number'] != 'None':
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
        elif kargs['assembled_dir'] is not None:
            pass
        else:
            click.secho('Error: Argument & Option 조합 오류.', fg='red')

    if (kargs['mode'] == 'length_filter') or (kargs['mode'] == 'all'):
        if kargs['length_min'] is None:
            click.secho('--length_min 옵션의 값이 필요합니다.', fg='magenta')
            exit()
        elif kargs['length_max'] is None:
            click.secho('--length_max 옵션의 값이 필요합니다.', fg='magenta')
            exit()
        else:
            pass
    end_time = time()
    _, run_time = length_filter_pipeline(kargs)
    print_cofi_legnth_filter_run_time(start_time, run_time, end_time)


@main.command('CD-HIT-OTU', short_help='OTU: Denovo 방식')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--pooled_sample', '-ps',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='Order Number 또는 Analysis_Number Argument에 None 을 기입.')
@click.option('--out_dir', '-o',
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='결과 출력 디렉터리. 입력하지 않으면 Pooled_Sample.fastq 파일을 기준으로 자동 생성.')
@click.option('--OTU_cutoff', '-c',
              default=OTU_CUTOFF,
              show_default=True,
              type=click.FloatRange(0.7, 1))
@click.option('--prefix_length_or_primer', '-p',
              required=True,
              help='cd-hit-otu의 기본값: 6. '
                   'If a prefix-length (a digit number) is provided get the consensus of prefix of the all reads '
                   'remove the reads without this consensus). '
                   'If a primers_sequence is provided, read primers from this file, '
                   'remove the reads don\'t match the primers.')
@click.option('--trim_cutoff', '-t',
              default=TRIM_CUTOFF,
              show_default=True,
              help='default 1.0 (means no trimming) '
                   'if cutoff is a integer number > 1 (like 200), the program will trim reads to this length '
                   'if cutoff is a fraction (like 0.8), the program will trim reads to its length times this fraction'
                   ''.format(TRIM_CUTOFF))
@click.option('--singleton', 'small_cluster',
              default='singleton',
              show_default=True,
              flag_value='singleton',
              help='기본적용. Singleton 포함 여부')
@click.option('--doubleton', 'small_cluster',
              flag_value='doubleton',
              help='Doubleton 포함 여부')
@click.option('--check_chimera', '-m',
              default=CHECK_CHIMERA,
              show_default=True,
              type=click.BOOL)
@click.option('--pcr_error', '-e',
              default=PCR_ERROR,
              show_default=True,
              help='per-base PCR error.'.format(PCR_ERROR))
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def cd_hit_otu(**kargs):
    """
    \b
    CD-HIT-OTU 프로그램을 이용하여 denovo방식의 Clustering & OTU Picking을 진행합니다.
    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI CD-HIT-OTU HN001XXXXX Analysis_1
    \b
    \033[1;32;40m<실행방법2>
    --pooled_sample 옵션 사용\033[m
    CoIF CD-HIT-OTU None None --pooled_sample HN001XXXXX.pooled.fastq
    """
    start_time = time()
    from CoFI.core import cd_hit_otu_pipeline
    from SpoON.run_time import print_cofi_cd_hit_run_time
    pprint(kargs)
    if ((kargs['order_number'] != 'None') or (kargs['analysis_number'] != 'None')) \
            and (kargs['pooled_sample'] is not None):
        click.secho('Error: ORDER_NUMBER & ANALYSIS_NUMBER Argument 와 '
                    '--pooled_sample 옵션을 같이 사용할 수 없습니다.', fg='red')
        exit()

    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if kargs['order_number'] != 'None':
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
        elif kargs['pooled_sample'] is not None:
            pass
        else:
            click.secho('Error: Argument & Option 조합 오류.', fg='red')

    cd_hit_runt_time = cd_hit_otu_pipeline(kargs)
    end_time = time()
    print_cofi_cd_hit_run_time(start_time, cd_hit_runt_time, end_time)


@main.command('FLASH_HIT', short_help='FLASH --> CD-HIT-OTU(denovo)')
@click.argument('order_number')
@click.option('--target_region', '-r',
              required=True,
              type=click.Choice(['Bakt_341F-805R', 'ITS_3F-4R', '27F-Bact338']))
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--target_dir_suffix', '-td',
              default=TARGET_DIR_SUFFIX,
              show_default=True)
@click.option('--sample_suffix', '-s',
              help='검출하고자 하는 시료의 접미사. ex) V3.4 ITS3.4')
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']))
@click.option('--read_length', '-l',
              default=READ_LENGTH,
              show_default=True,
              type=click.IntRange(51, 301))
@click.option('--no_use_worker', '-noW',
              is_flag=True,
              help='FLASH 실행시 Worker를 사용하지 않음. 시료가 많아 Error 발생하는 경우 사용.')
@click.option('--percent_combined_alarm', '-pca',
              default=PERCENT_COMBINED_ALARM,
              show_default=True,
              type=click.FloatRange(0, 100))
@click.option('--no_compress',
              is_flag=True,
              help='FLASH의 결과파일을 압축하지 않음.')
@click.option('--OTU_cutoff', '-c',
              default=OTU_CUTOFF,
              show_default=True,
              type=click.FloatRange(0.7, 1))
@click.option('--trim_cutoff', '-t',
              default=TRIM_CUTOFF,
              show_default=True,
              help='default 1.0 (means no trimming) '
                   'if cutoff is a integer number > 1 (like 200), the program will trim reads to this length '
                   'if cutoff is a fraction (like 0.8), the program will trim reads to its length times this fraction'
                   ''.format(TRIM_CUTOFF))
@click.option('--singleton', 'small_cluster',
              default='singleton',
              show_default=True,
              flag_value='singleton',
              help='기본적용. Singleton 포함 여부')
@click.option('--doubleton', 'small_cluster',
              flag_value='doubleton',
              help='Doubleton 포함 여부')
@click.option('--check_chimera', '-m',
              default=CHECK_CHIMERA,
              show_default=True,
              type=click.BOOL)
@click.option('--pcr_error', '-e',
              default=PCR_ERROR,
              show_default=True,
              help='per-base PCR error.'.format(PCR_ERROR))
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def flash_hit(**kargs):
    """
    \b
    16S(Bacteria)와 ITS(Fungi)의 기본 Primer Set에 의해 생산된 데이터를 분석하는 파이프라인.
    목표 영역에 적합한 설정값으로 Read Assembly 및 Clustering(OTU Picking)을 진행합니다.

    \b
    FLASH_HIT 작업 구성
    시료 목록 구성 --> 작업 디렉터리 생성 --> FLASH -->
    Length Filer & Header 변경 --> Pooled Sample File 생성 --> CD-HIT-OTU

    \b
    영역별 기본 Primer Set
    16S(Bacteria) V3V4 - Bakt(341F-805R)
    16S(Bacteria)(유산균) V1V2 - 27F-Bact338
    ITS(Fungi) ITS2 - ITS(3F-4R)
    해당 Primer Set에 대한 정보는 primer.yaml 에 등록된 정보입니다.
    """
    start_time = time()
    from CoFI.core import flash_hit_pipeline
    from SpoON.run_time import print_flash_hit_run_time

    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass

    pprint(kargs)
    if check_file_type(os.path.join(*(kargs['analysis_base_path'], kargs['order_number'])), 'exists'):
        dic_primer = parse_config('primer')
        if kargs['target_region'] == 'Bakt_341F-805R':
            dic_target_info = dic_primer['Bacteria_16S']['Bakt(341F-805R)']
        elif kargs['target_region'] == 'ITS_3F-4R':
            dic_target_info = dic_primer['Fungi_ITS']['ITS(3F-4R)']
        elif kargs['target_region'] == '27F-Bact338':
            dic_target_info = dic_primer['Bacteria_16S']['27F-Bact338']
        else:
            raise ValueError(click.style('예기치 않은 오류: {}'.format(kargs['target_region'])))
        flash_run_time, cd_hit_run_time = flash_hit_pipeline(kargs, dic_target_info)
        end_time = time()
        print_flash_hit_run_time(start_time, flash_run_time, cd_hit_run_time, end_time)
    else:
        click.secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')


@main.command('CLOSED', short_help='Closed-Reference OTU Picking(UCLUST)')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--mode', '-m',
              type=click.Choice(['pooled_sample', 'by_sample']),
              default='pooled_sample',
              show_default=True,
              required=True,
              help='Closed-Reference OTU Picking 실행 방법. '
                   'pooled_sample: 1개의 FASTA 파일에 모든 시료의 Assembled Read가 포함된 경우. '
                   'by_sample: 시료별로 Closed-Reference OTU Picking을 진행할 경우. 시료별로 디렉터리가 생성됨. '
                   'OTU Table(biom)을 하나로 합치는 작업이 필요.')
@click.option('--database', '-db',
              type=click.Choice(DATABASE_LIST),
              required=True,
              help='Closed-Reference OTU Picking 에 사용할 DB. '
                   'Taxonomy Assignment를 진행하는 경우 Assignment 에도 사용.')
@click.option('--similarity', '-s',
              default='0.97',
              show_default=True,
              type=click.FloatRange(0.90, 1),
              help='OTU Picking Cut-off. 서열 유사도 기준.')
@click.option('--pooled_sample', '-ps',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='pooled.fasta.')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='--by_sample 옵션 사용시')
@click.option('--target_dir_suffix', '-td',
              default=TARGET_DIR_SUFFIX,
              show_default=True,
              help=f'--by_sample 옵션 사용시')
@click.option('--assembly_suffix',
              default=ASSEMBLY_SUFFIX,
              show_default=True,
              help='--by_sample 옵션 사용시. Assembled Read 파일들의 접미사.')
@click.option('--jobs',
              default=CLOSED_JOBS_TO_START,
              show_default=True,
              type=click.INT,
              help='작업의 개수. 12개 이상부터 시간이 줄지 않음.')
@click.option('--no_parallel',
              is_flag=True,
              help='작업을 나눠서 병렬로 처리 안함.')
@click.option('--no_assign',
              is_flag=True,
              help='Taxonomy Assignment 진행 안함.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def closed_otu_picking(**kargs):
    """
    \b
    Closed-Reference 방법을 이용하여 OTU Picking을 진행합니다.
    입력되는 서열들은 \033[1;5;32;40mFASTA\033[m 형식이여야 합니다.
    CoFI.py FLASH 에서 --out_seq_type 옵션을 fasta로 설정하면 FASTA로 생성이 됩니다.

    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI CLOSE HN001XXXXX Analysis_1 --mode [pooled_sample | by_sample] --database []

    \b
    \033[1;32;40m<실행방법2>
    --i_filtered 옵션 사용: --mode pooled_sample 만 가능\033[m
    CoFI CLOSE None None --mode pooled_sample
                         --pooled_sample 디렉터리경로/HN001XXXXX.pooled.fasta
                         --database []

    \b
    작업을 나눠서 병렬로 처리할 경우
    --pooled_sample 옵션에 \033[1;5;32;40m절대경로\033[m로 파일을 전달해야 합니다.
    jobs의 개수를 12개 이상 설정하는 경우,
    처리시간이 크게 개선되지 않으며, 시간이 더 걸릴 수도 있습니다.
    by_sample 방식으로 실행하세요.

    \b
    실행 방식
    pooled_sample : 여러 개의 시료들을 하나의 파일로 만들어서 진행.
    by_sample     : 시료별로 진행.
    시료의 개수가 많고 빠른 처리가 필요할 경우 by_sample 방식으로 실행.
    (Table 병합 시간 테스트 안함)
    \033[1;32;40m<by_sample 방식으로 진행시>\033[m
      - FLASH --no_compress 옵션을 사용하여야 합니다.
      - FLASH --pooled_sample 옵션을 False 로 설정해야 합니다.
      ex) CoFI.py FLASH -ts 465 -lm 400 -lM 500 -ps False -st FASTA --no_compress
      \b
      - Table(biom)을 하나로 합치는 작업이 필요합니다.
        merge_otu_tables.py 사용.
    """
    start_time = time()
    from CoFI.core import closed_otu_pipeline
    from SpoON.run_time import print_cofi_closed_otu_run_time
    pprint(kargs)

    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    elif kargs['pooled_sample'] is not None:
        pass
    else:
        click.secho('Error: Argument & Option 조합 오류.', fg='red')

    if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
        pass
    else:
        if (kargs['mode'] == 'pooled_sample') and (kargs['pooled_sample'] is None):
            click.secho('Error: --pooled_sample 옵션을 사용하세요.', fg='red')
            exit()
        else:
            if os.path.splitext(kargs['pooled_sample']) == 'fastq':
                click.secho('Error: FASTQ 형식을 사용할 수 없습니다.', fg='red')
                exit()
    closed_otu_run_time = closed_otu_pipeline(kargs)
    end_time = time()
    print_cofi_closed_otu_run_time(start_time, closed_otu_run_time, end_time)


@main.command('Q2_DADA2', short_help='ASVs: Denoising by QIIME2')
@click.argument('order_number')
@click.option('--metadata', '-m',
              required=True,
              help='Metadata가 없을 경우 no 입력. 시료명으로 구성된 metadata file 생성.')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--target_dir_suffix', '-t',
              default=TARGET_DIR_SUFFIX,
              show_default=True)
@click.option('--sample_suffix', '-s',
              help='검출하고자 하는 시료의 접미사. ex) V3.4 ITS3.4')
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--adapter_trim_tool', '-at',
              default=ADAPTER_TRIM_TOOL,
              show_default=True,
              type=click.Choice(['fastp', 'SeqPurge', 'Scythe']))
@click.option('--random_read',
              default=RANDOM_READ,
              show_default=True,
              type=click.INT,
              help='Quality Score Plot 작성시 반영되는 Read의 개수. '
                   'QIIME2 demux summarize 수행시 사용.')
@click.option('--p_trunc_len_f',
              default=P_TRUNC_LEN_F,
              show_default=True,
              type=click.IntRange(0, 301),
              help='Position at which forward read sequences should be truncated due to decrease in quality. '
                   'This truncates the 3\' end of the of the input sequences, '
                   'which will be the bases that were sequenced in the last cycles. '
                   'Reads that are shorter than this value will be discarded. '
                   'After this parameter is applied there must still be at least a 20 nucleotide overlap '
                   'between the forward and reverse reads. '
                   'If 0 is provided, no truncation or length filtering will be performed'.format(P_TRUNC_LEN_F))
@click.option('--p_trunc_len_r',
              default=P_TRUNC_LEN_R,
              show_default=True,
              type=click.IntRange(0, 301),
              help='Position at which reverse read sequences; should be truncated due to decrease in; quality. '
                   'This truncates the 3\' end of the of; the input sequences, '
                   'which will be the bases that were sequenced in the last cycles. '
                   'Reads that are shorter than this value will be discarded. '
                   'After this parameter is applied there must still be at least a 20 nucleotide overlap '
                   'between the forward and; reverse reads. '
                   'If 0 is provided, no truncation or length filtering will be performed'.format(P_TRUNC_LEN_R))
@click.option('--p_trim_left_f',
              default=P_TRIM_LEFT_F,
              show_default=True,
              type=click.IntRange(0, 100),
              help='Position at which forward read sequences should be trimmed due to low quality. '
                   'This trims the 5\' end of the input sequences, '
                   'which will be the bases that were sequenced in the first cycles.'.format(P_TRIM_LEFT_F))
@click.option('--p_trim_left_r',
              default=P_TRIM_LEFT_R,
              show_default=True,
              type=click.IntRange(0, 100),
              help='Position at which reverse read sequences should be trimmed due to low quality. '
                   'This trims the 5\' end of the input sequences, '
                   'which will be the bases that were sequenced in the first cycles.'.format(P_TRIM_LEFT_R))
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def q2_dada2(**kargs):
    """
    \b
    DADA2를 통해 Denoising을 진행하며, ASVs(Amplicon Sequence Variants, 100% OTUs)를 생성합니다.
    - noisy sequences 제거, marginal sequences의 에러 교정, chimeric sequences 제거, singletons 제거
    - Read Assembly(join denoised paired-end reads), 대표서열 선정(dereplicate sequences)

    \b
    DADA2 PipeLine 작업 진행
    시료 목록 구성 --> 작업 디렉터리 생성 --> (Metadata File 생성) -->
    FASTQ Import(qza) --> FASTQ Summary(qzv) --> DADA2 -->
    DATA2 STATs(qav) --> Feature Table summary(qzv) --> Feature Table Sequences(qzv)

    \b
    \033[1;31;40m-주의-\033[m
    \033[1;31;40mITS 영역의 경우 사용하지 마세요.\033[m
    ITS 영역의 경우 InDel에 의한 Target Size의 변화로 인해 DADA2의 옵션값 설정이 어렵습니다.
    """
    from CoFI.core import dada2_pipeline
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    pprint(kargs)
    dada2_pipeline(kargs)


@main.command('R_DADA2', short_help='ASVs: Denoising by R')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--target_dir_suffix', '-pt',
              default=TARGET_DIR_SUFFIX,
              show_default=True,
              help='Analysis Run Dir')
@click.option('--R1_suffix', '-r1',
              default=R1_SUFFIX,
              show_default=True)
@click.option('--R2_suffix', '-r2',
              default=R2_SUFFIX,
              show_default=True)
@click.option('--analysis_dir_base_name', '-abn',
              default=ANALYSIS_DIR_BASE_NAME,
              show_default=True)
@click.option('--trim_left_F', '-tmf',
              required=True,
              type=click.INT,
              help='Forward Primer 서열의 길이. Read1에서 Primer 서열 제거. '
                   'Bakt_341F: 17. ITS_3F: 20.')
@click.option('--trim_left_R', '-tmr',
              required=True,
              type=click.INT,
              help='Reversed Primer 서열의 길이. Read2에서 Primer 서열 제거. '
                   'Bakt_805R: 21. ITS_4R: 20.')
@click.option('--trunc_Len_F', '-tcf',
              required=True,
              type=click.INT,
              help='Read1을 잘라내고 남기고 싶은 Read1의 길이. '
                   'Read1의 길이가 301bp일 경우, 250bp를 입력하면 301bp를 잘라서 250bp를 만든다.')
@click.option('--trunc_Len_R', '-tcr',
              required=True,
              type=click.INT,
              help='Read2을 잘라내고 남기고 싶은 Read2의 길이. '
                   'Read2의 길이가 301bp일 경우, 250bp를 입력하면 301bp를 잘라서 250bp를 만든다.')
@click.option('--maxN',
              default=0,
              show_default=True,
              type=click.INT,
              help='')
@click.option('--maxEE_F',
              default=5,
              show_default=True,
              type=click.INT,
              help='')
@click.option('--maxEE_R',
              default=5,
              show_default=True,
              type=click.INT,
              help='')
@click.option('--truncQ',
              default=2,
              show_default=True,
              type=click.INT,
              help='')
@click.option('--min_overlap',
              default=15,
              show_default=True,
              type=click.INT,
              help='')
@click.option('--chimera_method',
              default='consensus',
              show_default=True,
              type=click.Choice(['consensus', 'pooled', 'per-sample']),
              help='')
@click.option('--dada2_dir_name', '-dn',
              default='R_DADA2',
              show_default=True,
              help='DADA2 디렉터리 이름.')
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
    """
    \b
    DADA2 프로그램을 이용하여 ASVs(Amplicon Sequence Variants) 방식으로 분석을 진행한다.
    DADA2(R 사용)는 직접 작성한 스크립트를 통해 실행된다.
    \b
    \033[1;32m16S V3V4 : Bakt_341F-805R & 300bp X 2\033[m
    --trim_left_F 17
    --trim_left_R 21
    --trunc_Len_F 250
    --trunc_Len_R 250
    \b
    \033[1;31m** 주의: 옵션 설정값에 따라 결과의 차이가 심하게 발생한다.\033[m
    --trunc_Len_F & --trunc_Len_R 옵션은 FASTQ의 Quality Score에 따라 적절하게 선택한다.
    \b
    DADA2 옵션
    --trim_left_F    --trim_left_R
    --trunc_Len_F    --trunc_Len_R
    --maxN
    --maxEE_F        --maxEE_R
    --truncQ
    --min_overlap
    --chimera_method
    """
    from CoFI.core import r_dada2_pipeline
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    r_dada2_pipeline(kargs)


@main.command('MAFFT_FastTree', short_help='Alignment & Phylogeny')
@click.option('--i_sequences', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='FeatureData[Sequence]')
@click.option('--out_dir', '-o',
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='미입력시 --i_sequences 옵션에 전달된 경로를 참고하여 상위 디렉터리에 자동 생성.')
def mafft_fasttree(**kargs):
    from CoFI.core import maff_fasttree_pipeline
    pprint(kargs)
    maff_fasttree_pipeline(kargs)


@main.command('DIVERSITY', short_help='')
@click.option('--i_phylogeny', '-ip',
              required=True,
              help='Phylogeny[Rooted]. '
                   'Phylogenetic tree containing tip identifiers that correspond to the feature identifiers '
                   'in the table. This tree can contain tip ids that are not present in the table, '
                   'but all feature ids in the table must be present in this tree.')
@click.option('--i_table', '-it',
              required=True,
              help='FeatureTable[Frequency]. table.qza.'
                   'The feature table containing the samples over which diversity metrics should be computed.')
@click.option('--out_dir', '-o',
              help='미입력시 --i_phylogeny와 --i_table 옵션에 전달된 경로를 참고하여 상위 디렉터리에 자동 생성.'
                   'i_phylogeny와 --i_table의 상위 디렉터리 경로가 다를 경우  --i_table의 상위 디렉터리에 결과 생성.')
@click.option('--sample_depth', '-s',
              type=click.INT,
              help='미입력시 --i_table 옵션에 전달된 파일의 데이터에서 최소값을 추출.'
                   'The total frequency that each sample should be rarefied to prior to computing diversity metrics.')
@click.option('-metadata', '-m',
              required=True,
              multiple=True,
              help='')
def diversity(**kargs):
    from CoFI.core import diversity_pipeline
    pprint(kargs)
    diversity_pipeline(kargs)


@main.command('ALIGNMENT', short_help='QIIME1 muscle')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--i_otus_rep', '-i',
              required=False,
              help='otus_rep.fasta')
@click.option('--out_dir', '-o',
              help='미입릭시 --i_otus_rep 옵션에 전달된 경로를 참조하여 2단계 상위 디렉터리에 자동 생성.')
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 관련 파일 인식 경로가 달라짐.')
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='ALIGNMENT를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def alignment_muscle(**kargs):
    """
    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI ALIGNMENT HN001XXXXX Analysis_1

    \b
    \033[1;32;40m<실행방법2>
    --i_filtered 옵션 사용\033[m
    CoIF ALIGNMENT None None --i_otus_rep 디렉터리경로/otus_rep.fasta
    """
    from CoFI.core import alignment_muscle_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    elif kargs['i_otus_rep'] is not None:
        pass
    else:
        click.secho('Error: Argument & Option 조합 오류.', fg='red')
    alignment_muscle_pipeline(kargs)


@main.command('PHYLOGENY', short_help='QIIME1 Phylogeny')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 관련 파일 인식 경로가 달라짐.')
@click.option('--i_filtered', '-i',
              required=False,
              help='otus_rep_aligned_pfiltered.fasta')
@click.option('--out_dir', '-o',
              help='미입릭시 --i_otus_rep 옵션에 전달된 경로를 참조하여 2단계 상위 디렉터리에 자동 생성.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def phylogeny_qiime1(**kargs):
    """
    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI PHYLOGENY HN001XXXXX Analysis_1

    \b
    \033[1;32;40m<실행방법2>
    --i_filtered 옵션 사용\033[m
    CoIF PHYLOGENY None None --i_filtered 디렉터리경로/Alignment/filtered_alignment/otus_rep_aligned_pfiltered.fasta
    """
    from CoFI.core import phylogeny_qiime1_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    elif kargs['i_filtered'] is not None:
        pass
    else:
        click.secho('Error: Argument & Option 조합 오류.', fg='red')
    phylogeny_qiime1_pipeline(kargs)


@main.command('TAXONOMY', short_help='QIIME1 Taxonomy')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--i_otus_rep', '-i',
              required=False,
              help='otus_rep.fasta')
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 대표서열을 인식하는 경로가 달라짐.')
@click.option('--uclust_db', '-udb',
              type=click.Choice(UCLUST_DATABASE_LIST),
              multiple=True,
              help='RDP DB. 다중 선택 가능. ex) -udb [DB] -udb [DB]')
@click.option('--blast_db', '-bdb',
              type=click.Choice(BLAST_DATABASE_LIST),
              multiple=True,
              help='BLAST DB. 다중 선택 가능. ex) -dbd [DB] -bdb [DB]')
@click.option('--mode', '-m',
              default='all',
              show_default=True,
              type=click.Choice(['blast', 'parse', 'all']),
              help='BLAST DB 선택시 적용.')
@click.option('--remove_taxon', '-rt',
              type=click.Choice(['no_hit', 'filtered', 'uncultured', 'environmental', 'bacterium',
                                 'marine', 'unidentified', 'eukaryote', 'Cyanobacteria', 'Archaea',
                                 'all']),
              multiple=True,
              help='taxon 제거 조건. 다중 선택 가능. ex) -rm no_hit -rm filtered')
@click.option('--keep_taxon', '-kt',
              type=click.Choice(['Cyanobacteria', 'Archaea']),
              multiple=True,
              help='taxon 유지 조건. 다중 선택 가능. ex) -k Cyanobacteria')
@click.option('--query_coverage', '-qc',
              default=85,
              show_default=True,
              help='Query Coverage 필터 조건(%%). 기준값미만의  데이터의 Organism정보를 \"Unassignment\"로 표시')
@click.option('--identity_percentage', '-ip',
              default=85,
              show_default=True,
              help='Identity Percentage 필터 조건(%%). 기준값미만의  데이터의 Organism정보를 '
                   '\"Unassignment\"로 표시"')
@click.option('--read_per_job', '-rj',
              default=BLAST_READ_PER_JOB,
              show_default=True,
              help='Queue Job 당 OTU 대표 서열의 개수')
@click.option('--nt_max_job', '-nmj',
              default=BLAST_NT_MAX_JOB,
              show_default=True,
              help='NT DB일 경우 최대 생성할 수 있는 Queue Job의 개수. '
                   '분석 시간이 짧은 다른 작업이 진행되지 않는 문제 방지 목적. 변경시 문의 요망.')
@click.option('--queue', '-q',
              default=BLAST_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='작업을 등록할 Queue.')
@click.option('--out_dir', '-o',
              help='미입릭시 --i_otus_rep 옵션에 전달된 경로를 참조하여 2단계 상위 디렉터리에 자동 생성.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def taxonomy(**kargs):
    """
    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI TAXONOMY HN001XXXXX Analysis_1 --{blast_db | uclust_db}

    \b
    \033[1;32;40m<실행방법2>
    --i_filtered 옵션 사용\033[m
    CoIF TAXONOMY None None --i_otus_rep 디렉터리경로/Clustering/otus_rep.fasta --{blast_db | uclust_db}

    \b
    \033[1;32;40m<BLAST DB Option>\033[m
    --blast_db
    --mode
    --remove_taxon
    --query_coverage
    --identity_percentage
    --read_per_job
    --nt_max_job
    --queue
    """
    from CoFI.core import taxonomy_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    elif kargs['i_otus_rep'] is not None:
        pass
    else:
        click.secho('Error: Argument & Option 조합 오류.', fg='red')
    if all([kargs['remove_taxon'], kargs['keep_taxon']]):
        click.secho('Error: Argument & Option 조합 오류.', fg='red')
    taxonomy_pipeline(kargs)


@main.command('BIOM', short_help='QIIME1 BIOM')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--otu_map_fp', '-i',
              required=False,
              help='pick_otus.txt')
# @click.option('--out_dir_suffix', '-o',
#               help='BIOM 디렉터리뒤에 붙여질 접미사.')
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 관련 파일 인식 경로가 달라짐.')
@click.option('--biom', '-b',
              default='otu_table',
              show_default=True,
              help='생성할 biom 파일의 기본명.')
@click.option('--taxonomy', '-t',
              help='')
@click.option('--metadata', '-m',
              help='')
@click.option('--field', '-f',
              help='시료들을 통합할 metadata 필드 이름.')
@click.option('--collapse', '-c',
              type=click.Choice(['mean', 'sum', 'random', 'median', 'first']),
              help='metadata 파일에 기재된 각 그룹의 시료들을 지정된 방식으로 하나로 통합.')
@click.option('--sort_for_collapse', '-s',
               help='시료통합 후 그룹명(시료명)을 정렬할 순서가 담긴 파일.')
@click.option('--remove_sample', '-rs',
              multiple=True,
              help='BIOM 파일에서 특정 시료을 제거. ex) -rs A1 -rs A2. 동작 안 함.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def biom(**kargs):
    """
    Order Number & Analysis Number Argument 사용시, Taxonomy_assignment 디렉터리내의 모든 결과에 대해서 BIOM 파일을
    생성한다.

    \b
    \033[1;32;40m<실행방법1>
    Order Number & Analysis Number 사용\033[m
    CoFI BIOM HN001XXXXX Analysis_1

    \b
    \033[1;32;40m<실행방법2>
    --i_filtered 옵션 사용\033[m
    CoIF BIOM None None --otu_map_fp [] -

    \b
    \033[1;32;40m<시료 통합>
    각 그룹의 해당되는 시료들을 그룹명으로 통합\033[m
    CoFI BIOM HN001XXXXX Analysis_1 --field [] --collapse []
    """
    from CoFI.core import biom_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    elif kargs['otu_map_fp'] is not None:
        pass
    else:
        click.secho('Error: Argument & Option 조합 오류.', fg='red')
    biom_pipeline(kargs)


@main.command('AP', short_help='QIIME1 muscle & Phylogeny')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 관련 파일 인식 경로가 달라짐.')
@click.option('--queue', '-q',
              default=DEFAULT_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='ALIGNMENT를 실행하기 위한 Queue.')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def alignment_phylogeny(**kargs):
    from CoFI.core import alignment_muscle_pipeline, phylogeny_qiime1_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    else:
        click.secho('Error: Argument 오류.', fg='red')
    alignment_muscle_pipeline(kargs)
    phylogeny_qiime1_pipeline(kargs)


@main.command('APT', short_help='QIIME1 muscle & Phylogeny & Taxonomy')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--otu_method', '-om',
              type=click.Choice(['denovo', 'closed', 'r_dada2', 'auto']),
              default='auto',
              show_default=True,
              help='구성된 OTU를 생성한 방법. 값에 따라 관련 파일 인식 경로가 달라짐.')
@click.option('--uclust_db', '-udb',
              type=click.Choice(UCLUST_DATABASE_LIST),
              multiple=True,
              help='RDP DB')
@click.option('--blast_db', '-bdb',
              type=click.Choice(BLAST_DATABASE_LIST),
              multiple=True,
              help='BLAST DB. Probiotics DB 미설정(문의요망).')
@click.option('--mode', '-m',
              default='all',
              show_default=True,
              type=click.Choice(['blast', 'parse', 'all']),
              help='BLAST DB 선택시 적용.')
@click.option('--remove_taxon', '-rm',
              type=click.Choice(['no_hit', 'filtered', 'uncultured', 'environmental', 'bacterium',
                                 'marine', 'unidentified', 'eukaryote', 'all']),
              multiple=True,
              help='taxon 제거 조건. 다중 선택 가능. ex) -rm no_hit -rm filtered')
@click.option('--keep_taxon', '-k',
              type=click.Choice(['Cyanobacteria', 'Archaea']),
              multiple=True,
              help='taxon 유지 조건. 다중 선택 가능. ex) -k Cyanobacteria')
@click.option('--query_coverage', '-qc',
              default=85,
              show_default=True,
              help='Query Coverage 필터 조건(%%). 기준값미만의  데이터의 Organism정보를 \"Unassignment\"로 표시')
@click.option('--identity_percentage', '-ip',
              default=85,
              show_default=True,
              help='Identity Percentage 필터 조건(%%). 기준값미만의  데이터의 Organism정보를 '
                   '\"Unassignment\"로 표시"')
@click.option('--read_per_job', '-rj',
              default=BLAST_READ_PER_JOB,
              show_default=True,
              help='Queue Job 당 OTU 대표 서열의 개수.')
@click.option('--nt_max_job', '-nmj',
              default=BLAST_NT_MAX_JOB,
              show_default=True,
              help='NT DB일 경우 최대 생성할 수 있는 Queue Job의 개수. '
                   '분석 시간이 짧은 다른 작업이 진행되지 않는 문제 방지 목적. 변경시 문의 요망.')
@click.option('--queue', '-q',
              default=BLAST_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help=f'작업을 등록할 Queue')
@click.option('--no_queue',
              is_flag=True,
              help='Queue에 작업을 등록하지 않고 현재 서버에서 실행.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def alignment_phylogeny_taxonomy(**kargs):
    from CoFI.core import alignment_muscle_pipeline, phylogeny_qiime1_pipeline, taxonomy_pipeline
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['order_number'] != 'None':
        if kargs['no_order_number'] is True:
            pass
        else:
            if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
                pass
    else:
        click.secho('Error: Argument 오류.', fg='red')
    alignment_muscle_pipeline(kargs)
    phylogeny_qiime1_pipeline(kargs)
    taxonomy_pipeline(kargs)


if __name__ == '__main__':
    main()

