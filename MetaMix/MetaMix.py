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
__version__ = '1.1.2'

# --------------------------------------------
# 2020.05.20 - 1.1.2
# EXPRESSO

import click
import sys
import time
import os
from pprint import pprint
from SpoON.util import check_order_number_system, run_cmd, check_run_cmd, parse_config, glob_dir, check_file_type
from SpoON.run_time import print_metamix_run_time


# 기본값 설정
CONFIG = parse_config()
PREMA = CONFIG['MetaMix']['PreMA']
COFI = CONFIG['MetaMix']['CoFI']
THECUPS = CONFIG['MetaMix']['theCups']
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
QUEUE_LIST = CONFIG['Queue_List']
# TAXONOMY
BLAST_DATABASE_LIST = list(CONFIG['CoFI_TAXONOMY']['blast_db'].keys())
UCLUST_DATABASE_LIST = list(CONFIG['CoFI_TAXONOMY']['uclust_db'].keys())
S_DATABASE_LIST = set(BLAST_DATABASE_LIST + UCLUST_DATABASE_LIST)
BLAST_READ_PER_JOB = CONFIG['CoFI_TAXONOMY']['queue']['read_per_job']
BLAST_NT_MAX_JOB = CONFIG['CoFI_TAXONOMY']['queue']['nt_max_job']
BLAST_QUEUE = CONFIG['CoFI_TAXONOMY']['queue']['name']
try:
    S_DATABASE_LIST.remove('iBOL')
except KeyError:
    pass
except ValueError:
    pass
L_16S_DATABASE_LIST = list()
for ele in S_DATABASE_LIST:
    if ele == 'UNITE':
        continue
    elif 'UNITE' in ele:
        continue
    else:
        L_16S_DATABASE_LIST.append(ele)

L_ITS_DATABASE_LIST = list(S_DATABASE_LIST)
for ele in ['NCBI_16S', 'RDP', 'NCBI_Probiotics']:
    try:
        L_ITS_DATABASE_LIST.remove(ele)
    except KeyError:
        pass
    except ValueError:
        pass


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def metamix(**kargs):
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

        \b
             \033[1;32;40m<Amplicon Metagenome Analysis Pipeline>\033[m
        -------------------------------------------------------------
        Command      RawData -->    Analysis    -->  Report(별도 실행)
        -------------------------------------------------------------
        VANILLA      PreMA   --> CoFI FLASH_HIT -->  theCups
                                      ALIGNMENT
                                      PHYLOGENY
                                      TAXONOMY
                                      BIOM
        DoubleShot   PreMA   --> CoFI CLOSED    -->  theCups
        LongShot                      ALIGNMENT
                                      PHYLOGENY
                                      BIOM
        EXPRESSO    PreMA    --> CoFI R_DADA2   -->  teCups
                                      ALIGNMENT
                                      PHYLOGENY
                                      TAXONOMY
                                      BIOM
        -------------------------------------------------------------
        \033[1;32;40m***\033[m: 경우에 따라 직접 실행해야 되는 경우.
    """
    pass


@metamix.command('VANILLA', short_help='De Novo method OTU Picking: FLASH_HIT')
@click.version_option(version=__version__)
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--custom_sample_name', '-csn',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='변경할 시료명과 그룹 정보가 담긴 파일명. 각 행이 "#"으로 시작하면 주석 처리. '
                   '구조: #SampleID   New_name   Group1   Group2. 자세한 방법은 개발자에게 문의.')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_prema',
              is_flag=True,
              help='데이터 복사 없이(PreMA 실행없이, 이미 복사되었을 경우) 이후 단계를 진행할 때 사용.')
@click.option('--V1V2',
              help='V1V2(27F-Bact338) 시료을 구분하는 접미사(.) 주의. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--V3V4',
              help='V3V4(Bakt_341F-805R 시료을 구분하는 접미사(.) 주의. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--ITS34',
              help='ITS3-4(ITS_3F-4R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--target_region', '-r',
              required=True,
              type=click.Choice(['Bakt_341F-805R', 'ITS_3F-4R', '27F-Bact338', 'both']),
              help='Target Region 및 Primer Set.')
@click.option('--V1V2_DB',
              type=click.Choice(L_16S_DATABASE_LIST),
              multiple=True,
              help='BLAST: NCBI_16S, NCBI_NT, NCBI_Probiotics. UCLUST: RDP. '
                   '중복 사용 가능. ex) --V1V2_DB NCBI_16S --V1V2_db NCBI_Probiotics')
@click.option('--V3V4_DB',
              type=click.Choice(L_16S_DATABASE_LIST),
              multiple=True,
              help='BLAST: NCBI_16S, NCBI_NT, Probiotics. UCLUST: RDP. '
                   '중복 사용 가능. ex) --V3V4_DB NCBI_16S --V3V4_DB NCBI_Probiotics')
@click.option('--ITS34_DB',
              type=click.Choice(L_ITS_DATABASE_LIST),
              multiple=True,
              help='BLAST: NCBI_NT. UCLUST: UNITE. '
                   '중복 사용 가능. ex) --ITS34_DB NCBI_NT --ITS34_DB UNITE')
@click.option('--remove_taxon', '-rm',
              type=click.Choice(['no_hit', 'filtered', 'uncultured', 'environmental', 'bacterium',
                                 'marine', 'unidentified', 'eukaryote', 'all']),
              multiple=True,
              help='taxon 제거 조건. 다중 선택 가능. ex) -rm no_hit -rm filtered')
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
              help=f'Queue Job 당 OTU 대표 서열의 개수.')
@click.option('--nt_max_job', '-nmj',
              default=BLAST_NT_MAX_JOB,
              show_default=True,
              help='NT DB일 경우 최대 생성할 수 있는 Queue Job의 개수. '
                   '분석 시간이 짧은 다른 작업이 진행되지 않는 문제 방지 목적. 변경시 문의 요망.')
@click.option('--queue', '-q',
              default=BLAST_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='작업을 진행할 Queue.')
@click.option('--no_use_FLASH_worker', '-noW',
              is_flag=True,
              help='FLASH 실행시 Worker를 사용하지 않음. 시료가 많아 Error 발생하는 경우 사용.')
@click.option('--theCups_cmd', '-cups',
              default='MiSeq_V1',
              show_default=True,
              help='theCups command.')
@click.option('--no_theCups', '-noCups',
              is_flag=True,
              help='theCups 실행 안 함.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def vanilla(**kargs):
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

    \b
                      _ _ _
    __   ____ _ _ __ (_) | | __ _
    \ \ / / _` | '_ \| | | |/ _` |
     \ V / (_| | | | | | | | (_| |
      \_/ \__,_|_| |_|_|_|_|\__,_|

    \b
    \033[1;32;40m<Denovo OTU Picking PipeLine>\033[m
    PreMA CORE --> CoFI FLASH_HIT -->  theCups
                        ALIGNMENT
                        PHYLOGENY
                        TAXONOMY
                        BIOM
    """
    pprint(kargs)
    start_time = time.time()
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass

    if any([kargs['v1v2_db'], kargs['v3v4_db'], kargs['its34_db']]) is False:
        click.secho('Error: --V1V2_DB, --V3V4_DB 또는 --ITS34_DB 옵션을 입력하세요.', fg='red')
        exit()

    # PreMA - Data Copy & QC
    if kargs['integrate']:
        integrate_text = ''
        for ele in kargs['integrate']:
            integrate_text += '-i {} '.format(ele)
    else:
        integrate_text = ''

    if kargs['custom_sample_name']:
        custom_sample_name_text = f'-csn {kargs["custom_sample_name"]}'
    else:
        custom_sample_name_text = ''

    if kargs['no_prema'] is False:
        global PREMA
        prema_cmd = \
            'python3 {prema} ' \
            'CORE ' \
            '{integrate} ' \
            '{custom_sample_name} ' \
            '{copy} ' \
            '{no_order_number}' \
            '{order_number} '.format(
                prema=PREMA,
                integrate=integrate_text,
                custom_sample_name=custom_sample_name_text,
                copy='-c' if kargs['copy_rawdata'] else '',
                no_order_number='--no_order_number' if kargs['no_order_number'] else '',
                order_number=kargs['order_number'])
        run = run_cmd(prema_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'PreMA CORE 완료',
                'false_meg': 'PreMA CORE',
            }, p_exit=True)
        del prema_cmd, run
    prema_end_time = time.time()

    # Read Assembly & CD-HIT-OTU : FLASH_HIT 등
    if kargs['copy_rawdata'] or kargs['no_prema']:
        # 다중 옵션 처리 : remove_taxon
        if len(kargs['remove_taxon']) != 0:
            remove_taxon = ' '.join([f'-rm {x}' for x in kargs['remove_taxon']])
        else:
            remove_taxon = None
        # 여러 영역을 진행하지 않은 경우
        if kargs['target_region'] != 'both':
            l_db_name_tool = get_db_and_tool(kargs['target_region'], kargs['v1v2_db'],
                                             kargs['v3v4_db'], kargs['its34_db'])
            # Suffix 미입력.
            if (kargs['v1v2'] is None) and (kargs['v3v4'] is None) and (kargs['its34'] is None):
                pipeline_time = run_pipeline(
                    {
                        'base_path': CONFIG['base_path']['analysis_base_path'],
                        'order_number': kargs['order_number'],
                        'target_region': kargs['target_region'],
                        'l_db_name_tool': l_db_name_tool,
                        'sample_suffix': '',
                        'no_use_flash_worker': kargs['no_use_flash_worker'],
                        'custom_sample_name': kargs['custom_sample_name'],
                        'remove_taxon': remove_taxon,
                        'query_coverage': kargs['query_coverage'],
                        'identity_percentage': kargs['identity_percentage'],
                        'read_per_job': kargs['read_per_job'],
                        'nt_max_job': kargs['nt_max_job'],
                        'queue': kargs['queue'],
                        'theCups_cmd': kargs['thecups_cmd'],
                        'no_theCups': kargs['no_thecups'],
                        'no_order_number': kargs['no_order_number']
                    })
            else:
                if kargs['v1v2'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["v1v2"]}'
                elif kargs['v3v4'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["v3v4"]}'
                elif kargs['its34'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["its34"]}'
                else:
                    raise ValueError(f'알 수 없는 에러가 발생하였습니다.')
                pipeline_time = run_pipeline(
                    {
                        'base_path': CONFIG['base_path']['analysis_base_path'],
                        'order_number': kargs['order_number'],
                        'target_region': kargs['target_region'],
                        'l_db_name_tool': l_db_name_tool,
                        'sample_suffix': sample_suffix,
                        'no_use_flash_worker': kargs['no_use_flash_worker'],
                        'custom_sample_name': kargs['custom_sample_name'],
                        'remove_taxon': remove_taxon,
                        'query_coverage': kargs['query_coverage'],
                        'identity_percentage': kargs['identity_percentage'],
                        'read_per_job': kargs['read_per_job'],
                        'nt_max_job': kargs['nt_max_job'],
                        'queue': kargs['queue'],
                        'theCups_cmd': kargs['thecups_cmd'],
                        'no_theCups': kargs['no_thecups'],
                        'no_order_number': kargs['no_order_number']
                    })

        elif kargs['target_region'] == 'both':
            # V1V2 V3V4 영역을 동시 진행한 경우
            if (kargs['v1v2'] is not None) and (kargs['v3v4'] is not None) and (kargs['its34'] is None):
                pipeline_time = list()
                for region in ['v1v2', 'v3v4']:
                    if region == 'v1v2':
                        target_region = '27F-Bact338'
                    elif region == 'v3v4':
                        target_region = 'Bakt_341F-805R'
                    else:
                        raise ValueError()

                    l_db_name_tool = get_db_and_tool(target_region, kargs['v1v2_db'],
                                                     kargs['v3v4_db'], kargs['its34_db'])
                    pipeline_time.append(
                        run_pipeline({
                            'base_path': CONFIG['base_path']['analysis_base_path'],
                            'order_number': kargs['order_number'],
                            'target_region': target_region,
                            'l_db_name_tool': l_db_name_tool,
                            'sample_suffix': f'--sample_suffix {kargs[region]}',
                            'no_use_flash_worker': kargs['no_use_flash_worker'],
                            'custom_sample_name': kargs['custom_sample_name'],
                            'remove_taxon': remove_taxon,
                            'query_coverage': kargs['query_coverage'],
                            'identity_percentage': kargs['identity_percentage'],
                            'read_per_job': kargs['read_per_job'],
                            'nt_max_job': kargs['nt_max_job'],
                            'queue': kargs['queue'],
                            'theCups_cmd': kargs['thecups_cmd'],
                            'no_theCups': kargs['no_thecups'],
                            'no_order_number': kargs['no_order_number']
                        })
                    )

            # V3V4 ITS3-4 역영을 동시 진행한 경우
            elif (kargs['v1v2'] is None) and (kargs['v3v4'] is not None) and (kargs['its34'] is not None):
                pipeline_time = list()
                for region in ['v3v4', 'its34']:
                    if region == 'v3v4':
                        target_region = 'Bakt_341F-805R'
                    elif region == 'its34':
                        target_region = 'ITS_3F-4R'
                    else:
                        raise ValueError()

                    l_db_name_tool = get_db_and_tool(target_region, kargs['v1v2_db'],
                                                     kargs['v3v4_db'], kargs['its34_db'])
                    pipeline_time.append(
                        run_pipeline({
                            'base_path': CONFIG['base_path']['analysis_base_path'],
                            'order_number': kargs['order_number'],
                            'target_region': target_region,
                            'l_db_name_tool': l_db_name_tool,
                            'sample_suffix': f'--sample_suffix {kargs[region]}',
                            'no_use_flash_worker': kargs['no_use_flash_worker'],
                            'custom_sample_name': kargs['custom_sample_name'],
                            'remove_taxon': remove_taxon,
                            'query_coverage': kargs['query_coverage'],
                            'identity_percentage': kargs['identity_percentage'],
                            'read_per_job': kargs['read_per_job'],
                            'nt_max_job': kargs['nt_max_job'],
                            'queue': kargs['queue'],
                            'theCups_cmd': kargs['thecups_cmd'],
                            'no_theCups': kargs['no_thecups'],
                            'no_order_number': kargs['no_order_number']
                        })
                    )

        cofi_end_time = time.time()
        print_metamix_run_time(start_time, prema_end_time, pipeline_time, cofi_end_time)


@metamix.command('ESPRESSO - 개발 중', short_help='ASVs Method Pipleline: R_DADA2')
@click.version_option(version=__version__)
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--custom_sample_name', '-csn',
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='변경할 시료명과 그룹 정보가 담긴 파일명. 각 행이 "#"으로 시작하면 주석 처리. '
                   '구조: #SampleID   New_name   Group1   Group2. 자세한 방법은 개발자에게 문의.')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_prema',
              is_flag=True,
              help='데이터 복사 없이(PreMA 실행없이, 이미 복사되었을 경우) 이후 단계를 진행할 때 사용.')
# @click.option('--V1V2',
#               help='V1V2(27F-Bact338) 시료을 구분하는 접미사(.) 주의. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--V3V4',
              help='V3V4(Bakt_341F-805R 시료을 구분하는 접미사(.) 주의. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
# @click.option('--ITS34',
#               help='ITS3-4(ITS_3F-4R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--target_region', '-r',
              required=True,
              type=click.Choice(['Bakt_341F-805R']),
              help='Target Region 및 Primer Set.')
# @click.option('--V1V2_DB',
#               type=click.Choice(L_16S_DATABASE_LIST),
#               multiple=True,
#               help='BLAST: NCBI_16S, NCBI_NT, NCBI_Probiotics. UCLUST: RDP. '
#                    '중복 사용 가능. ex) --V1V2_DB NCBI_16S --V1V2_db NCBI_Probiotics')
@click.option('--V3V4_DB',
              type=click.Choice(L_16S_DATABASE_LIST),
              multiple=True,
              help='BLAST: NCBI_16S, NCBI_NT, Probiotics. UCLUST: RDP. '
                   '중복 사용 가능. ex) --V3V4_DB NCBI_16S --V3V4_DB NCBI_Probiotics')
# @click.option('--ITS34_DB',
#               type=click.Choice(L_ITS_DATABASE_LIST),
#               multiple=True,
#               help='BLAST: NCBI_NT. UCLUST: UNITE. '
#                    '중복 사용 가능. ex) --ITS34_DB NCBI_NT --ITS34_DB UNITE')
@click.option('--remove_taxon', '-rm',
              type=click.Choice(['no_hit', 'filtered', 'uncultured', 'environmental', 'bacterium',
                                 'marine', 'unidentified', 'eukaryote', 'all']),
              multiple=True,
              help='taxon 제거 조건. 다중 선택 가능. ex) -rm no_hit -rm filtered')
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
              help=f'Queue Job 당 OTU 대표 서열의 개수.')
@click.option('--nt_max_job', '-nmj',
              default=BLAST_NT_MAX_JOB,
              show_default=True,
              help='NT DB일 경우 최대 생성할 수 있는 Queue Job의 개수. '
                   '분석 시간이 짧은 다른 작업이 진행되지 않는 문제 방지 목적. 변경시 문의 요망.')
@click.option('--queue', '-q',
              default=BLAST_QUEUE,
              show_default=True,
              type=click.Choice(QUEUE_LIST),
              help='작업을 진행할 Queue.')
# @click.option('--theCups_cmd', '-cups',
#               default='MiSeq_V1',
#               show_default=True,
#               help='theCups command.')
# @click.option('--no_theCups', '-noCups',
#               is_flag=True,
#               help='theCups 실행 안 함.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def espresso(**kargs):
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

    \b
      ___  ___ _ __  _ __ ___  ___ ___  ___
     / _ \\/ __| '_ \\| '__/ _ \\/ __/ __|/ _ \\
    |  __/\\__ \\ |_) | | |  __/\\__ \\__ \\ (_) |
     \\___||___/ .__/|_|  \\___||___/___/\\___/
              | |
              |_|

    \b
    \033[1;32;40m<ASVs PipeLine>\033[m
    PreMA CORE --> CoFI R_DADA2    -->  theCups
                        ALIGNMENT
                        PHYLOGENY
                        TAXONOMY
                        BIOM
    """
    pprint(kargs)
    start_time = time.time()
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass

    if any([kargs['v1v2_db'], kargs['v3v4_db'], kargs['its34_db']]) is False:
        click.secho('Error: --V1V2_DB, --V3V4_DB 또는 --ITS34_DB 옵션을 입력하세요.', fg='red')
        exit()

    # PreMA - Data Copy & QC
    if kargs['integrate']:
        integrate_text = ''
        for ele in kargs['integrate']:
            integrate_text += '-i {} '.format(ele)
    else:
        integrate_text = ''

    if kargs['custom_sample_name']:
        custom_sample_name_text = f'-csn {kargs["custom_sample_name"]}'
    else:
        custom_sample_name_text = ''

    if kargs['no_prema'] is False:
        global PREMA
        prema_cmd = \
            'python3 {prema} ' \
            'CORE ' \
            '{integrate} ' \
            '{custom_sample_name} ' \
            '{copy} ' \
            '{no_order_number}' \
            '{order_number} '.format(
                prema=PREMA,
                integrate=integrate_text,
                custom_sample_name=custom_sample_name_text,
                copy='-c' if kargs['copy_rawdata'] else '',
                no_order_number='--no_order_number' if kargs['no_order_number'] else '',
                order_number=kargs['order_number'])
        run = run_cmd(prema_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'PreMA CORE 완료',
                'false_meg': 'PreMA CORE',
            }, p_exit=True)
        del prema_cmd, run
    prema_end_time = time.time()

    # Read Assembly & CD-HIT-OTU : FLASH_HIT 등
    if kargs['copy_rawdata'] or kargs['no_prema']:
        # 다중 옵션 처리 : remove_taxon
        if len(kargs['remove_taxon']) != 0:
            remove_taxon = ' '.join([f'-rm {x}' for x in kargs['remove_taxon']])
        else:
            remove_taxon = None
        # 여러 영역을 진행하지 않은 경우
        if kargs['target_region'] != 'both':
            l_db_name_tool = get_db_and_tool(kargs['target_region'], kargs['v1v2_db'],
                                             kargs['v3v4_db'], kargs['its34_db'])
            # Suffix 미입력.
            if (kargs['v1v2'] is None) and (kargs['v3v4'] is None) and (kargs['its34'] is None):
                pipeline_time = run_pipeline(
                    {
                        'base_path': CONFIG['base_path']['analysis_base_path'],
                        'order_number': kargs['order_number'],
                        'target_region': kargs['target_region'],
                        'l_db_name_tool': l_db_name_tool,
                        'sample_suffix': '',
                        'no_use_flash_worker': kargs['no_use_flash_worker'],
                        'custom_sample_name': kargs['custom_sample_name'],
                        'remove_taxon': remove_taxon,
                        'query_coverage': kargs['query_coverage'],
                        'identity_percentage': kargs['identity_percentage'],
                        'read_per_job': kargs['read_per_job'],
                        'nt_max_job': kargs['nt_max_job'],
                        'queue': kargs['queue'],
                        'theCups_cmd': kargs['thecups_cmd'],
                        'no_theCups': kargs['no_thecups'],
                        'no_order_number': kargs['no_order_number']
                    })
            else:
                if kargs['v1v2'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["v1v2"]}'
                elif kargs['v3v4'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["v3v4"]}'
                elif kargs['its34'] is not None:
                    sample_suffix = f'--sample_suffix {kargs["its34"]}'
                else:
                    raise ValueError(f'알 수 없는 에러가 발생하였습니다.')
                pipeline_time = run_pipeline(
                    {
                        'base_path': CONFIG['base_path']['analysis_base_path'],
                        'order_number': kargs['order_number'],
                        'target_region': kargs['target_region'],
                        'l_db_name_tool': l_db_name_tool,
                        'sample_suffix': sample_suffix,
                        'no_use_flash_worker': kargs['no_use_flash_worker'],
                        'custom_sample_name': kargs['custom_sample_name'],
                        'remove_taxon': remove_taxon,
                        'query_coverage': kargs['query_coverage'],
                        'identity_percentage': kargs['identity_percentage'],
                        'read_per_job': kargs['read_per_job'],
                        'nt_max_job': kargs['nt_max_job'],
                        'queue': kargs['queue'],
                        'theCups_cmd': kargs['thecups_cmd'],
                        'no_theCups': kargs['no_thecups'],
                    })

        elif kargs['target_region'] == 'both':
            # V1V2 V3V4 영역을 동시 진행한 경우
            if (kargs['v1v2'] is not None) and (kargs['v3v4'] is not None) and (kargs['its34'] is None):
                pipeline_time = list()
                for region in ['v1v2', 'v3v4']:
                    if region == 'v1v2':
                        target_region = '27F-Bact338'
                    elif region == 'v3v4':
                        target_region = 'Bakt_341F-805R'
                    else:
                        raise ValueError()

                    l_db_name_tool = get_db_and_tool(target_region, kargs['v1v2_db'],
                                                     kargs['v3v4_db'], kargs['its34_db'])
                    pipeline_time.append(
                        run_pipeline({
                            'base_path': CONFIG['base_path']['analysis_base_path'],
                            'order_number': kargs['order_number'],
                            'target_region': target_region,
                            'l_db_name_tool': l_db_name_tool,
                            'sample_suffix': f'--sample_suffix {kargs[region]}',
                            'no_use_flash_worker': kargs['no_use_flash_worker'],
                            'custom_sample_name': kargs['custom_sample_name'],
                            'remove_taxon': remove_taxon,
                            'query_coverage': kargs['query_coverage'],
                            'identity_percentage': kargs['identity_percentage'],
                            'read_per_job': kargs['read_per_job'],
                            'nt_max_job': kargs['nt_max_job'],
                            'queue': kargs['queue'],
                            'theCups_cmd': kargs['thecups_cmd'],
                            'no_theCups': kargs['no_thecups'],
                        })
                    )

            # V3V4 ITS3-4 역영을 동시 진행한 경우
            elif (kargs['v1v2'] is None) and (kargs['v3v4'] is not None) and (kargs['its34'] is not None):
                pipeline_time = list()
                for region in ['v3v4', 'its34']:
                    if region == 'v3v4':
                        target_region = 'Bakt_341F-805R'
                    elif region == 'its34':
                        target_region = 'ITS_3F-4R'
                    else:
                        raise ValueError()

                    l_db_name_tool = get_db_and_tool(target_region, kargs['v1v2_db'],
                                                     kargs['v3v4_db'], kargs['its34_db'])
                    pipeline_time.append(
                        run_pipeline({
                            'base_path': CONFIG['base_path']['analysis_base_path'],
                            'order_number': kargs['order_number'],
                            'target_region': target_region,
                            'l_db_name_tool': l_db_name_tool,
                            'sample_suffix': f'--sample_suffix {kargs[region]}',
                            'no_use_flash_worker': kargs['no_use_flash_worker'],
                            'custom_sample_name': kargs['custom_sample_name'],
                            'remove_taxon': remove_taxon,
                            'query_coverage': kargs['query_coverage'],
                            'identity_percentage': kargs['identity_percentage'],
                            'read_per_job': kargs['read_per_job'],
                            'nt_max_job': kargs['nt_max_job'],
                            'queue': kargs['queue'],
                            'theCups_cmd': kargs['thecups_cmd'],
                            'no_theCups': kargs['no_thecups'],
                        })
                    )

        cofi_end_time = time.time()
        print_metamix_run_time(start_time, prema_end_time, pipeline_time, cofi_end_time)


@metamix.command('DoubleShot', short_help='Closed-Ref. method OTU Picking: HiSeq & MiSeq')
@click.version_option(version=__version__)
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_prema',
              is_flag=True,
              help='데이터 복사 없이(PreMA 실행없이, 이미 복사되었을 경우) 이후 단계를 진행할 때 사용.')
@click.option('--V1V2',
              help='V1V2(27F-Bact338) 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--V3V4',
              help='V3V4(Bakt_341F-805R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--ITS34',
              help='ITS3-4(ITS_3F-4R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--target_region', '-r',
              required=True,
              type=click.Choice(['Bakt_341F-805R', 'ITS_3F-4R', '27F-Bact338', 'both']),
              help='Target Region 및 Primer Set.')
@click.option('--blast_db', '-bdb',
              is_flag=True,
              help='BLAST DB. 영역에 따라 자동 선택. 현재 사용 불가능.')
@click.option('--uclust_db', '-udb',
              is_flag=True,
              help='UCLUST DB 사용. 영역에 따라 자동 선택. RDP: 16S, UNITE: ITS')
@click.option('--no_use_FLASH_worker', '-noW',
              is_flag=True,
              help='FLASH 실행시 Worker를 사용하지 않음. 시료가 많아 Error 발생하는 경우 사용.')
def double_shot(**kargs):
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

    \b
     ____              _     _      ____  _           _
    |  _ \  ___  _   _| |__ | | ___/ ___|| |__   ___ | |_
    | | | |/ _ \| | | | '_ \| |/ _ \___ \| '_ \ / _ \| __|
    | |_| | (_) | |_| | |_) | |  __/___) | | | | (_) | |_
    |____/ \___/ \__,_|_.__/|_|\___|____/|_| |_|\___/ \__|

    \b
    \033[1;32;40m<Closed-Reference OTU Picking PipeLine>\033[m
    HiSeq & MiSeq
    PreMA CORE --> CoFI CLOSED
    """
    pprint(kargs)


@metamix.command('LongShot', short_help='Closed-Ref. method OTU Picking: PacBio RSII & Sequel')
@click.version_option(version=__version__)
@click.argument('order_number')
@click.option('--integrate', '-i',
              multiple=True,
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              help='통합분석할 수주번호. ex) -i 수주번호 -i 수주번호...')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_prema',
              is_flag=True,
              help='데이터 복사 없이(PreMA 실행없이, 이미 복사되었을 경우) 이후 단계를 진행할 때 사용.')
@click.option('--V1V2',
              help='V1V2(27F-Bact338) 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--V3V4',
              help='V3V4(Bakt_341F-805R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--ITS34',
              help='ITS3-4(ITS_3F-4R 시료을 구분하는 접미사. 한 수주에 여러 영역을 진행하는 경우에만 적용.')
@click.option('--target_region', '-r',
              required=True,
              type=click.Choice(['Bakt_341F-805R', 'ITS_3F-4R', '27F-Bact338', 'both']),
              help='Target Region 및 Primer Set.')
@click.option('--blast_db', '-bdb',
              is_flag=True,
              help='BLAST DB. 영역에 따라 자동 선택. 현재 사용 불가능.')
@click.option('--uclust_db', '-udb',
              is_flag=True,
              help='UCLUST DB 사용. 영역에 따라 자동 선택. RDP: 16S, UNITE: ITS')
@click.option('--no_use_FLASH_worker', '-noW',
              is_flag=True,
              help='FLASH 실행시 Worker를 사용하지 않음. 시료가 많아 Error 발생하는 경우 사용.')
def long_shot(**kargs):
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

    \b
     _                      ____  _           _
    | |    ___  _ __   __ _/ ___|| |__   ___ | |_
    | |   / _ \| '_ \ / _` \___ \| '_ \ / _ \| __|
    | |__| (_) | | | | (_| |___) | | | | (_) | |_
    |_____\___/|_| |_|\__, |____/|_| |_|\___/ \__|
                      |___/
    \b
    \033[1;32;40m<Closed-Reference OTU Picking PipeLine>\033[m
    PacBio RSII & Sequel
    PreMA CORE --> CoFI CLOSED
    """
    pprint(kargs)


def run_pipeline(kargs):
    """
    Amplicon Metagenome Analysis 를 분석 파이프라인 절차에 따라 진행한다.

    :param kargs:
                base_path:
                order_number:
                target_region:
                sample_suffix:
                l_db_name_tool: tuple을 원소로 가지는 리스트. [(db_name, db_tool), (db_name_ db_tool)]
                # db_tool:
                # db_name:
                no_use_flash_worker:
                custom_sample_name: 시료명 변경과 그룹 정보. 파일 지정시 Biom 생성 및 Report 생성 중지.
                remove_taxon:
                query_coverage:
                identity_percentage:
                read_per_job:
                nt_max_job:
                queue:
                theCups_cmd: MiSeq_V1, MiSeq_V2
                no_theCups:
                no_order_number
    :return:
    """
    # FLASH & CD-HIT-OTU
    start_time = time.time()
    run_cofi_flash_hit({
        'order_number': kargs['order_number'],
        'target_region': kargs['target_region'],
        'sample_suffix': kargs['sample_suffix'],
        'no_use_flash_worker': kargs['no_use_flash_worker'],
        'no_order_number': kargs['no_order_number'],
    })
    flash_hit_end_time = time.time()

    # Multiple Alignment
    order_path = os.path.join(CONFIG['base_path']['analysis_base_path'], kargs['order_number'])
    analysis_dir_path = get_work_path(order_path, 'Analysis_?')
    analysis_dir_name = os.path.basename(analysis_dir_path)
    # cd_hit_dir = os.path.join(analysis_dir_path, 'CD-HIT-OTU')
    # cluster_dir = os.path.join(cd_hit_dir, f'{kargs["order_number"]}*')
    run_cofi_alignment(kargs['order_number'], analysis_dir_name, kargs['no_order_number'])
    alignment_end_time = time.time()

    # Phylogeny
    run_cofi_phylogeny(kargs['order_number'], analysis_dir_name, kargs['no_order_number'])
    phylogeny_end_time = time.time()

    # Taxonomy Assignment
    for db_name, db_tool in kargs['l_db_name_tool']:
        run_cofi_taxonomy(
            {
                'order_number': kargs['order_number'],
                'analysis_number': analysis_dir_name,
                'db_tool': db_tool,
                'db_name': db_name,
                'remove_taxon': kargs['remove_taxon'],
                'query_coverage': kargs['query_coverage'],
                'identity_percentage': kargs['identity_percentage'],
                'read_per_job': kargs['read_per_job'],
                'nt_max_job': kargs['nt_max_job'],
                'queue': kargs['queue'],
                'no_order_number': kargs['no_order_number']
            }
        )
    taxonomy_end_time = time.time()

    # BIOM
    run_cofi_biom(kargs['order_number'], analysis_dir_name, kargs['no_order_number'])
    biom_end_time = time.time()

    # theCups
    if kargs['no_theCups'] is True:
        thecups_end_time = time.time()
    else:
        run_thecups(kargs['theCups_cmd'], kargs['order_number'], analysis_dir_name, kargs['no_order_number'])
        thecups_end_time = time.time()
    return {
                'flash_hit_start_time': start_time,
                'flash_hit_end_time': flash_hit_end_time,
                'alignment_end_time': alignment_end_time,
                'phylogeny_end_time': phylogeny_end_time,
                'taxonomy_end_time': taxonomy_end_time,
                'biom_end_time': biom_end_time,
                'theCups_end_time': thecups_end_time,
            }


def get_work_path(p_path, p_pattern):
    l_analysis_dir, _ = check_file_type(glob_dir(p_path, p_pattern, p_mode='many', p_verbose=False), 'isdir')
    if len(l_analysis_dir) == 1:
        analysis_dir = l_analysis_dir[0]
    else:
        l_analysis_dir.sort()
        analysis_dir = l_analysis_dir[-1]
    return analysis_dir


def get_db_and_tool(p_region, p_v1v2_db, p_v3v4_db, p_its34_db):
    l_db_name_tool = list()
    if p_region == 'Bakt_341F-805R':
        if p_v3v4_db:
            for db_name in p_v3v4_db:
                if 'NCBI' in db_name:
                    db_tool = 'blast'
                else:
                    db_tool = 'uclust'
                l_db_name_tool.append((db_name, db_tool))
        else:
            click.secho('Error: 목표영역에 DB가 설정되지 않았습니다.', fg='red')
            click.echo(f'--V3V4_DB: {p_v3v4_db}')
            exit()
    elif p_region == 'ITS_3F-4R':
        if p_its34_db:
            for db_name in p_its34_db:
                if 'NCBI' in db_name:
                    db_tool = 'blast'
                else:
                    db_tool = 'uclust'
                l_db_name_tool.append((db_name, db_tool))
        else:
            click.secho('Error: 목표영역에 DB가 설정되지 않았습니다.', fg='red')
            click.echo(f'--ITS34_DB: {p_v3v4_db}')
            exit()
    elif p_region == '27F-Bact338':
        if p_v1v2_db:
            for db_name in p_v1v2_db:
                if 'NCBI' in db_name:
                    db_tool = 'blast'
                else:
                    db_tool = 'uclust'
                l_db_name_tool.append((db_name, db_tool))
        else:
            click.secho('Error: 목표영역에 DB가 설정되지 않았습니다.', fg='red')
            click.echo(f'--V1V2_DB: {p_v3v4_db}')
            exit()
    return l_db_name_tool


def get_cofi_flash_hit_cmd(kargs):
    """
    CoFI FLASH_HIT 실행을 위한 명령어를 반환한다.

    :param kargs:
                target_region:
                sample_suffix:
                no_use_flash_worker:
                order_number:
    :return:
    """
    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'FLASH_HIT ' \
        '--target_region {target} ' \
        '{sample_suffix}' \
        '{no_order_number}' \
        '{no_worker}' \
        '{order_number}'.format(
            cofi=COFI,
            target=kargs['target_region'],
            sample_suffix=kargs['sample_suffix'],
            no_order_number='--no_order_number ' if kargs['no_order_number'] else ' ',
            no_worker='--no_use_worker ' if kargs['no_use_flash_worker'] else ' ',
            order_number=kargs['order_number'])
    return cofi_cmd


def run_cofi_flash_hit(kargs):
    """
    CoFI FLASH_HIT 실행을 실행한다.

    :param kargs:
                target_region:
                sample_suffix:
                no_use_flash_worker:
                order_number:
                no_order_number:

    :return:
    """
    cofi_cmd = get_cofi_flash_hit_cmd(kargs)
    run = run_cmd(cofi_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'CoFI FLASH_HIT 완료',
            'false_meg': 'CoFI FLASH_HIT',
        }, p_exit=True)


def get_cofi_flash_hit_cmd_ori(kargs, p_region):
    # Not Use
    if p_region is not None:
        if p_region == 'v3v4':
            target_region = 'Bakt_341F-805R'
        elif p_region == 'its34':
            target_region = 'ITS_3F-4R'
        elif p_region == 'v1v2':
            target_region = '27F-Bact338'
        else:
            raise ValueError(f'target_region: {p_region}')
        sample_suffix = f'--sample_suffix {kargs[p_region]}'
    else:
        target_region = kargs['target_region']
        sample_suffix = ''

    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'FLASH_HIT ' \
        '--target_region {target} ' \
        '{sample_suffix}' \
        '{no_worker} ' \
        '{order_number}'.format(
            cofi=COFI,
            target=target_region,
            sample_suffix=sample_suffix,
            no_worker='--no_use_worker' if kargs['no_use_flash_worker'] else '',
            order_number=kargs['order_number'])
    return cofi_cmd


def get_cofi_alignment_cmd(p_order_number, p_analysis_number, p_no_order_number):
    """
    CoFI ALIGNMENT 실행을 위한 명령어를 반환한다.

    :param p_order_number:
    :param p_analysis_number:
    :param p_no_order_number:
    :return:
    """
    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'ALIGNMENT ' \
        '{no_order_number}' \
        '{order_number} ' \
        '{analysis_number}'.format(
            cofi=COFI,
            no_order_number='--no_order_number ' if p_no_order_number else ' ',
            order_number=p_order_number,
            analysis_number=p_analysis_number,
        )
    return cofi_cmd


def run_cofi_alignment(p_order_number, p_analysis_number, p_no_order_number):
    """
    CoFI ALIGNMENT 을 실행한다.

    :param p_order_number:
    :param p_analysis_number:
    :param p_no_order_number:
    :return:
    """
    global OTUS_REP_FASTA
    cofi_alignment_cmd = get_cofi_alignment_cmd(p_order_number, p_analysis_number, p_no_order_number)
    run = run_cmd(cofi_alignment_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'CoFI ALIGNMENT 완료',
            'false_meg': 'CoFI ALIGNMENT',
        }, p_exit=True)


def get_cofi_phylogeny_cmd(p_order_number, p_analysis_number, p_no_order_number):
    """
    CoFI PHYLOGENY 실행을 위한 명령어를 반환한다.

    :param p_order_number:
    :param p_analysis_number:
    :param p_no_order_number:
    :return:
    """
    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'PHYLOGENY ' \
        '{no_order_number}' \
        '{order_number} ' \
        '{analysis_number}'.format(
            cofi=COFI,
            no_order_number='--no_order_number ' if p_no_order_number else ' ',
            order_number=p_order_number,
            analysis_number=p_analysis_number,
        )
    return cofi_cmd


def run_cofi_phylogeny(p_order_number, p_analysis_number, p_no_order_number):
    """
    CoFI PHYLOGENY 를 실행한다.

    :param p_order_number:
    :param p_analysis_number:
    :param p_no_order_number:
    :return:
    """
    cofi_phylogeny_cmd = get_cofi_phylogeny_cmd(p_order_number, p_analysis_number, p_no_order_number)
    run = run_cmd(cofi_phylogeny_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'CoFI PHYLOGENY 완료',
            'false_meg': 'CoFI PHYLOGENY',
        }, p_exit=True)


def get_cofi_taxonomy_cmd(p_kargs):
    """
    CoFI TAXONOMY 실행을 위한 명령어를 반환한다.

    :param p_kargs: 다음을 key로 가지는 딕션너리
        order_number:
        analysis_number:
        db_tool:
        db_name:
        remove_taxon:
        query_coverage:
        identity_percentage:
        read_per_job:
        nt_max_job:
        queue:
        no_order_number:
    :return:
    """
    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'TAXONOMY ' \
        '{order_number} ' \
        '{analysis_number} ' \
        '{no_order_number}' \
        '--{tool}_db {db_name} ' \
        '{remove_taxon} ' \
        '-qc {query_coverage} ' \
        '-ip {identity_percentage} ' \
        '-rj {read_per_job} ' \
        '-nmj {nt_max_job} ' \
        '-q {queue}'.format(
            cofi=COFI,
            order_number=p_kargs['order_number'],
            analysis_number=p_kargs['analysis_number'],
            no_order_number='--no_order_number ' if p_kargs['no_order_number'] else ' ',
            tool=p_kargs['db_tool'],
            db_name=p_kargs['db_name'],
            remove_taxon='' if p_kargs['remove_taxon'] is None else p_kargs['remove_taxon'],
            query_coverage=p_kargs['query_coverage'],
            identity_percentage=p_kargs['identity_percentage'],
            read_per_job=p_kargs['read_per_job'],
            nt_max_job=p_kargs['nt_max_job'],
            queue=p_kargs['queue'])
    return cofi_cmd


def run_cofi_taxonomy(p_kargs):
    """
    CoFI TAXONOMY 을 실행한다.

     :param p_kargs: 다음을 key로 가지는 딕션너리
        order_number:
        analysis_number:
        db_tool:
        db_name:
        remove_taxon:
        query_coverage:
        identity_percentage:
        read_per_job:
        nt_max_job:
        queue:
        no_order_number:
    :return:
    """
    cofi_taxonomy_cmd = get_cofi_taxonomy_cmd(p_kargs)
    run = run_cmd(cofi_taxonomy_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': f'CoFI TAXONOMY 완료({p_kargs["db_tool"]}, {p_kargs["db_name"]})',
            'false_meg': 'CoFI TAXONOMY'
        }, p_exit=True)


def get_cofi_biom_cmd(p_order_number, p_analysis_number, p_no_order_number):
    global COFI
    cofi_cmd = \
        'python3 {cofi} ' \
        'BIOM ' \
        '{no_order_number}' \
        '{order_number} ' \
        '{analysis_number}'.format(
            cofi=COFI,
            no_order_number='--no_order_number ' if p_no_order_number else ' ',
            order_number=p_order_number,
            analysis_number=p_analysis_number)
    return cofi_cmd


def run_cofi_biom(p_order_number, p_analysis_number, p_no_order_number):
    cofi_biom_cmd = get_cofi_biom_cmd(p_order_number, p_analysis_number, p_no_order_number)
    run = run_cmd(cofi_biom_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'CoFI BIOM 완료',
            'false_meg': 'CoFI BIOM'
        }, p_exit=True)


def get_thecups_cmd(p_thecups_cmd, p_order_number, p_anslysis_number, p_no_order_number):
    global THECUPS
    if p_no_order_number is True:
        no_order_number = '--no_order_number '
    else:
        no_order_number = ''
    thecups_cmd = \
        f'python3 {THECUPS} ' \
        f'{p_thecups_cmd} ' \
        f'{p_order_number} ' \
        f'{p_anslysis_number} ' \
        f'{no_order_number}' \
        f'--drink'
    return thecups_cmd


def run_thecups(p_thecups_cmd, p_order_number, p_analsyis_number, p_no_order_number):
    thecups_cmd = get_thecups_cmd(p_thecups_cmd, p_order_number, p_analsyis_number, p_no_order_number)
    run = run_cmd(thecups_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
    check_run_cmd(
        {
            'run': run,
            'true_meg': f'theCups {p_thecups_cmd} 완료',
            'false_meg': f'theCups {p_thecups_cmd}',
        }, p_exit=True)


if __name__ == '__main__':
    metamix()
