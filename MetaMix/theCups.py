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
#  _   _       _____
# | |_| |_ ___|     |_ _ ___ ___
# |  _|   | -_|   --| | | . |_ -|
# |_| |_|_|___|_____|___|  _|___|
#                       |_|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.1.0'

# ------------------------------------
# 1.0.2 : 2020.04.24
# 2D_PCoA : y축 제목 짤림 현상 해결
#   - make_2d_plot_size.py 생성
#   - theCups 옵션 추가
# ------------------------------------
# Ver. 1.0.3 : 2020.04.28
# Probiotics21 DB 추가 --> 보고서 추가
# ------------------------------------
# Ver. 1.0.4 : 2020.05.19
# --no_order_number 추가
# --no_3d_pcoa 추가
# ------------------------------------
# Ver. 1.1.0 : 2020.06.04
# ASVs 보고서 추가

import click
from pprint import pprint
from SpoON.util import parse_config, check_order_number_system


# 기본값 설정
CONFIG = parse_config()
ANALYSIS_BASE_PATH = CONFIG['base_path']['analysis_base_path']
TABLE_PATH_SUFFIX_19 = CONFIG['theCups']['Probiotics19']['table_path_suffix']
TABLE_PATH_SUFFIX_21 = CONFIG['theCups']['Probiotics21']['table_path_suffix']
V1V2_REPORT_DIR = CONFIG['theCups']['Probiotics']['V1V2_report_dir']
V1V2_REPORT_NUM = CONFIG['theCups']['Probiotics']['V1V2_report_num']
V3V4_REPORT_DIR = CONFIG['theCups']['Probiotics']['V3V4_report_dir']
V3V4_REPORT_NUM = CONFIG['theCups']['Probiotics']['V3V4_report_num']

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
     _   _       _____
    | |_| |_ ___|     |_ _ ___ ___
    |  _|   | -_|   --| | | . |_ -|
    |_| |_|_|___|_____|___|  _|___|
                          |_|

    작성자 : 박정원(JungWon Park, KOSST)
    """
    pass


@main.command('MiSeq_V1', short_help='MiSeq Report Version 1. (QIIME 1.9)')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--report_number', '-rn',
              help='Report_? 디렉터리명')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--drink',
              is_flag=True,
              help='모든 결과에 대한 분석 보고서를 작성.')
@click.option('--read_assembly',
              is_flag=True,
              help='STAT.html')
@click.option('--summary',
              is_flag=True,
              help='OTU Picking Summary. Summary.html')
@click.option('--alpha_rarefaction',
              is_flag=True,
              help='')
@click.option('--max_rare_depth',
              default=None,
              type=click.INT,
              help='--alpha_rarefaction 옵션과 함께 사용. the upper limit of rarefaction depths')
@click.option('--beta_diversity_3d',
              is_flag=True,
              help='')
@click.option('--no_3d_pcoa',
              is_flag=True,
              help='3D PCoA 생성 안 함. 3D PCoA 생성이 안 되는 경우인데 2D PCoA를 생성하거나 2D PCoA 결과가 있을 경우 '
                   'PC Value 페이지 생성이 필요할 경우 사용.')
@click.option('--pcoa_2d',
              is_flag=True,
              help='beta diversity 데이터가 있어야 결과가 생성됩니다.')
@click.option('--pcoa_2d_x_size',
              type=click.FLOAT,
              help='2D PCoA Y축 제목이 짤릴 경우 사용. 그래프의 X축 크기 조절(기본값 4.5).')
@click.option('--pcoa_2d_y_size',
              type=click.FLOAT,
              help='2D PCoA Y축 제목이 짤릴 경우 사용. 그래프의 Y축 크기 조절(기본값 4.5).')
@click.option('--diversity_index',
              is_flag=True,
              help='')
@click.option('--taxonomy_assignment',
              is_flag=True,
              help='가장 긴 시료명 길이에 따라 그래프의 y축 크기 결정(max, max+1 길이의 그래프 생성). '
                   'max(시료명 길이) < 6 : 6. 시료의 개수 > 30: 6.')
@click.option('--taxonomy_y_size', '-y',
              type=click.INT,
              multiple=True,
              help='Taxonomy Assignment의 결과 중 Bar & Area Plot의 y축의 크기 조절. 시료명이 잘리는 경우 사용. '
                   '다중 선택 가능. ex) -y 6 -y 7')
@click.option('--upgma_tree',
              is_flag=True,
              help='')
@click.option('--main',
              is_flag=True,
              help='OTUanalysis.html 생성.')
@click.option('--taste',
              is_flag=True,
              help='생성된 분석 보고서의 이상 여부를 확인.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def miseq_report_1(**kargs):
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
     _   _       _____
    | |_| |_ ___|     |_ _ ___ ___
    |  _|   | -_|   --| | | . |_ -|
    |_| |_|_|___|_____|___|  _|___|
                          |_|

    작성자 : 박정원(JungWon Park, KOSST)
    """
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['drink'] is True) and (kargs['report_number'] is not None):
        click.secho('Error: --drink 와 --report_number 옵션을 동시에 사용하지 못합니다.', fg='red')
        exit()
    from theCups.core import miseq_report_pipeline_v1
    miseq_report_pipeline_v1(kargs)


@main.command('ASVs_V1', short_help='ASVs Report Version 1. (DADA2 + QIIME 1.9)')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--report_number', '-rn',
              help='Report_? 디렉터리명')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--drink',
              is_flag=True,
              help='모든 결과에 대한 분석 보고서를 작성.')
@click.option('--summary',
              is_flag=True,
              help='OTU Picking Summary. Summary.html')
@click.option('--alpha_rarefaction',
              is_flag=True,
              help='')
@click.option('--max_rare_depth',
              default=None,
              type=click.INT,
              help='--alpha_rarefaction 옵션과 함께 사용. the upper limit of rarefaction depths')
@click.option('--beta_diversity_3d',
              is_flag=True,
              help='')
@click.option('--no_3d_pcoa',
              is_flag=True,
              help='3D PCoA 생성 안 함. 3D PCoA 생성이 안 되는 경우인데 2D PCoA를 생성하거나 2D PCoA 결과가 있을 경우 '
                   'PC Value 페이지 생성이 필요할 경우 사용.')
@click.option('--pcoa_2d',
              is_flag=True,
              help='beta diversity 데이터가 있어야 결과가 생성됩니다.')
@click.option('--pcoa_2d_x_size',
              type=click.FLOAT,
              help='2D PCoA Y축 제목이 짤릴 경우 사용. 그래프의 X축 크기 조절(기본값 4.5).')
@click.option('--pcoa_2d_y_size',
              type=click.FLOAT,
              help='2D PCoA Y축 제목이 짤릴 경우 사용. 그래프의 Y축 크기 조절(기본값 4.5).')
@click.option('--diversity_index',
              is_flag=True,
              help='')
@click.option('--taxonomy_assignment',
              is_flag=True,
              help='가장 긴 시료명 길이에 따라 그래프의 y축 크기 결정(max, max+1 길이의 그래프 생성). '
                   'max(시료명 길이) < 6 : 6. 시료의 개수 > 30: 6.')
@click.option('--taxonomy_y_size', '-y',
              type=click.INT,
              multiple=True,
              help='Taxonomy Assignment의 결과 중 Bar & Area Plot의 y축의 크기 조절. 시료명이 잘리는 경우 사용. '
                   '다중 선택 가능. ex) -y 6 -y 7')
@click.option('--upgma_tree',
              is_flag=True,
              help='')
@click.option('--main',
              is_flag=True,
              help='OTUanalysis.html 생성.')
@click.option('--taste',
              is_flag=True,
              help='생성된 분석 보고서의 이상 여부를 확인.')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def asvs_report_1(**kargs):
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if (kargs['drink'] is True) and (kargs['report_number'] is not None):
        click.secho('Error: --drink 와 --report_number 옵션을 동시에 사용하지 못합니다.', fg='red')
        exit()
    from theCups.core import asvs_report_pipeline_v1
    asvs_report_pipeline_v1(kargs)


@main.command('MiSeq_V2', short_help='MiSeq Report Version 2. (QIIME 1.9 + Alpha)')
@click.argument('order_number')
@click.argument('analysis_number')
@click.option('--report_number', '-rn',
              help='Report_? 디렉터리명')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True),
              help='기본값: {}'.format(ANALYSIS_BASE_PATH))
@click.option('--drink',
              is_flag=True,
              help='모든 결과에 대한 분석 보고서를 작성.')
@click.option('--read_assembly',
              is_flag=True,
              help='STAT.html')
@click.option('--length_filter',
              is_flag=True,
              help='Length_Filter.html')
@click.option('--clustering',
              is_flag=True,
              help='Clustering_OTU_Picking.html')
def miseq_report_2(**kargs):
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
     _   _       _____
    | |_| |_ ___|     |_ _ ___ ___
    |  _|   | -_|   --| | | . |_ -|
    |_| |_|_|___|_____|___|  _|___|
                          |_|

    작성자 : 박정원(JungWon Park, KOSST)
    """
    # 수주번호 형식 확인
    if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
        pass
    if (kargs['drink'] is True) and (kargs['report_number'] is not None):
        click.secho('Error: --drink 와 --report_number 옵션을 동시에 사용하지 못합니다.', fg='red')
        exit()
    from theCups.core import miseq_report_pipeline_v2
    miseq_report_pipeline_v2(kargs)


@main.command('PROBIOTICS', short_help='식약처 제출용 보고서 19종')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--v1v2_dir', '-12',
              default=V1V2_REPORT_DIR,
              show_default=True,
              help='V1V2 영역의 보고서의 디렉터리명.')
@click.option('--v1v2_num', '-12n',
              default=V1V2_REPORT_NUM,
              show_default=True,
              help='V1V2 영역의 보고서 번호.')
@click.option('--v3v4_dir', '-34',
              default=V3V4_REPORT_DIR,
              show_default=True,
              help='V3V4 영역의 보고서 디렉터리명.')
@click.option('--v3v4_num', '-34n',
              default=V3V4_REPORT_NUM,
              show_default=True,
              help='V3V4 영역의 보고서 번호.')
@click.option('--table_path_suffix',
              default=TABLE_PATH_SUFFIX_19,
              show_default=True)
@click.option('--mixed_species', '-m',
              help='투입된 유산균 종의 목록 파일. '
                   '해당 옵션을 사용하지 않으면 모든 시료에 동일한 목록 적용(프로그램 실행시 입력 대화창 실행).'
                   '투입된 유산균 종의 목록 파일은 탭으로 구분. '
                   '제품명이 없을 경우 \'None\'으로 표기. '
                   '\'#\'문자로 시작하면 해당 열은 제외.'
                   '시로명 입력시 v12, V1V2, v34, v3v4 등 목표 영역을 나타내는 접미사는 제거하고 입력한다. '
                   'ex1) A1 유산균제품A 1 2 10 ex2) B1 None 1 2')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def probiotics(**kargs):
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
     _   _       _____
    | |_| |_ ___|     |_ _ ___ ___
    |  _|   | -_|   --| | | . |_ -|
    |_| |_|_|___|_____|___|  _|___|
                          |_|

    작성자 : 박정원(JungWon Park, KOSST)

    \b
    식약처 제출용 Probiotics 엑셀 분석 보고서를 생성합니다.
    생성된 엑세파일을 이용하여 PDF 보고서를 생성할 수 있습니다.
    식약처 제출용 보고서는 목표 영역이 V1V2, V3V4 또는 V1V2 V3V4인 경우에만 생성할 수 있습니다.
    info.txt 파일은 V1V2 또는 V3V4 데이터 분석시 사용한 파일을 하나만 지정하면 됩니다.
    info.txt 파일에 기재된 정보를 바탕으로 고객 정보를 출력합니다.
    \b
    제품에 투입된 유산균 종 목록을 입력하는 2가지 방법
    1. 프로그램 실행시 선택 창에서 입력 : 모든 시료에 대하여 동일하게 적용.
    2. --mixed_species(-m) 옵션을 이용하여 파일 형태로 입력 : 시료별로 입력 가능.
    """
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if kargs['v1v2_dir'] == 'None':
        kargs['v1v2_dir'] = None
    if kargs['v1v2_num'] == 'None':
        kargs['v1v2_num'] = None
    if kargs['v3v4_dir'] == 'None':
        kargs['v3v4_dir'] = None
    if kargs['v3v4_num'] == 'None':
        kargs['v3v4_num'] = None
    from theCups.core import probiotics_report
    probiotics_report(kargs, 'Probiotics19')


@main.command('PROBIOTICS21', short_help='식약처 제출용 보고서 21종(19종+2종)')
@click.argument('order_number')
@click.option('--analysis_base_path', '-ap',
              default=ANALYSIS_BASE_PATH,
              show_default=True,
              type=click.Path(exists=True))
@click.option('--v1v2_dir', '-12',
              default=V1V2_REPORT_DIR,
              show_default=True,
              help='V1V2 영역의 보고서의 디렉터리명.')
@click.option('--v1v2_num', '-12n',
              default=V1V2_REPORT_NUM,
              show_default=True,
              help='V1V2 영역의 보고서 번호.')
@click.option('--v3v4_dir', '-34',
              default=V3V4_REPORT_DIR,
              show_default=True,
              help='V3V4 영역의 보고서 디렉터리명.')
@click.option('--v3v4_num', '-34n',
              default=V3V4_REPORT_NUM,
              show_default=True,
              help='V3V4 영역의 보고서 번호.')
@click.option('--table_path_suffix',
              default=TABLE_PATH_SUFFIX_21,
              show_default=True)
@click.option('--mixed_species', '-m',
              help='투입된 유산균 종의 목록 파일. '
                   '해당 옵션을 사용하지 않으면 모든 시료에 동일한 목록 적용(프로그램 실행시 입력 대화창 실행).'
                   '투입된 유산균 종의 목록 파일은 탭으로 구분. '
                   '제품명이 없을 경우 \'None\'으로 표기. '
                   '\'#\'문자로 시작하면 해당 열은 제외.'
                   '시로명 입력시 v12, V1V2, v34, v3v4 등 목표 영역을 나타내는 접미사는 제거하고 입력한다. '
                   'ex1) A1 유산균제품A 1 2 10 ex2) B1 None 1 2')
@click.option('--no_order_number',
              is_flag=True,
              help='수주번호가 아닌 디렉터리를 실행할 경우.')
def probiotics21(**kargs):
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
     _   _       _____
    | |_| |_ ___|     |_ _ ___ ___
    |  _|   | -_|   --| | | . |_ -|
    |_| |_|_|___|_____|___|  _|___|
                          |_|

    작성자 : 박정원(JungWon Park, KOSST)

    \b
    식약처 제출용 Probiotics21(19종+2종) 엑셀 분석 보고서를 생성합니다.
    생성된 엑세파일을 이용하여 PDF 보고서를 생성할 수 있습니다.
    식약처 제출용 보고서는 목표 영역이 V1V2, V3V4 또는 V1V2 V3V4인 경우에만 생성할 수 있습니다.
    info.txt 파일은 V1V2 또는 V3V4 데이터 분석시 사용한 파일을 하나만 지정하면 됩니다.
    info.txt 파일에 기재된 정보를 바탕으로 고객 정보를 출력합니다.
    \b
    제품에 투입된 유산균 종 목록을 입력하는 2가지 방법
    1. 프로그램 실행시 선택 창에서 입력 : 모든 시료에 대하여 동일하게 적용.
    2. --mixed_species(-m) 옵션을 이용하여 파일 형태로 입력 : 시료별로 입력 가능.
    """
    pprint(kargs)
    # 수주번호 형식 확인
    if kargs['no_order_number'] is True:
        pass
    else:
        if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
            pass
    if kargs['v1v2_dir'] == 'None':
        kargs['v1v2_dir'] = None
    if kargs['v1v2_num'] == 'None':
        kargs['v1v2_num'] = None
    if kargs['v3v4_dir'] == 'None':
        kargs['v3v4_dir'] = None
    if kargs['v3v4_num'] == 'None':
        kargs['v3v4_num'] = None
    from theCups.core import probiotics_report
    probiotics_report(kargs, 'Probiotics21')


@main.command('SEQUEL', short_help='Sequel Report Version 1. (QIIME 1.9 + Alpha)')
@click.argument('order_number')
@click.argument('analysis_number')
def sequel(**kargs):
    pass


if __name__ == '__main__':
    main()
