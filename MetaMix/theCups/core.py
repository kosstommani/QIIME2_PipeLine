#!/crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
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
__version__ = '1.0.3'

# -----------------------------------
# 1.0.2 - 2020.04.24
# 2D PCoA Y축 제목 짤림 현상 해결 - x, y 크기 조절 옵션 추가
# -----------------------------------
# Ver. 1.0.3 - 2020.04.28
# Probiotics21 DB 및 보고서 추가
# -----------------------------------


# import os
# from time import time
from click import secho, echo, style
from SpoON.util import parse_config


# 기본값 설정
CONFIG = parse_config()


def miseq_report_pipeline_v1(kargs):
    from theCups.report_miseq_v1 import MiSeqReportV1
    report = MiSeqReportV1(kargs)
    report.check_analysis_data()
    if any(
            [
                kargs['drink'], kargs['read_assembly'], kargs['main'],
                kargs['summary'], kargs['alpha_rarefaction'], kargs['summary'],
                kargs['alpha_rarefaction'], kargs['beta_diversity_3d'], kargs['pcoa_2d'],
                kargs['diversity_index'], kargs['taxonomy_assignment'], kargs['upgma_tree'],
            ]
    ) is False:
        if all([kargs['taste'], kargs['report_number']]) is True:
            pass
        else:
            exit()

    report.check_and_set_for_report_path()
    sample_count = report.check_sample_count()
    if kargs['report_number'] is None:
        report.make_dir(report.report_paths.report_base_path)
        report.make_dir(report.report_paths.report_dir_path)
        report.make_dir(report.report_paths.otu_results_path)
    report.set_report_template()
    if (kargs['drink'] is True) or (kargs['read_assembly'] is True):
        # OTU Picking 방식이 CLOSED(UCLUST)일 경우, STAT.txt 생성을 위해 Filtering 전의 FASTQ를 사용.
        # CLOSED 방식일 경우 FASTA 파일이 필요. LENGTH_FILER 작업시 Filtering 과 FASTA 변환를 동시에 진행.
        # 따라서 STAT 생성에 필요한 Quality Score 정보가 없으므로 STAT.txt 파일을 생성할 수 없음.
        # if report.analysis_data.cd_hit_otu is not None:
        #     otu_picking_method = 'DENOVO'
        # elif report.analysis_data.closed is not None:
        #     otu_picking_method = 'CLOSED'
        report.make_read_assembly_page()
    if (kargs['drink'] is True) or (kargs['summary'] is True):
        # OTU Picking 방식 확인: Denovo(CD-HIT-OTU), CLOSED(UCLUST)
        if report.analysis_data.cd_hit_otu is not None:
            report.make_summary_page_for_cd_hit_otu()
        elif report.analysis_data.closed is not None:
            report.make_summary_page_for_closed()
    if (kargs['drink'] is True) or (kargs['diversity_index'] is True):
        report.make_alpha_diversity_page()
    if (kargs['drink'] is True) or (kargs['alpha_rarefaction'] is True):
        report.make_alpha_rarefaction_page()
    if kargs['drink'] is True:
        report.make_beta_diversity_2d_3d_plot_page(sample_count)

    if kargs['beta_diversity_3d'] is True:
        if sample_count <= 2:
            secho('Error: 시료의 개수가 적어(<=2) Beta Diversity(2D & 3D PCoA)를 생성할 수 없습니다.', fg='red')
            exit()
        elif sample_count == 3:
            secho(f'Warning: 시료의 개수가 적어(==3) 3D PCoA의 생성을 생략합니다.', fg='yellow')
            report.make_beta_diversity_page(p_no_3d=True)
        else:
            report.make_beta_diversity_page(p_no_3d=False)
    if kargs['pcoa_2d'] is True:
        if sample_count <= 2:
            secho('Error: 시료의 개수가 적어(<=2) 2D PCoA를 생성할 수 없습니다.', fg='red')
            exit()
        else:
            if (kargs['pcoa_2d_x_size'] is None) and (kargs['pcoa_2d_y_size'] is None):
                report.make_2d_pcoa_plot_page(p_unifrac='weighted')
            else:
                report.make_2d_pcoa_plot_page(p_unifrac='weighted',
                                              p_x_size=kargs['pcoa_2d_x_size'],
                                              p_y_size=kargs['pcoa_2d_y_size'])

    if (kargs['drink'] is True) or (kargs['upgma_tree'] is True):
        if sample_count <= 2:
            secho(f'Warning: 시료의 개수가 적어(<=2) UPGMA Tree(Beta Diversity)의 생성을 생략합니다.', fg='yellow')
            report.skipped_page.append('upgma_tree')
        else:
            report.make_upgma_tree_page(p_unifrac='weighted')
    if (kargs['drink'] is True) or (kargs['taxonomy_assignment'] is True):
        if sample_count == 1:
            secho(f'Warning: 시료의 개수가 적어(==1) Area Plot의 생성을 생략합니다.', fg='yellow')
            report.skipped_page.append('area')
            report.make_taxonomy_assignment_page('bar')
        else:
            report.make_taxonomy_assignment_page('bar,area')

    # Main HTML 작성
    if kargs['drink'] is True:
        report.make_main_page()
        report.copy_src()

    elif kargs['main'] is True:
        if kargs['report_number'] is not None:
            report.taste_cofi()
            report.convert_checked_results_to_page_status()
            report.make_main_page()
            report.copy_src()
        else:
            report.make_main_page()
            report.copy_src()

    # 생성된 HTML 확인
    # TODO: 함수화
    echo()  # 빈줄 추가
    try:
        max_length = max([len(x) for x in report.page_status.keys()])
    except ValueError:  # max() arg is an empty sequence
        max_length = 0
    l_echo_text = list()
    for key in report.page_status.keys():
        status = report.page_status[key]
        status_color = 'cyan' if status else 'red'
        status_text = style(f'{status}', fg=status_color)
        l_echo_text.append(f'  {key:{max_length}}\t[ {status_text} ]')
    secho('- 보고서 페이지 생성 여부 확인 -'.center(max_length+6), fg='white', bg='cyan', bold=True)
    if l_echo_text:
        secho('\n'.join(l_echo_text))
    else:
        secho('None'.center(33), fg='red', bold=True)
        echo()  # 빈 줄
    if report.skipped_page:
        secho(f'생략 페이지: {", ".join(report.skipped_page)}', fg='yellow')
    # TODO: 보고서 생성이 중단된 디렉터리 및 파일에 대한 삭제 코드 생성

    if kargs['taste'] is True:
        if kargs['report_number'] is None:
            if any(
                [
                    kargs['drink'], kargs['read_assembly'], kargs['taste'],
                    kargs['summary'], kargs['alpha_rarefaction'], kargs['summary'],
                    kargs['alpha_rarefaction'], kargs['beta_diversity_3d'], kargs['pcoa_2d'],
                    kargs['diversity_index'], kargs['taxonomy_assignment'], kargs['upgma_tree'],
                ]
            ) is False:
                secho('Warning: --report_number 옵션을 사용하세요.', fg='')
                echo('--taste 옵션을 무시합니다.')
            else:
                report.taste_cofi()
        else:
            report.taste_cofi()


def miseq_report_pipeline_v2(kargs):
    from theCups.report_miseq_v2 import MiSeqReportV2
    report = MiSeqReportV2(kargs)
    report.check_analysis_data()

    report.check_and_set_for_report_path()
    sample_count = report.check_sample_count()
    if kargs['report_number'] is None:
        report.make_dir(report.report_paths.report_base_path)
        report.make_dir(report.report_paths.report_dir_path)
        report.make_dir(report.report_paths.otu_results_path)

    report.set_report_template()
    if (kargs['drink'] is True) or (kargs['read_assembly'] is True):
        report.make_read_assembly_page()
    if (kargs['drink'] is True) or (kargs['length_filter'] is True):
        report.make_length_filter_page()
    if (kargs['drink'] is True) or (kargs['clustering'] is True):
        report.make_cd_hit_otu_page()

    # 생성된 HTML 확인
    # TODO: 함수화
    echo()  # 빈줄 추가
    try:
        max_length = max([len(x) for x in report.page_status.keys()])
    except ValueError:  # max() arg is an empty sequence
        max_length = 0
    l_echo_text = list()
    for key in report.page_status.keys():
        status = report.page_status[key]
        status_color = 'cyan' if status else 'red'
        status_text = style(f'{status}', fg=status_color)
        l_echo_text.append(f'  {key:{max_length}}\t[ {status_text} ]')
    secho('- 보고서 페이지 생성 여부 확인 -'.center(max_length+6), fg='white', bg='cyan', bold=True)
    if l_echo_text:
        secho('\n'.join(l_echo_text))
    else:
        secho('None'.center(33), fg='red', bold=True)
        echo()  # 빈 줄
    if report.skipped_page:
        secho(f'생략 페이지: {", ".join(report.skipped_page)}', fg='yellow')
    # TODO: 보고서 생성이 중단된 디렉터리 및 파일에 대한 삭제 코드 생성


def probiotics_report(kargs, p_type):
    if p_type == 'Probiotics19':
        from theCups.probiotics import Data, check_sample_name, Report19
        Report = Report19
    elif p_type == 'Probiotics21':
        from theCups.probiotics import Data, check_sample_name, Report21
        Report = Report21
    else:
        raise ValueError(f'p_type 파라미터의 값이 올바르지 않습니다. p_type: {p_type}')
    i_data = Data(kargs, p_type)
    i_data.check_info_file()
    i_data.read_info_file()
    i_data.check_table_file()

    # 투입종 입력 방법 선택
    i_data.read_species_list()
    if kargs['mixed_species'] is None:
        i_data.check_mixed_species()
    else:
        i_data.read_mixed_species_list()
    i_data.read_table_file()

    # 유산균 제품 이름 입력 - 직접 입력
    if kargs['mixed_species'] is None:
        i_data.input_product_name()

    # 시료별 보고서 생성
    v1v2_status = True if (kargs['v1v2_dir'] is not None) and (kargs['v1v2_num'] is not None) else False
    v3v4_status = True if (kargs['v3v4_dir'] is not None) and (kargs['v3v4_num'] is not None) else False
    if v1v2_status and v3v4_status:  # V1V2영역 + V3V4영역인 경우
        for v1v2_sample in i_data.l_v1v2_sample_table_list:
            v1v2_sample_name = check_sample_name(v1v2_sample.sample_name, 'v1v2')
            if kargs['mixed_species'] is None:
                pass
            elif v1v2_sample_name not in i_data.d_mixed_species_from_file.keys():  # 보고서 생성 제외 시료 확인
                continue
            last_sample_indicator = len(i_data.l_v3v4_sample_table_list)
            for count, v3v4_sample in enumerate(i_data.l_v3v4_sample_table_list, 1):  # V1V2의 시료명과 V3V4의 시료명 확인
                v3v4_sample_name = check_sample_name(v3v4_sample.sample_name, 'v3v4')
                if v1v2_sample_name == v3v4_sample_name:  # 두 영역의(V1V2, V3V4)의 시료명이 일치하는 경우
                    if kargs['mixed_species'] is None:
                        mixed_species = i_data.mixed_species
                    else:
                        if v3v4_sample_name not in i_data.d_mixed_species_from_file:  # 보고서 생성 제외 시료 확인
                            continue
                        mixed_species = i_data.d_mixed_species_from_file[v1v2_sample_name]

                    product_name = i_data.d_sample_product_name[v1v2_sample_name]
                    report = Report(
                        {
                            'o_v1v2_sample_table': v1v2_sample,
                            'o_v3v4_sample_table': v3v4_sample,
                            'product_name': product_name,
                            'analysis_base_path': i_data.analysis_base_path,
                            'order_number': i_data.order_number,
                            'info': i_data.info,
                            'd_species_list_index': i_data.d_species_list_index,
                            'mixed_species': mixed_species,
                        })
                    report.make_report()
                    break
                else:  # 두 영역의 시료명이 일치하지 않는 경우
                    if count == last_sample_indicator:  # V1V2영역의 시료가 V3V4영역에 없는 경우 처리
                        secho('Error : V1V2영역의 시료가 V3V4영역의 데이터에 없습니다!', fg='red', blink=True)
                        secho('각 영역(V1V2, V3V4)의 시료 구성을 확인하세요. (개발자 문의)', fg='magenta')
                        exit()
                    else:
                        continue
            del last_sample_indicator, v3v4_sample
        del v1v2_sample
    else:  # 1개 영역만 진행한 경우
        selected_region_sample_table_list = None
        selected_target_region = None
        if (v1v2_status is True) and (v3v4_status is False):  # V1V2영역만 있는 경우
            selected_region_sample_table_list = i_data.l_v1v2_sample_table_list
            selected_target_region = 'v1v2'
        elif (v1v2_status is False) and (v3v4_status is True):  # V3V4영역만 있는 경우
            selected_region_sample_table_list = i_data.l_v3v4_sample_table_list
            selected_target_region = 'v3v4'
        for selected_sample in selected_region_sample_table_list:
            selected_sample_name = check_sample_name(selected_sample.sample_name, selected_target_region)
            if kargs['mixed_species'] is None:  # 대화창으로 입력
                mixed_species = i_data.mixed_species
            else:  # 파일로 입력
                if selected_sample_name not in i_data.d_mixed_species_from_file:  # 보고서 생성 제외 시료 확인
                    continue
                mixed_species = i_data.d_mixed_species_from_file[selected_sample_name]

            product_name = i_data.d_sample_product_name[selected_sample_name]
            if selected_target_region == 'v1v2':
                report = Report(
                    {
                        'o_v1v2_sample_table': selected_sample,
                        'o_v3v4_sample_table': None,
                        'product_name': product_name,
                        'analysis_base_path': i_data.analysis_base_path,
                        'order_number': i_data.order_number,
                        'info': i_data.info,
                        'd_species_list_index': i_data.d_species_list_index,
                        'mixed_species': mixed_species,
                    })
            elif selected_target_region == 'v3v4':
                report = Report(
                    {
                        'o_v1v2_sample_table': None,
                        'o_v3v4_sample_table': selected_sample,
                        'product_name': product_name,
                        'analysis_base_path': i_data.analysis_base_path,
                        'order_number': i_data.order_number,
                        'info': i_data.info,
                        'd_species_list_index': i_data.d_species_list_index,
                        'mixed_species': mixed_species,
                    })
            report.make_report()
