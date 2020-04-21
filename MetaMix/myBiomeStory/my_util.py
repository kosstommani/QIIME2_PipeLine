# ----------------------------------------------------------------------------------------------------------------------
#                        888b     d888          888             888b     d888 d8b
#                        8888b   d8888          888             8888b   d8888 Y8P
#                        88888b.d88888          888             88888b.d88888
# 88888b.d88b.  888  888 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b 888  888 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 Y88b 888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888  "Y88888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#                    888
#               Y8b d88P
#                "Y88P"
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '0.2.3'

import os
import time
import csv
from pprint import pprint
import string
from click import echo, secho
from SpoON.run_time import compute_run_time
from SpoON.util import check_file_type, glob_dir, get_order_json
from CoFI.core import make_sample_list


def print_run_time(p_start, p_prema, p_cofi, p_closed, p_diversity, p_summarize, p_make_story, p_insert):
    date_fmt = '%Y-%m-%d %a %H:%M:%S'
    start_time = time.strftime(date_fmt, time.localtime(p_start))
    prema_run_time = compute_run_time(p_start, p_prema)
    cofi_run_time = compute_run_time(p_prema, p_cofi)
    closed_run_time = compute_run_time(p_cofi, p_closed)
    diversity_run_time = compute_run_time(p_closed, p_diversity)
    summarize_run_time = compute_run_time(p_diversity, p_summarize)
    make_story_run_time = compute_run_time(p_summarize, p_make_story)
    insert_run_time = compute_run_time(p_make_story, p_insert)
    total_run_time = compute_run_time(p_start, p_insert)
    end_time = time.strftime(date_fmt, time.localtime(p_insert))
    msg = string.Template('''\
               ,        ,              주문시간: $start 
              /(        )`              - myMetaMix 제조 시간 - 
              \\ \\___   / |                  PreMA(Copy & fastp)  
              /- _  `-/  '                       : $prema
             (/\\/ \\ \\   /\\                  CoFI(FLASH & Length Filtering)     
             / /   | `    \\                      : $cofi  
             O O   ) /    |                 Closed - Reference OTU Picking
             `-^--'`<     '                      : $closed
            (_.)  _  )   /                  Alpha-Diversity
             `.___/`    /                        : $diversity
               `-----' /                    Taxa Summary
  <----.     __ / __   \\                         : $summary
  <----|====O)))==) \\) /====                Make StoryIndex:
  <----'    `--' `.__,' \\                        : $story
               |        |                   Insert to DB:
                \\       /                        : $insert
           ______( (_  / \\______            Total
         ,'  ,-----'   |        \\                : $total
         `--{__________)        \\/     완성시간: $end
      악마의 유혹 - MetaMix 카페
''')
    echo(msg.substitute({
        'start': start_time,
        'prema': prema_run_time,
        'cofi': cofi_run_time,
        'closed': closed_run_time,
        'diversity': diversity_run_time,
        'summary': summarize_run_time,
        'story': make_story_run_time,
        'insert': insert_run_time,
        'total': total_run_time,
        'end': end_time
    }))


def check_sample_list(p_base_path, p_order_number, p_dir_suffix, p_no_check):
    """
    분석을 진행하고자 하는 수주의 시료들을 목록화한다.
    시료들의 목록은 해당 수주의 작업 디렉터리내에 존재하는 RawData 디렉터리에서 시료들의 디렉터리 이름을 기준으로 생성한다.
      base_path --- order_number --- RawData*/*_*_? --- 시료명

    :type p_base_path: str
    :param p_base_path: 분석 경로
    :type p_order_number: str
    :param p_order_number: 수주번호
    :type p_dir_suffix: str
    :param p_dir_suffix: 목표 디렉터리 접미사. 와일드 카드 사용 가능.
    :type p_no_check: bool
    :param p_no_check: 분석 작업 디렉터리(Analysis_?)가 여러 개인지 확인여부.
                       True: 작업 디렉터리가 여러 개 존재할 경우 목록 확인 후 선택.
                       False: 작업 디렉터리가 여러 개 존재할 경우 작업 종료.
    :return: l_sample, analysis_path
             l_sample: 시료의 목록. 시료명을 원소로 가지는 리스트.
             analysis_path: 분석 작업 디렉터리(~/Anaysis_?)의 경로.
    """
    order_path = os.path.join(*(p_base_path, p_order_number))
    if check_file_type(order_path, 'exists'):
        l_nt_run_list = make_sample_list(p_base_path, p_order_number, p_dir_suffix, None)
        l_samples = list()
        for run in l_nt_run_list:
            for sample in run.samples:
                l_samples.append(sample.name)
        del l_nt_run_list

        if p_no_check is True:
            analysis_path = glob_dir(order_path, 'Analysis_*', 'only')
        elif p_no_check is False:
            analysis_path = glob_dir(order_path, 'Analysis_*', 'many')
            if len(analysis_path) > 1:
                secho('작업 디렉터리가 여러개 존재합니다.', fg='red', err=True)
                echo(analysis_path, err=True)
                exit(1)
            else:
                analysis_path = analysis_path[0]
        else:
            raise ValueError('no_check : {}'.format(p_no_check))
    else:
        secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red', err=True)
        exit(1)
    return l_samples, analysis_path


def filter_sample(p_l_list, p_l_exclude):
    """

    :type p_l_list: list
    :param p_l_list:
    :type p_l_exclude: list
    :param p_l_exclude:
    :return:
    """
    l_final = list()
    for index in range(len(p_l_list)):
        sample = p_l_list.pop()
        if sample in p_l_exclude:
            continue
        else:
            l_final.append(sample)
    return l_final


def set_sample_list(p_l_list, p_l_include, p_l_exclude):
    """

    :type p_l_list: list
    :param p_l_list:
    :type p_l_include: list
    :param p_l_include:
    :type p_l_exclude: list
    :param p_l_exclude:
    :return:
    """
    l_final_samples = None
    if p_l_include is not None:
        l_final_samples = p_l_include
    elif p_l_exclude is not None:
        l_final_samples = filter_sample(p_l_list, p_l_exclude)
    secho('>>> 분석 대상 시료 설정 완료', fg='cyan')
    secho(f'시료개수: {len(l_final_samples)}', fg='yellow', blink=True, bold=True)
    echo(l_final_samples)
    return l_final_samples


def get_sample_info(p_order_number, p_output):
    order_num_type, json_data = get_order_json(p_order_number, p_service='myBiomeStory')
    if json_data['resultCode'] == 200:
        if p_output.endswith('.csv'):
            output_file = p_output
        else:
            output_file = f'{p_output}.csv'
        with open(output_file, 'w') as o_output:
            fieldnames = ['KitId', 'client_name', 'Sex', 'Birth', 'sample_source',
                          'sample_kit', 'sample_type', 'date_received', 'date_collected', 'code_reseller',
                          'sample_note', 'client_id_number']
            csv_writer = csv.DictWriter(o_output, fieldnames=fieldnames)
            csv_writer.writeheader()
            csv_writer.writerows(json_data['resultData'])
        secho('>>> SampleInfo.csv 생성 완료', fg='cyan')
        return True

    elif json_data['resultCode'] == -1:
        secho('Error: 고객정보 JSON 에 문제가 있습니다.', fg='red')
        echo(f'resultCode: {json_data["resultCode"]}')
        echo(f'resultMsg: {json_data["resultMsg"]}')
        return False
    else:
        secho('Error: 고객정보 JSON 에 알 수 없는 문제가 발생했습니다.', fg='red', blink=True)
        pprint(json_data)
        return False


def read_csv_to_dict(p_file, p_verbose=True):
    """
    CSV 파일을 읽고 데이터를 딕션너리로 반환한다.

    :param p_file:
    :type p_verbose: bool
    :param p_verbose: 안내메시지 출력 여부
    :rtype list
    :return csv_data - 딕션너리를 원소로 가지는 리스트
    """
    with open(p_file, 'r') as o_csv:
        o_csv_data = csv.DictReader(o_csv)
        csv_data = [x for x in o_csv_data]
    if p_verbose:
        secho('>>> SampleInfo.csv 읽기 완료', fg='cyan')
        echo(p_file)
    return csv_data


def read_csv(p_file, p_verbose=True):
    """
    CSV 파일을 읽고 데이터를 리스트로 반환한다.

    :param p_file:
    :param p_verbose: 안내메시지 출력 여부
    :return: csv_data
    """
    with open(p_file, 'r') as o_csv:
        o_csv_data = csv.reader(o_csv)
        csv_data = [x for x in o_csv_data]
    if p_verbose:
        secho('>>> csv 파일 읽기 완료', fg='cyan')
        echo(p_file)
    return csv_data


def check_db_version(p_list):
    l_version = [x[-1] for x in p_list]
    s_version = set(l_version)
    if len(s_version) == 1:
        return s_version.pop()
    else:
        secho(f'DB version: {s_version}')
        return False
