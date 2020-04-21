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
__version__ = '1.0.0'

from click import secho, echo, style, progressbar
import os
import subprocess
from humanize import intcomma
from collections import defaultdict
import shutil
from SpoON.util import glob_dir, check_file_type, read_data, run_cmd, check_run_cmd, parse_config, \
                       get_order_json, launcher_cmd
from PreMA.data_structure import RawDataPath, SampleList
from pprint import pprint

# from multiprocessing import Pool, current_process
# import datetime
# import pdb


CONFIG = parse_config()
FASTQC = CONFIG['PreMA_RawData_Handler']['FastQC']
MULTIQC = CONFIG['PreMA_RawData_Handler']['MultiQC']
POOL_WORKER = CONFIG['PreMA_PoolWorker']


def find_rawdata(p_path, p_target_dir_suffix):
    """
    RawData 보관디렉터리내 해당 수주번호의 Run 단위 데이터를 검색한다.
    분석에 사용할 Run 디렉터리를 선택하고, 선택되어진 Run 디렉터리내의 시료 디렉터리를 목록화한다.

    :type p_path: str
    :param p_path: 목표 경로
           ex) /crystal/Analysis/Project/1808AMI-0008

    :type p_target_dir_suffix: str
    :param p_target_dir_suffix: Run 디렉터리 검색 조건(와일드 카드 사용 가능)
           ex) 180828_BVGWT_1

    :rtype: object
    :return: nt_rawdata_path
             named tuple인 RawDataPath
             - RawDataPath.rundata_path : run 디렉터리. ex) /crystal/Analysis/Project/1810KMI-0007/181101_C2BNW_1
                          .sample_list : run 디렉터리내의 시료명 디렉터리의 목록.
                                  ex) ['/crystal/Analysis/Project/1810KMI-0007/181101_C2BNW_1/Anaerobic',
                                      '/crystal/Analysis/Project/1810KMI-0007/181101_C2BNW_1/Aerobic1']
                          .stat_file : stat.txt 파일. ex) /crystal/Analysis/Project/1810KMI-0007/181101_C2BNW_1/stat.txt
                          .nread_file : NreadPercent.csv 파일. ex) /crystal/Analysis/Project/HN00102528
                                                            /190221_C92TC_1/NreadPercent.csv
    """
    target_dir_path = glob_dir(p_path, p_target_dir_suffix)
    dir_list = glob_dir(target_dir_path, '*', 'many')
    sample_dir_list, not_dir_list = check_file_type(dir_list, 'isdir')
    stat_file = '{0}/stat.txt'.format(target_dir_path)
    nread_file = '{0}/NreadPercent.csv'.format(target_dir_path)
    if stat_file in not_dir_list:
        pass
    else:
        secho('Warning: {stat} 파일이 없습니다.'.format(stat=stat_file), fg='yellow')
        secho('--> MAGEBI에서 init 실행!', fg='cyan')

    if nread_file in not_dir_list:
        pass
    else:
        secho('Warning: {nread} 파일이 없습니다.'.format(nread=nread_file), fg='yellow')
        secho('--> N Read %을 확인하세요.', fg='cyan')

    # 제외할 디렉터리 또는 데이터가 없는 샘플 제외
    no_sample = ['report_dir', 'tmp_sh', 'md5', 'sqs', 'fastqc']
    exclusion_index = list()
    for index, element in enumerate(sample_dir_list):
        if element.strip().split('/')[-1] in no_sample:
            exclusion_index.append(index)
    if exclusion_index:
        secho('>>> 복사 제외 목록', fg='cyan')
        exclusion_index.sort()
        exclusion_index.reverse()
        for index in exclusion_index:
            echo(sample_dir_list.pop(index))
    nt_rawdata_path = RawDataPath(rundata_path=target_dir_path,
                                  sample_list=sample_dir_list,
                                  stat_file=stat_file,
                                  nread_file=nread_file)
    return nt_rawdata_path


# TODO: 삭제
# Header 변경을 FlowCell ID을 이용한 변경이 아니라 Header 전체를 변경할 것이므로
def check_flowcell_id(p_list, p_platform='MiSeq'):
    """

    :param p_list:
    :param p_platform:
    :return:
    """
    dic_flowcell_id = dict()
    for sample_dir in p_list:
        sample_name = os.path.split(sample_dir)[-1]
        completed_process = subprocess.run('zcat -q {file} | head -n 1'.format(
                file=os.path.join(*(sample_dir, sample_name + '_1.fastq.gz'))),
            shell=True, check=True, stdout=subprocess.PIPE, encoding='utf-8')
        fastq_header = completed_process.stdout
        if p_platform == 'HiSeq':
            # '@HISEQ:425:AV9JG:1:1101:16867:1483 1:N:0:ATGCGCAGCTATTAAG'
            flowcell_id = fastq_header.sttrip().split(':')[2]
        elif p_platform == 'MiSeq':
            # '@MISEQ:425:000000000-AV9JG:1:1101:16867:1483 1:N:0:ATGCGCAGCTATTAAG'
            flowcell_id = fastq_header.strip().split('000000000-')[1].split(':')[0]
        else:
            raise ValueError(p_platform)
        dic_flowcell_id[sample_name] = flowcell_id
    echo('>>> 시료별 FlowCell ID 읽기 완료')
    return dic_flowcell_id


def check_stat(p_file, p_nbase, p_q30, p_sample_read):
    """
    RawData 디렉터리에 있는 stat.txt 파일을 읽어 출력한다.

    :type p_file: str
    :param p_file: stat.txt 파일의 경로

    :type p_nbase: list or tuple
    :param p_nbase: [float, float] N base의 알람 기준치. 첫 번째 요소(>): 주의, 두 번째 요소(>): 경고.

    :type p_q30: list or tuple
    :param p_q30: [int, int] Q30의 알람 기준치. 첫 번째 요소(<): 경고, 두 번째 요소(<): 주의.

    :type p_sample_read: int
    :param p_sample_read: 시료별 생산해야 될 Read의 개수(기준치).

    :return: 없음. STDOUT 출력.
    """
    stat_data = read_data(p_file)
    stat_data.sort()
    max_name = max([len(x[0]) for x in stat_data])
    max_base = max([len(intcomma(x[1])) for x in stat_data])
    max_read = max([len(intcomma(x[2])) for x in stat_data])
    header = '\t'.join(['{sample:>{digit1}}', '{base:>{digit2}}', '{read:>{digit3}}',
                        'Length', 'N_base(%)', 'GC(%)', '{q20:>5}', '{q30:>5}']).format(
                         digit1=max_name,
                         digit2=max_base,
                         digit3=max_read,
                         sample='Sample',
                         base='Base',
                         read='Read',
                         q20='Q20',
                         q30='Q30',
                        )
    size = int(len(header) / 2)
    secho('-' * size + 'STAT 확인' + '-' * size, fg='cyan')
    echo(header)
    for stat in stat_data:
        read_count = int(stat[2])
        if read_count < p_sample_read:
            read_text = '{read:>{digit}}'.format(read=intcomma(stat[2]), digit=max_read)
            checked_read = style(read_text, fg='red')
        else:
            checked_read = '{read:>{digit}}'.format(read=intcomma(stat[2]), digit=max_read)

        nbase = float(stat[3])
        if nbase >= p_nbase[1]:
            stat[3] = style('{:>9}'.format(nbase), fg='red')
        elif nbase >= p_nbase[0]:
            stat[3] = style('{:>9}'.format(nbase), fg='yellow')
        else:
            stat[3] = '{:>9}'.format(nbase)

        q30 = float(stat[6])
        if q30 < p_q30[0]:
            stat[6] = style('{:.2f}'.format(q30), fg='red')
        elif q30 < p_q30[1]:
            stat[6] = style('{:.2f}'.format(q30), fg='yellow')
        else:
            stat[6] = '{:.2f}'.format(q30)
        echo('\t'.join(['{sample:>{digit1}}', '{base:>{digit2}}', '{read}', '{length:>6.1f}',
                        '{n_base}', '{gc:.2f}', '{q20:.2f}', '{q30}']).format(
            digit1=max_name,
            digit2=max_base,
            sample=stat[0],
            base=intcomma(stat[1]),
            read=checked_read,
            length=float(stat[1]) / float(stat[2]),
            n_base=stat[3],
            gc=float(stat[4]),
            q20=float(stat[5]),
            q30=stat[6],
        ))
    secho('-' * size + '-' * len('STAT 확인'.encode('utf-8')) + '-' * size, fg='cyan')


def check_nread(p_file, p_nread):
    try:
        nread_data = read_data(p_file, ',')
    except FileNotFoundError:
        secho('Warning: {file}을 열 수 없습니다.'.format(file=p_file), fg='yellow')
        secho('--> N Read % 을 확인하세요.', fg='magenta')
        return
    else:
        nread_data.sort()
        max_name = max([len(x[0]) for x in nread_data])

        header = '\t'.join(['{sample:>{digit1}}', '{n_read1}', '{n_read2}']).format(
            sample='Sample', digit1=max_name,
            n_read1='N_Read_1(%)', n_read2='N_Read_2(%)')
        size = int(len(header) / 2)
        secho('-' * size + 'N Read % 확인' + '-' * size, fg='cyan')
        echo(header)

        d_nread_data = defaultdict(dict)
        for stat in nread_data:
            d_nread_data[stat[0]][stat[1]] = stat[2]

        for key in sorted(d_nread_data.keys()):
            nread_1 = float(d_nread_data[key]['_1'])
            nread_2 = float(d_nread_data[key]['_2'])
            if nread_1 >= p_nread:
                nread_1 = style('{:>10.2f}'.format(nread_1), fg='red')
            else:
                nread_1 = '{:>10.2f}'.format(nread_1)

            if nread_2 >= p_nread:
                nread_2 = style('{:>10.2f}'.format(nread_2), fg='red')
            else:
                nread_2 = '{:>10.2f}'.format(nread_2)
            echo('\t'.join(['{sample:>{digit1}}', '{n_read1:>5}', '{n_read2:>5}']).format(
                digit1=max_name,
                sample=key,
                n_read1=nread_1,
                n_read2=nread_2,
            ))
        secho('-' * size + '-' * len('N Read % 확인'.encode('utf-8')) + '-' * size, fg='cyan')


# TODO : 추후 삭제 예정. LIMS2 에서 사용.
def get_index_kit_info(p_file):
    """
    directory_path_kit.txt 파일을 읽어 라이브러리 제작에 사용한 index kit 정보를 확인한다.

    :type p_file: str
    :param p_file: index kit 정보가 담긴 파일(directory_path_kit.txt)

    :rtype: str
    :return: 라이브러리 제작에 사용한 index kit 의 이름.
    """
    completed_process = subprocess.run('cut -f 2 -d , {file} | sort | uniq'.format(file=p_file),
                                       shell=True, check=True, stdout=subprocess.PIPE, encoding='utf-8')
    return completed_process.stdout.strip()


# TODO : 추후 삭제 예정. LIMS2 에서 사용.
def check_index_kit_info(p_file, p_run_path, p_library_kit):
    """
    라이브러리 제작시 사용한 index kit의 정보를 확인한다.
    해당 run 디렉터리에서 index kit 정보를 확인하지 못할 경우(kit file이 없을 경우), 상위 디렉터리에서 검색한다.

    :type p_file: str
    :param p_file: directory_path_kit.txt 파일명.

    :type p_run_path: str
    :param p_run_path: RawData의 run 디렉터리 경로.

    :type p_library_kit: str
    :param p_library_kit: auto, Nextera, TruSeq.
                auto일 경우 RawData 정보에서 확인.
                Nextera, TruSeq일 경우 RawData 정보에서 확인은 하지만 강제 실행.

    :rtype: str
    :return: index_kit_from_rawdata.
             index Kit 정보.
    """
    try:
        kit_info = get_index_kit_info(os.path.join(*(p_run_path, p_file)))
        if kit_info.lower() == 'nextera':
            secho('>>> Detected Library Index Kit : {}'.format(kit_info), fg='cyan')
            index_kit_from_rawdata = 'nextera'
        elif kit_info == '':
            echo('>>> Library Index Kit 정보 찾기 시작')
            order_path = os.path.split(p_run_path)[0]
            globed_dir, globed_file = check_file_type(glob_dir(order_path, '*', 'many', False), 'isdir')
            l_target_kit_info = list()
            with progressbar(globed_dir, label='directory_path_kit.txt 검색 중') as bar__globed_dir:
                for dir_path in bar__globed_dir:
                    target_file = glob_dir(dir_path, 'p_file', p_verbose=False)
                    if target_file is not None:
                        target_kit_info = get_index_kit_info(target_file)
                        l_target_kit_info.append((target_file, target_kit_info))

            for element in globed_file:
                if p_file in element:
                    target_kit_info = get_index_kit_info(element)
                    l_target_kit_info.append((element, target_kit_info))

            for target_file, target_kit_info in l_target_kit_info:
                echo('{file} : {kit}'.format(file=target_file, kit=target_kit_info))

            unique_target_kit_info = set([x[1] for x in l_target_kit_info])
            if len(unique_target_kit_info) == 1:
                secho('>>> Detected Library Index Kit : {}'.format(list(unique_target_kit_info)), fg='cyan')
                kit = unique_target_kit_info.pop()
                if kit.lower() == 'nextera':
                    index_kit_from_rawdata = 'nextera'
                else:
                    index_kit_from_rawdata = kit
            else:
                echo(unique_target_kit_info)
                secho('>>> 파일이 없거나 여러 개의 Library Kit이 확인됩니다.!', fg='red', blink=True)
                if p_library_kit == 'auto':
                    echo('특정 Library Index Kit으로 진행할려면 -l(--library_kit)을 사용!', err=True)
                    exit()
                else:
                    index_kit_from_rawdata = p_library_kit
                    secho('>>> Library Index Kit: {} 적용'.format(index_kit_from_rawdata), fg='yellow')
        else:
            secho('>>> Library Index Kit : {}'.format(kit_info), fg='red', blink=True)
            index_kit_from_rawdata = kit_info
    except subprocess.CalledProcessError as err:
        secho('>>>  Library Index Kit 확인 필요!', fg='red', blink=True)
        echo('명령어 오류 : {}'.format(err))
    return index_kit_from_rawdata


def check_index_kit_info_using_json(p_json, p_index_kit):
    """
    LIMS의 JSON 파일에서 Library Index Kit 정보를 확인한 후 요청된 Index Kit 정보와 비교한다.
    확인된 Index Kit 정보와 요청된 Index Kit 정보와 다를 경우에는 요청된 Index Kit 정보가 적용된다.

    :type p_json: dict
    :param p_json: LIMS에 등록된 수주의 JSON Data
    :param p_index_kit: 요청된 Index Kit 정보(--index_kit)
    :return: index_kit_from_rawdata
    """
    d_kit = dict()
    if p_json is None:
        secho('Warning: JSON 파일에 데이터가 없습니다.(None)', fg='yellow')
        s_kit = list()
    else:
        try:
            for sample in p_json['demux']:
                for run in sample['run']:
                    d_kit[sample['sample_id']] = run.get('process')[0].get('Library Kit')
        except KeyError:
            secho('Warning: JSON 파일에서 특정 항목(demux or run)을 찾을 수 없습니다.(KeyError)', fg='yellow')
            s_kit = list()
        else:
            s_kit = set(d_kit.values())

    if len(s_kit) == 1:
        kit_text = list(s_kit)[0]
        if 'nextera' in kit_text.lower():
            index_kit_from_rawdata = 'Nextera'
        elif 'truseq' in kit_text.lower():
            index_kit_from_rawdata = 'TruSeq'
        else:
            secho('Warning: LIMS - Library Index Kit : {}'.format(kit_text), fg='red')
            index_kit_from_rawdata = kit_text

        if p_index_kit != 'auto':
            if index_kit_from_rawdata.lower() != p_index_kit.lower():
                secho('>>> 확인된 Library Index Kit : {}'.format(index_kit_from_rawdata), fg='yellow', blink=True)
                echo(kit_text)
                secho('>>> 요청된 Library Index Kit : {}'.format(p_index_kit), fg='cyan')
                return p_index_kit
            else:
                return index_kit_from_rawdata
        else:
            secho('>>> 확인된 Library Index Kit : {}'.format(index_kit_from_rawdata), fg='cyan')
            echo(kit_text)
            return index_kit_from_rawdata
    elif len(s_kit) > 1:
        secho('>>> 확인된 Library Index Kit : {}'.format(len(s_kit)), fg='yellow', blink=True)
        echo(list(s_kit))
        pprint(d_kit)
        if p_index_kit == 'auto':
            secho('---> 특정 Library Index Kit 으로 진행할려면 -ik(--index_kit)을 사용하세요.', fg='magenta')
            exit()
        else:
            secho('>>> 요청된 Library Index Kit : {}'.format(p_index_kit), fg='cyan')
            index_kit_from_rawdata = p_index_kit
            return index_kit_from_rawdata
    else:
        secho('Warning: LIMS의 JSON 파일에서 Index Kit 정보를 확인할 수 없습니다.', fg='yellow', blink=True)
        secho('---> LIMS에서 Library Index Kit 정보를 확인하세요!', fg='magenta')
        if p_index_kit == 'auto':
            secho('---> 특정 Library Index Kit 으로 진행할려면 -ik(--index_kit)을 사용하세요.', fg='magenta')
            exit()
        else:
            secho('>>> 요청된 Library Index Kit : {}'.format(p_index_kit), fg='cyan')
            index_kit_from_rawdata = p_index_kit
            return index_kit_from_rawdata


def check_output_file_requirements_from_json(p_json):
    output_file = dict()
    try:
        for sample in p_json.get('analysis'):
            output_file[sample.get('description')] = sample.get('output_file_requirements')
    except KeyError:
        pass
    except AttributeError:  # p_json 이 None인 경우
        pass

    if output_file is False:
        secho('Warning: JSON 데이터에서 \'output_file_requirements\' 정보를 확인할 수 없습니다.', fg='yellow')
    else:
        secho('>>> 데이터 전달 방법', fg='cyan')
        for key in output_file.keys():
            echo('{key}: {value}'.format(key=key, value=output_file[key]))


def add_zero(p_max, p_str):
    """
    최대시료명을 기준으로 시료명에 0을 추가한다.
    ex) p_max: 300, p_str: 5 --> 005

    :type p_max: str
    :param p_max: 시료명이 문자열 숫자인 것 중에서 길이가 가장 큰 시료명.

    :type p_str: str
    :param p_str: 변경할 시료명(문자열 숫자).

    :rtype: str
    :return: new_str
             str - 0을 추가한 시료명.
    """
    if p_str.isdigit():
        if len(p_str) < len(str(p_max)):
            new_str = '0'*(len(str(p_max))-len(p_str)) + p_str
        else:
            new_str = p_str
    else:
        new_str = p_str
    return new_str


def check_sample_name(p_list, p_custom_name, p_mode='all'):
    """
    분석에 사용가능한 시료명으로 변환한다.
    시료명을 새로운 이름으로 변경한다.

    :type p_list: list
    :param p_list: 경로를 포함하는 시료 디렉터리명을 원소로 가지는 리스트

    :type p_custom_name: dict
    :param p_custom_name:  변경할 이름을 담고 있는 딕션너리. Key: 현재 시료명, Value: 변경할 시료명.
    :param p_mode: all, all_no, no_zero, no_mid
                    all - zero, middle zero 모두 추가.
                    all_no - no_zero, no_mid 모두 적용.
                    no_zero - 문자형 숫자에 0을 추가하지 않음.
                    no_mid - '.'를 구분자로 하는 문자열에서 구분자를 기준으로 문자형 숫자가 있는 경우 0을 추가하지 않음.
    :rtype: list
    :return: l_nt_new_sample_list
             SampleList 를 원소로 가지는 리스트
             - SampleList - path cur_name new_name 의 이름을 가지는 named tuple
    """
    l_nt_temp_sample_list = list()
    l_nt_new_sample_list = list()
    l_digit = list()
    l_not_digit = list()
    p_list.sort()
    for path_dir in p_list:
        path, name = os.path.split(path_dir)
        if p_custom_name is not None:
            custom_name = p_custom_name.get(name)
            if custom_name:
                new_name = custom_name.replace('-', '.').replace('_', '.')
            else:
                new_name = name.replace('-', '.').replace('_', '.')
        else:
            new_name = name.replace('-', '.').replace('_', '.')
        l_nt_temp_sample_list.append(SampleList(path, name, new_name))
        if new_name.isdigit():
            l_digit.append(new_name)
        else:
            l_not_digit.append(new_name)
    del path_dir, path, name, new_name

    if len(l_digit) != 0:
        max_digit = max([int(x) for x in l_digit])

    if len(l_not_digit) != 0:
        l_status = list()
        for name in l_not_digit:
            pieces = name.split('.')
            status = [int(x) if x.isdigit() else 0 for x in pieces]
            l_status.append(status)
        del name, pieces, status
        d_max = defaultdict(int)
        for ele in l_status:
            for index, i in enumerate(ele):
                value = d_max[index]
                if value == 0:
                    d_max[index] = i
                elif value < i:
                    d_max[index] = i
                else:
                    continue
        del ele, l_status

    if p_mode == 'all_no':
        l_nt_new_sample_list = l_nt_temp_sample_list
    else:
        for sample in l_nt_temp_sample_list:
            if sample.new_name in l_digit:
                if p_mode != 'no_zero':
                    l_nt_new_sample_list.append(
                        SampleList(
                            sample.path,
                            sample.cur_name,
                            add_zero(max_digit, sample.new_name)
                        )
                    )
                elif p_mode == 'no_zero':
                    l_nt_new_sample_list.append(sample)
            elif sample.new_name in l_not_digit:
                l_elements = sample.new_name.split('.')
                temp_name = list()
                for index, i in enumerate(l_elements):
                    if i.isdigit():
                        temp_name.append(add_zero(d_max[index], i))
                    else:
                        temp_name.append(i)
                if p_mode != 'no_mid':
                    l_nt_new_sample_list.append(SampleList(sample.path, sample.cur_name, '.'.join(temp_name)))
                elif p_mode == 'no_mid':
                    l_nt_new_sample_list.append(sample)
            else:
                raise RuntimeError(secho('논리적 오류', fg='red', blink=True))

    secho('>> 시료명 변경({num})'.format(num=len(l_nt_new_sample_list)), fg='cyan')
    cur_max = max([len(x.cur_name) for x in l_nt_new_sample_list])
    new_max = max([len(x.new_name) for x in l_nt_new_sample_list])
    for ele in l_nt_new_sample_list:
        echo('{cur:>{di1}} --> {new:>{di2}}'.format(cur=ele.cur_name, new=ele.new_name, di1=cur_max, di2=new_max))
    return l_nt_new_sample_list


# TODO: check_flowcell_id 함수 실행 코드 삭제
def check_rawdata(kargs):
    """
    시료별 FlowCell ID를 목록화하고(삭제예정), FASTQ STAT 에 대해서 기준 충족여부를 확인하여 출력한다.
    라이브러리 제작에 사용한 index kit 정보를 확인한다.
    결과 전달 방법을 확인한다.
    시료명을 분석에 적합하게 변경한다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                sample_list(list): 경로를 포함하는 시료 디렉터리명을 원소로 가지는 리스트
                custom_name(dict): 변경할 이름을 담고 있는 딕션너리. Key: 현재 시료명, Value: 변경할 시료명.
                stat_file(str): stat.txt 파일 위치.
                nread_file(str): NreadPercent.csv 파일 위치.
                n_base(float, float): N base의 기준. 첫번째 인자(>): 주의, 두 번째 인자(>): 경고.
                n_read(float): N Read 비율의 기준.
                sample_read(int): 시료별 생산되어야 할 Read의 개수.
                q30[int, int]: Q30의 알람 기준치. 첫 번째 요소(<): 경고, 두 번째 요소(<): 주의.
                order_number: 수주번호.
                index_kit: auto, Nextera, TruSeq.
                sample_name_mode(str): 시료명 변경 모드(all, all_no, no_zero, no_mid).
                (삭제예정)run_path(str): RawData에 대한 Run 디렉터리 경로.
                (삭제예정)kit_file(str): index kit 정보가 있는 directory_path_kit.txt 파일 위치.

    :return: index_kit_info_from_rawdata, l_nt_sample_list
             - index_kit_info_from_rawdata(str): index kit 정보
             - l_nt_sample_list: SampleList를 원소로 가지는 리스트.
                    SampleList - path cur_name new_name 의 이름을 가지는 named tuple.

            # - dic_flowcell_id: 시료명이 Key이며, FlowCell ID를 Value로 가지는 딕션너리.
    """
    # dic_flowcell_id = check_flowcell_id(kargs['sample_list'])
    check_stat(kargs['stat_file'], kargs['n_base'], kargs['q30'], kargs['sample_read'])
    check_nread(kargs['nread_file'], kargs['n_read'])
    # index_kit_info_from_rawdata = check_index_kit_info(kargs['kit_file'], kargs['run_path'], kargs['library_kit'])
    order_num_type, json_data = get_order_json(kargs['order_number'])
    index_kit_info_from_rawdata = check_index_kit_info_using_json(json_data, kargs['index_kit'])
    check_output_file_requirements_from_json(json_data)
    l_nt_sample_list = check_sample_name(kargs['sample_list'], kargs['custom_name'], kargs['sample_name_mode'])
    # return dic_flowcell_id, index_kit_info_from_rawdata, l_nt_sample_list
    return index_kit_info_from_rawdata, l_nt_sample_list


def copy_rawdata(kargs, p_mode='main'):
    """
    fastq 파일들과 FastQC 디렉터리 또는 Zip 파일을 복사한다.

    :type p_mode: str
    :param p_mode: 통합분석 수주인지 확인. 생성되는 디렉터리 이름의 규칙이 다름.
                - main: 주 수주
                - integrate: 통합 수주

    :type kargs: dict
    :param kargs: 다음을 키(key)로 가지는 딕션너리(Dictionary)
                analysis_base_path(str): 분석 디렉터리 경로
                order_number(str): 수주번호
                run_path(str): RawData에 대한 Run 디렉터리 경로.
                sample_list(list): SampleList를 원소로 가지는 리스트.
                             SampleList - path cur_name new_name 의 이름을 가지는 named tuple.
                R1_suffix(str): Read1에 대한 fastq 파일명의 접미사.
                R2_suffix(str): Read2에 대한 fastq 파일명의 접미사.
                copy_fastqc(boolen): FastQC 데이터 복사 여부

    :return: l_nt_sample_list_in_analysis, l_nt_sample_list_in_analysis_for_fastqc
             namedtuple인 SampleList 를 원소로 가지는 리스트.
             SampleList.rawdata_path_in_analysis
                       .cur_name
                       .new_name
    """
    # 분석디렉터리에 해당 수주번호의 디렉터리 존재 유무 확인
    if p_mode == 'main':
        analysis_workdir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'])
        if check_file_type(analysis_workdir_path, 'exists'):
            secho('Warning: 분석디렉터리에 해당 수주번호 존재!', fg='yellow')
            echo(analysis_workdir_path)
        else:
            os.mkdir(analysis_workdir_path)
            os.mkdir(os.path.join(analysis_workdir_path, 'RawData'))
            echo('>>> 수주 분석 디렉터리 생성')
            echo(analysis_workdir_path + '\n')
        rawdata_path_in_analysis = os.path.join(analysis_workdir_path, 'RawData', os.path.split(kargs['run_path'])[1])
        os.mkdir(rawdata_path_in_analysis)
    elif p_mode == 'integrate':
        analysis_workdir_path = os.path.join(kargs['analysis_base_path'], 'RawData_' + kargs['order_number'])
        if check_file_type(analysis_workdir_path, 'exists'):
            secho('Warning: 분석디렉터리에 해당 수주번호의 RawData 디렉터리(통합분석) 존재!', fg='yellow')
            echo(analysis_workdir_path)
        else:
            os.mkdir(analysis_workdir_path)
            echo('>>> 수주 분석 디렉터리 생성')
            echo(analysis_workdir_path + '\n')
        rawdata_path_in_analysis = os.path.join(analysis_workdir_path, os.path.split(kargs['run_path'])[1])
        os.mkdir(rawdata_path_in_analysis)
    else:
        raise ValueError(style('p_mode: {}'.format(p_mode)))

    # TODO 데이터 복사 병렬처리
    l_nt_sample_list_in_analysis = list()
    l_nt_sample_list_in_analysis_for_fastqc = list()
    bar_label = '데이터{fastqc}복사'.format(fastqc=' & FastQC ' if kargs['copy_fastqc'] else '')
    with progressbar(kargs['sample_list'], label=bar_label,
                     fill_char=style('#', fg='green')) as bar_l_nt_sample_list:
        for nt_sample in bar_l_nt_sample_list:
            src_path = os.path.join(nt_sample.path, nt_sample.cur_name, nt_sample.cur_name)
            des_path = os.path.join(rawdata_path_in_analysis, nt_sample.new_name, nt_sample.new_name)
            os.mkdir(os.path.split(des_path)[0])
            shutil.copyfile(src_path + kargs['R1_suffix'], des_path + kargs['R1_suffix'])
            shutil.copyfile(src_path + kargs['R2_suffix'], des_path + kargs['R2_suffix'])
            if kargs['copy_fastqc']:
                r1_status = copy_fastqc(src_path, des_path, 1)
                r2_status = copy_fastqc(src_path, des_path, 2)
                if (r1_status is False) or (r2_status is False):
                    l_nt_sample_list_in_analysis_for_fastqc.append(
                        SampleList(
                            rawdata_path_in_analysis,
                            nt_sample.cur_name,
                            nt_sample.new_name,
                        )
                    )

            l_nt_sample_list_in_analysis.append(
                SampleList(
                    rawdata_path_in_analysis,
                    nt_sample.cur_name,
                    nt_sample.new_name,
                )
            )
    return l_nt_sample_list_in_analysis, l_nt_sample_list_in_analysis_for_fastqc


def copy_fastqc(p_src, p_des, p_read):
    """
    RawData 보관 디렉터리에 있는 FastQC 디렉터리 또는 Zip파일을 복사한다.

    :type p_src: str
    :param p_src: 원본 경로

    :type p_des: str
    :param p_des: 복사 경로

    :type p_read: int
    :param p_read: Paired End의 Read 방향 표시 (1 or 2)

    :return: True, False
    """
    dir_suffix = '_{read}_fastqc'.format(read='1' if p_read == 1 else '2')
    file_suffix = '_{read}_fastq.zip'.format(read='1' if p_read == 1 else '2')
    src_dir = p_src + dir_suffix
    src_file = p_src + file_suffix
    des_dir = p_des + dir_suffix
    des_file = p_des + file_suffix
    if check_file_type(src_dir, 'exists') and check_file_type(src_dir, 'isdir'):
        shutil.copytree(src_dir, des_dir)
        return True
    elif check_file_type(src_file, 'exists') and check_file_type(src_dir, 'isfile'):
        shutil.copyfile(src_file, des_file)
        return True
    else:
        return False


def get_fastqc_cmd(p_output_path, p_fastq):
    """

    :type p_output_path: str
    :param p_output_path: 결과 출력 경로.

    :type p_fastq: str
    :param p_fastq: FASTQ 파일.

    :return: fastqc_cmd
    """
    # TODO : RawData --> v0.11.6 사용
    # FastQC - v0.11.8
    # """
    #             FastQC - A high throughput sequence QC analysis tool
    #
    # SYNOPSIS
    #
    # 	fastqc seqfile1 seqfile2 .. seqfileN
    #
    #     fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
    #            [-c contaminant file] seqfile1 .. seqfileN
    #
    # DESCRIPTION
    #
    #     FastQC reads a set of sequence files and produces from each one a quality
    #     control report consisting of a number of different modules, each one of
    #     which will help to identify a different potential type of problem in your
    #     data.
    #
    #     If no files to process are specified on the command line then the program
    #     will start as an interactive graphical application.  If files are provided
    #     on the command line then the program will run with no user interaction
    #     required.  In this mode it is suitable for inclusion into a standardised
    #     analysis pipeline.
    #
    #     The options for the program as as follows:
    #
    #     -h --help       Print this help file and exit
    #
    #     -v --version    Print the version of the program and exit
    #
    #     -o --outdir     Create all output files in the specified output directory.
    #                     Please note that this directory must exist as the program
    #                     will not create it.  If this option is not set then the
    #                     output file for each sequence file is created in the same
    #                     directory as the sequence file which was processed.
    #
    #     --casava        Files come from raw casava output. Files in the same sample
    #                     group (differing only by the group number) will be analysed
    #                     as a set rather than individually. Sequences with the filter
    #                     flag set in the header will be excluded from the analysis.
    #                     Files must have the same names given to them by casava
    #                     (including being gzipped and ending with .gz) otherwise they
    #                     won't be grouped together correctly.
    #
    #     --nano          Files come from nanopore sequences and are in fast5 format. In
    #                     this mode you can pass in directories to process and the program
    #                     will take in all fast5 files within those directories and produce
    #                     a single output file from the sequences found in all files.
    #
    #     --nofilter      If running with --casava then don't remove read flagged by
    #                     casava as poor quality when performing the QC analysis.
    #
    #     --extract       If set then the zipped output file will be uncompressed in
    #                     the same directory after it has been created.  By default
    #                     this option will be set if fastqc is run in non-interactive
    #                     mode.
    #
    #     -j --java       Provides the full path to the java binary you want to use to
    #                     launch fastqc. If not supplied then java is assumed to be in
    #                     your path.
    #
    #     --noextract     Do not uncompress the output file after creating it.  You
    #                     should set this option if you do not wish to uncompress
    #                     the output when running in non-interactive mode.
    #
    #     --nogroup       Disable grouping of bases for reads >50bp. All reports will
    #                     show data for every base in the read.  WARNING: Using this
    #                     option will cause fastqc to crash and burn if you use it on
    #                     really long reads, and your plots may end up a ridiculous size.
    #                     You have been warned!
    #
    #     --min_length    Sets an artificial lower limit on the length of the sequence
    #                     to be shown in the report.  As long as you set this to a value
    #                     greater or equal to your longest read length then this will be
    #                     the sequence length used to create your read groups.  This can
    #                     be useful for making directly comaparable statistics from
    #                     datasets with somewhat variable read lengths.
    #
    #     -f --format     Bypasses the normal sequence file format detection and
    #                     forces the program to use the specified format.  Valid
    #                     formats are bam,sam,bam_mapped,sam_mapped and fastq
    #
    #     -t --threads    Specifies the number of files which can be processed
    #                     simultaneously.  Each thread will be allocated 250MB of
    #                     memory so you shouldn't run more threads than your
    #                     available memory will cope with, and not more than
    #                     6 threads on a 32 bit machine
    #
    #     -c              Specifies a non-default file which contains the list of
    #     --contaminants  contaminants to screen overrepresented sequences against.
    #                     The file must contain sets of named contaminants in the
    #                     form name[tab]sequence.  Lines prefixed with a hash will
    #                     be ignored.
    #
    #     -a              Specifies a non-default file which contains the list of
    #     --adapters      adapter sequences which will be explicity searched against
    #                     the library. The file must contain sets of named adapters
    #                     in the form name[tab]sequence.  Lines prefixed with a hash
    #                     will be ignored.
    #
    #     -l              Specifies a non-default file which contains a set of criteria
    #     --limits        which will be used to determine the warn/error limits for the
    #                     various modules.  This file can also be used to selectively
    #                     remove some modules from the output all together.  The format
    #                     needs to mirror the default limits.txt file found in the
    #                     Configuration folder.
    #
    #    -k --kmers       Specifies the length of Kmer to look for in the Kmer content
    #                     module. Specified Kmer length must be between 2 and 10. Default
    #                     length is 7 if not specified.
    #
    #    -q --quiet       Supress all progress messages on stdout and only report errors.
    #
    #    -d --dir         Selects a directory to be used for temporary files written when
    #                     generating report images. Defaults to system temp directory if
    #                     not specified.
    #
    # """
    fastqc_cmd = \
        '{fastqc} ' \
        '--extract ' \
        '--nogroup ' \
        '-t 10 ' \
        '-o {output_path} {fastq}'.format(
            fastqc=FASTQC,
            output_path=p_output_path,
            fastq=p_fastq
        )
    return fastqc_cmd


def run_fastqc(p_list, p_r1_suffix, p_r2_suffix):
    """
    FastQC 실행.

    :type p_list: list
    :param p_list: namedtuple인 SampleList 를 원소로 가지는 리스트.
             SampleList.rawdata_path_in_analysis
                       .cur_name
                       .new_name

    :type p_r1_suffix: str
    :param p_r1_suffix:

    :type p_r2_suffix: str
    :param p_r2_suffix:

    :return: None
    """
    echo('>>> FastQC 시작')
    echo([x.new_name for x in p_list])
    l_cmd = list()
    for nt_sample in p_list:
        src_path = os.path.join(*(nt_sample.path, nt_sample.new_name))
        r1_file = '{path}/{name}{suffix}'.format(
            path=src_path,
            name=nt_sample.new_name,
            suffix=p_r1_suffix)
        r2_file = '{path}/{name}{suffix}'.format(
            path=src_path,
            name=nt_sample.new_name,
            suffix=p_r2_suffix)
        for file in r1_file, r2_file:
            cmd = get_fastqc_cmd(src_path, file)
            l_cmd.append(cmd)

    launcher_cmd(l_cmd, 'FastQC', POOL_WORKER, True)


def get_multiqc_cmd(kargs):
    """
    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
    
    :return: multiqc_cmd
    """
    # MultiQC 1.6
    # """
    # Usage: multiqc [OPTIONS] <analysis directory>
    #
    #   MultiQC aggregates results from bioinformatics analyses across many
    #   samples into a single report.
    #
    #   It searches a given directory for analysis logs and compiles a HTML
    #   report. It's a general use tool, perfect for summarising the output from
    #   numerous bioinformatics tools.
    #
    #   To run, supply with one or more directory to scan for analysis results. To
    #   run here, use 'multiqc .'
    #
    #   See http://multiqc.info for more details.
    #
    #   Author: Phil Ewels (http://phil.ewels.co.uk)
    #
    # Options:
    #   -f, --force                     Overwrite any existing reports
    #   -d, --dirs                      Prepend directory to sample names
    #   -dd, --dirs-depth INTEGER       Prepend [INT] directories to sample names.
    #                                   Negative number to take from start of path.
    #   -s, --fullnames                 Do not clean the sample names (leave as full
    #                                   file name)
    #   -i, --title TEXT                Report title. Printed as page header, used
    #                                   for filename if not otherwise specified.
    #   -b, --comment TEXT              Custom comment, will be printed at the top
    #                                   of the report.
    #   -n, --filename TEXT             Report filename. Use 'stdout' to print to
    #                                   standard out.
    #   -o, --outdir TEXT               Create report in the specified output
    #                                   directory.
    #   -t, --template [default|default_dev|geo|sections|simple]
    #                                   Report template to use.
    #   --tag TEXT                      Use only modules which tagged with this
    #                                   keyword, eg. RNA
    #   --view-tags, --view_tags        View the available tags and which modules
    #                                   they load
    #   -x, --ignore TEXT               Ignore analysis files (glob expression)
    #   --ignore-samples TEXT           Ignore sample names (glob expression)
    #   --ignore-symlinks               Ignore symlinked directories and files
    #   --sample-names PATH             File containing alternative sample names
    #   -l, --file-list                 Supply a file containing a list of file
    #                                   paths to be searched, one per row
    #   -e, --exclude [module name]     Do not use this module. Can specify multiple
    #                                   times.
    #   -m, --module [module name]      Use only this module. Can specify multiple
    #                                   times.
    #   --data-dir                      Force the parsed data directory to be
    #                                   created.
    #   --no-data-dir                   Prevent the parsed data directory from being
    #                                   created.
    #   -k, --data-format [tsv|json|yaml]
    #                                   Output parsed data in a different format.
    #                                   Default: tsv
    #   -z, --zip-data-dir              Compress the data directory.
    #   -p, --export                    Export plots as static images in addition to
    #                                   the report
    #   -fp, --flat                     Use only flat plots (static images)
    #   -ip, --interactive              Use only interactive plots (HighCharts
    #                                   Javascript)
    #   --lint                          Use strict linting (validation) to help code
    #                                   development
    #   --pdf                           Creates PDF report with 'simple' template.
    #                                   Requires Pandoc to be installed.
    #   --no-megaqc-upload              Don't upload generated report to MegaQC,
    #                                   even if MegaQC options are found
    #   -c, --config PATH               Specific config file to load, after those in
    #                                   MultiQC dir / home dir / working dir.
    #   --cl-config, --cl_config TEXT   Specify MultiQC config YAML on the command
    #                                   line
    #   -v, --verbose                   Increase output verbosity.
    #   -q, --quiet                     Only show log warnings
    #   --version                       Show the version and exit.
    #   -h, --help                      Show this message and exit.
    # """
    multiqc_cmd = \
        '{multiqc} ' \
        '-d {target_dir} ' \
        '-dd {depth} ' \
        '-i {order_number} ' \
        '-o {output_path} ' \
        '-m fastp ' \
        '-m fastqc ' \
        '-b "내부용입니다. 보안(서버내 경로)상 고객에게 전달되면 안됩니다." ' \
        '--no-data-dir'.format(
            multiqc=MULTIQC,
            target_dir=kargs['target_dir'],
            depth=kargs['depth'],
            order_number=kargs['order_number'],
            output_path=kargs['output_path'],
        )
    return multiqc_cmd


def run_multiqc(kargs):
    """

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                target_dir: 대상 디렉터리명, 또는 보고서에 출력할 시료명의 접두어.
                order_number: 수주번호.
                depth: 시료명에 포함할 접두어의 깊이.
                output_path: 출력 결과 경로.
    :return: None
    """
    cmd = get_multiqc_cmd(kargs)
    echo('>>> MultiQC 시작')
    run = run_cmd(cmd, p_stdout=None, p_stderr=None)
    run.stderr = '상단에 출력'
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'MultiQC 생성 완료',
            'false_meg': 'MultiQC 작업 오류'
        }
    )


def run_qc():
    #TODO: QC 실행 코드 작성
    pass

