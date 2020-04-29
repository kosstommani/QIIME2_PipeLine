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
# |   __|___ ___|     |   | |
# |__   | . | . |  |  | | | |
# |_____|  _|___|_____|_|___|
#       |_|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.1.2'  # 2020.04.24

from click import secho, echo, style
import glob
import subprocess
import os
from shutil import rmtree
from multiprocessing import current_process, Pool


def glob_dir(p_path, p_pattern, p_mode='only', p_verbose=True):
    """
    주어진 경로에서 패턴에 일치하는 파일이나 디렉터리명을 찾음.
    mode가 'only'인 경우 패턴과 일치하는 1개의 대상만 선택.

    :param p_path: 디렉터리 경로
    :param p_pattern: 찾고자하는 파일 또는 디렉터리의 패턴
    :param p_mode: 패턴과 일치하는 대상 선택. only - 1개(str), many - 다수(list)
    :param p_verbose: 미검출 여부 출력(only mode 에서)
    :return:
        default(raw) : 패턴 일치 -  디렉터리명 또는 파일명 1개(문자열)
                       패턴 불일치 - None
        경로입력시 : 패턴 일치 -  디렉터리명 또는 파일명의 리스트(list)
                    패턴 불일치 - None
    """
    path = '{path}/{pattern}'.format(path=p_path, pattern=p_pattern)
    if p_mode == 'only':
        glob_list = glob.glob(path)
        glob_list.sort()
        if len(glob_list) == 1:
            if p_verbose:
                secho('>>> 패턴 일치', fg='blue')
                echo('{0}'.format(glob_list[0]))
            return glob_list[0]
        elif len(glob_list) > 1:
            secho('>>> {0} 패턴과 일치 목록'.format(p_pattern), fg='yellow')
            for num, element in enumerate(glob_list, 1):
                secho('{0}. {1}'.format(num, element))
            glob_list = glob_list[int(input('분석에 사용할 디렉터리 선택? ')) - 1]
            echo('>>> 선택완료')
            echo(glob_list + "\n")
            return glob_list
        else:
            if p_verbose:
                secho('Error : {0} 패턴과 일치하는 목록이 없음.'.format(p_pattern), fg='red', blink=True)
            return None

    elif p_mode == 'many':
        glob_list = glob.glob(path)
        if glob_list:
            return glob_list
        else:
            return None
    else:
        secho('Error : 잘못된 파라미터!', fg='red', blink=True, err=True)
        exit(1)


def check_file_type(p_list, p_type):
    """
    주어진 목록(경로포함)의 대상들이 특정 유형인지 확인.

    :type p_list: list or str
    :param p_list: 경로가 포함된 목록(리스트, 문자열)
    :param p_type: 확인하고자 하는 유형.
                    isabs : 절대경로?
                    isfile : 파일?
                    isdir : 디렉터리?
                    islink : 링크파일?
                    ismount : 마운트포인트?
                    exists : 존재?
                    lexists : 링크존재?
    :return: p_list가 list 또는 tuple인 경우 지정한 유형의 대상를 목록과 대상이 아닌 목록을 반환.
                반환 타입 : list
                type_list, not_type_list
             p_list가 문자열인 경우 True, False 반환
    """
    type_list = list()
    not_type_list = list()
    if type(p_list) is list:
        for i in p_list:
            type_state = eval('os.path.{type}("{dir}")'.format(type=p_type, dir=i))
            if type_state:
                type_list.append(i)
            else:
                not_type_list.append(i)
        return type_list, not_type_list
    else:
        type_state = eval('os.path.{type}("{dir}")'.format(type=p_type, dir=p_list))
        return type_state


def check_order_number_system(p_order_number):
    """

    :param p_order_number:
    :return:
    """
    import re
    lims2_re = re.compile(r'^\d{4}\S{3}-\d{4}$')
    lims3_re = re.compile(r'^\S{2}\d{8}$')
    lims2_result = lims2_re.findall(p_order_number)
    lims3_result = lims3_re.findall(p_order_number)
    lims2_num = len(lims2_result)
    lims3_num = len(lims3_result)
    if (lims2_num == 1) and (lims3_num == 0):
        return 'LIMS2'
    elif (lims2_num == 0) and (lims3_num == 1):
        return 'LIMS3'
    else:
        secho('Error: 수주번호에 문제가 있습니다.', fg='red', err=True)
        echo('입력된 수주번호 : {}'.format(p_order_number), err=True)
        echo('수주번호 예시', err=True)
        echo('\tLIMS2(ex: 1812KMI-0001)', err=True)
        echo('\tLIMS3(ex: HN00100001)', err=True)
        exit(1)


def get_order_json_from_json_server(p_url):
    """

    :param p_url:
    :return:
    """
    import socket

    max_size = 1024
    address = ('172.19.85.10', 8877)
    echo('>>> json server 연결')
    client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        client.connect(address)
    except (ConnectionRefusedError, OSError) as err:
        secho('Warning: json server와 연결이 되지 않습니다.', fg='magenta', blink=True)
        echo(err)
        json_text = 'ConnectionRefusedError'
    else:
        client.sendall(p_url)
        total_text = b''
        while True:
            data = client.recv(max_size)
            if data:
                total_text += data
                continue
            else:
                break
        echo('json 정보를 받았습니다.')
        json_text = total_text.strip()
    finally:
        client.close()
        echo('>>> json server 연결 해제')
    return json_text


def get_order_json(p_order_number, p_service='NGS'):
    """
    LIMS에서 해당 수주번호에 대한 JSON 데이터를 가져옴.
    인터넷 연결 필수!

    :param p_order_number: 수주번호
    :param p_service: [NSG | myBiomeStory] 서비스 종류.
    :rtype: tuple
    :return: (order_num_type, json_data)
             order_num_type: LIMS2 or LIMS3
             json_data
                - JSON 파일이 있는 경우 - json data의 Dictionary
                - JSON 파일이 없는 경우 - None
    """
    from urllib import request
    import ssl
    import json
    import platform

    order_num_type = check_order_number_system(p_order_number)
    if p_service == 'NGS':
        if order_num_type == 'LIMS2':
            json_prefix = 'https://lims.macrogen.com/admin/ngs2/analysis/ngsDataMod2.jsp?on='
        elif order_num_type == 'LIMS3':
            json_prefix = 'http://lims02.macrogen.com/LIMSMOD/ngsDataMod2.jsp?on='
    elif p_service == 'myBiomeStory':
        json_prefix = 'https://lims3.macrogen.com/ngs/openapi/demulti/retrieveMiBioCsvInfo.do?cmpnyCd=1000&ordNo='
    else:
        raise ValueError(f'p_service: {p_service}')

    # 실행 서버의 인터넷 연결 여부에 따라 LIMS JSON 파일 획득 방법 다름.
    # 인터넷 연결 여부와 상관 없음. - 2020.03.10
    if platform.node() in ['denovo06', 'denovo07', 'denovo08', 'denovo09', 'denovo10']:
        response = request.urlopen(json_prefix + p_order_number, context=ssl._create_unverified_context())
        if response:
            try:
                json_data = json.load(response)
            except json.decoder.JSONDecodeError as err:
                secho('Error: JSON DecodeError 발생', fg='red')
                echo(f'내용: {err}')
                echo('JSON 파일에 문제가 있어 데이터를 변환할 수 없습니다.')
                secho(' --> JSON 파일 파싱을 생략합니다.', fg='magenta')
                return order_num_type, None
            else:
                return order_num_type, json_data
        else:
            secho('Warning : 수주에 대한 JSON 파일을 읽어 올 수 없습니다.', fg='magenta')
            return order_num_type, None
    else:
        response = get_order_json_from_json_server(json_prefix.encode('utf-8') + p_order_number.encode('utf-8'))
        if response and response != 'ConnectionRefusedError':
            try:
                json_data = json.loads(response.decode('utf-8'))
            except UnicodeDecodeError as err:
                json_data = json.loads(response.decode('utf-16'))
            return order_num_type, json_data
        else:
            secho('Warning : 수주에 대한 JSON 파일을 읽어 올 수 없습니다.', fg='magenta', blink=True)
            return order_num_type, None


# TODO : demux key 정보 파싱 - 재런 또는 재디멀티일 경우 데이터 구조가 어떻게 되는지 확인 필요.
def make_info_file(p_order_number, p_output, p_sample, p_target_region):
    """
    OTU분석보고서(HTML)에 기입될 고객정보를 담고 있는 info.txt 파일을 생성한다.

    :type p_order_number: str
    :param p_order_number: 수주번호

    :type p_output: str
    :param p_output: info.txt 파일이 생성될 경로.

    :type p_sample: int
    :param p_sample: OTU분석에 사용되는 시료의 개수.

    :type p_target_region: str
    :param p_target_region: region & primer set 정보(있는 경우)
                     ex) (V3V4(Bakt_341F-805R)
    :return: None
    """
    order_type, data = get_order_json(p_order_number)
    if data is None:
        secho('Warning: info.txt 파일을 생성할 수 없습니다.', fg='magenta', blink=True)
        return
    customer = data['customer']
    info_format = """organization: {organization}
name: {name}
order_number: {order_number}
sample_count: {sample_count}
target_region: {target_region}
e_mail: {e_mail}
""".format(organization=customer['organization'] if order_type == 'LIMS3' else customer['organiztion'],
           name=customer['name'],
           order_number=p_order_number,
           sample_count='{}Samples'.format(p_sample) if p_sample > 1 else '{}Sample'.format(p_sample),
           target_region=p_target_region,
           e_mail='ngskr@macrogen.com' if 'korea' in customer['country'].lower() else 'ngs@macrogen.com')
    with open(os.path.join(p_output, 'info.txt'), 'w') as o_output:
        o_output.write(info_format)
    secho('>>> info.txt 생성 완료', fg='cyan')
    return


def read_data(p_file, p_sep=None):
    """
    행단위로 Text 파일을 읽음.

    :param p_file: 입력 파일
    :param p_sep: 행단위 텍스트를 나눌 구분자.
    :return: 행단위 데이터를 구분자로 구분한 원소를 리스트로 가지고, 해당 리스트를 원소로 가지는 리스트 (중첩 리스트)
    """
    data = list()
    with open(p_file, 'r') as o_input:
        for text in o_input:
            data.append(text.strip().split(p_sep))
    echo('>>> 파일 읽기 완료')
    echo(p_file)
    return data


def read_all_data(p_file, p_verbose=True):
    with open(p_file, 'r') as o_input:
        data = o_input.read()
    if p_verbose:
        echo('>>> 파일 읽기 완료')
        echo(p_file)
    return data


def typeWrite(p_worksheet, p_row, p_col, p_string, p_int_format=None):
    """
    xlsxwriter

    :param p_worksheet:
    :param p_row:
    :param p_col:
    :param p_string:
    :param p_int_format:
    :return:
    """
    try:
        data = int(p_string)
        p_worksheet.write_number(p_row, p_col, data, p_int_format)
    except ValueError:
        try:
            data = float(p_string)
            p_worksheet.write_number(p_row, p_col, data)
        except ValueError:
            p_worksheet.write_string(p_row, p_col, p_string)


def typeWrite_float(p_worksheet, p_row, p_col, p_string, p_int_format=None):
    try:
        data = float(p_string)
        p_worksheet.write_number(p_row, p_col, data, p_int_format)
    except ValueError:
        p_worksheet.write_string(p_row, p_col, p_string)


def cast_num(p_text):
    """
    문자형태의 숫자를 적절한 형태의 숫자형 자료로 변환한다.
    text -> int -> float -> str

    :param p_text:
    :return:
    """
    try:
        num = int(p_text)
        num_type = int
    except ValueError:
        try:
            num = float(p_text)
            num_type = float
        except ValueError:
            num = p_text
            num_type = str
    return num_type, num


def run_cmd(p_cmd, p_stdout=subprocess.PIPE, p_stderr=subprocess.PIPE, p_shell=True):
    """
    :type p_cmd: str
    :param p_cmd: 실행 명령어

    :param p_stdout:
    :param p_stderr:
    :param p_shell:

    :return: completed_run
    """
    # python 3.6+
    completed_run = subprocess.run(p_cmd, shell=p_shell, stdout=p_stdout, stderr=p_stderr, encoding='utf-8')
    # python 3.5
    # completed_run = subprocess.run(p_cmd, shell=True, stdout=p_stdout, stderr=p_stderr)
    return completed_run


def launcher_cmd(p_l_cmd, p_meg, p_pool_worker, p_exit):
    """
    외부 명령어들을 병렬로 실행시키고 완료 여부를 확인하는 일련의 과정들을 진행한다.

    :type p_l_cmd: list
    :param p_l_cmd: 실행 명령어(str)를 원소로 가지는 리스트.
    :type p_meg: str
    :param p_meg: 출력 메세지.
    :type p_pool_worker: int
    :param p_pool_worker: 프로세스 개수.
    :type p_exit: bool
    :param p_exit: Eror 발생시 종료 여부.
    :rtype: bool
    :return: 0 - 실행된 모든 명령어가 완료됨.
             1 - 실행된 명령 중 실패한 명령어 존재함.
    """
    process = p_pool_worker
    if len(p_l_cmd) < p_pool_worker:
        process = len(p_l_cmd)
    echo(f'>>> {p_meg} 시작')
    with Pool(process, initializer=start_process) as pool:
        pool_outputs = pool.map(run_cmd, p_l_cmd)
        pool.close()
        pool.join()

    done_count = len([1 for run in pool_outputs if run.returncode == 0])
    error_count = len([1 for run in pool_outputs if run.returncode != 0])
    done_text = style('완료: {}'.format(done_count), fg='cyan')
    error_text = style('에러: {}'.format(error_count),
                       fg='red' if (error_count != 0) else 'cyan',
                       blink=True if (error_count != 0) else False)
    echo('{msg} {done}, {error}'.format(
        msg=style(f'>>> {p_meg}', fg='cyan'),
        done=done_text,
        error=error_text))
    if error_count != 0:
        for run in pool_outputs:
            if run.returncode == 0:
                continue
            else:
                check_run_cmd({
                    'run': run,
                    'true_meg': None,
                    'false_meg': f'{p_meg} 실행'
                }, p_exit=False)
        if p_exit is True:
            echo('Error: check stdout Message', err=True)
            exit(1)
        else:
            return 1
    return 0


def run_move_and_cmd(p_path, p_cmd, p_stdout=subprocess.PIPE, p_stderr=subprocess.PIPE):
    """
    실행할 명령어(스크립트)가 위치하는 경로로 이동한 후, 해당 명령어(스크립트를) 실행한다.
    실행 후, 이동하기 전의 경로로 돌아온다.

    :type p_path: str
    :param p_path: 이동할 경로(명령어를 실행할 경로).

    :type p_cmd: str
    :param p_cmd: 실행할 명령어. p_path 경로에 해당 명령어 또는 관련 파일들이 존재해야 됨.
    :param p_stdout:
    :param p_stderr:
    :return: completed_run
    """
    cur_dir = os.getcwd()
    os.chdir(p_path)
    completed_run = subprocess.run(p_cmd, shell=True, stdout=p_stdout, stderr=p_stderr)
    os.chdir(cur_dir)
    return completed_run


def check_run_cmd(p_args, p_exit=True, p_stdout=True, p_stderr=True, p_t_color='cyan', p_f_color='red'):
    """
    실행된 명령어가 정상적으로 완료되었는지를 확인하고, 관련 메시지를 출력한다.
    실행된 명령어가 비정상적으로 완료되었을 경우(에러 발생)에는 관련 메시지 및 실행 명령어가 출력된다.

    :param p_args: 다음을 key로 가지는 Dictionary
                run:
                true_meg: 완료 메시지. None이면 미출력
                false_meg: 실패 메시지. None이면 미출력
    :param p_exit: 프로그램 종료 여부.
    :param p_stdout: stdout 출력 여부
    :param p_stderr: stderr 출력 여부
    :param p_t_color: 완료 메시지의 색상 
    :param p_f_color: 실패 메시지의 색상
    :return:
    """
    run = p_args['run']
    if run.returncode == 0:
        if p_args['true_meg'] is not None:
            secho('>>> {meg}'.format(meg=p_args['true_meg']), fg=p_t_color)
            if p_stdout:
                echo(run.stdout)
    elif run.returncode != 0:
        if p_args['false_meg'] is not None:
            secho('Error: {meg}'.format(meg=p_args['false_meg']), fg=p_f_color, blink=True, err=True)
            echo('------ Error Message ------', err=True)
            if p_stderr:
                echo(run.stderr, err=True)
            echo('----------- CMD -----------', err=True)
            echo(run.args, err=True)
            if p_exit:
                echo('Error: check stdout Message', err=True)
                exit(1)
    else:
        raise RuntimeError('예기치 않은 에러 발생! 개발자에게 문의하세요.')


def parse_config(p_mode='config'):
    """

    :param p_mode:
    :return:
    """
    import yaml
    if p_mode == 'config':
        if check_file_type(os.path.expanduser('~/MetaMix_config.yaml'), 'exists'):
            config_yaml = os.path.expanduser('~/MetaMix_config.yaml')
        elif check_file_type(os.path.expanduser('~/config.yaml'), 'exists'):
            config_yaml = os.path.expanduser('~/config.yaml')
        else:
            config_yaml = '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/config.yaml'
    elif p_mode == 'primer':
        if check_file_type(os.path.expanduser('~/MetaMix_primer.yaml'), 'exists'):
            config_yaml = os.path.expanduser('~/MetaMix_primer.yaml')
        elif check_file_type(os.path.expanduser('~/primer.yaml'), 'exists'):
            config_yaml = os.path.expanduser('~/primer.yaml')
        else:
            config_yaml = '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/primer.yaml'

    with open(config_yaml, 'r') as o_config_yaml:
        config = yaml.safe_load(o_config_yaml)
    return config


def run_dir_tree(p_path, p_dir, p_level=3):
    """

    :param p_path:
    :param p_dir:
    :param p_level:
    :return:
    """
    TREE = '/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/tree'
    target_path = os.path.join(p_path, p_dir)
    cmd = '{tree} -L {level} -C {path}'.format(tree=TREE, level=p_level, path=target_path)
    run = run_cmd(cmd)
    if run.returncode == 0:
        text = run.stdout
        if 'error opening dir' in text:
            secho('Error: {meg}'.format(meg='목표를 찾을 수 없습니다.'), fg='red', blink=True, err=True)
            echo('------ Error Message ------', err=True)
            echo(run.stdout, err=True)
            echo('----------- CMD -----------', err=True)
            echo(run.args, err=True)
            exit(1)
        else:
            secho('>>> {meg}'.format(meg='목표위치의 디렉토리 구조'))
            echo(run.stdout)
    elif run.returncode == 1:
        secho('Error: {meg}'.format(meg='명령어 오류.'), fg='red', blink=True, err=True)
        echo('------ Error Message ------', err=True)
        echo(run.stderr, err=True)
        echo('----------- CMD -----------', err=True)
        echo(run.args, err=True)
        exit(1)
    else:
        raise RuntimeError('예기치 않은 에러 발생! 개발자에게 문의하세요.')


def delete_dir(p_path, p_dir):
    """

    :param p_path:
    :param p_dir:
    :return:
    """
    target_path = os.path.join(p_path, p_dir)
    rmtree(target_path)
    secho('>>> 목표 디렉터리 삭제 완료', fg='cyan')
    echo(target_path)


def get_file_copy_cmd(p_src, p_des):
    cmd = 'cp {} {}'.format(p_src, p_des)
    return cmd


def start_process():
    echo('Starting - {}'.format(current_process().name))


def get_upper_path(p_path, p_n):
    """
    경로 또는 경로를 포함하는 파일명이 주어졌을 때, 특정한 위치의 상위 경로를 반환한다.

    :type p_path: str
    :param p_path: 경로 또는 경로를 포함하는 파일명
    :type p_n: int
    :param p_n: 현재위치 기준에서 얻고자하는 경로(상위)의 상대적 위치거리. ex) /A/B/C. p_n=2이면, /A
    :rtype: str
    :return: 제시된 위치거리에 해당되는 경로
    """
    if p_n == 0:
        return p_path
    else:
        dir_name = os.path.dirname(p_path)
        return get_upper_path(dir_name, p_n - 1)


def make_dir_using_input_file(p_target_file, p_new_dir, p_upper_n, p_check=True):
    """
    주어진 파일명을 이용하여 절대경로를 획득하고, 원하는 위치에 새로운 디렉터리를 생성한다.
    생성된 디렉터리명을 포함하는 경로를 반환한다.
    생성하고자 하는 디렉터리가 있는 경우 p_check Warning 메시지를 출력한다.
    ex) make_dir_using_input_file('/A/B/C/D/test.txt', 'New', 3)
        --> /A/B/New 디렉터리 생성. 경로 반환.

    :type p_target_file: str
    :param p_target_file: 파일명 또는 경로.
    :type p_new_dir: str
    :param p_new_dir: 생성하고자하는 디렉터리명
    :type p_upper_n: int
    :param p_upper_n: 새로운 디렉터리를 생성할 상위 디렉터리와의 위치거리.
    :type p_check: bool
    :param p_check: 생성할 디렉터리 존재 유무 확인. True: 존재시 Error 발생. False: Warning 메시지 출력.
    :rtype: str
    :return: 생성된 디렉터리명을 포함하는 경로
    """
    if os.path.isabs(p_target_file):
        target = p_target_file
    else:
        target = os.path.abspath(p_target_file)
    path = get_upper_path(target, p_upper_n)
    new_path = os.path.join(path, p_new_dir)
    try:
        os.mkdir(new_path)
    except FileExistsError as err:
        if p_check is True:
            secho(f'Error: {p_new_dir} 디렉터리가 이미 있습니다.', fg='red', err=True, blink=True)
            echo(err, err=True)
            exit(1)
        else:
            secho(f'Warning: {p_new_dir} 디렉터리가 이미 있습니다.', fg='yellow')
            echo(new_path)
    else:
        echo(f'>>> {p_new_dir} 디렉터리 생성')
        echo(new_path)
    return new_path


def parse_html(p_html):
    from bs4 import BeautifulSoup
    with open(p_html, 'r') as o_html:
        html_text = o_html.read()
    html = BeautifulSoup(html_text, 'html5lib')
    return html


def read_metadata_for_sort(p_metadata):
    """
    고객맞춤정렬을 위해 metadata.txt 파일을 읽어서, 시료명을 리스트로 반환한다.

    :type p_metadata: str
    :param p_metadata: metadata.txt
    :return: 시료명을 원소로하는 리스트
    """
    metadata = read_data(p_metadata)
    metadata_order = [x[0] for x in metadata[1:]]  # Header 제외
    return metadata_order


def sort_data_by_custom(p_data, p_order, p_key=0):
    """
    입력받은 데이터의 목록을 정해진 순서대로 정렬한 후, 정렬된 데이터를 반환한다.(고객맞춤정렬)
    Data: [[Key2, Any, Any], [Key1, Any, Any]]
    Order: [Key1, Key2]
    Out: [[Key1, Any, Any], [Key2, Any, Any]]

    :type p_data: list
    :param p_data: 정렬을 원하는 데이터
    :type p_order: list
    :param p_order: 정렬 순서
    :type p_key: int
    :param p_key: 정렬을 위한 키 위치.
    :return: 맞춤정렬된 데이터(리스트)를 포함하는 리스트
    """
    d_data = dict()
    sorted_data = list()
    for ele in p_data:
        d_data[ele[p_key]] = ele
    for ele in p_order:
        sorted_data.append(d_data[ele])
    return sorted_data


def get_target_dir_number(p_dir, p_target):
    """
    경로에 존재하는 특정 대상을 검출하여 번호를 매긴 후 반환한다.
    대상이 없을 경우 _1 이 추가되고, 대상이 있을 경우 +1 적용된 이름이 반환된다.
    ex) p_target: Report --> Report_[0-9]* --> Report_1 or Report_2 등

    :type p_dir: str
    :param p_dir: 검출 대상 경로
    :type p_target: str
    :param p_target: 검출 대상
    :return
    """
    l_dir_name = glob_dir(p_dir, f'{p_target}_[0-9]*', 'many', False)
    if l_dir_name is None:
        return f'{p_target}_1'
    else:
        new_num = max([int(x.split(f'{p_target}_')[-1]) for x in l_dir_name]) + 1
        return f'{p_target}_{new_num}'


def check_queue_node() -> tuple:
    """
    작업 중인 서버가 Queue 사용이 가능한 서버인지 확인한다.

    :return: bool, str
    """
    import platform
    node = platform.node()
    if node in ['denovo06', 'denovo10', 'cm456']:
        return True, node
    else:
        return False, node


def run_cmd_qsub(kargs: dict):
    """

    :param kargs: dict
            stream: [y(es)|n(o)]merge stdout and stderr stream of job
            queue: [bi3m.q|meta.q]
            interpreter: command interpreter to be used
            cmd: command or file
            out_path: qsub.run, qsub.log, qsub.recipe 저장 위치
    :return:
    """
    qsub_out = os.path.join(kargs['out_path'], 'qsub.run')
    qsub_log = os.path.join(kargs['out_path'], 'qsub.log')
    cmd = 'qsub ' \
          '-V ' \
          f'-j {kargs["stream"]} ' \
          f'-q {kargs["queue"]} ' \
          f'-S {kargs["interpreter"]} ' \
          f'-o {qsub_out} ' \
          f'{kargs["cmd"]} >> {qsub_log}'
    queue_state, node = check_queue_node()
    if queue_state is True:
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'Queue 작업 등록 완료',
                'false_meg': 'run_cmd_qsub',
            }
        )
    else:
        secho('Error: Queue 사용할 수 없는 서버입니다.', fg='red', blink=True, err=True)
        echo(f'서버: {node}', err=True)
    qsub_recipe_file = os.path.join(kargs['out_path'], 'qsub.recipe')
    with open(qsub_recipe_file, 'w') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')


def read_qsub_log(p_log):
    with open(p_log, 'r') as o_log:
        # log_text : Your job 5125107 ("Job_MBS191230NE001") has been submitted
        job_id = o_log.read().strip().replace('Your job ', '').split()[0]
    return int(job_id)


def parse_job_id(p_path, p_pattern, glob_type):
    qsub_log = glob_dir(p_path, p_pattern, glob_type)
    if qsub_log is None:
        secho('Error: qsub.log 파일이 없습니다.', fg='red', blink=True, err=True)
        echo(p_path)
    if glob_type == 'many':
        l_job_id = list()
        for log in qsub_log:
            job_id = read_qsub_log(log)
            l_job_id.append(int(job_id))
        return l_job_id
    elif glob_type == 'only':
        return read_qsub_log(qsub_log)


def check_jobs_in_queue(l_job_id: list):
    from sgeparse import get_jobs
    from time import time, sleep
    from SpoON.run_time import compute_run_time
    start_time = time()
    done_job_count = 0  # init
    status_check = 1
    while done_job_count != len(l_job_id):
        s_running_job = set()
        s_wait_job = set()
        done_job_count = 0
        l_parsed_jobs = get_jobs()
        l_parsed_jobs_id = [x.get('job_number') for x in l_parsed_jobs]
        for job in l_job_id:
            if job in l_parsed_jobs_id:
                index = l_parsed_jobs_id.index(job)
                if l_parsed_jobs[index].get('state') == 'r':
                    s_running_job.add(job)
                elif l_parsed_jobs[index].get('state') == 'qw':
                    s_wait_job.add(job)
            else:
                done_job_count += 1
        wait_job_text = style(f'waiting: {len(s_wait_job)}', fg='magenta')
        done_job_text = style(f'done: {done_job_count}', fg='cyan')
        check_time = time()
        run_time = compute_run_time(start_time, check_time)
        status_check += 1
        echo(f'\rrunning: {len(s_running_job)},  {wait_job_text}, {done_job_text}'
             f' -- {run_time} 경과, 상태확인: {status_check}회', nl=False)
        sleep(1)


def read_yaml(p_file):
    """
    yaml 파일을 읽어 딕션너리 객체로 반환한다.

    :param p_file: yaml 파일
    :return: data
    """
    import yaml
    with open(p_file, 'r') as o_yaml:
        data = yaml.safe_load(o_yaml)
    return data

