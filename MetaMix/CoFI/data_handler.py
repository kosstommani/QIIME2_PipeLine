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
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.3'

from collections import Counter
import os
from click import secho, echo
from SpoON.util import glob_dir, check_file_type, read_all_data
from CoFI.data_structure import RunData, SampleList


def find_rawdata(p_path, p_target_suffix='*', p_sample_suffix=None):
    """

    ex) RunData.path: '/crystal/Analysis/BI/AmpliconMetagenomics/PipeLine_Test/1808KMI-0057/RawData/180910_BVPP3_1'
        RunData.run_info: '180910_BVPP3_1'
        RunData.samples: [SampleList, ..., ...]
            SampleList.path: '/crystal/Analysis/BI/AmpliconMetagenomics/PipeLine_Test/
                             1808KMI-0057/RawData/180910_BVPP3_1/YD.080'
            SampleLIst.name:'YD.080'

    :type p_path: str
    :param p_path: rawdata가 있는 경로

    :type p_target_suffix: str
    :param p_target_suffix:

    :type p_sample_suffix: str
    :param p_sample_suffix: 선택하고자 하는 시료의 접미사. 값이 None일 경우 검색하지 않음.

    :return: named tuple인 RunData 를 원소로 가지는 리스트
             RunData.path
                    .run_info
                    .samples: named tuple인 SampleList을 원소로 가지는 리스트
                        SampleList.path: 경로
                                  .name: 시료명
    """
    temp_dir_list = glob_dir(p_path, p_target_suffix, 'many')
    run_dir_list, not_run_dir_list = check_file_type(temp_dir_list, 'isdir')
    del temp_dir_list
    l_nt_run_list = list()
    for path in run_dir_list:
        if os.path.basename(path).lower() in ['multiqc_data']:
            secho('>>> 시료 목록 작성 제외 디렉터리', fg='magenta')
            echo(path)
            continue
        dir_list = glob_dir(path, '*', 'many')
        sample_run_dir, not_sample_run_dir = check_file_type(dir_list, 'isdir')
        if p_sample_suffix is None:
            l_nt_sample_list = [SampleList(path=path, name=os.path.split(path)[1]) for path in sample_run_dir]
            l_nt_sample_list.sort()
            l_nt_run_list.append(RunData(path=path, run_info=os.path.split(path)[1], samples=l_nt_sample_list))
        else:
            l_nt_sample_list = list()
            for run_path in sample_run_dir:
                name = os.path.split(run_path)[1]
                if p_sample_suffix == name[-(len(p_sample_suffix)):]:
                    l_nt_sample_list.append(SampleList(path=run_path, name=name))
                else:
                    continue
            l_nt_sample_list.sort()
            l_nt_run_list.append(RunData(path=path, run_info=os.path.split(path)[1], samples=l_nt_sample_list))
    secho('>>> Run Directory 확인', fg='cyan')
    _ = [echo(run.path) for run in l_nt_run_list]
    return l_nt_run_list


def check_duplication(p_run):
    l_sample = list()
    for run in p_run:
        l_sample.extend([sample.name for sample in run.samples])
    count = Counter(l_sample)
    duplicated_sample = [x for x in count if count[x] > 1]
    if duplicated_sample:
        secho('Error: 시료 중복', fg='red', blink=True, err=True)
        for run in p_run:
            name = [sample.name for sample in run.samples if sample.name in duplicated_sample]
            if name:
                echo('{run}: {name}'.format(run=run.run_info, name=name))
        # TODO: 중복 시료 처리
        exit(1)
    else:
        secho('>>> 시료 중복 확인 완료', fg='cyan')
        secho('시료 개수: {}'.format(len(l_sample)), fg='yellow')
        echo('시료 목록: {}'.format(l_sample))
        return 0


def read_metadata_group_file(p_file):
    d_group = dict()
    with open(p_file, 'r') as o_file:
        for i in o_file:
            l_text = i.strip().split('\t')
            d_group[l_text[0]] = i.strip()
    return d_group


def read_custom_order_for_metadata(p_order_file):
    """
    고객맞춤정렬을 위해 시료 순서가 정리된 파일을 읽고, 시료의 순서를 Tuple 로 반환한다.
    해당 파일은 수주번호 디렉터리에 있으며 파일명은 다음과 같다.
    수주번호_customer_sort.txt

    :param p_order_file: 수주번호_customer_sort.txt
    :rtype: tuple
    :return: 시료의 순서가 반영된 시료명을 원소로 가지는 Tuple.
    """
    if check_file_type(p_order_file, 'exists'):
        customer_sort_order = read_all_data(p_order_file)
        l_customer_sort_order = customer_sort_order.split()
        return tuple(l_customer_sort_order)
    else:
        secho('Warning: Custom Order File 이 없습니다.', fg='yellow', blink=True)
        echo(f'{p_order_file}')
        secho('\t---> 고객맞춤정렬이 적용되지 않습니다.', fg='magenta')
        return None


def make_metadata_file(p_l_nt_run, p_path, p_d_group=None, p_t_custom_order=None):
    """
    시료명으로만 구성된 Metadata File을 생성한다.

    :type  p_l_nt_run: list
    :param p_l_nt_run: named tuple인 RunData 를 원소로 가지는 리스트
                RunData.path
                       .run_info
                       .samples: named tuple인 SampleList을 원소로 가지는 리스트
                            SampleList.path: 경로
                                      .name: 시료명

    :type p_path: str
    :param p_path: 파일 생성 경로.
    :type p_d_group: dict
    :param p_d_group: 시료에 대한 그룹 정보.
                      Key: 시료명
                      value: 시료명을 포함한 그룹정보의 문자열. ex) Test1\tGroup1\tGroup2
    :type p_t_custom_order: tuple
    :param p_t_custom_order: 고객맞춤정렬을 위한 시료 순서가 반영된 tuple.
    :return: metadata_file
    """
    l_sample_for_analysis = list()
    for nt_run in p_l_nt_run:
        for nt_sample in nt_run.samples:
            l_sample_for_analysis.append(nt_sample.name)

    metadata_file = os.path.join(p_path, 'metadata.txt')
    # metadata_file = os.path.join(*(p_path, 'mappingfile.txt'))
    with open(metadata_file, 'w') as metadata:
        if p_d_group is None:  # 그룹 정보가 없는 경우
            # metadata.write('ID\n')
            metadata.write('#SampleID\n')
            if p_t_custom_order is None:
                for nt_run in p_l_nt_run:
                    for nt_sample in sorted(nt_run.samples):
                        metadata.write(nt_sample.name)
                        metadata.write('\n')
            else:  # 고객맞춤정렬
                l_no_samples = list()
                for order_ele in p_t_custom_order:
                    if order_ele in l_sample_for_analysis:
                        index = l_sample_for_analysis.index(order_ele)
                        sample_name = l_sample_for_analysis.pop(index)
                        if order_ele == sample_name:
                            metadata.write(sample_name)
                            metadata.write('\n')
                    else:
                        l_no_samples.append(order_ele)
                s_sample_for_analysis = set(l_sample_for_analysis)
                s_no_samples = set(l_no_samples)
                difference = s_sample_for_analysis - s_no_samples
                if len(s_sample_for_analysis) == len(difference):
                    pass
                else:
                    secho('Warning: Metadata.txt 파일에 분석 대상 시료가 모두 기재되지 않았습니다.', fg='yellow', blink=True)
                    echo(f'분석대상시료: {l_sample_for_analysis}')
                    echo(f'Metadata 기재 제외 시료: {l_no_samples}')
                    secho('\t---> Metadata.txt 파일을 확인하세요.', fg='magenta')
        else:  # 그룹 정보가 있는 경우
            text = p_d_group.get('#SampleID')
            if text is None:
                secho('Error: metadata_group 파일에 Header 정보가 없습니다.', fg='red', blink=True)
                echo(p_d_group)
                exit()
            else:
                metadata.write(p_d_group.get('#SampleID'))
                metadata.write('\n')
            if p_t_custom_order is None:
                for nt_run in p_l_nt_run:
                    for nt_sample in sorted(nt_run.samples):
                        text = p_d_group.get(nt_sample.name)
                        if text is None:
                            secho('Error: metadata_group 파일에 해당 시료의 정보가 없습니다.', fg='red', blink=True)
                            echo(f'미검출 시료명: {nt_sample.name}')
                            echo(p_d_group)
                            exit()
                        else:
                            metadata.write(text)
                            metadata.write('\n')
            else:  # p_t_custom_order 가 있는 경우 --> 고객맞춤정렬
                l_no_samples = list()
                for order_ele in p_t_custom_order:
                    if order_ele in l_sample_for_analysis:
                        index = l_sample_for_analysis.index(order_ele)
                        sample_name = l_sample_for_analysis.pop(index)
                        if order_ele == sample_name:
                            text = p_d_group.get(sample_name)
                            if text is None:
                                secho('Error: metadata_group 파일에 해당 시료의 정보가 없습니다.', fg='red', blink=True)
                                echo(f'미검출 시료명: {nt_sample.name}')
                                echo(p_d_group)
                                exit()
                            else:
                                metadata.write(text)
                                metadata.write('\n')
                        else:
                            secho('Error: 고객맞춤정렬이 적용된 Metadata 파일생성에 문제가 있습니다.', fg='red', blink=True)
                            secho('\t---> 개발자에게 문의하세요.', fg='magenta')
                            echo(f'order_ele  : {order_ele}')
                            echo(f'sample_name: {sample_name}')
                            echo(f'index : {index}')
                            echo(f'l_sample_for_analysis: {l_sample_for_analysis}')
                            exit()
                    else:  # 고객맞춤정렬 목록의 시료가 분석 대상에 없는 경우.
                        l_no_samples.append(order_ele)
                # Metadata.txt 파일에 분석 대상 시료들이 모두 포함되었는지 확인.
                s_sample_for_analysis = set(l_sample_for_analysis)
                s_no_samples = set(l_no_samples)
                difference = s_sample_for_analysis - s_no_samples
                if len(s_sample_for_analysis) == len(difference):
                    pass
                else:
                    secho('Warning: Metadata.txt 파일에 분석 대상 시료가 모두 기재되지 않았습니다.', fg='yellow', blink=True)
                    echo(f'분석대상시료: {l_sample_for_analysis}')
                    echo(f'Metadata 기재 제외 시료: {l_no_samples}')
                    secho('\t---> Metadata.txt 파일을 확인하세요.', fg='magenta')

    secho('>>> Metadata File 생성완료', fg='cyan')
    echo(metadata_file)
    return metadata_file


def export_artifact(p_qza, p_out_dir):
    from SpoON.util import parse_config, run_cmd, check_run_cmd
    echo('>>> Export Artifact 시작')
    CONFIG = parse_config()
    QIIME2 = CONFIG['CoFI_QIIME2']
    cmd = '{qiime2} tools export ' \
        '--input-path {qza} ' \
        '--output-path {dir} '.format(
            qiime2=QIIME2,
            qza=p_qza,
            dir=p_out_dir
        )
    run = run_cmd(cmd)
    check_run_cmd({
        'run': run,
        'true_meg': 'Export Artifact 완료',
        'false_meg': 'QIIME2 - export',
    })
