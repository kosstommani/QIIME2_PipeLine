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
__version__ = '1.0.5'

# Ver. 1.0.2
# 2020.04.17
# TAXONOMY - blast parser 변경
# ----------------------------
# Ver. 1.0.3
# 2020.04.28
# TAXONOMY DB 다중 선택 기능 추가
# ----------------------------
# Ver. 1.0.4
# R_DADA2 명령어 추가
# ----------------------------
# Ver. 1.0.5
# R_DADA2 : Summary Files 통합
# R_DADA2 : metadata.txt, info.txt 생성


import os
from time import time
from click import secho, echo, style, progressbar
from CoFI.data_handler import find_rawdata, check_duplication, make_metadata_file, export_artifact, \
    read_metadata_group_file, read_custom_order_for_metadata
from CoFI.read_assembly import run_flash, make_excel_for_flash
from CoFI.OTU import run_cd_hit_otu, make_pick_otus_rep_file, run_closed_reference_otu
from SpoON.util import glob_dir, check_file_type, make_info_file, read_all_data, make_dir_using_input_file, \
    parse_config, run_cmd, check_run_cmd, get_target_dir_number, detect_otu_method, read_data
from SpoON.run_time import compute_run_time
from SpoON.fastq_handler import filter_length_all, exist_assembled_fastq
from CoFI.ASVs import run_import_fastq, run_demux_summarize, run_dada2, run_dada2_stats_qzv, \
    run_feature_table_summarize, run_feature_table_tabulate_seqs
from CoFI.phylogeny import run_align_to_tree_mafft_fasttree, run_phylogeny_qiime1
from CoFI.alignment import run_alignment_muscle, run_filter_alignment
from CoFI.taxonomy import run_taxonomy_qiime1, run_standalone_blast2
from CoFI.biom import run_make_otu_table, run_collapse_samples
from MicrobeAndMe.data_structure import Microbe, ASVsMerger, FilesCollector
# from CoFI.diversity import run_core_metrics_phylogenetic
# from pprint import pprint

# 기본값 설정
CONFIG = parse_config()
ANALYSIS_DIR_NAME = CONFIG['Analysis_Dir_Name']


def make_sample_list(p_path, p_order_number, p_target_suffix, p_sample_suffix):
    """
    OTU분석을 진행할 시료들의 목록을 정리한다.

    :type p_path: str
    :param p_path: 경로

    :type p_order_number: str
    :param p_order_number: 수주번호

    :type p_target_suffix: str
    :param p_target_suffix: 시료 검출을 위한 대상 디렉터리의 접미사

    :type p_sample_suffix: str
    :param p_sample_suffix: 검출하고자하는 시료의 접미사. None일 경우 검색하지 않음.

    :rtype: list
    :return: l_nt_run_list - named tuple인 RunData 를 원소로 가지는 리스트
                 RunData.path
                        .run_info
                        .samples: named tuple인 SampleList을 원소로 가지는 리스트
                            SampleList.path: 경로
                                      .name: 시료명
    """
    order_path = os.path.join(p_path, p_order_number)
    l_nt_run_list = find_rawdata(order_path, p_target_suffix, p_sample_suffix)
    check_duplication(l_nt_run_list)
    return l_nt_run_list


def make_analysis_dir(p_path, p_order_number, p_dir_base_name):
    """
    OTU분석을 위한 디렉터리를 생성한다.

    :param p_path: str
    :param p_order_number: str
    :param p_dir_base_name: str
    :return:
    """
    num = 1
    order_path = os.path.join(p_path, p_order_number)
    while glob_dir(order_path, p_dir_base_name.format(num), p_verbose=False):
        num += 1
    analysis_path = os.path.join(order_path, p_dir_base_name.format(num))
    os.mkdir(analysis_path)
    secho('>>> 분석 디렉터리 생성', fg='cyan')
    echo(analysis_path)
    return analysis_path


def save_metadata_file(kargs):
    """
    metadata 파일을 생성하기 위한 데이터 처리 과정을 수행하고 metadata파일을 저장한다.
    
    :param kargs:
        analysis_base_path
        order_number
        target_dir_suffix
        l_nt_run_list
        analysis_path
    :return:
    """
    work_order_number_dir = os.path.join(kargs['analysis_base_path'], kargs['order_number'])
    l_group_file = glob_dir(work_order_number_dir, f"{kargs['target_dir_suffix']}/metadata_group.txt", p_mode='many')
    custom_sort_name = f'{kargs["order_number"]}_custom_sort.txt'
    custom_sort_file = os.path.join(work_order_number_dir, custom_sort_name)
    if l_group_file is not None:
        l_group_dict = list()
        for group_file in l_group_file:
            l_group_dict.append(read_metadata_group_file(group_file))
        d_sample_group = dict()
        for group_dict in l_group_dict:
            header_value = d_sample_group.get('#SampleID')
            if header_value is not None:
                if header_value == group_dict.get('#SampleID'):
                    del group_dict['#SampleID']
                else:
                    l_keep_header = header_value.split('\t')
                    l_new_header = group_dict.get('#SampleID').split('\t')
                    if len(l_keep_header) == len(l_new_header):
                        for keep, new in zip(l_keep_header, l_new_header):
                            if keep == new:
                                continue
                            else:
                                secho('Error: metadata_group 파일들의 Header의 이름이 다릅니다.', fg='red', blink=True)
                                echo(f'l_keep_header: {l_keep_header}')
                                echo(f'l_new_header: {l_new_header}')
                                echo(f'Kept Header: {keep}')
                                echo(f'New Header : {new}')
                                exit()
                        else:
                            del group_dict['#SampleID']
                    else:
                        secho('Error: metadata_group 파일들의 Header의 개수가 다릅니다.', fg='red', blink=True)
                        echo(f'l_keep_header: {l_keep_header}')
                        echo(f'l_new_header: {l_new_header}')
                        exit()
            d_sample_group.update(group_dict)
        else:
            t_custom_order = read_custom_order_for_metadata(custom_sort_file)
            if t_custom_order is None:
                make_metadata_file(kargs['l_nt_run_list'], kargs['analysis_path'], d_sample_group, None)
            else:
                make_metadata_file(kargs['l_nt_run_list'], kargs['analysis_path'], d_sample_group, t_custom_order)
    else:  # l_group_file 파일이 없는 경우 --> (그룹 정보 없음)
        t_custom_order = read_custom_order_for_metadata(custom_sort_file)
        if t_custom_order is None:
            make_metadata_file(kargs['l_nt_run_list'], kargs['analysis_path'])
        else:  # 고객맞춤정렬
            make_metadata_file(kargs['l_nt_run_list'], kargs['analysis_path'], None, t_custom_order)


def save_info_file(kargs):
    """
    info.txt 파일을 생성하기 위해 필요한 데이터 처리를 진행하고 info.txt 파일을 저장한다.
    :param kargs:
            target_info
            target_region
            l_nt_run_list
            analysis_path
            order_number
    :return:
    """
    if kargs['target_info'] is not None:
        if kargs['target_region'] == 'Bakt_341F-805R':
            target_region = 'V3V4(Bakt_341F-805R)'
        elif kargs['target_region'] == 'ITS_3F-4R':
            target_region = 'ITS(ITS_3F-4R)'
        elif kargs['target_region'] == '27F-Bact338':
            target_region = 'V1V2(27F-Bact338)'
        else:
            raise ValueError('입력된 값이 등록된 목록에 없습니다.\n'
                             '입력값 : {}'.format(kargs['target_region']))
    else:
        target_region = None
    sample_count = 0
    for run in kargs['l_nt_run_list']:
        sample_count += len(run.samples)
    make_info_file(kargs['order_number'], kargs['analysis_path'], sample_count, target_region)


def flash_pipeline(kargs, p_target_info=None):
    """
    - Read Assembly(FLASH)
    Assembled Read에 대한 Length Filter, FASTQ Header 변경, Pooled Sample File 생성의 작업들을 수행한다.
    
    - 기타
    metadata.txt 파일을 생성한다.
    info.txt 파일을 생성한다.

    :type kargs: dict
    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            order_number: str - 수주번호.
            analysis_base_path: str - 분석 기본 경로.
            target_dir_suffix: str - 대상 디렉터리 접미사.
            analysis_dir_base_name: str - 작업 디렉터리명.
            sample_suffix: str - 검출 대상 시료의 접미사.
            r1_suffix: str - R1 접미사.
            r2_suffix: str - R2 접미사.
            adapter_trim_tool: str - [fastp | SeqPurge | Scythe]
            read_length: int
            no_use_worker: bool
            target_size: int
            min_overlap: [int | None] - None --> 자동 계산
            max_overlap: [int | None] - None --> 자동 계산
            percent_combined_alarm: int
            no_compress: bool
            no_csv: bool
            no_excel: bool
            no_Hcsv: bool
            length_filter: bool
            length_min: [int | False] - False --> 미적용
            length_max: [int | False] - False --> 미적용
            modify_header: bool
            pooled_sample: bool
            out_seq_type: str - [FASTQ | FASTA]

    :type p_target_info: None or dict
    :param p_target_info: Primer 정보(primer.yaml 형식 참조)

    :rtype: None or str
    :return: pooled_file
    """
    # FLASH_HIT COMMAND 사용시 설정 추가
    if p_target_info is not None:
        kargs['target_size'] = p_target_info['target_size']
        kargs['min_overlap'] = None
        kargs['max_overlap'] = None
        kargs['no_excel'] = False
        kargs['no_hcsv'] = False
        kargs['no_csv'] = False
        kargs['length_filter'] = True
        kargs['length_min'], kargs['length_max'] = p_target_info['length_filter']
        kargs['stat'] = True

    # 시료 목록 작성
    l_nt_run_list = make_sample_list(
        kargs['analysis_base_path'],
        kargs['order_number'],
        kargs['target_dir_suffix'],
        kargs['sample_suffix'])

    # 작업 디렉터리 생성
    analysis_path = make_analysis_dir(
        kargs['analysis_base_path'],
        kargs['order_number'],
        kargs['analysis_dir_base_name'])

    # FLASH - Read Assembly
    flash_start = time()
    l_nt_flash_list = run_flash(
        {
            'run_list': l_nt_run_list,
            'R1_suffix': kargs['r1_suffix'].format(trim=kargs['adapter_trim_tool']),
            'R2_suffix': kargs['r2_suffix'].format(trim=kargs['adapter_trim_tool']),
            'output_path': analysis_path,
            'order_number': kargs['order_number'],
            'read_length': kargs['read_length'],
            'target_size': kargs['target_size'],
            'min_overlap': kargs['min_overlap'],
            'max_overlap': kargs['max_overlap'],
            'percent_combined_alarm': kargs['percent_combined_alarm'],
            'no_compress': kargs['no_compress'],
            'no_csv': kargs['no_csv'],
            'no_excel': kargs['no_excel'],
            'no_hcsv': kargs['no_hcsv'],
            'no_use_worker': kargs['no_use_worker'],
        }
    )
    flash_end = time()
    flash_run_time = compute_run_time(flash_start, flash_end)

    # metadata.txt 파일 생성
    save_metadata_file(
        {
            'analysis_base_path': kargs['analysis_base_path'],
            'order_number': kargs['order_number'],
            'target_dir_suffix': kargs['target_dir_suffix'],
            'l_nt_run_list': l_nt_run_list,
            'analysis_path': analysis_path
        }
    )
    # work_order_number_dir = os.path.join(kargs['analysis_base_path'], kargs['order_number'])
    # l_group_file = glob_dir(work_order_number_dir, f"{kargs['target_dir_suffix']}/metadata_group.txt", p_mode='many')
    # custom_sort_name = f'{kargs["order_number"]}_custom_sort.txt'
    # custom_sort_file = os.path.join(work_order_number_dir, custom_sort_name)
    # if l_group_file is not None:
    #     l_group_dict = list()
    #     for group_file in l_group_file:
    #         l_group_dict.append(read_metadata_group_file(group_file))
    #     d_sample_group = dict()
    #     for group_dict in l_group_dict:
    #         header_value = d_sample_group.get('#SampleID')
    #         if header_value is not None:
    #             if header_value == group_dict.get('#SampleID'):
    #                 del group_dict['#SampleID']
    #             else:
    #                 l_keep_header = header_value.split('\t')
    #                 l_new_header = group_dict.get('#SampleID').split('\t')
    #                 if len(l_keep_header) == len(l_new_header):
    #                     for keep, new in zip(l_keep_header, l_new_header):
    #                         if keep == new:
    #                             continue
    #                         else:
    #                             secho('Error: metadata_group 파일들의 Header의 이름이 다릅니다.', fg='red', blink=True)
    #                             echo(f'l_keep_header: {l_keep_header}')
    #                             echo(f'l_new_header: {l_new_header}')
    #                             echo(f'Kept Header: {keep}')
    #                             echo(f'New Header : {new}')
    #                             exit()
    #                     else:
    #                         del group_dict['#SampleID']
    #                 else:
    #                     secho('Error: metadata_group 파일들의 Header의 개수가 다릅니다.', fg='red', blink=True)
    #                     echo(f'l_keep_header: {l_keep_header}')
    #                     echo(f'l_new_header: {l_new_header}')
    #                     exit()
    #         d_sample_group.update(group_dict)
    #     else:
    #         t_custom_order = read_custom_order_for_metadata(custom_sort_file)
    #         if t_custom_order is None:
    #             make_metadata_file(l_nt_run_list, analysis_path, d_sample_group, None)
    #         else:
    #             make_metadata_file(l_nt_run_list, analysis_path, d_sample_group, t_custom_order)
    # else:  # l_group_file 파일이 없는 경우 --> (그룹 정보 없음)
    #     t_custom_order = read_custom_order_for_metadata(custom_sort_file)
    #     if t_custom_order is None:
    #         make_metadata_file(l_nt_run_list, analysis_path)
    #     else:  # 고객맞춤정렬
    #         make_metadata_file(l_nt_run_list, analysis_path, None, t_custom_order)

    # info.txt 파일 생성
    save_info_file(
        {
            'target_info': p_target_info,
            'target_region': kargs['target_region'],
            'l_nt_run_list': l_nt_run_list,
            'analysis_path': analysis_path,
            'order_number': kargs['order_number'],
        }
    )
    # if p_target_info is not None:
    #     if kargs['target_region'] == 'Bakt_341F-805R':
    #         target_region = 'V3V4(Bakt_341F-805R)'
    #     elif kargs['target_region'] == 'ITS_3F-4R':
    #         target_region = 'ITS(ITS_3F-4R)'
    #     elif kargs['target_region'] == '27F-Bact338':
    #         target_region = 'V1V2(27F-Bact338)'
    #     else:
    #         raise ValueError('입력된 값이 등록된 목록에 없습니다.\n'
    #                          '입력값 : {}'.format(kargs['target_region']))
    # else:
    #     target_region = None
    # sample_count = 0
    # for run in l_nt_run_list:
    #     sample_count += len(run.samples)
    # make_info_file(kargs['order_number'], analysis_path, sample_count, target_region)

    # Length Filtering & Filtered fastq(Assembled) STAT 파일 생성 & Pooled Sample File 생성
    filter_start = time()
    if kargs['length_filter']:
        l_flash_outputs = [os.path.join(x.path, x.output) for x in l_nt_flash_list]
        pooled_file = filter_length_all(
            l_flash_outputs, kargs['length_min'], kargs['length_max'],
            l_nt_flash_list[0].path, kargs['order_number'], p_modify_header=kargs['modify_header'],
            p_stat=kargs['stat'], p_pooled=kargs['pooled_sample'], p_out_seq_type=kargs['out_seq_type'])
        end_time = time()
        filter_run_time = compute_run_time(filter_start, end_time)
        return pooled_file, {'flash_run_time': flash_run_time, 'filter_run_time': filter_run_time}
    elif kargs['length_filter'] is False:
        l_flash_outputs = [os.path.join(x.path, x.output) for x in l_nt_flash_list]
        pooled_file = filter_length_all(
            l_flash_outputs, 'no', 'no',
            l_nt_flash_list[0].path, kargs['order_number'], p_modify_header=kargs['modify_header'],
            p_stat=kargs['stat'], p_pooled=kargs['pooled_sample'], p_out_seq_type=kargs['out_seq_type'])
        end_time = time()
        filter_run_time = compute_run_time(filter_start, end_time)
        return pooled_file, {'flash_run_time': flash_run_time, 'filter_run_time': filter_run_time}
    end_time = time()
    filter_run_time = compute_run_time(filter_start, end_time)
    return None, {'flash_run_time': flash_run_time, 'filter_run_time': filter_run_time}


def length_filter_pipeline(kargs):
    if (kargs['analysis_number'] != 'None') and \
            ((kargs['assembled_dir'] == 'None') or (kargs['assembled_dir'] is None)):
        assembly_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'],
                                         kargs['analysis_number'], ANALYSIS_DIR_NAME['read_assembly'])
    elif (kargs['assembled_dir'] != 'None') or (kargs['assembled_dir'] is not None):
        if os.path.isabs(kargs['assembled_dir']):
            assembly_dir_path = kargs['assembled_dir']
        else:
            assembly_dir_path = os.path.abspath(kargs['assembled_dir'])
    else:
        raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')

    if (kargs['mode'] == 'excel') or (kargs['mode'] == 'all'):
        excel_start = time()
        l_glob_log = glob_dir(assembly_dir_path, kargs['flash_log_suffix'], p_mode='many')
        secho(f'>>> FLASH Log 파일 확인(개수: {len(l_glob_log)}).', fg='cyan')

        # 기존 엑셀 파일 삭제
        l_excel = [f'{kargs["order_number"]}.FLASH.Hcsv', f'{kargs["order_number"]}.FLASH.csv',
                   f'{kargs["order_number"]}.FLASH.xlsx', f'{kargs["order_number"]}.hist.csv']
        for excel in l_excel:
            target = os.path.join(assembly_dir_path, excel)
            cmd = f'rm {target}'
            run = run_cmd(cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': f'{excel} 삭제 완료',
                    'false_meg': None,
                }, p_exit=False, p_stdout=False,
            )

        l_log_texts = list()
        with progressbar(l_glob_log, label='FLASH Log', fill_char=style('#', fg='green')) as bar_l_glob_log:
            for log in bar_l_glob_log:
                l_log_texts.append(read_all_data(log, p_verbose=False))
        make_excel_for_flash(
            {
                'log_texts': l_log_texts,
                'save_path': assembly_dir_path,
                'output_name': kargs['order_number'],
                'R1_suffix': kargs['r1_suffix'].format(trim=kargs['adapter_trim_tool']),
                'R2_suffix': kargs['r2_suffix'].format(trim=kargs['adapter_trim_tool']),
                'target_size': kargs['target_size'],
                'percent_combined_alarm': kargs['percent_combined_alarm']
            },
            p_csv=True, p_excel=True, p_Hcsv=True, p_no_compress=False if kargs['compression'] else True
        )
        del l_log_texts, l_glob_log, log
        excel_end = time()
        excel_time = compute_run_time(excel_start, excel_end)

    # Length Filtering & Filtered fastq(Assembled) STAT 파일 생성 & Pooled Sample File 생성
    if (kargs['mode'] == 'length_filter') or (kargs['mode'] == 'all'):
        filter_start = time()
        # 기존 파일 삭제
        l_will_del = [f'{kargs["order_number"]}.pooled.fastq', f'{kargs["order_number"]}.F_length.csv',
                      'F_*extendedFrags.fastq.gz', '*extendedFrags.fastq.gz.log',
                      'STAT.bash', 'STAT.txt', 'soft_link.bash']
        for file in l_will_del:
            target = os.path.join(assembly_dir_path, file)
            cmd = f'rm {target}'
            run = run_cmd(cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': f'{file} 삭제 완료',
                    'false_meg': None,
                }, p_exit=False, p_stdout=False,
            )
        cmd = f'find {assembly_dir_path} -type l | xargs rm'  # _1.fastq.gz 링크 파일 삭제
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'FASTQ 링크 파일 삭제 완료',
                'false_meg': None,
            }, p_exit=False, p_stdout=False,
        )

        l_glob_file = glob_dir(assembly_dir_path, '*.extendedFrags.fastq*[!log]', p_mode='many')
        secho(f'>>> Assembled FASTQ 파일 확인(개수: {len(l_glob_file)}).', fg='cyan')
        pooled_file = filter_length_all(
            l_glob_file, kargs['length_min'], kargs['length_max'], assembly_dir_path,
            kargs['order_number'], p_modify_header=kargs['modify_header'], p_stat=kargs['stat'],
            p_pooled=kargs['pooled_sample'], p_out_seq_type=kargs['out_seq_type'])
        filter_end = time()
        filter_run_time = compute_run_time(filter_start, filter_end)
    if kargs['mode'] == 'excel':
        return None, {'flash_excel_run_time': excel_time, 'filter_run_time': None}
    elif kargs['mode'] == 'length_filter':
        return pooled_file, {'flash_excel_run_time': None, 'filter_run_time': filter_run_time}
    elif kargs['mode'] == 'all':
        return None, {'flash_excel_run_time': excel_time, 'filter_run_time': filter_run_time}
    else:
        raise ValueError(f'mode : {kargs["mode"]}')


def cd_hit_otu_pipeline(kargs):
    """
    CD-HIT-OTU 프로그램을 실행한다.
    QIIME 프로그램에서 BIOM 파일을 생성하기 위해 필요한 pick_otus.txt 파일을 생성한다.
    OTU 대표 서열에 대한 fasta 파일(otus_rep.fasta)을 생성한다.
    STAT.txt 파일을 복사한다.

    :type kargs: dict
    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            pooled_sample: str
            out_dir: None or str
            otu_cutoff: float
            prefix_length_or_primer: int or str
            trim_cutoff: int or float
            small_cluster: str
            check_chimera: bool
            pcr_error: float
            analysis_number: str
    :return:
    """
    # CD-HIT-OTU
    start_time = time()
    out_dir = run_cd_hit_otu(kargs)
    cd_hit_end_time = time()
    cd_hit_run_time = compute_run_time(start_time, cd_hit_end_time)
    # pick_otus.txt & outs_rep.fasta
    make_pick_otus_rep_file(out_dir)
    otus_file_end_time = time()
    otus_file_run_time = compute_run_time(cd_hit_end_time, otus_file_end_time)
    return {'cd-hit_run_time': cd_hit_run_time, 'otus_file_run_time': otus_file_run_time}


def flash_hit_pipeline(kargs, p_dic_target_info):
    """
    Bakt(341F-805R)와 ITS(3F-4R) 영역, 27F-Bact338(유산균 V1V2)에 대해서 기본 설정으로 아래의 분석을 진행한다.
    Read Assembly(FLASH) --> Length Filter & Header 변경 --> Pooled Sample File --> Clustering(CD-HIT-OTU)

    :type kargs: dict
    :param kargs: flash_pipeline 참조

    :type p_dic_target_info: dict
    :param p_dic_target_info: primer.yaml 에 등록된 Primer 정보.
    :return:
    """
    secho('>>> FLASH-HIT 시작: {}'.format(style(kargs['target_region'], fg='yellow')))
    kargs['modify_header'] = True
    kargs['pooled_sample'] = True
    kargs['out_seq_type'] = 'FASTQ'
    kargs['pooled_sample'], flash_run_time = flash_pipeline(kargs, p_dic_target_info)
    kargs['analysis_number'] = 'None'
    kargs['out_dir'] = None
    if kargs['target_region'] == 'Bakt_341F-805R':
        kargs['prefix_length_or_primer'] = p_dic_target_info['Bakt_341F']
    elif kargs['target_region'] == 'ITS_3F-4R':
        kargs['prefix_length_or_primer'] = p_dic_target_info['3F']
    elif kargs['target_region'] == '27F-Bact338':
        kargs['prefix_length_or_primer'] = p_dic_target_info['27F']
    else:
        raise ValueError(style('예기치 않은 오류: {}'.format(kargs['target_region']), fg='red'))
    cd_hit_run_time = cd_hit_otu_pipeline(kargs)
    return flash_run_time, cd_hit_run_time


def closed_otu_pipeline(kargs):
    global ANALYSIS_DIR_NAME
    start_time = time()
    if kargs['mode'] == 'pooled_sample':
        if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
            assembly_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'],
                                             kargs['analysis_number'], ANALYSIS_DIR_NAME['read_assembly'])
            pooled_fasta_name = f'{kargs["order_number"]}.pooled.fasta'
            pooled_fasta_file = os.path.join(assembly_dir_path, pooled_fasta_name)
        elif kargs['pooled_sample'] is not None:
            if os.path.isabs(kargs['pooled_sample']):
                pooled_fasta_file = kargs['pooled_sample']
            else:
                pooled_fasta_file = os.path.abspath(kargs['pooled_sample'])
        _ = make_dir_using_input_file(pooled_fasta_file, ANALYSIS_DIR_NAME['closed'], 2, p_check=False)
        run_closed_reference_otu(
            {
                'pooled_sample': pooled_fasta_file,
                'no_parallel': kargs['no_parallel'],
                'no_assign': kargs['no_assign'],
                'jobs': kargs['jobs'],
                'database': kargs['database'],
                'similarity': kargs['similarity'],
            }, p_mode='pool')
    elif kargs['mode'] == 'by_sample':
        order_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'])
        if check_file_type(order_path, 'exists'):
            l_nt_run_list = make_sample_list(
                kargs['analysis_base_path'],
                kargs['order_number'],
                kargs['target_dir_suffix'],
                None,
            )
            l_samples = list()
            for run in l_nt_run_list:
                for sample in run.samples:
                    l_samples.append(sample.name)
            del l_nt_run_list

            analysis_path = os.path.join(order_path, kargs['analysis_number'])
            assembly_path = glob_dir(analysis_path, ANALYSIS_DIR_NAME['read_assembly'], 'only')
            exist_assembled_fastq(assembly_path,
                                  ['F_{}'.format(sample) for sample in l_samples],
                                  kargs['assembly_suffix'])
        else:
            secho('해당 수주번호의 디렉터리가 존재하지 않습니다.', fg='red')
        closed_dir_path = make_dir_using_input_file(analysis_path, ANALYSIS_DIR_NAME['closed'], 0, p_check=False)
        run_closed_reference_otu(
            {
                'sample_list': l_samples,
                'closed_path': closed_dir_path,
                'assembly_path': assembly_path,
                'assembly_suffix': kargs['assembly_suffix'],
                'no_parallel': kargs['no_parallel'],
                'no_assign': kargs['no_assign'],
                'jobs': kargs['jobs'],
                'database': kargs['database'],
                'similarity': kargs['similarity'],
            }, p_mode='sample')
    end_time = time()
    closed_otu_run_time = compute_run_time(start_time, end_time)
    return closed_otu_run_time


def dada2_pipeline(kargs):
    """
    QIIME2 프로그램내의 DADA2 명령어를 실행하기 위한 전반적인 작업을 진행한다.

    DADA2 PipeLine 작업 진행
    시료 목록 구성 --> 작업 디렉터리 생성 --> (Metadata File 생성) --> FASTQ Import(qza) --> FASTQ Summary(qzv) -->
    DADA2 --> DATA2 STATs(qav) --> Feature Table summary(qzv) --> Feature Table Sequences(qzv)

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                order_number: str
                metadata: str
                analysis_base_path: str
                target_dir_suffix: str
                analysis_dir_base_name: str
                r1_suffix: str
                r2_suffix: str
                adapter_trim_tool: str
                random_read: int
                p_trunc_len_f: int
                p_trunc_len_r: int
                p_trim_left_f: int
                p_trim_left_r: int

    :return:
    """
    # 시료 목록 작성
    l_nt_run_list = make_sample_list(
        kargs['analysis_base_path'],
        kargs['order_number'],
        kargs['target_dir_suffix'],
        kargs['sample_suffix'])

    # 작업 디렉터리 생성
    analysis_path = make_analysis_dir(
        kargs['analysis_base_path'],
        kargs['order_number'],
        kargs['analysis_dir_base_name'])

    # MetaData File - 기본 형식 생성
    if kargs['metadata'].lower() == 'no':
        kargs['metadata'] = make_metadata_file(l_nt_run_list, analysis_path)

    # FASTQ import --> qza
    fastq_artifact = run_import_fastq(
        {
            'run_list': l_nt_run_list,
            'R1_suffix': kargs['r1_suffix'].format(trim=kargs['adapter_trim_tool']),
            'R2_suffix': kargs['r2_suffix'].format(trim=kargs['adapter_trim_tool']),
            'output_path': analysis_path,
            'order_number': kargs['order_number']
        }
    )

    # FASTQ Summarize --> qzv
    summary_qzv = os.path.join(os.path.dirname(fastq_artifact), kargs['order_number'] + '.Summary.qzv')
    run_demux_summarize(fastq_artifact, summary_qzv, kargs['random_read'])

    # DADA2 실행
    run_dada2(
        {
            'qza': fastq_artifact,
            'output_path': analysis_path,
            'p_trunc_len_f': kargs['p_trunc_len_f'],
            'p_trunc_len_r': kargs['p_trunc_len_r'],
            'p_trim_left_f': kargs['p_trim_left_f'],
            'p_trim_left_r': kargs['p_trim_left_r'],
        }
    )

    # DADA2 STATs --> qzv
    dada2_path_dir = os.path.join(analysis_path, 'DADA2')
    run_dada2_stats_qzv(
        {
            'qza': os.path.join(dada2_path_dir, 'denoising_stats.qza'),
            'qzv': os.path.join(dada2_path_dir, 'denoising_stats.qzv'),
        }
    )

    # feature table summarize --> qzv
    run_feature_table_summarize(
        {
            'qza': os.path.join(dada2_path_dir, 'table.qza'),
            'qzv': os.path.join(dada2_path_dir, 'table.qzv'),
            'metadata': kargs['metadata'],
        }
    )

    # feature table seqs --> qzv
    run_feature_table_tabulate_seqs(
        {
            'qza': os.path.join(dada2_path_dir, 'representative_sequences.qza'),
            'qzv': os.path.join(dada2_path_dir, 'representative_sequences.qzv'),
        }
    )


def r_dada2_pipeline(kargs):
    kargs['cofi_target_dir_suffix'] = kargs['target_dir_suffix']
    r_dada2 = Microbe(kargs, mode='r_dada2', analysis_type='NGS')
    r_dada2.make_analysis_number_dir()
    r_dada2.set_path()
    r_dada2.make_sample_list()
    r_dada2.run_dada2_using_r()
    r_dada2.check_jobs()
    # Make Order Number File: ex) HN00124490   Analysis_1
    order_number_name = f'{r_dada2.order_number}_ASVs.txt'
    order_number_file = os.path.join(r_dada2.i_path.analysis_number_path, order_number_name)
    with open(order_number_file, 'w') as o_file:
        o_file.write(f'{r_dada2.order_number}\t{r_dada2.analysis_number}\n')

    # metadata.txt 파일 생성
    save_metadata_file(
        {
            'analysis_base_path': kargs['analysis_base_path'],
            'order_number': r_dada2.order_number,
            'target_dir_suffix': kargs['target_dir_suffix'],
            'l_nt_run_list': r_dada2.l_nt_run_list,
            'analysis_path': r_dada2.i_path.analysis_number_path,
        }
    )

    # info.txt 파일 생성
    save_info_file(
        {
            'target_info': None,
            'target_region': None,
            'l_nt_run_list': r_dada2.l_nt_run_list,
            'analysis_path': r_dada2.i_path.analysis_number_path,
            'order_number': r_dada2.order_number,
        }
    )

    # Merge ASVs
    kargs['order_number_file'] = order_number_file
    kargs['out_dir'] = r_dada2.i_path.analysis_number_path
    asv_merge = ASVsMerger(kargs)
    asv_merge.read_order_number_file()
    asv_merge.glob_samples()
    asv_merge.check_rds_files()
    all_rds_dir_path = asv_merge.copy_rds_files('rds_files')
    asv_merge.run_merge_asvs(all_rds_dir_path)

    # Merge DADA2 Summary Files
    kargs['file_name'] = 'dada2.summary'
    summary_collector = FilesCollector(kargs)
    summary_collector.read_order_number_file()
    summary_collector.glob_samples()
    summary_collector.check_files()
    summary_dir_path = summary_collector.copy_files('R_DADA2_Summary')
    # Save all summary files
    all_summary_file = os.path.join(summary_dir_path, 'all_dada2.summary')
    with open(all_summary_file, 'w') as o_all:
        first_sample = True
        for nt_run in summary_collector.l_nt_analysis_list:
            for sample in nt_run.samples:
                sample_summary_file = os.path.join(summary_dir_path, f'{sample.name}_dada2.summary')
                summary_data = read_data(sample_summary_file)
                if first_sample is True:
                    o_all.write('\t'.join(summary_data[0]))
                    o_all.write('\n')
                    o_all.write('\t'.join(summary_data[1]))
                    o_all.write('\n')
                    first_sample = False
                else:
                    o_all.write('\t'.join(summary_data[1]))
                    o_all.write('\n')
    secho('>>> all_data2.summary 파일 생성 완료', fg='cyan')


def maff_fasttree_pipeline(kargs):
    """
    This pipeline will start by creating a sequence alignment using MAFFT, after which any alignment columns
    that are phylogenetically uninformative or ambiguously aligned will be removed (masked).
    The resulting masked alignment will be used to infer a phylogenetic tree and then subsequently rooted at
    its midpoint. Output files from each step of the pipeline will be saved. This includes both the unmasked
    and masked MAFFT alignment from q2-alignment methods, and both the rooted and unrooted phylogenies from
    q2-phylogeny methods.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                i_sequences: str - FeatureData[Sequence]
                             The sequences to be used for creating a fasttree based rooted phylogenetic tree.
                out_dir: None or str
    :return:
    """
    out_dir_name = 'Alignment_Phylogeny'
    if kargs['out_dir'] is None:
        path = os.path.dirname(os.path.dirname(kargs['i_sequences']))  # 상위 디렉터리
        if path == '':
            kargs['out_dir'] = out_dir_name
        else:
            kargs['out_dir'] = os.path.join(path, out_dir_name)
    run_align_to_tree_mafft_fasttree(kargs)


def diversity_pipeline(kargs):
    out_dir_name = 'Diversity'
    if kargs['out_dir'] is None:
        path1 = os.path.dirname(os.path.dirname(kargs['i_phylogeny']))
        path2 = os.path.dirname(os.path.dirname(kargs['i_table']))
        if (path1 == '') and (path2 == ''):
            kargs['out_dir'] = out_dir_name
        else:
            kargs['out_dir'] = os.path.join(path2, out_dir_name)
    if kargs['sample_depth'] is None:
        from tempfile import mkdtemp
        path = os.path.dirname(kargs['out_dir'])
        temp_dir = mkdtemp('_Div', 'biom_', path)
        export_artifact(kargs['i_table'], temp_dir)

    # run_core_metrics_phylogenetic(kargs)


def alignment_muscle_pipeline(kargs):
    global ANALYSIS_DIR_NAME
    out_dir_name = ANALYSIS_DIR_NAME['alignment']
    
    # order_number & analysis_number 지정 방식
    if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
        analysis_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'], kargs['analysis_number'])
        if kargs['otu_method'] == 'auto':
            otu_method = detect_otu_method(analysis_dir_path)
        elif kargs['otu_method'] == 'denovo':
            otu_method = 'denovo'
        elif kargs['otu_method'] == 'closed':
            otu_method = 'closed'
        elif kargs['otu_method'] == 'r_dada2':
            otu_method = 'r_dada2'
        else:
            secho(f'Error: --otu_method의 값이 잘못되었습니다.', fg='read', blink=True)
            echo(f'kargs["otu_method"]: {kargs["otu_method"]}')
            exit()

        # CD-HIT-OTU
        if otu_method == 'denovo':
            cd_hit_otu_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['cd_hit_otu'])
            clustering_dir_path = glob_dir(cd_hit_otu_dir_path, f'{kargs["order_number"]}*[!.log]',
                                           p_mode='only', p_verbose=False)
            otus_rep_file = os.path.join(clustering_dir_path, 'otus_rep.fasta')
        # CLOSED
        elif otu_method == 'closed':
            closed_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['closed'])
            rep_set_dir_path = os.path.join(closed_dir_path, 'rep_set')
            otus_rep_file = os.path.join(rep_set_dir_path, 'otus_rep.fasta')
        # ASVs
        elif otu_method == 'r_dada2':
            r_dada2_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['r_dada2'])
            otus_rep_file = os.path.join(r_dada2_dir_path, 'all_ASVs.fasta')
        align_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, p_check=False)
    # 파일 지정 방식
    elif (kargs['i_otus_rep'] != 'None') or (kargs['i_otus_rep'] is not None):
        if os.path.isabs(kargs['i_otus_rep']):
            otus_rep_file = kargs['i_otus_rep']
        else:
            otus_rep_file = os.path.abspath(kargs['i_otus_rep'])
        if kargs['out_dir'] is None:
            # 분석디렉터리경로/수주번호/Analysis_?
            align_path = make_dir_using_input_file(kargs['i_otus_rep'], out_dir_name, 3)
        else:
            align_path = kargs['out_dir']
    else:
        raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')
    muscle_path = run_alignment_muscle(otus_rep_file, align_path, kargs['queue'], kargs['no_queue'])
    aligned_file = '{0}_aligned{1}'.format(*(os.path.splitext(os.path.basename(otus_rep_file))))
    run_filter_alignment(os.path.join(muscle_path, aligned_file), align_path)


def phylogeny_qiime1_pipeline(kargs):
    global ANALYSIS_DIR_NAME
    out_dir_name = ANALYSIS_DIR_NAME['phylogeny']
    if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
        analysis_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'], kargs['analysis_number'])
        if kargs['otu_method'] == 'auto':
            otu_method = detect_otu_method(analysis_dir_path)
        elif kargs['otu_method'] == 'denovo':
            otu_method = 'denovo'
        elif kargs['otu_method'] == 'closed':
            otu_method = 'closed'
        elif kargs['otu_method'] == 'r_dada2':
            otu_method = 'r_dada2'
        else:
            secho(f'Error: --otu_method의 값이 잘못되었습니다.', fg='read', blink=True)
            echo(f'kargs["otu_method"]: {kargs["otu_method"]}')
            exit()
        if otu_method == 'r_dada2':
            aligned_fasta = 'all_ASVs_aligned_pfiltered.fasta'
        else:
            aligned_fasta = 'otus_rep_aligned_pfiltered.fasta'
        filtered_aligned_fasta_file = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['alignment'],
                                                   'filtered_alignment', aligned_fasta)
        phylo_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, p_check=False)
    elif (kargs['i_filtered'] != 'None') or (kargs['i_filtered'] is not None):
        filtered_aligned_fasta_file = kargs['i_filtered']
        if kargs['out_dir'] is None:
            # 분석디렉터리경로/수주번호/Analysis_?
            phylo_path = make_dir_using_input_file(kargs['i_filtered'], out_dir_name, 3)
        else:
            phylo_path = kargs['out_dir']
    else:
        raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')
    run_phylogeny_qiime1(filtered_aligned_fasta_file, phylo_path)


def taxonomy_pipeline(kargs):
    if (not kargs['blast_db']) and (not kargs['uclust_db']):
        secho('Error: --blast_db 또는 --rdp_db 옵션의 설정이 필요합니다.', fg='red', err=True)
        exit()

    global ANALYSIS_DIR_NAME
    out_dir_name = ANALYSIS_DIR_NAME['taxonomy_assignment']
    if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
        analysis_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'], kargs['analysis_number'])
        if kargs['otu_method'] == 'auto':
            otu_method = detect_otu_method(analysis_dir_path)
        elif kargs['otu_method'] == 'denovo':
            otu_method = 'denovo'
        elif kargs['otu_method'] == 'closed':
            otu_method = 'closed'
        elif kargs['otu_method'] == 'r_dada2':
            otu_method = 'r_dada2'
        else:
            secho(f'Error: --otu_method의 값이 잘못되었습니다.', fg='read', blink=True)
            echo(f'kargs["otu_method"]: {kargs["otu_method"]}')
            exit()
        if otu_method == 'denovo':
            otu_dir = ANALYSIS_DIR_NAME['cd_hit_otu']
            sub_otu_dir = f'{kargs["order_number"]}*[!.log]'
        elif otu_method == 'closed':
            otu_dir = ANALYSIS_DIR_NAME['closed']
            sub_otu_dir = 'rep_set'
        elif otu_method == 'r_dada2':
            otu_dir = ANALYSIS_DIR_NAME['r_dada2']
            sub_otu_dir = None
        else:
            raise ValueError(f'otu_method: {otu_method} 알맞지 않은 값입니다. ')
        otu_dir_path = os.path.join(analysis_dir_path, otu_dir)
        if otu_method == 'r_dada2':
            otus_rep_file = os.path.join(otu_dir_path, 'all_ASVs.fasta')
        else:
            clustering_dir_path = glob_dir(otu_dir_path, sub_otu_dir, p_mode='only', p_verbose=False)
            otus_rep_file = os.path.join(clustering_dir_path, 'otus_rep.fasta')
        taxa_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, p_check=False)
    elif (kargs['i_otus_rep'] != 'None') or (kargs['i_otus_rep'] is not None):
        otus_rep_file = kargs['i_otus_rep']
        if kargs['out_dir'] is None:
            # 분석디렉터리경로/수주번호/Analysis_?
            taxa_path = make_dir_using_input_file(kargs['i_otus_rep'], out_dir_name, 3)
        else:
            taxa_path = kargs['out_dir']
    else:
        raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')

    if kargs['blast_db']:
        db_tool = 'blast'
        blast_method = 'blastn'
        for db_name in kargs['blast_db']:
            # db_name = kargs['blast_db']
            secho('»-(¯`·.·´¯)->'.center(80), fg='red', bold=True, blink=True)
            secho(f'++++++ ======= {db_name} ======= +++++'.center(80), fg='yellow')
            if len(kargs['remove_taxon']) != 0:
                remove_taxon = ' '.join([f'-rt {x}' for x in kargs['remove_taxon']])
            else:
                remove_taxon = None
            if len(kargs['keep_taxon']) != 0:
                keep_taxon = ' '.join([f'-kt {x}' for x in kargs['keep_taxon']])
            else:
                keep_taxon = None
            if (kargs['mode'] == 'blast') or (kargs['mode'] == 'all'):
                blast_work_path = make_dir_using_input_file(taxa_path, f'{db_tool}_{db_name}', 0, p_check=True)
            else:
                blast_work_path = os.path.join(taxa_path, f'{db_tool}_{db_name}')
            run_standalone_blast2(
                {
                    'fasta': otus_rep_file,
                    'tool': blast_method,
                    'db': db_name,
                    'out_path': blast_work_path,
                    'remove_taxon': remove_taxon,
                    'keep_taxon': keep_taxon,
                    'query_coverage': kargs['query_coverage'],
                    'identity_percentage': kargs['identity_percentage'],
                    'read_per_job': kargs['read_per_job'],
                    'nt_max_job': kargs['nt_max_job'],
                    'queue': kargs['queue'],
                    'mode': kargs['mode'],
                },
            )
            # secho('Error: blast 는 현재 지원되지 않습니다.', fg='red', err=True)
    elif kargs['uclust_db']:
        db_tool = 'uclust'
        for db_name in kargs['uclust_db']:
            # db_name = kargs['uclust_db']
            secho('»-(¯`·.·´¯)->'.center(80), fg='red', bold=True, blink=True)
            secho(f'++++++ ======= {db_name} ======= +++++'.center(80), fg='yellow')
            run_taxonomy_qiime1(otus_rep_file, db_tool, db_name, taxa_path)


def biom_pipeline(kargs):
    global ANALYSIS_DIR_NAME
    analysis_dir_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'], kargs['analysis_number'])

    if kargs['collapse'] and kargs['remove_sample']:
        secho('Error: --collapse 와 --remove_sample 같이 사용할 수 없습니다.', fg='red', err=True)
        secho('\t--> 현재 지원되지 않는 기능입니다. 개발자와 상의하세요.', fg='magenta', err=True)
        exit()
    elif kargs['collapse']:
        out_dir_name_base = f'{ANALYSIS_DIR_NAME["biom"]}_{kargs["collapse"].upper()}'
        out_dir_name = get_target_dir_number(analysis_dir_path, out_dir_name_base)
        out_dir_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, False)
        # metadata.txt 확인
        from theCups.report import Report
        i_report = Report(kargs, 'BIOM_Collapsed_')
        i_report.make_analysis_data_table({'set_style': 'box', 'column_headers': ['DIR', 'FILE']})
        i_report.check_metadata_file()
        metadata_status, metadata_file = i_report.confirm_metadata(f'BIOM - Collapsed_Sample: {kargs["collapse"]}')
        if metadata_status is False:
            exit()
        # biom 파일 확인
        i_report.check_taxonomy_assignment_dir()
        i_report.check_biom_dir()
        # Taxonomy Assignment 결과 확인
        d_analysis_dir_name = i_report.analysis_paths.analysis_dir_name
        l_taxonomy_data = list()
        for dir_name in i_report.analysis_data.taxonomy_assignment.keys():
            if dir_name == d_analysis_dir_name['taxonomy_assignment']:
                continue
            else:
                d_status = i_report.analysis_data.taxonomy_assignment.get(dir_name)
                if d_status.get('dir') is True:
                    if all(d_status.get('file').values()) is True:
                        l_taxonomy_data.append(dir_name)
                    else:
                        secho(f'--- {dir_name} 제외', fg='yellow')
                        for file_name in d_status.get('file').keys():
                            if d_status.get('file').get(file_name) is False:
                                secho(f'\t{file_name} 가 없습니다.', fg='red', blink=True)
                else:
                    secho(f'--- {dir_name} 제외', fg='yellow')
        try:
            del dir_name, d_status, file_name
        except UnboundLocalError:
            pass

        echo('+++ BASE BIOM')
        biom_name = f'otu_table.biom'
        collapsed_biom_base = f'otu_table.collapsed'
        collapsed_metadata_name = f'metadata.collapsed.txt'
        biom_file = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['biom'], biom_name)
        collapsed_metadata_file = os.path.join(analysis_dir_path, collapsed_metadata_name)
        run_collapse_samples({
            'biom': biom_file,
            'metadata': metadata_file,
            'field': kargs['field'],
            'mode': kargs['collapse'],
            'out_path': out_dir_path,
            'out_biom_base': collapsed_biom_base,
            'out_metadata': collapsed_metadata_file,
        })
        del biom_name, collapsed_biom_base, biom_file
        for db_name in l_taxonomy_data:
            biom_name = f'otu_table.{db_name}.biom'
            if i_report.analysis_data.biom['BIOM']['file'][biom_name] is True:
                echo(f'+++ {db_name}')
                collapsed_biom_base = f'otu_table.{db_name}.collapsed'
                biom_file = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['biom'], biom_name)
                # collapsed_metadata_file = os.path.join(analysis_dir_path, collapsed_metadata_name)
                run_collapse_samples({
                    'biom': biom_file,
                    'metadata': metadata_file,
                    'field': kargs['field'],
                    'mode': kargs['collapse'],
                    'out_path': out_dir_path,
                    'out_biom_base': collapsed_biom_base,
                    'out_metadata': collapsed_metadata_file,
                })

    elif kargs['remove_sample']:
        secho('지원하지 않습니다.', fg='red')
        exit()
        out_dir_name_base = f'{ANALYSIS_DIR_NAME["biom"]}_RM'
        out_dir_name = get_target_dir_number(out_dir_name_base)
        biom_dir_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, False)
    else:  # BIOM
        out_dir_name_base = None
        out_dir_name = ANALYSIS_DIR_NAME['biom']
        if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
            if kargs['otu_method'] == 'auto':
                otu_method = detect_otu_method(analysis_dir_path)
            elif kargs['otu_method'] == 'denovo':
                otu_method = 'denovo'
            elif kargs['otu_method'] == 'closed':
                otu_method = 'closed'
            elif kargs['otu_method'] == 'r_dada2':
                otu_method = 'r_dada2'
            else:
                secho(f'Error: --otu_method의 값이 잘못되었습니다.', fg='read', blink=True)
                echo(f'kargs["otu_method"]: {kargs["otu_method"]}')
                exit()
            taxonomy_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['taxonomy_assignment'])
            cd_hit_otu_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['cd_hit_otu'])
            closed_dir_path = os.path.join(analysis_dir_path, ANALYSIS_DIR_NAME['closed'])
            biom_dir_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, False)
            metadata_file = glob_dir(analysis_dir_path, 'metadata.txt', p_mode='only', p_verbose=True)
            # Taxonomy_Assignment 디렉터리 목록 확인
            l_glob_taxonomy = glob_dir(taxonomy_dir_path, '*', p_mode='many', p_verbose=False)
            if l_glob_taxonomy is not None:
                l_assignment_dir, not_dir = check_file_type(l_glob_taxonomy, 'isdir')
            # CD-HIT-OTU
            if otu_method == 'denovo':
                clustering_dir_path = glob_dir(cd_hit_otu_dir_path, f'{kargs["order_number"]}*[!.log]',
                                               p_mode='only', p_verbose=False)
                pick_otus_file = os.path.join(clustering_dir_path, 'pick_otus.txt')
                echo('+++ BASE BIOM')
                run_make_otu_table(biom_dir_path, pick_otus_file, kargs['biom'], None, metadata_file)
                for db_dir_path in l_assignment_dir:
                    echo(f'+++ CD-HIT-OTU: {os.path.basename(db_dir_path)}')
                    biom_name = f'{kargs["biom"]}.{os.path.basename(db_dir_path)}'
                    taxonomy_file = os.path.join(db_dir_path, 'otus_rep_tax_assignments.txt')
                    run_make_otu_table(biom_dir_path, pick_otus_file, biom_name, taxonomy_file, metadata_file)
            # CLOSED
            elif otu_method == 'closed':
                pick_otus_name = f'{kargs["order_number"]}.pooled_otus.txt'
                pick_otus_file = os.path.join(closed_dir_path, 'uclust_ref_picked_otus', pick_otus_name)
                echo('+++ BASE BIOM')
                run_make_otu_table(biom_dir_path, pick_otus_file, kargs['biom'], None, metadata_file)
                for db_dir_path in l_assignment_dir:
                    echo(f'+++ CLOSED: {os.path.basename(db_dir_path)}')
                    biom_name = f'{kargs["biom"]}.{os.path.basename(db_dir_path)}'
                    taxonomy_file = os.path.join(db_dir_path, 'otus_rep_tax_assignments.txt')
                    run_make_otu_table(biom_dir_path, pick_otus_file, biom_name, taxonomy_file, metadata_file)
            # ASVs
            elif otu_method == 'r_dada2':
                l_db = list()
                for db_dir_path in l_assignment_dir:
                    db_dir_name = os.path.basename(db_dir_path)
                    if db_dir_name.startswith('blast'):
                        db_name = os.path.basename(db_dir_name).replace('blast_', '')
                        l_db.append(db_name)
                    elif db_dir_name.startswith('uclust'):
                        secho('Error: 개발 중입니다.', fg='red')
                        exit()
                kargs['database'] = l_db
                kargs['metadata'] = metadata_file
                i_biom = Microbe(kargs, mode='biom', analysis_type='NGS')
                i_biom.set_path()
                i_biom.make_biom()

        elif (kargs['otu_map_fp'] != 'None') or (kargs['otu_map_fp'] is not None):
            if out_dir_name_base is None:
                # 분석디렉터리경로/수주번호/Analysis_?
                biom_dir_path = make_dir_using_input_file(kargs['otu_map_fp'], out_dir_name, 3, False)
            else:
                biom_dir_path = make_dir_using_input_file(analysis_dir_path, out_dir_name, 0, False)

            if kargs['taxonomy'] is None:
                run_make_otu_table(biom_dir_path, kargs['otu_map_fp'], kargs['biom'], None, kargs['metadata'])
            else:
                if os.path.isabs(kargs['taxonomy']):
                    taxonomy_file_path = kargs['taxonomy']
                else:
                    taxonomy_file_path = os.path.abspath(kargs['taxonomy'])
                file_path_split = taxonomy_file_path.split('/')
                db_dir_path = file_path_split[-2]
                biom_name = f'{kargs["biom"]}.{db_dir_path}'
                run_make_otu_table(biom_dir_path, kargs['otu_map_fp'], biom_name, kargs['taxonomy'], kargs['metadata'])
        else:
            raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')


# def check_clustering_dir(p_cd_hit_otu_dir_path, p_closed_dir_path):
#     """
#     분석 디렉터리에 CD-HIT-OTU 디렉터리와 CLOSED 디렉터리의 존재 유무를 확인한다.
#     두 개의 디렉터리 모두 존재할 경우, Error 메시지를 출력하고 프로그램을 종료한다.
#     1개만 존재할 경우, 각 디렉터리의 존재유무를 bool 로 반환한다.
#
#     :param p_cd_hit_otu_dir_path: CD-HIT-OTU 디렉터리 경로
#     :param p_closed_dir_path: CLOSED 디렉터리 경로
#     :return: cd_hit_otu_dir_status, closed_dir_status
#     """
#     cd_hit_otu_dir_status = check_file_type(p_cd_hit_otu_dir_path, 'exists')
#     closed_dir_status = check_file_type(p_closed_dir_path, 'exists')
#     if (cd_hit_otu_dir_status is True) and (closed_dir_status is True):
#         secho('Error: CD-HIT-OTU 디렉터리와 CLOSED 디렉터리가 모두 존재합니다.', fg='red')
#         secho('       분석 파이프라인 설계상 같이 존재할 수 없습니다.', fg='red')
#         secho('\t --> Analysis_? 디렉터리를 생성하여, Clustering 디렉터리를 분리하세요.\n'
#               '       참고용일 경우, 디렉터리 이름을 변경하세요.', fg='magenta')
#         exit()
#     else:
#         return cd_hit_otu_dir_status, closed_dir_status
