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

import os
from time import time
from click import secho, echo, progressbar, style
from PreMA.data_structure import PreMA_Pack
from PreMA.rawdata_handler import find_rawdata, check_rawdata, copy_rawdata, check_sample_name, run_multiqc, run_fastqc
from PreMA.adapter_trim import run_adapter_trim, make_excel
from SpoON.util import check_file_type
from SpoON.run_time import compute_run_time


# TODO : 주석 처리된 매개변수 삭제
def core_pipeline(p_pack, p_mode='main'):
    target_data_path = os.path.join(p_pack.rawdata_base_path, p_pack.order_number)
    echo()
    p_pack.nt_rundata_path = find_rawdata(target_data_path, p_pack.target_dir_suffix)
    print(f'p_pack.index_kit: {p_pack.index_kit}')
    p_pack.index_kit_info_from_rawdata, p_pack.l_nt_sample_list = check_rawdata(
        {
            'sample_list': p_pack.nt_rundata_path.sample_list,
            'custom_name': p_pack.dic_custom_sample_name,
            'stat_file': p_pack.nt_rundata_path.stat_file,
            'nread_file': p_pack.nt_rundata_path.nread_file,
            'n_base': p_pack.n_base,
            'n_read': p_pack.n_read,
            'sample_read': p_pack.sample_read,
            'q30': p_pack.q30,
            'order_number': p_pack.order_number,
            'index_kit': p_pack.index_kit,
            'sample_name_mode': p_pack.sample_name_mode,
        }
    )
    if p_pack.copy_rawdata:
        copy_start = time()
        p_pack.l_nt_sample_list_in_analysis, p_pack.l_nt_sample_list_in_analysis_for_fastqc = copy_rawdata(
            {
                'analysis_base_path': p_pack.analysis_base_path,
                'order_number': p_pack.order_number,
                'run_path': p_pack.nt_rundata_path.rundata_path,
                'sample_list': p_pack.l_nt_sample_list,
                'R1_suffix': p_pack.r1_suffix,
                'R2_suffix': p_pack.r2_suffix,
                'copy_fastqc': p_pack.copy_fastqc,
            }, p_mode=p_mode,
        )
        p_pack.save_metadata_group_file(p_pack.l_nt_sample_list_in_analysis[0].path)
        copy_end = time()
        copy_run_time = compute_run_time(copy_start, copy_end)

        if p_pack.copy_fastqc and p_pack.l_nt_sample_list_in_analysis_for_fastqc:
            fastqc_start = time()
            run_fastqc(p_pack.l_nt_sample_list_in_analysis_for_fastqc, p_pack.r1_suffix, p_pack.r2_suffix)
            fastqc_end = time()
            fastqc_run_time = compute_run_time(fastqc_start, fastqc_end)
        else:
            fastqc_run_time = 'None'

        trim_start = time()
        p_pack.l_summary_files = run_adapter_trim(
            {
                'sample_list': p_pack.l_nt_sample_list_in_analysis,
                'R1_suffix': p_pack.r1_suffix,
                'R2_suffix': p_pack.r2_suffix,
                'index_kit': p_pack.index_kit_info_from_rawdata,
                'my_story': p_pack.my_story,
                'microbe_and_me': p_pack.microbe_and_me,
                'trim_tail': p_pack.trim_tail,
            },
            p_pack.adapter_trim_tool
        )
        trim_end = time()
        trim_run_time = compute_run_time(trim_start, trim_end)

        excel_start = time()
        make_excel(
            {
                'files': p_pack.l_summary_files,
                'output_path': p_pack.l_nt_sample_list_in_analysis[0].path,
                'order_number': p_pack.order_number,
            },
            p_pack.adapter_trim_tool
        )
        excel_end = time()
        excel_run_time = compute_run_time(excel_start, excel_end)
    else:
        copy_run_time = 'None'
        fastqc_run_time = 'None'
        trim_run_time = 'None'
        excel_run_time = 'None'
    d_run_time = {'copy_run_time': copy_run_time, 'fastqc_run_time': fastqc_run_time,
                  'trim_run_time': trim_run_time, 'excel_run_time': excel_run_time}
    return d_run_time


def check_integrated_order(kargs):
    l_stats = list()
    for target in kargs['integrate']:
        target_data_path = os.path.join(kargs['rawdata_base_path'], target)
        l_stats.append(check_file_type(target_data_path, 'exists'))
    del target
    l_error_target = list()
    for (target, status) in zip(kargs['integrate'], l_stats):
        if status is False:
            l_error_target.append(target)
    if all(l_stats) is False:
        secho('Error: RawData 보관 디렉터리에 해당 수주가 존재하지 않습니다.', fg='red', blink=True, err=True)
        secho('{}'.format(', '.join(l_error_target, )), fg='bright_magenta', err=True)
        secho('** RawData 를 요청하세요! **', fg='green', reverse=True, err=True)
        exit(1)


def run_core(kargs):
    pack = PreMA_Pack(kargs)
    if kargs['integrate']:  # 미지정시 --> () 빈튜플
        check_integrated_order(kargs)

    if kargs['custom_sample_name'] is not None:
        pack.parse_custom_sample_name(kargs['custom_sample_name'])
        print(pack.dic_custom_sample_name)

    main_run_time = core_pipeline(pack, p_mode='main')
    # if (kargs['custom_sample_name'] is not None) and (kargs['copy_rawdata'] is True):
    #     pack.copy_custom_sample_name_file(kargs['custom_sample_name'])
    #     pack.save_custom_sort_file(kargs['custom_sample_name'])

    if kargs['integrate']:
        integ_run_time = list()
        for integ_order_num in kargs['integrate']:
            new_kargs = dict()
            new_kargs.update(kargs)
            new_kargs['analysis_base_path'] = os.path.join(pack.analysis_base_path, pack.order_number)
            new_kargs['order_number'] = integ_order_num
            integ_pack = PreMA_Pack(new_kargs)
            if kargs['custom_sample_name'] is not None:
                integ_pack.parse_custom_sample_name(kargs['custom_sample_name'])
            integ_run_time.append(core_pipeline(integ_pack, p_mode='integrate'))
            # 메인 수주에 통합수주의 시료 정보 삽입
            pack.l_nt_sample_list += integ_pack.l_nt_sample_list
    else:
        integ_run_time = 'None'

    if (kargs['custom_sample_name'] is not None) and (kargs['copy_rawdata'] is True):
        pack.copy_custom_sample_name_file(kargs['custom_sample_name'])
        pack.save_custom_sort_file(kargs['custom_sample_name'])

    if kargs['copy_rawdata'] and (kargs['no_multiqc'] is False):
        multiqc_start = time()
        run_multiqc(
            {
                'target_dir': os.path.join(pack.analysis_base_path, pack.order_number, 'RawData*'),
                'order_number': pack.order_number,
                'depth': 4,
                'output_path': os.path.join(pack.analysis_base_path, pack.order_number),
            }
        )
        multiqc_end = time()
        multiqc_run_time = compute_run_time(multiqc_start, multiqc_end)
    else:
        multiqc_run_time = 'None'
    return main_run_time, integ_run_time, multiqc_run_time


def check_custom_name(kargs):
    l_mini_pack = list()

    l_mini_pack.append(PreMA_Pack(kargs, 'mini'))
    if kargs['integrate']:
        for integ_order_num in kargs['integrate']:
            new_kargs = dict()
            new_kargs.update(kargs)
            new_kargs['order_number'] = integ_order_num
            l_mini_pack.append(PreMA_Pack(new_kargs, 'mini'))
        del integ_order_num
    if kargs['custom_sample_name'] is not None:
        for mini_pack in l_mini_pack:
            mini_pack.parse_custom_sample_name(kargs['custom_sample_name'])
        del mini_pack
    for mini_pack in l_mini_pack:
        target_data_path = os.path.join(mini_pack.rawdata_base_path, mini_pack.order_number)
        mini_pack.nt_rundata_path = find_rawdata(target_data_path, mini_pack.target_dir_suffix)
        mini_pack.l_nt_sample_list = check_sample_name(
            mini_pack.nt_rundata_path.sample_list,
            mini_pack.dic_custom_sample_name,
            mini_pack.sample_name_mode
        )
        echo("=" * 30)
        mini_pack.check_custom_name()
        echo("=" * 30)


def run_length_trim(p_l_nt_run_list, p_length, p_r1_suffix, p_r2_suffix):
    out_r1_suffix = p_r1_suffix.replace('.fastq', f'.{p_length}bp.fastq')
    out_r2_suffix = p_r2_suffix.replace('.fastq', f'.{p_length}bp.fastq')
    from SpoON.fastq_handler import get_file_object, length_trim
    for run in p_l_nt_run_list:
        with progressbar(run.samples, label=f'Length Trimming({p_length}bp)',
                         fill_char=style('#', fg='green')) as bar_samples:
            for sample in bar_samples:
                output_r1 = os.path.join(sample.path, sample.name+out_r1_suffix)
                output_r2 = os.path.join(sample.path, sample.name+out_r2_suffix)
                o_in_r1 = get_file_object(os.path.join(sample.path, sample.name+p_r1_suffix))
                o_in_r2 = get_file_object(os.path.join(sample.path, sample.name+p_r2_suffix))
                o_out_r1 = get_file_object(output_r1, p_mode='write')
                o_out_r2 = get_file_object(output_r2, p_mode='write')
                length_trim(o_in_r1, o_out_r1, p_length)
                length_trim(o_in_r2, o_out_r2, p_length)
