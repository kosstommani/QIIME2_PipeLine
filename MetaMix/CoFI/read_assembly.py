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
__version__ = '1.0.1'

import os
from multiprocessing import Pool
from click import secho, echo, style, progressbar
from SpoON.util import typeWrite, typeWrite_float, run_cmd, check_run_cmd, parse_config, start_process
from CoFI.data_structure import FLASH_List

# 기본값 설정
CONFIG = parse_config()
FLASH = CONFIG['CoFI_ReadAssembly']['FLASH']
COMPRESS = CONFIG['CoFI_ReadAssembly']['compress']
COMPRESS_SUFFIX = CONFIG['CoFI_ReadAssembly']['compress_suffix']
FLASH_THREADS = CONFIG['CoFI_ReadAssembly']['FLASH_threads']
POOL_WORKER = CONFIG['CoFI_PoolWorker']
FLASH_POOL_WORKER = CONFIG['CoFI_FLASH_PoolWorker']
ANALYSIS_DIR_NAME = CONFIG['Analysis_Dir_Name']


# FLASH -t, --threads=NTHREADS 조절 완료.
# PoolWorker: 20, -t Default(number of processors) ==> BlockingIOError: [Errno 11] Resource temporarily unavailable
# PoolWorker: 10, -t Default ==> OK.
def get_flash1_cmd(kargs):
    """

    :type kargs: dict
    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                in1:
                in2:
                out_path:
                out_name:
                min_overlap:
                max_overlap:
                no_compress:
                no_use_worker:
    :return:
    """
    global FLASH, COMPRESS, COMPRESS_SUFFIX, FLASH_THREADS
    max_overlap = kargs['read_length'] * 2 - kargs['target_size'] + 10
    mix_overlap = max_overlap - 20
    if kargs['min_overlap'] is not None:
        mix_overlap = kargs['min_overlap']
    if kargs['max_overlap'] is not None:
        max_overlap = kargs['max_overlap']
    if kargs['no_compress']:
        compress = ''
        compress_suffix = ''
    else:
        compress = COMPRESS
        compress_suffix = COMPRESS_SUFFIX
    flash_cmd = \
        '{flash} ' \
        '{input1} ' \
        '{input2} ' \
        '-d {out_path} ' \
        '-o {out_name} ' \
        '-M {max} ' \
        '-m {min} ' \
        '{compress} ' \
        '{compress_suffix} ' \
        '--threads {threads} ' \
        '2>&1 | tee ' \
        '{log}'.format(
            flash=FLASH,
            input1=kargs['in1'],
            input2=kargs['in2'],
            out_path=kargs['out_path'],
            out_name=kargs['out_name'],
            max=max_overlap,
            min=mix_overlap,
            compress=compress,
            compress_suffix=compress_suffix,
            threads=FLASH_THREADS,
            log=os.path.join(*(kargs['out_path'], kargs['out_name'] + '.FLASH.log'))
        )
    return flash_cmd


def run_flash(kargs):
    """

    :type kargs: dict
    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                run_list:
                R1_suffix:
                R2_suffix:
                output_path:
                order_number:
                read_length:
                target_size:
                min_overlap:
                max_overlap:
                no_compress:
                no_csv:
                no_excel:
                no_hcsv:
                no_use_worker:
    :return: named tuple 인 FLASH_List 을 원소로 가지는 리스트 반환.
             FLASH_List = namedtuple('FLASH_File_List', 'path name output')
    """
    global POOL_WORKER, FLASH_POOL_WORKER, ANALYSIS_DIR_NAME
    l_cmd = list()
    l_nt_flash_list = list()
    if kargs['no_compress']:
        output_suffix = '{}.extendedFrags.fastq'
    else:
        output_suffix = '{}.extendedFrags.fastq.' + COMPRESS_SUFFIX.strip().split()[1]
    flash_output_path = os.path.join(*(kargs['output_path'], ANALYSIS_DIR_NAME['read_assembly']))

    o_flash_cmd_log = open(os.path.join(*(kargs['output_path'], 'FLASH_CMD.recipe')), 'w')
    for run in kargs['run_list']:
        for sample in run.samples:
            cmd = get_flash1_cmd(
                {
                    'in1': os.path.join(*(sample.path, sample.name + kargs['R1_suffix'])),
                    'in2': os.path.join(*(sample.path, sample.name + kargs['R2_suffix'])),
                    'out_path': flash_output_path,
                    'out_name': sample.name,
                    'read_length': kargs['read_length'],
                    'target_size': kargs['target_size'],
                    'min_overlap': kargs['min_overlap'],
                    'max_overlap': kargs['max_overlap'],
                    'no_compress': kargs['no_compress'],
                    'no_use_worker': kargs['no_use_worker'],
                }
            )
            l_cmd.append(cmd)
            l_nt_flash_list.append(FLASH_List(path=flash_output_path,
                                              name=sample.name,
                                              output=output_suffix.format(sample.name)
                                              )
                                   )
            echo(cmd, file=o_flash_cmd_log)
            echo(file=o_flash_cmd_log)
    o_flash_cmd_log.close()
    del cmd

    os.mkdir(flash_output_path)
    secho('>>> Read Assembly 디렉터리 생성', fg='cyan')
    echo(flash_output_path)

    # 시료가 많은 경우 FLASH 실행이 Error 나는 경우 발생.
    if kargs['no_use_worker']:
        echo('>>> Read Assembly 시작 : FLASH')
        pool_outputs = list()
        with progressbar(l_cmd, label='FLASH', fill_char=style('#', fg='green')) as bar_l_cmd:
            for cmd in bar_l_cmd:
                pool_outputs.append(run_cmd(cmd))
    else:
        # BlockingIOError 해결
        process = FLASH_POOL_WORKER
        if len(l_cmd) < FLASH_POOL_WORKER:
            process = len(l_cmd)
        echo('>>> Read Assembly 시작 : FLASH')
        with Pool(process, initializer=start_process) as pool:
            pool_outputs = pool.map(run_cmd, l_cmd)
            pool.close()
            pool.join()

    done_count = len([1 for run in pool_outputs
                      if (run.returncode == 0) and ('FLASH v1.2.11 complete!' in run.stdout)])
    # FLASH 에서 Error 발생시 Error 에 대한 return code 을 출력하지 않아 출력 메시지로 확인.
    error_count = len([1 for run in pool_outputs
                       if (run.returncode != 0) or
                       ('FLASH did not complete successfully' in run.stdout)])
    done_text = style('완료: {}'.format(done_count), fg='cyan')
    error_text = style('에러: {}'.format(error_count),
                       fg='red' if (error_count != 0) else 'cyan',
                       blink=True if (error_count != 0) else False)
    echo('{msg} {done}, {error}'.format(
        msg=style('>>> Read Assembly', fg='cyan'),
        done=done_text,
        error=error_text))

    if error_count != 0:
        for run in pool_outputs:
            if run.returncode != 0 or 'FLASH did not complete successfully' in run.stdout:
                run.returncode = 1  # 에러 메시지로 확인한 경우
                if run.stderr:
                    pass
                else:
                    run.stderr = run.stdout
                check_run_cmd({
                    'run': run,
                    'true_meg': None,
                    'false_meg': 'FLASH 실행'
                }, p_exit=False)
        exit(1)

    l_flash_logs = [run.stdout for run in pool_outputs if run.returncode == 0]
    csv = False if kargs['no_csv'] else True
    excel = False if kargs['no_excel'] else True
    hcsv = False if kargs['no_hcsv'] else True
    make_excel_for_flash(
        {
            'log_texts': l_flash_logs,
            'save_path': flash_output_path,
            'output_name': kargs['order_number'],
            'R1_suffix': kargs['R1_suffix'],
            'R2_suffix': kargs['R2_suffix'],
            'target_size': kargs['target_size'],
            'percent_combined_alarm': kargs['percent_combined_alarm']
        },
        p_csv=csv, p_excel=excel, p_Hcsv=hcsv,
        p_no_compress=kargs['no_compress']
    )
    return l_nt_flash_list


def read_hist(p_file):
    """
    FLASH 프로그램에서 출력된 .hist 파일을 읽어 딕션너리 구조로 반환한다.
    
    :param p_file: FLASH.hist 파일
    :return: 딕션너리
            key: 길이(int)
            value: 개수(int)
    """
    d_data = dict()
    with open(p_file, 'r') as o_hist:
        lines = o_hist.readlines()
        for line in lines:
            key, value = line.strip().split('\t')
            d_data[int(key)] = int(value)
    return d_data


def make_excel_for_flash(kargs, p_csv=True, p_excel=True, p_Hcsv=True, p_no_compress=False):
    """
    FLASH 프로그램에서 출력하는 메세지를 캡쳐하여 정리한 후, 정리한 데이터를 다양한 형식(csv, Hcsv, xlsx)으로 나타낸다.
    FLASH.hist 파일을 이용하여 길이 분포 그래프를 엑셀에 생성한다.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                log_texts: [list] - FLASH 로그를 하나의 문자열로 읽은 원소를 여러개 가지는 리스트. [log, [log]...]
                save_path: [str] - 결과물 저장 위치
                output_name: [str] - 결과파일명
                R1_suffix: [str] - Read1의 파일명 접미사
                R2_suffix: [str] - Read2의 파일명 접미사
                target_size: [int] - 목표영역의 크기
                percent_combined_alarm: [float] - percent combined 의 비율이 해당값 미만일 경우 강조. 기본값: 40.0(%)
    :param p_csv: [True(default), False] - csv 파일 생성 여부
    :param p_excel: [True(default), False] - 엑셀 파일 생성 여부
    :param p_Hcsv: [True(default), False] - 서식 적용 csv 파일 생성 여부
    :param p_no_compress: FLASH 출력물의 압축여부
    :return:
    """
    l_d_hist_data = list()
    sample_name = ["Sample"]
    min_overlap = ["min_overlap"]
    max_overlap = ["max_overlap"]
    max_mismatch_density = ["max_mismatch_density"]
    allow_outie_pairs = ["allow_outie_pairs"]
    cap_mismatch_quals = ["cap_mismatch_quals"]
    target_size = ["target_size"]
    total_pairs = ["total_pairs"]
    combined_pairs = ["combined_pairs"]
    uncombined_pairs = ["uncombined_pairs"]
    percent_combined = ["percent_combined"]
    if p_no_compress:
        output_suffix = '.extendedFrags.fastq'
    else:
        output_suffix = '.extendedFrags.fastq.' + COMPRESS_SUFFIX.strip().split()[1]
    for sample_text in kargs['log_texts']:
        l_text = sample_text.split('\n')
        while True:
            try:
                text = l_text.pop(0)
            except IndexError:
                break

            if 'Input files' in text:
                input_file1 = l_text.pop(0).split()[1].split('/')[-1].replace(kargs['R1_suffix'], '')
                input_file2 = l_text.pop(0).split()[1].split('/')[-1].replace(kargs['R2_suffix'], '')
                l_text.pop(0)  # 공백 제거
                l_text.pop(0)  # 제거: [FLASH] Output files:
                output_file = l_text.pop(0).split()[1].split('/')[-1].replace(output_suffix, '')
                l_text.pop(0)  # 제거: notCombined_1
                l_text.pop(0)  # 제거: notCombined_2
                hist_file = l_text.pop(0).split()[1]
                l_text.pop(0)  # 제거: .histogram
                if input_file1 == input_file2 == output_file:
                    sample_name.append(output_file)
                    l_d_hist_data.append(read_hist(hist_file))
                    continue
                else:
                    msg = style('Error: FLASH 시료명 다름.\n' +
                                'input: {0}, {1}, output: {2}'.format(
                                    input_file1,
                                    input_file2,
                                    output_file), fg='red'
                                )
                    raise RuntimeError(msg)
            elif 'Parameters' in text:
                for ele in [min_overlap, max_overlap, max_mismatch_density, allow_outie_pairs, cap_mismatch_quals]:
                    ele.append(l_text.pop(0).split(':')[1].strip())
                target_size.append(kargs['target_size'])
                continue
            elif 'Read combination statistics' in text:
                for ele in [total_pairs, combined_pairs, uncombined_pairs, percent_combined]:
                    ele.append(l_text.pop(0).split(':')[1].strip())
                continue
    
    # hist 데이터 정리
    l_s_length = [set(x.keys()) for x in l_d_hist_data]
    s_length_list = set()
    _ = [s_length_list.update(length) for length in l_s_length]

    l_data_vars = [sample_name, total_pairs, combined_pairs, uncombined_pairs, percent_combined,
                   target_size, min_overlap, max_overlap, max_mismatch_density, allow_outie_pairs,
                   cap_mismatch_quals]

    # FLASH_log.csv 생성 #
    if p_csv:
        csv_file = os.path.join(kargs['save_path'], kargs['output_name'] + ".FLASH.csv")
        with open(csv_file, "w") as csv:
            for csv_data in l_data_vars:
                csv.write(",".join([str(x) for x in csv_data]))
                csv.write("\n")
        echo('>>> FLASH_log.csv 생성 완료')

    # FLASH_log.excel 생성 & FLASH_hist.csv 생성
    if p_excel:
        import xlsxwriter
        percent_combined_remove_mark = [X.replace("%", "") for X in percent_combined]  # % 기호 제거
        excel_file = os.path.join(kargs['save_path'], kargs['output_name'] + ".FLASH.xlsx")
        hist_csv_file = os.path.join(kargs['save_path'], kargs['output_name'] + '.hist.csv')
        workbook = xlsxwriter.Workbook(excel_file)
        o_hist_csv = open(hist_csv_file, 'w')

        number_type_format = workbook.add_format({'num_format': '#,##0'})

        sheet_name = 'STAT'
        stat_worksheet = workbook.add_worksheet(sheet_name)
        header_list = []
        for i in sample_name:
            header_list.append({"header": i})
        stat_worksheet.add_table(0, 0, 4, len(sample_name) - 1,
                                 {'columns': header_list, 'style': 'Table Style Light 8', 'autofilter': False})
        for col_num, cell in enumerate(total_pairs):
            typeWrite(stat_worksheet, 1, col_num, cell, number_type_format)

        for col_num, cell in enumerate(combined_pairs):
            typeWrite(stat_worksheet, 2, col_num, cell, number_type_format)

        for col_num, cell in enumerate(uncombined_pairs):
            typeWrite(stat_worksheet, 3, col_num, cell, number_type_format)

        for col_num, cell in enumerate(percent_combined_remove_mark):
            typeWrite(stat_worksheet, 4, col_num, cell, number_type_format)

        read_count_chart = workbook.add_chart({'type': 'column'})
        percent_combined_chart = workbook.add_chart({'type': 'column'})

        # total pairs chart #
        read_count_chart.add_series(
            {
                'values': [sheet_name, 1, 1, 1, len(sample_name) - 1],
                'categories': [sheet_name, 0, 1, 0, len(sample_name) - 1],
                'name': [sheet_name, 1, 0],
            })

        # combined pairs chart #
        read_count_chart.add_series(
            {
                'values': [sheet_name, 2, 1, 2, len(sample_name) - 1],
                'categories': [sheet_name, 0, 1, 0, len(sample_name) - 1],
                'name': [sheet_name, 2, 0],
            })

        # uncombined pairs chart #
        read_count_chart.add_series(
            {
                'values': [sheet_name, 3, 1, 3, len(sample_name) - 1],
                'categories': [sheet_name, 0, 1, 0, len(sample_name) - 1],
                'name': [sheet_name, 3, 0],
            })

        # percent combined chart #
        percent_combined_chart.add_series(
            {
                'values': [sheet_name, 4, 1, 4, len(sample_name) - 1],
                'categories': [sheet_name, 0, 1, 0, len(sample_name) - 1],
                'name': [sheet_name, 4, 0],
            })

        read_count_chart.set_x_axis({'name': 'Sample Name'})
        read_count_chart.set_title({'name': 'FLASH - Read count by states'})
        percent_combined_chart.set_x_axis({'name': 'Sample Name'})
        percent_combined_chart.set_title({'name': 'FLASH - Percentage of combined reads'})

        read_count_chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
        percent_combined_chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
        stat_worksheet.insert_chart('A8', read_count_chart)
        stat_worksheet.insert_chart('A25', percent_combined_chart)

        # worksheet : hist_count, hist_percentage #
        hist_count_sheet_name = 'HIST_Count'
        hist_percentage_sheet_name = 'HIST_Percentage'
        worksheet_hist_count = workbook.add_worksheet(hist_count_sheet_name)
        worksheet_hist_percentage = workbook.add_worksheet(hist_percentage_sheet_name)
        l_length_worksheet = [worksheet_hist_count, worksheet_hist_percentage]
        hist_header = sample_name[:]
        hist_header[0] = 'Length'
        last_col_num = len(hist_header) - 1
        for col_num, cell in enumerate(hist_header):
            worksheet_hist_count.write(0, col_num, cell)
            worksheet_hist_percentage.write(0, col_num, cell)
            o_hist_csv.write(str(cell))
            if col_num < last_col_num:
                o_hist_csv.write(',')
        else:
            o_hist_csv.write('\n')
            del last_col_num
        for row_num, length in enumerate(sorted(list(s_length_list)), 1):
            typeWrite(worksheet_hist_count, row_num, 0, length)
            typeWrite_float(worksheet_hist_percentage, row_num, 0, length)
            o_hist_csv.write(str(length) + ',')
            last_col_num = len(l_d_hist_data)
            for col_num, d_hist in enumerate(l_d_hist_data, 1):
                total_read_count = sum(d_hist.values())
                count = d_hist.get(length)
                if count is None:
                    count = 0
                typeWrite(worksheet_hist_count, row_num, col_num, count)
                typeWrite_float(worksheet_hist_percentage, row_num, col_num, (count/total_read_count)*100)
                o_hist_csv.write(str(count))
                if col_num < last_col_num:
                    o_hist_csv.write(',')
            else:
                o_hist_csv.write('\n')
                del last_col_num
        else:
            o_hist_csv.close()
            echo('>>> FLASH_hist.csv 생성 완료')

        last_row_num = row_num
        del col_num, row_num
        for col_num, _ in enumerate(list(s_length_list), 1):
            for length_worksheet in l_length_worksheet:
                length_worksheet.conditional_format(1, 1, last_row_num, col_num,
                                                    {'type': 'data_bar', 'bar_solid': True, 'bar_color': 'green'})

        hist_header_list = header_list[:]
        hist_header_list[0] = {'header': 'Length'}
        for length_worksheet in l_length_worksheet:
            length_worksheet.add_table(0, 0, last_row_num, len(hist_header) - 1,
                                       {
                                           'columns': hist_header_list,
                                           'style': 'Table Style Light 8',
                                           'autofilter': False
                                       })

        length_count_chart = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
        length_percentage_chart = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
        for col_num, _ in enumerate(list(l_d_hist_data), 1):
            length_count_chart.add_series(
                {
                    'values': [hist_count_sheet_name, 1, col_num, last_row_num, col_num],
                    'categories': [hist_count_sheet_name, 1, 0, last_row_num, 0],
                    'name': [hist_count_sheet_name, 0, col_num]
                }
            )
        for col_num, _ in enumerate(list(l_d_hist_data), 1):
            length_percentage_chart.add_series(
                {
                    'values': [hist_percentage_sheet_name, 1, col_num, last_row_num, col_num],
                    'categories': [hist_percentage_sheet_name, 1, 0, last_row_num, 0],
                    'name': [hist_percentage_sheet_name, 0, col_num]
                }
            )
        for chart in [length_count_chart, length_percentage_chart]:
            chart.set_x_axis({'name': 'Assembled Read Length(bp)'})
            chart.set_title({'name': 'FLASH - Length Distribution of Assembled Read'})
            chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
        length_count_chart.set_y_axis({'name': 'Count'})
        length_percentage_chart.set_y_axis({'name': 'Percentage(%)'})
        stat_worksheet.insert_chart('A42', length_count_chart)
        stat_worksheet.insert_chart('A59', length_percentage_chart)
        workbook.close()
        echo('>>> FLASH_log.xlsx 생성 완료')

    # FLASH_log.Hcsv 생성 #
    if p_Hcsv:
        Hcsv_file = os.path.join(*(kargs['save_path'], kargs['output_name'] + ".FLASH.Hcsv"))
        with open(Hcsv_file, "w") as Hcsv:
            temp_index_number = range(1, len(sample_name) + 1)
            index_number = [str(X) for X in temp_index_number]
            del temp_index_number

            for index in range(1, len(percent_combined[1:]) + 1):
                if float(percent_combined[index].replace("%", "")) < kargs['percent_combined_alarm']:
                    l_data_vars.insert(0, index_number)
                    for excel_data in l_data_vars:
                        try:
                            excel_data[index] = style(excel_data[index], fg='red')
                        except TypeError:
                            excel_data[index] = style(str(excel_data[index]), fg='red')

            for Hcsv_data in l_data_vars:
                Hcsv.write(','.join([str(x) for x in Hcsv_data]))
                Hcsv.write('\n')
        echo('>>> FLASH_log.Hcsv 생성 완료')
