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
__version__ = '1.0.0'

import os
from click import style, secho, echo
from multiprocessing import Pool
from SpoON.util import run_cmd, run_move_and_cmd, check_run_cmd, parse_config, check_file_type, start_process, \
    typeWrite
from CoFI.data_structure import FilteredFastq  # 사용


CONFIG = parse_config()
POOL_WORKER = CONFIG['CoFI_PoolWorker']


def get_file_object(p_file, p_mode='read'):
    """
    압축파일의 확장자명에 따라 적합한 파일 객체를 반환한다.
    압축파일이 아닐 경우, 기본 파일 객체를 반환한다.
    가능한 압축 형식: gz, bz, xz

    :type p_file: str
    :param p_file: 파일명.

    :type p_mode: str
    :param p_mode: read or write
    :return: 파일 객체
    """
    if p_mode == 'read':
        if p_file.endswith('gz'):
            import gzip
            return gzip.open(p_file, 'rt', encoding='utf-8')
        elif p_file.endswith('bz'):
            import bz2
            return bz2.open(p_file, 'rt', encoding='utf-8')
        elif p_file.endswith('xz'):
            import lzma
            return lzma.open(p_file, 'rt', encoding='utf-8')
        else:
            return open(p_file, 'r')
    elif p_mode == 'write':
        if p_file.endswith('gz'):
            import gzip
            return gzip.open(p_file, 'wt', encoding='utf-8')
        elif p_file.endswith('bz'):
            import bz2
            return bz2.open(p_file, 'wt', encoding='utf-8')
        elif p_file.endswith('xz'):
            import lzma
            return lzma.open(p_file, 'wt', encoding='utf-8')
        else:
            return open(p_file, 'w')
    else:
        raise ValueError('read or write: 입력값 --> {}'.format(p_mode))


def filter_length(p_file, p_lower, p_upper, p_modify_header, p_seq_type):
    """
    Assembled Read 들의 서열 길이가 제시된 필터 조건을 충족하지 않는 서열들을 제거한다.
    또한, FASTQ Header 정보를 '>시료명_SeqNum'으로 변경한다.
    필터 미적용(p_lower=p_upper=no)일 때, p_modify_header 가 True 이면 'F_파일명' 으로 데이터 생산.
    필터 미적용(p_lower=p_upper=no)일 때, p_modify_header 가 False 이면 미적용.
    필터 미적용(p_lower=p_upper=no)이면서 헤더 미변경(p_modify_header=False)시 p_seq_type 은 적용되지 않는다.
    필터 미적용(p_lower=p_upper=no)이면서 헤더 변경(p_modify_header=True)시 p_seq_type 은 적용된다.

    :type p_file: str
    :param p_file: Assembled Read File

    :type p_lower: no or int
    :param p_lower: Assembled Read 의 최소 필터 길이. 
                False: 미적용.
                숫자: 기준값미만 모두 제거.

    :type p_upper: no or int
    :param p_upper: Assembled Read 의 최대 필터 길이.
                False: 미적용.
                숫자: 기준값초과 모두 제거.
                
    :type p_modify_header: bool
    :param p_modify_header: FASTQ Header 정보 변경.

    :type p_seq_type: str
    :param p_seq_type: 서열정보의 파일 양식.('FASTQ' or 'FASTA')
                       필터 미적용 & 헤더 미변경일 경우 적용 안됨.
    :return: named tuple - FilteredFastq
                FilteredFastq(
                    path= 경로,
                    src= 경로를 포함한 Assembled Read 파일명. 시료명.extendedFrags.fastq[.gz|.bz|.xz],
                    des= 경로를 포함한 Filtered 파일명. F_시료명.extendedFrags.fastq[.gz|.bz|.xz],
                    log= 경로를 포함한 로그파일명,
                    stat=
                    {
                        'name': 시료명,
                        'total_read': total_read,
                        'removed_read': removed_read,
                        'lower_read': lower_read,
                        'upper_read': upper_read,
                        'final_read': final_read,
                        'removed_rate': removed_rate,
                        'lower_bp': p_lower,
                        'upper_bp': p_upper,
                    })
    """
    # TODO: 속도 느림. C 코드 변환. 또는 C로 작성.
    input_file = get_file_object(p_file, 'read')
    path = os.path.split(p_file)

    if (p_lower == 'no') and (p_upper == 'no') and (p_modify_header is False):
        secho('>>> No Filter 적용', fg='yellow')
        echo(p_file)
        object_text = 'FilteredFastq(path="{path}", src="{file}", des="{des}", log="{log}", stat={stat})'.format(
            path=path[0], file=p_file, des=p_file, log=None, stat=None)
        return object_text
    else:
        if p_seq_type.upper() == 'FASTQ':
            if 'fastq' in path[1]:
                output_file = os.path.join(*(path[0], "F_" + path[1]))
            else:
                secho('Error: 해당 파일은 FASTQ 형식으로 만들 수 없습니다.', fg='red', err=True)
                echo(path[1], err=True)
                echo('옵션을 변경하세요!', err=True)
                exit(1)
        elif p_seq_type.upper() == 'FASTA':
            if 'fastq' in path[1]:
                output_file = os.path.join(*(path[0], "F_" + path[1].replace('fastq', 'fasta')))
            elif 'fasta' in path[1]:  # FLASH - FASTQ 만 가능.
                output_file = os.path.join(*(path[0], "F_" + path[1]))
            else:
                secho('>>> 지원하지 않는 형식입니다.(가능 형식: fastq or fasta)', fg='red', err=True)
                echo(path[1], err=True)
                exit(1)
        else:
            raise ValueError('p_seq_type: {}. FASTQ or FASTA 만 가능.'.format(p_seq_type))

    if p_modify_header:
        sample_name = path[1][:path[1].find('.extendedFrags.fastq')]  # FLASH - FASTQ만 가능.

    with get_file_object(output_file, 'write') as output:
        total_read = 0
        removed_read = 0
        lower_read = 0
        upper_read = 0
        final_read = 0
        line = 0
        error = False
        while True:
            header = input_file.readline().strip()  # header
            line += 1
            if header:
                pass
            elif (line == 1) and header is False:
                secho('Error: FASTQ Format - Header', fg='red')
                echo('Line: {}'.format(line))
                echo('File: {}'.format(p_file))
                echo('Header: {}'.format(header))
                error = True
                break  # FASTQ Format Error
            else:
                break  # EOF
            sequence = input_file.readline().strip()  # sequence
            line += 1
            plus = input_file.readline().strip()  # +
            line += 1
            if plus == "+":
                pass
            else:
                secho('Error: FASTQ Format - (+)', fg='red')
                echo('Line: {}'.format(line))
                echo('File: {}'.format(p_file))
                echo('Separator: {}'.format(plus))
                error = True
                break  # FASTQ Format Error
            quality_score = input_file.readline().strip()  # quality score
            line += 1

            if (p_lower == 'no') and (p_upper == 'no') and p_modify_header:
                final_read += 1
                output.write('>{name}_{count}\n'.format(name=sample_name, count=final_read))
                output.write(sequence + '\n')
                if p_seq_type.upper() == 'FASTQ':
                    output.write(plus + '\n')
                    output.write(quality_score + '\n')

            elif (p_lower != 'no') and (p_upper != 'no'):
                if int(p_lower) <= len(sequence) <= int(p_upper):
                    final_read += 1
                    if p_modify_header:
                        output.write('>{name}_{count}\n'.format(name=sample_name, count=final_read))
                    else:
                        output.write(header + '\n')
                    output.write(sequence + '\n')
                    if p_seq_type.upper() == 'FASTQ':
                        output.write(plus + '\n')
                        output.write(quality_score + '\n')
                else:
                    removed_read += 1
                    if len(sequence) <= int(p_lower):
                        lower_read += 1
                    if len(sequence) >= int(p_upper):
                        upper_read += 1
                    del header, sequence, plus, quality_score

            elif (p_lower != 'no') and (p_upper == 'no'):
                if len(sequence) >= int(p_lower):
                    final_read += 1
                    if p_modify_header:
                        output.write('>{name}_{count}\n'.format(name=sample_name, count=final_read))
                    else:
                        output.write(header + '\n')
                    output.write(sequence + '\n')
                    if p_seq_type.upper() == 'FASTQ':
                        output.write(plus + '\n')
                        output.write(quality_score + '\n')

                else:
                    removed_read += 1
                    if len(sequence) <= int(p_lower):
                        lower_read += 1
                    del header, sequence, plus, quality_score

            elif (p_lower == 'no') and (p_upper != 'no'):
                if len(sequence) <= int(p_upper):
                    final_read += 1
                    if p_modify_header:
                        output.write('>{name}_{count}\n'.format(name=sample_name, count=final_read))
                    else:
                        output.write(header + '\n')
                    output.write(sequence + '\n')
                    if p_seq_type.upper() == 'FASTQ':
                        output.write(plus + '\n')
                        output.write(quality_score + '\n')
                else:
                    removed_read += 1
                    if len(sequence) >= int(p_upper):
                        upper_read += 1
                    del header, sequence, plus, quality_score
            total_read += 1
    input_file.close()

    if error:
        return 1

    log_file = os.path.join(*(path[0], path[1] + '.log'))
    with open(log_file, 'w') as log_output:
        removed_rate = 0 if float(removed_read) == 0 else float(removed_read) / total_read
        log_output.write("Total read   : {}\n".format(total_read))
        log_output.write("Removed read : {}\n".format(removed_read))
        log_output.write("lower        : {}\n".format(lower_read))
        log_output.write("upper        : {}\n".format(upper_read))
        log_output.write("Final read   : {}\n".format(final_read))
        log_output.write("Removed rate : {}\n".format(removed_rate))
        log_output.write("Removed < bp : {}\n".format(p_lower))
        log_output.write("Removed > bp : {}\n".format(p_upper))

    stat = "'name': '{name}', 'total_read': {total_read}, 'removed_read': {removed_read}, " \
           "'lower_read': {lower_read}, 'upper_read': {upper_read},'final_read': {final_read}, " \
           "'removed_rate': {removed_rate}, 'lower_bp': '{lower}', 'upper_bp': '{upper}'".format(
            name=path[1], total_read=total_read, removed_read=removed_read,
            lower_read=lower_read, upper_read=upper_read, final_read=final_read,
            removed_rate=removed_rate, lower=p_lower, upper=p_upper)
    dict_stat_text = '{' + stat + '}'
    object_text = f'FilteredFastq(path="{path[0]}", src="{p_file}", ' \
                  f'des="{output_file}", log="{log_file}", stat={dict_stat_text})'
    return object_text


def filer_length_launcher(p_args):
    return p_args[0](p_args[1], p_args[2], p_args[3], p_args[4], p_args[5])


def filter_length_all(p_list, p_lower, p_upper, p_summary_path, p_order_number,
                      p_modify_header=True, p_stat=True, p_pooled=True, p_out_seq_type='FASTQ'):
    """
    Assembled Read 의 Length Filter, Length Filter STAT 파일 생성(CSV), FASTQ Header 변경, Pooled Sample File 생성을
    전달받은 목록에 모두 적용한다.

    :type p_list: list
    :param p_list: 경로를 포함한 Assembled Read 파일. 시료명.extendedFrags.fastq[.gz|.bz|.xz].

    :type p_lower: int or
    :param p_lower: 최소 길이.

    :type p_upper: int or False
    :param p_upper: 최대 길이.

    :type p_summary_path: str
    :param p_summary_path: Summary file 생성 경로.

    :type p_order_number: str
    :param p_order_number: 수주번호.

    :type p_modify_header: bool
    :param p_modify_header: Header 변경 여부. ex) >시료명_seqNum --> >Sample1_1

    :type p_stat: bool
    :param p_stat : F_extendedFrags.fastq.gz 파일에 대한 STAT 파일 생성 여부.

    :type p_pooled: bool
    :param p_pooled: 여러 개의 FASTQ 파일들을 하나로 만듬.

    :type p_out_seq_type: str
    :param p_out_seq_type: Filtered & Pooled Sample 파일의 형식을 지정. FASTQ or FASTA

    :return: pooled_file - p_pooled = True 이면, 경로를 포함한 파일명. 아닐 경우, None
    """
    global POOL_WORKER
    filter_condition = style('{min} <= Good Sequence <= {max}'.format(min=p_lower, max=p_upper), fg='yellow')
    echo('>>> Length Filter 시작: {}'.format(filter_condition))
    if p_modify_header:
        echo('+++ Header 변경 적용')
    else:
        secho('--- Header 변경 미적용', fg='yellow')
    if p_pooled:
        echo('+++ Pooled Sample File 생성 적용')
    else:
        secho('--- Pooled Sample File 미생성', fg='yellow')
    if p_out_seq_type.upper() == 'FASTQ':
        secho('+++ Filtered & Pooled : FASTQ 형식 적용', fg='yellow')
    elif p_out_seq_type.upper() == 'FASTA':
        secho('+++ Filtered & Pooled : FASTA 형식 적용', fg='yellow')
    else:
        secho('--- out_seq_type 형식 지정 오류', fg='red', blink=True)

    l_cmd = list()
    for file in p_list:
        l_cmd.append((filter_length, file, p_lower, p_upper, p_modify_header, p_out_seq_type))

    process = POOL_WORKER
    if len(p_list) < POOL_WORKER:
        process = len(p_list)
    with Pool(process, initializer=start_process) as pool:
        pool_outputs = pool.map(filer_length_launcher, l_cmd)
        pool.close()
        pool.join()

    error_count = 0
    done_count = 0
    l_nt_filtered = list()
    for result in pool_outputs:
        if result == 1:
            error_count += 1
        else:
            l_nt_filtered.append(eval(result))
            done_count += 1
    done_text = style('완료: {}'.format(done_count), fg='cyan')
    error_text = style('에러: {}'.format(error_count),
                       fg='red' if (error_count != 0) else 'cyan',
                       blink=True if (error_count != 0) else False)
    echo('{msg} {done}, {error}'.format(
        msg=style('>>> Length Filter: ', fg='cyan'),
        done=done_text,
        error=error_text))

    # 직렬 코드
    # l_nt_filtered = list()
    # echo('>>> Length Filter 시작: {min} <= Sequence <= {max}'.format(min=p_lower, max=p_upper))
    # with progressbar(p_list, label='Length Filter') as bar_list:
    #     for file in bar_list:
    #         l_nt_filtered.append(filter_length(file, p_lower, p_upper))

    if (p_lower is False) and (p_upper is False) and (p_modify_header is False):
        pass
    else:
        name = [data.stat.get('name') for data in l_nt_filtered]
        total_read = [str(data.stat.get('total_read')) for data in l_nt_filtered]
        removed_read = [str(data.stat.get('removed_read')) for data in l_nt_filtered]
        lower_read = [str(data.stat.get('lower_read')) for data in l_nt_filtered]
        upper_read = [str(data.stat.get('upper_read')) for data in l_nt_filtered]
        final_read = [str(data.stat.get('final_read')) for data in l_nt_filtered]
        removed_rate = [str(data.stat.get('removed_rate')) for data in l_nt_filtered]
        lower_bp = [str(data.stat.get('lower_bp')) for data in l_nt_filtered]
        upper_bp = [str(data.stat.get('upper_bp')) for data in l_nt_filtered]
        csv_file = os.path.join(p_summary_path, p_order_number + '.F_length.csv')
        excel_file = os.path.join(p_summary_path, p_order_number + '.F_length.xlsx')
        l_log_stat = \
            [
                ('Name', name),
                ('Total_Read', total_read), ('Removed_Read', removed_read),
                ('Removed_Read(Lower)', lower_read), ('Removed_Read(Upper)', upper_read),
                ('Final_Read', final_read), ('Removed_Rate', removed_rate),
                ('lower_bp', lower_bp), ('upper_bp', upper_bp)
            ]
        save_csv_for_length_filter_log(csv_file, l_log_stat)
        save_excel_for_length_filter_log(excel_file, l_log_stat)

    if p_stat:
        filtered_file = [data.des for data in l_nt_filtered]
        soft_link = list()
        if p_out_seq_type.upper() == 'FASTA':
            secho('>>> Warning: STAT 생성 대상 변경', fg='yellow', blink=True, bold=True)
            secho('Filtered Assembled Read --> Assembled Read', fg='magenta', bold=True)
        for src in filtered_file:
            src_file = os.path.basename(src)
            if p_out_seq_type.upper() == 'FASTQ':
                if src.endswith('.extendedFrags.fastq.gz'):
                    des_file = src_file.replace('F_', '').replace('.extendedFrags.fastq.gz', '_1.fastq.gz')
                elif src.endswith('.extendedFrags.fastq'):
                    des_file = src_file.replace('F_', '').replace('.extendedFrags.fastq', '_1.fastq')
                else:
                    secho('>>> Warning: STAT 생성 불가능 파일 형식입니다.', fg='magenta', blink=True)
                    echo(src)
                    continue
                soft_link.append((src_file, des_file))
            elif p_out_seq_type.upper() == 'FASTA':
                if src.endswith('.extendedFrags.fasta.gz'):
                    new_src_file = src_file.replace('F_', '').replace('.fasta.gz', '.fastq.gz')
                    des_file = new_src_file.replace('.extendedFrags.fastq.gz', '_1.fastq.gz')
                elif src.endswith('.extendedFrags.fasta'):
                    new_src_file = src_file.replace('F_', '').replace('.fasta', '.fastq')
                    des_file = new_src_file.replace('.extendedFrags.fastq', '_1.fastq')
                else:
                    secho('>>> Warning: STAT 생성 불가능 파일 형식입니다.', fg='magenta', blink=True)
                    echo(src)
                    continue
                soft_link.append((new_src_file, des_file))

        # FastQStatGz.sh 실행시 파일명에 경로가 있으면 STAT 결과의 파일명에도 경로가 포함된다.
        if len(filtered_file) == 1:  # 리스트에 원소가 1개일 경우 처리.
            link_script_path = os.path.dirname(filtered_file[0])
            if link_script_path.endswith('Read_Assembly'):
                pass
            else:
                secho('Warning: Assembled Read가 보관된 디렉터리명이 Read_Assembly가 아닙니다.',
                      fg='magenta', blink=True)
        else:
            link_script_path = os.path.commonpath(filtered_file)
        if link_script_path == '/':
            secho('Error: Filter FASTQ 파일들의 경로가 동일하지 않습니다.', fg='red')
            echo(filtered_file)
            exit()
        link_script = 'soft_link.bash'
        make_soft_link_script(soft_link, os.path.join(*(link_script_path, link_script)))
        run = run_move_and_cmd(link_script_path, 'bash {}'.format(link_script))
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'soft_link.bash 실행 완료',
                'false_meg': 'soft_link.bash',
            }, p_exit=False)
        del run

        stat_script = 'STAT.bash'
        l_linked_files = [des for src, des in soft_link]
        make_fastq_stat_script(l_linked_files, os.path.join(link_script_path, stat_script))
        run = run_move_and_cmd(link_script_path, 'bash {}'.format(stat_script))
        if p_out_seq_type.upper() == 'FASTQ':
            true_meg = 'Filtered FASTQ STAT 생성 완료'
            false_meg = 'Filtered FASTQ STAT'
        elif p_out_seq_type.upper() == 'FASTA':
            true_meg = 'Assembled FASTQ STAT 생성 완료'
            false_meg = 'Assembled FASTQ STAT'
        check_run_cmd(
            {
                'run': run,
                'true_meg': true_meg,
                'false_meg': false_meg,
            }, p_exit=False)
        del run
        run = run_move_and_cmd(link_script_path, 'head -n 1 -q *.stat > STAT.txt')
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'STAT.txt 생성 완료',
                'false_meg': 'STAT.txt',
            }, p_exit=False)

    if p_pooled:
        if p_out_seq_type.upper() == 'FASTQ':
            pooled_file = os.path.join(p_summary_path, p_order_number + '.pooled.fastq')
        elif p_out_seq_type.upper() == 'FASTA':
            pooled_file = os.path.join(p_summary_path, p_order_number + '.pooled.fasta')
        sample_pool_cmd = get_sample_pool_cmd([nt_filtered.des for nt_filtered in l_nt_filtered], pooled_file)
        run = run_cmd(sample_pool_cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'Pooled Sample File 생성 완료',
                'false_meg': 'Pooled Sample File'
            }, p_exit=True)
        echo(pooled_file)
        return pooled_file
    else:
        secho('!!! Pooled Sample 미적용 !!!', fg='yellow')
        return None


def save_csv_for_length_filter_log(p_csv, p_data):
    """
    Legnth Filter에 대한 STAT를 CSV 파일 형태로 저장한다.

    :type p_csv: str
    :param p_csv:
    :type p_data: list
    :param p_data: 튜플을 원소로 가지는 리스트
    :return:
    """
    with open(p_csv, 'w') as o_csv:
        for header, data in p_data:
            o_csv.write(f'{header},')
            o_csv.write(','.join(data))
            o_csv.write('\n')
    echo('>>> Length Filter Log : csv 생성 완료')
    echo(p_csv)


def save_excel_for_length_filter_log(p_excel, p_data):
    import xlsxwriter
    workbook = xlsxwriter.Workbook(p_excel)
    worksheet = workbook.add_worksheet()
    percent_format = workbook.add_format({'num_format': '0.00%'})
    number_format = workbook.add_format({'num_format': '#,##0'})

    # make header dictionary
    l_header = list()
    header, data = p_data[0]
    l_header.append({'header': header})
    for header in data:
        del_index = header.find('.extendedFrags')
        l_header.append({'header': header[:del_index]})
    del header, data

    # Set column widths
    worksheet.set_column(0, 0, 19)

    # Write data
    for row, (header, data) in enumerate(p_data[1:], 1):
        worksheet.write_string(row, 0, header)
        for col, cell in enumerate(data, 1):
            if row == 6:
                worksheet.write_number(row, col, float(cell), percent_format)
            else:
                typeWrite(worksheet, row, col, cell, number_format)

    # Add table
    worksheet.add_table(0, 0, row, len(data),
                        {
                            'columns': l_header,
                            'autofilter': False,
                            'style': 'Table Style Light 8',
                        })

    # Read Count Chart
    read_count_chart = workbook.add_chart({'type': 'column'})
    read_count_chart.add_series(
        {
            'categories': ['Sheet1', 0, 1, 0, len(data)],
            'values': ['Sheet1', 1, 1, 1, len(data)],
            'name': ['Sheet1', 1, 0]
        }
    )
    read_count_chart.add_series(
        {
            'categories': ['Sheet1', 0, 1, 0, len(data)],
            'values': ['Sheet1', 5, 1, 5, len(data)],
            'name': ['Sheet1', 5, 0]
        }
    )
    read_count_chart.set_x_axis({'name': 'Sample Name'})
    read_count_chart.set_y_axis({'name': 'Read Count'})
    read_count_chart.set_title({'name': 'Length Filtering: Read Count'})

    # Removed Rate Chart
    removed_rate_chart = workbook.add_chart({'type': 'column'})
    removed_rate_chart.add_series(
        {
            'categories': ['Sheet1', 0, 1, 0, len(data)],
            'values': ['Sheet1', 6, 1, 6, len(data)],
            'name': ['Sheet1', 6, 0]
        }
    )
    removed_rate_chart.set_x_axis({'name': 'Sample Name'})
    removed_rate_chart.set_y_axis({'name': 'Removed Percentage'})
    removed_rate_chart.set_title({'name': 'Length Filtering: Removed Read'})

    # Set Chart Style
    read_count_chart.set_style(34)
    removed_rate_chart.set_style(34)
    if len(data) >= 20:
        x_scale = 2
    else:
        x_scale = 1.5
    read_count_chart.set_size({'x_scale': x_scale})
    removed_rate_chart.set_size({'x_scale': x_scale})

    worksheet.insert_chart('A11', read_count_chart)
    worksheet.insert_chart('A27', removed_rate_chart)

    workbook.close()
    echo('>>> Length Filter Log : Excel 생성 완료')
    echo(p_excel)


def make_soft_link_script(p_list, p_script):
    """
    입력받은 목록들에 대해 특정 파일명으로 소프트 링크를 걸 수 있는 bash 스크립트를 생성한다.

    :type p_list: list
    :param p_list: 원소가 튜플인 리스트. ex) [(src1, des1), (src2, des2)]
    :param p_script: 생성될 스크립트의 파일(경로 포함).
    :return: None
    """
    with open(p_script, 'w') as script:
        for src, des in p_list:
            cmd = 'ln -s {} {}\n'.format(src, des)
            script.write(cmd)
    secho('>>> soft link script 생성 완료', fg='cyan')
    echo(p_script)


def make_fastq_stat_script(p_list, p_script):
    """
    Assembled FASTQ 파일들에 대한 STAT 정보를 생성하는 bash 스크립트를 생성한다.

    :type p_list: list
    :param p_list: 파일명을 원소로 가지는 리스트(경로 포함시 시료명에 파일경로가 포함됨. 경로 제거)
    :type p_script: str
    :param p_script: bash 스크립트 파일명.
    :return:
    """
    with open(p_script, 'w') as script:
        for name in p_list:
            cmd = get_fastq_stat_cmd(name)
            if cmd is None:
                continue
            script.write(cmd)
            script.write('\n')
    secho('>>> fastq stat script 생성 완료', fg='cyan')
    echo(p_script)
    

def get_fastq_stat_cmd(p_name):
    """
    FASTQ 파일에 대한 STAT 정보를 생성하는 명령어를 반환한다.

    :param p_name: 파일명(경로 포함시 시료명에 파일경로가 포함됨. 경로 제거.)
    :rtype: str or None
    :return: cmd
    """
    if p_name.endswith('_1.fastq.gz'):
        name = p_name.replace('_1.fastq.gz', '')
        cmd = 'FastQStatGz.sh se gz {name} {name}.stat'.format(name=name)
    elif p_name.endswith('_1.fastq'):
        name = p_name.replace('_1.fastq', '')
        cmd = 'FastQStatGz.sh se no {name} {name}.stat'.format(name=name)
    else:
        echo('FASTQ STAT를 생성할 수 없는 파일 형식입니다.')
        return None
    return cmd


def get_sample_pool_cmd(p_list, p_out_file):
    """
    여러 시료들의 FASTQ파일들을 하나의 파일로 통합한다.

    :type p_list: list
    :param p_list: FASTQ 파일들.

    :param p_out_file: 출력파일명.

    :return: pool_cmd - 실행 명령어.
    """
    compressed_type = set([x.split('.')[-1] for x in p_list])
    if len(compressed_type) == 1:
        type_text = list(compressed_type)[0]
        if type_text == 'gz':
            cat_mode = 'zcat'
        elif type_text == 'bz':
            cat_mode = 'bzcat'
        elif type_text == 'xz':
            cat_mode = 'xzcat '
        elif type_text == 'fastq' or 'fasta':
            cat_mode = 'cat'
        else:
            raise RuntimeError(style('Error: 압축 형식(gz,bz,xz,fastq,fasta) --> {}'.format(compressed_type), fg='red'))
    elif len(compressed_type) > 1:
        raise RuntimeError(
            style('Error: 압축형식이 다른 파일이 존재합니다.\n'
                  '검출된 압축 형식: {}'.format(compressed_type), fg='red'))
    else:
        raise RuntimeError(
            style('Error: 예기치 않은 오류입니다(압축 형식 검출 문제 등).\n'
                  '검출된 압축 형식: {format}'
                  '파일 목록: {file}'.format(format=compressed_type, file=p_list), fg='red'))
    pool_cmd = '{cmd} {file} > {pooled_file}'.format(cmd=cat_mode, file=' '.join(p_list), pooled_file=p_out_file)
    return pool_cmd


def exist_fastq(p_l_nt_run_list, p_kargs):
    """

    :param p_l_nt_run_list:
    :param p_kargs:
                r1_suffix:
                r2_suffix:
                adatapter_trim_tool
    :return:
    """
    echo('>>> Read1 & Read2 FASTQ 파일 존재 여부 확인')
    l_error_status = list()
    for run in p_l_nt_run_list:
        l_no_exists = list()
        l_exits = list()
        for sample in run.samples:
            r1_suffix = p_kargs['r1_suffix'].format(trim=p_kargs['adapter_trim_tool'])
            r2_suffix = p_kargs['r2_suffix'].format(trim=p_kargs['adapter_trim_tool'])
            r1_target = os.path.join(*(sample.path, sample.name + r1_suffix))
            r2_target = os.path.join(*(sample.path, sample.name + r2_suffix))
            for target in r1_target, r2_target:
                if check_file_type(target, 'exists') is False:
                    l_no_exists.append(target)
                else:
                    l_exits.append(target)
        done_text = style('완료: {}'.format(len(l_exits)), fg='cyan')
        error_text = style('에러: {}'.format(len(l_no_exists)),
                           fg='red' if len(l_no_exists) != 0 else 'cyan',
                           blink=True if len(l_no_exists) != 0 else False)
        echo('확인 디렉터리: {dir1}/{dir2} --> {done}, {error}'.format(
            dir1=os.path.basename(os.path.dirname(run.path)),
            dir2=style(run.run_info, fg='yellow'),
            done=done_text,
            error=error_text))
        if len(l_no_exists) != 0:
            l_error_status.append(True)
            echo('------ 미확인 파일 목록 ------')
            _ = [echo(name) for name in l_no_exists]
    if any(l_error_status):
        exit(1)


def exist_something(p_path, p_list, p_suffix, p_meg):
    """

    :type p_path: str
    :param p_path: 경로
    :type p_list: list
    :param p_list: 검출 대상 목록
    :type p_suffix: str
    :param p_suffix: 검출 대상을 한정하기 위한 접미사. 와일드 카드 사용 가능.
    :param p_meg: 안내 메시지.
    :return:
    """
    echo(p_meg)
    echo(f'검출 접미사: {p_suffix}')
    l_error_status = list()
    l_no_exists = list()
    l_exits = list()
    for sample in p_list:
        target = os.path.join(p_path, sample + p_suffix)
        if check_file_type(target, 'exists') is False:
            l_no_exists.append(target)
        else:
            l_exits.append(target)
    done_text = style('완료: {}'.format(len(l_exits)), fg='cyan')
    error_text = style('에러: {}'.format(len(l_no_exists)),
                       fg='red' if len(l_no_exists) != 0 else 'cyan',
                       blink=True if len(l_no_exists) != 0 else False)
    echo('확인 디렉터리: {msg} \n\t--> {done}, {error}'.format(
        msg=style(p_path, fg='yellow'),
        done=done_text,
        error=error_text))
    if len(l_no_exists) != 0:
        l_error_status.append(True)
        echo('------ 미확인 파일 목록 ------')
        _ = [echo(name) for name in l_no_exists]
    if any(l_error_status):
        exit(1)


def exist_assembled_fastq(*args):
    path, list_data, suffix = args
    return exist_something(path, list_data, suffix, '>>> Assembled FASTQ(A) 파일 존재 여부 확인')


def length_trim(p_o_infile, p_o_outfile, p_base):
    for num, line in enumerate(p_o_infile, 1):
        num_index = num % 4
        if num_index == 1:  # Header
            p_o_outfile.write(line)
        elif num_index == 2:  # Sequence
            p_o_outfile.write(line.strip()[:p_base] + '\n')
        elif num_index == 3:  # Quality Indicator(+)
            p_o_outfile.write(line)
        elif num_index == 0:  # Quality Score
            p_o_outfile.write(line.strip()[:p_base] + '\n')
    p_o_infile.close()
    p_o_outfile.close()


def change_header_to_name(p_src, p_dec, p_name=None):
    # TODO 정리 필요
    if p_name is None:
        name = os.split(p_src)[1].remove('.fastq', '')
    else:
        name = p_name
