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
__version__ = '1.0.1'

import os
from SpoON.util import run_cmd, check_run_cmd, parse_config, read_metadata_for_sort, sort_data_by_custom


def get_alpha_diversity_cmd(p_biom, p_out, p_tre):
    # PD_whole_tree 매트릭스 사용시 -t 옵션 추가
    if p_tre is None:
        option = '-m observed_species,chao1,shannon,simpson,goods_coverage '
    else:
        option = f'-m observed_species,PD_whole_tree,chao1,shannon,simpson,goods_coverage -tre {p_tre}'
    cmd = \
        'alpha_diversity.py ' \
        f'-i {p_biom} ' \
        f'{option} ' \
        f'-o {p_out}'
    return cmd


def read_adiv(p_file):
    adiv = list()
    with open(p_file, 'r') as o_adiv:
        for i in o_adiv:
            if 'observed_species' in i:
                continue
            else:
                data = i.strip().split()
                data[1] = float(data[1])  # observed_species
                data[2] = float(data[2])  # chao1
                data[3] = float(data[3])  # shannon
                data[4] = float(data[4])  # simpson (Inverse)
                data[5] = float(data[5])  # goods_coverage
                adiv.append(data)
    return adiv


def run_alpha_diversity(p_biom, p_outpath, p_tre, p_metadata):
    output_file = os.path.join(p_outpath, 'adiv.txt')
    cmd = get_alpha_diversity_cmd(p_biom, output_file, p_tre)
    run = run_cmd(cmd)
    check_run_cmd({
        'run': run,
        'true_meg': 'Alpha Diversity Index 완료',
        'false_meg': 'alpha_diversity',
    }, p_exit=False, p_stdout=False)
    if run.returncode == 0:
        adiv = read_adiv(output_file)
        metadata_order = read_metadata_for_sort(p_metadata)
        sorted_adiv = sort_data_by_custom(adiv, metadata_order)
        adiv_excel_file = os.path.join(p_outpath, 'adiv.xlsm')
        make_excel(sorted_adiv, adiv_excel_file)
        return sorted_adiv
    else:
        return False


def make_excel(p_adiv, p_out_file):
    from copy import deepcopy
    l_adiv = deepcopy(p_adiv)
    import xlsxwriter
    config = parse_config()
    workbook = xlsxwriter.Workbook(p_out_file)
    worksheet = workbook.add_worksheet()
    workbook.add_vba_project(config['Excel_VBA']['vba_adiv'])
    worksheet.insert_button(
        'H1',
        {
            'macro': 'Chart2GIF',
            'caption': 'Save chart to GIF',
            'width': 100,
            'height': 25,
        })

    # setting format
    first_row_format = workbook.add_format({
        'bold': True,
        'border': 6,
        'align': 'center',
        'valign': 'vcenter',
        'font_size': 22,
        'font_name': 'Arial',
    })

    second_row_format = workbook.add_format({
        'font_size': 10,
        'font_name': 'Arial',
        'bg_color': '4F81BD',
        'font_color': '#FFFFFF',
        'align': 'center',
        'border': 1,
        'border_color': '#FFFFFF',
        'bold': True,
    })

    common_format = workbook.add_format({
        'font_name': 'Arial',
        'font_size': 9,
        'valign': 'vcenter',
    })

    table_header_format = workbook.add_format({
        'font_name': 'Arial',
        'font_size': 10,
        'align': 'center',
        'valign': 'vcenter',
    })

    # writing data and table
    worksheet.merge_range("A1:F1", "Community richness & diversity", first_row_format)
    worksheet.merge_range("A3:F3", "", second_row_format)
    worksheet.write_string('J1', 'Set PATH:')
    worksheet.write_string('K1', 'C:\\')
    start_table_row = 3
    table_header = ['SampleName', 'OTUs', 'Chao1', 'Shannon', 'Inverse Simpson', 'Good\'s Coverage']
    worksheet.add_table(start_table_row, 0, start_table_row + len(l_adiv), len(l_adiv[0]) - 1,
                        {
                            'columns':
                                [
                                    {'header': table_header[0]},
                                    {'header': table_header[1]},
                                    {'header': table_header[2]},
                                    {'header': table_header[3]},
                                    {'header': table_header[4]},
                                    {'header': table_header[5]},
                                ],
                            'autofilter': False,
                            'style': 'Table Style Light 16',
                        })

    # Find long length sample
    sample_name_header_width = table_header[0]
    for i in l_adiv:
        if len(i[0]) > len(sample_name_header_width):
            sample_name_header_width = i[0]
        else:
            continue

    # Set column widths
    table_header_copy = table_header[:]
    table_header_copy[0] = sample_name_header_width
    for i in range(len(table_header_copy)):
        worksheet.set_column(i, i, len(table_header_copy[i]) * 1.2)
        # print "%d: %s %f" % (i, table_header_copy[i], len(table_header_copy[i])*1.2)

    worksheet.set_row(3, 16.5, table_header_format)

    # 시료의 개수가 1개일 경우, Diversity 꺾은선 그래프 생성을 위해 데이터 추가.
    if len(l_adiv) == 1:
        l_adiv.append(l_adiv[0])

    for row_num, row in enumerate(l_adiv, start_table_row + 1):
        for col_num, cell in enumerate(row):
            if col_num == 0:
                out_cell = str(cell)
                worksheet.write_string(row_num, col_num, out_cell, common_format)
            else:
                try:
                    out_cell = float(cell)
                except ValueError:  # n/a value
                    out_cell = cell
                worksheet.write(row_num, col_num, out_cell, common_format)

    # OTU - column chart
    otu_chart = workbook.add_chart({'type': 'column'})
    otu_chart.add_series({
        'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(l_adiv), 0],
        'values': ['Sheet1', start_table_row + 1, 1, start_table_row + len(l_adiv), 1],
        'name': ['Sheet1', start_table_row, 1],
    })

    otu_chart.set_x_axis({'name': 'sample name'})
    otu_chart.set_y_axis({'name': 'count'})
    otu_chart.set_legend({'position': 'top'})

    # community diversity(shannon, inverse simpson)
    community_diversity_chart = workbook.add_chart({'type': 'line'})
    community_diversity_chart.add_series({
        'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(l_adiv), 0],
        'values': ['Sheet1', start_table_row + 1, 3, start_table_row + len(l_adiv), 3],
        'name': ['Sheet1', start_table_row, 3],
    })

    community_diversity_chart.add_series({
        'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(l_adiv), 0],
        'values': ['Sheet1', start_table_row + 1, 4, start_table_row + len(l_adiv), 4],
        'name': ['Sheet1', start_table_row, 4],
        'y2_axis': True,
    })
    community_diversity_chart.set_x_axis({'name': 'sample name'})
    community_diversity_chart.set_y_axis({
        'name': 'shannon',
        'num_font': {'color': 'blue'}
    })
    community_diversity_chart.set_y2_axis({
        'name': 'inverse simpson',
        'num_font': {'color': 'red'}
    })
    community_diversity_chart.set_title({'name': 'Community Diversity'})
    community_diversity_chart.set_legend({'position': 'top'})

    # Chart Size - using the number of sample
    correction_value = 1 + (0.08 * len(l_adiv) / 10)
    otu_chart.set_size({'x_scale': correction_value, 'y_scale': correction_value})
    community_diversity_chart.set_size({'x_scale': correction_value, 'y_scale': correction_value})

    # Chart Style
    otu_chart.set_style(34)
    community_diversity_chart.set_style(34)

    # inserting OTU chart and community diversity chart
    worksheet.insert_chart('H3', otu_chart)
    worksheet.insert_chart('H' + str(int(20 * correction_value) + 1),
                           community_diversity_chart)  # position correction according to the chart size
    workbook.close()
