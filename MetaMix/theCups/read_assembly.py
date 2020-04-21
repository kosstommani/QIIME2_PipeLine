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
__version__ = '1.0.0'

import os
from SpoON.util import parse_config, sort_data_by_custom


def read_assembly_stat(p_file):
    """
    MiSeq_V1
    Read Assembly --> Length Filter --> STAT.txt 생성
    
    :param p_file: 
    :return: 
    """
    stat = list()
    with open(p_file, 'r') as o_stat:
        for i in o_stat:
            if i.startswith('==>') or i.startswith('\n'):
                continue
            else:
                sample_stat = i.strip().split()
                sample_stat[1] = int(sample_stat[1])    # Base Count
                sample_stat[2] = int(sample_stat[2])    # Read Count
                sample_stat[3] = float(sample_stat[3])  # N Bases %
                sample_stat[4] = float(sample_stat[4])  # GC %
                sample_stat[5] = float(sample_stat[5])  # Q20 %
                sample_stat[6] = float(sample_stat[6])  # Q30 %
                stat.append(sample_stat)
    return stat


class StatExcel:
    def __init__(self, kargs):
        """
        MiSeq_V1
        Read Assembly --> Length Filter --> STAT.txt 생성

        :param kargs:
                    stat_file:
                    metadata_order:
                    out_path:
        """
        self.stat_file = kargs['stat_file']
        self.metadata_order = kargs['metadata_order']
        self.out_path = kargs['out_path']
        self.l_stat = None

    def read_assembly_stat(self):
        stat = read_assembly_stat(self.stat_file)
        sorted_stat = sort_data_by_custom(stat, self.metadata_order)
        self.l_stat = sorted_stat

    def make_excel(self):
        import xlsxwriter
        config = parse_config()
        stat_name = os.path.split(self.stat_file)[-1]
        excel_name = f'{os.path.splitext(stat_name)[0]}.xlsm'
        excel_file = os.path.join(self.out_path, excel_name)
        workbook = xlsxwriter.Workbook(excel_file)
        worksheet = workbook.add_worksheet()
        workbook.add_vba_project(config['Excel_VBA']['vba_stat'])
        worksheet.insert_button(
            'I1',
            {
                'macro': 'Chart2GIF',
                'caption': 'Save chart to GIF',
                'width': 100,
                'height': 25,
            })

        # Set the cell format
        first_row_rich_format1 = workbook.add_format({
            'font_size': 22,
            'bold': True,
            'font_name': 'Arial',
        })

        first_row_rich_format2 = workbook.add_format({
            'font_size': 10,
            'bold': True,
            'font_name': 'Arial',
        })

        first_row_rich_cell_format = workbook.add_format({
            'border': 6,
            'align': 'center',
            'valign': 'vcenter',
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

        count_data_format = workbook.add_format({
            'font_name': 'Arial',
            'font_size': 9,
            'valign': 'vcenter',
            'num_format': '#,##0',
        })

        # Write the data and table
        worksheet.merge_range("A1:G1", "", first_row_rich_cell_format)
        worksheet.write_rich_string("A1",
                                    first_row_rich_format1, "Result of Assembly",
                                    first_row_rich_format2, " (by FLASH)",
                                    first_row_rich_cell_format)

        worksheet.merge_range("A3:G3", "FLASH", second_row_format)
        worksheet.write_string('K1', 'Set PATH:')
        worksheet.write_string('L1', 'C:\\')

        table_header = ['SampleName', 'Total Bases', 'Read Count', 'N (%)', 'GC (%)', 'Q20 (%)', 'Q30 (%)']
        start_table_row = 3
        worksheet.add_table(start_table_row, 0, start_table_row + len(self.l_stat), len(self.l_stat[0]) - 1,
                            {'columns':
                                [
                                    {'header': table_header[0]},
                                    {'header': table_header[1]},
                                    {'header': table_header[2]},
                                    {'header': table_header[3]},
                                    {'header': table_header[4]},
                                    {'header': table_header[5]},
                                    {'header': table_header[6]},
                                ],
                                'autofilter': False,
                                'style': 'Table Style Light 16',
                            })

        # Find long length sample
        sample_name_header_width = table_header[0]
        for i in self.l_stat:
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

        for row_num, row in enumerate(self.l_stat, start_table_row + 1):
            for col_num, cell in enumerate(row):
                if col_num == 0:
                    out_cell = str(cell)
                    worksheet.write_string(row_num, col_num, out_cell, common_format)
                elif (col_num == 1) or (col_num == 2):
                    out_cell = int(cell)
                    worksheet.write_number(row_num, col_num, out_cell, count_data_format)
                else:
                    try:
                        out_cell = float(cell)
                    except ValueError:  # n/a value
                        out_cell = cell
                    worksheet.write(row_num, col_num, out_cell, common_format)

        # Set conditional format
        # len(self.l_stat) + start_table_row + 1 --> the number of sample + start_table_row + header(1)
        worksheet.conditional_format(4, 2, len(self.l_stat) + start_table_row + 1, 2,
                                     {'type': 'data_bar', 'bar_color': '#FFC000'})  # Read Count
        worksheet.conditional_format(4, 5, len(self.l_stat) + start_table_row + 1, 5,
                                     {'type': 'data_bar', 'bar_color': '#9ACD32'})  # Q20
        worksheet.conditional_format(4, 6, len(self.l_stat) + start_table_row + 1, 6,
                                     {'type': 'data_bar', 'bar_color': '#DC1687'})  # Q30

        # Read Count - column chart
        read_count_chart = workbook.add_chart({'type': 'column'})
        read_count_chart.add_series({
            'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(self.l_stat), 0],
            'values': ['Sheet1', start_table_row + 1, 2, start_table_row + len(self.l_stat), 2],
            'name': ['Sheet1', start_table_row, 2],
        })

        read_count_chart.set_x_axis({'name': 'sample name'})
        read_count_chart.set_y_axis({'name': 'read'})

        # Quality - Q20, Q30 - column chart
        quality_chart = workbook.add_chart({'type': 'column'})
        quality_chart.add_series({
            'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(self.l_stat), 0],
            'values': ['Sheet1', start_table_row + 1, 5, start_table_row + len(self.l_stat), 5],
            'name': ['Sheet1', start_table_row, 5],
        })

        quality_chart.add_series({
            'categories': ['Sheet1', start_table_row + 1, 0, start_table_row + len(self.l_stat), 0],
            'values': ['Sheet1', start_table_row + 1, 6, start_table_row + len(self.l_stat), 6],
            'name': ['Sheet1', start_table_row, 6],
            'overlap': 50,
        })
        quality_chart.set_x_axis({'name': 'sample name'})
        quality_chart.set_y_axis({'name': '%', 'min': 0, 'max': 100})
        quality_chart.set_title({'name': 'Quality'})

        # Chart Size - using the number of sample
        correction_value = 1 + (0.08 * len(self.l_stat) / 10)
        read_count_chart.set_size({'x_scale': correction_value, 'y_scale': correction_value})
        quality_chart.set_size({'x_scale': correction_value, 'y_scale': correction_value})

        # Chart Style
        read_count_chart.set_style(34)
        quality_chart.set_style(34)

        # Insert OTU chart and community diviersity chart
        worksheet.insert_chart('I3', read_count_chart)
        worksheet.insert_chart('I' + str(int(20 * correction_value)),
                               quality_chart)  # position correction according to the chart size
        workbook.close()


def transpose_data(p_l_data):
    row_count = len(p_l_data)
    col_count = len(p_l_data[0])
    transposed_data = [[0 for x in range(row_count)] for y in range(col_count)]
    for row in range(row_count):
        for col in range(col_count):
            transposed_data[col][row] = p_l_data[row][col]
    return transposed_data
