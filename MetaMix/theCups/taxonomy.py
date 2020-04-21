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

from SpoON.util import run_cmd, check_run_cmd, glob_dir, sort_data_by_custom, parse_config, run_move_and_cmd
from bs4 import BeautifulSoup
import os
from click import secho, echo


def get_summarize_taxa_cmd(p_biom, p_out_dir, p_level, p_absolute=False):
    """

    :param p_biom:
    :param p_out_dir:
    :type  p_level: str
    :param p_level: Taxonomy Level. ex) 2,3,4,5,6,7
    :param p_absolute: abundance 표현 방식.
                    True: absolute abundance.
                    False: relative abundance.
    :return: cmd
    """
    if p_absolute is True:
        absolute_abundance = '-a'
    else:
        absolute_abundance = ''
    cmd = \
        'summarize_taxa.py ' \
        '-i {biom} ' \
        '-o {out_dir} ' \
        '-L {level} ' \
        '--suppress_biom_table_output ' \
        '{a_a}'.format(
            biom=p_biom,
            out_dir=p_out_dir,
            level=p_level,
            a_a=absolute_abundance,
        )
    return cmd


def run_summarize_taxa(p_biom, p_out_dir, p_level, p_absolute):
    if p_absolute:
        method = 'Absolute'
    else:
        method = 'Relative'
    cmd = get_summarize_taxa_cmd(p_biom, p_out_dir, p_level, p_absolute)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': f'Summarize Taxa - {method} 완료',
            'false_meg': f'summarize_taxa - {method}',
        }, p_exit=False, p_stdout=False,
    )
    if run.returncode == 0:
        return True
    else:
        return False


def get_plot_taxa_summary_cmd(p_table, p_level, p_out_dir, p_y_size, p_chart):
    """

    :type p_table: str
    :param p_table:
    :type p_level: str
    :param p_level: lineage level name. ex) Phylum,Class,Order,Family,Genus,Species
    :type p_out_dir: str
    :param p_out_dir: 출력파일명
    :type p_y_size: int
    :param p_y_size: y축 크기
    :type p_chart: str
    :param p_chart: bar, area, pie. ex) 'bar,area'
    :return:
    """
    cmd = \
        'plot_taxa_summary.py ' \
        '-i {otu_table} ' \
        '-l {level} ' \
        '-c {chart} ' \
        '-o {out_dir} ' \
        '-y {size}'.format(
            otu_table=p_table,
            level=p_level,
            chart=p_chart,
            out_dir=p_out_dir,
            size=p_y_size,
        )
    return cmd


def run_plot_taxa_summary(p_table, p_level, p_out_file, p_y_size, p_chart):
    """

    :param p_table:
    :param p_level:
    :param p_out_file:
    :param p_y_size:
    :param p_chart: bar, area, pie. ex) 'bar,area'. 시료의 개수가 1개이면 생성 불가.
    :return:
    """
    cmd = get_plot_taxa_summary_cmd(p_table, p_level, p_out_file, p_y_size, p_chart)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': f'Plot Taxa summary: Size {p_y_size} 완료',
            'false_meg': f'plot_taxa_summary: Size {p_y_size}',
        }, p_exit=False, p_stdout=False,
    )
    if run.returncode == 0:
        return True
    else:
        return False


class TaxonHTML:
    # -----------------------------------------------------------------------------------------------------------------------
    # Ver 1.0
    # QIIME의 plot_taxa_summary.py에 의해서 생성되는 bar & area chart의 html를 보고서 양식으로 변경.
    # 그래프를 레벨별로 접근할 수 있도록 탭기능을 추가.
    # 탭기능은 Bootstrap Front-End Framework를 사용하여 구현.
    # Amplicon MetaGenome 파이프라인의 Modify_html_for_Taxonomy_assignment_Species.py 를 대체함.
    # -----------------------------------------------------------------------------------------------------------------------
    # 참고
    # QIIME CSS와 보고서 CSS(default.css)에서 border-collapse의 설정이 다름.
    # 보고서의 CSS(default.css)파일의 table에 대한 서식 지정을 주석처리.
    # -----------------------------------------------------------------------------------------------------------------------

    def __init__(self, p_ori_html, p_new_html):
        self.ori_html = p_ori_html
        self.new_html = p_new_html
        self.html_template = '''\
<html>
<head>
    <meta charset="utf-8">
    <link rel="stylesheet" href="../../src/bootstrap.css">
    <link rel="stylesheet" href="../../src/default.css">
    <link rel="stylesheet" href="./css/qiime_style.css" type="text/css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
    <script type="text/javascript" src="./js/overlib.js"></script>

    <script type="text/javascript">
        <!-- Begin
        function set_target(new_target)
        {
            sf = document.getElementById("search_form");
            sf.target = new_target;
        }
        function gg(targetq)
        {
                window.open("http://www.google.com/search?q=" + targetq, 'searchwin');
        }

        //  End -->
    </script>
    <title>Taxa Summaries</title>
</head>

<body>
    <div id="wrapper">
        <div id="header">
            <div class="title-sec">
                <div class="title-inner-sec">
                    <a href="http://www.macrogen.com/" target="_blank"><img src="../../src/img/logo.png" alt="MACROGEN" /></a>
                    <h1>NGS Analysis Report - OTU</h1>
                </div>
            </div>
        </div><!-- // header END -->
        <div class="clear"></div>
        <div id="container" class="paddingT20">
            <h2 class="title">Taxonomic assignment</h2>
            <li>The following chart shows the taxonomic composition for each sample from phylum to genus levels (x-axis: sample name;y-axis: OTU proportions).<br/>
            You can mouse over the plot to see which taxa are contributing to the percentage shown</li>
        </div><!-- // container END -->
        <div id="container" class="paddingT20">
            <div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;"></div>
            <br/><br/><br/>
            <table width=971>
                <tr>
                    <td>
                        <ul class="nav nav-tabs nav-pills nav-justified">
                            <li class="active"><a data-toggle="pill" href="#phylum">Phylum</a></li>
                            <li><a data-toggle="pill" href="#Mclass">Class</a></li>
                            <li><a data-toggle="pill" href="#order">Order</a></li>
                            <li><a data-toggle="pill" href="#family">Family</a></li>
                            <li><a data-toggle="pill" href="#genus">Genus</a></li>
                            <li><a data-toggle="pill" href="#species">Species</a></li>
                        </ul>
                    </td>
                </tr>
            </table>

            <div class="tab-content">
                <div id="phylum" class="tab-pane fade in active">
                </div><!-- phylum tab -->

                <div id="Mclass" class="tab-pane fade">
                </div><!-- class tab -->

                <div id="order" class="tab-pane fade">
                </div><!-- order tab -->

                <div id="family" class="tab-pane fade">
                </div><!-- family tab -->

                <div id="genus" class="tab-pane fade">
                </div><!-- genus tab -->

                <div id="species" class="tab-pane fade">
                </div><!-- species tab -->
            </div> <!-- class="tab-content" -->
        </div><!-- container paddingT20 END -->
        <div class="clear"></div>
        <div style="height:100px"></div><!-- space-->
    </div><!-- wrapper END -->
</body>
</html>'''

        with open(self.ori_html, 'r') as o_ori_html:
            ori_html_data = o_ori_html.read()
            self.ori_html_data = ori_html_data.replace('onmouseout="return nd();">',
                                                       'onmouseout="return nd();"></AREA>')
        self.soup = BeautifulSoup(self.ori_html_data, 'html.parser')
        self.new_soup = BeautifulSoup(self.html_template, 'html.parser')
        self.level_list = ['phylum', 'Mclass', 'order', 'family', 'genus', 'species']

    @staticmethod
    def _make_level_title(p_level):
        level_list = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        tag_text = '<tr><td class="title-mid" colspan=2 style="background-color:#eff5f8">' \
                   '<font size=5>{0} Level</font></td></tr>'.format(level_list[p_level])
        soup = BeautifulSoup(tag_text, 'html.parser')
        return soup

    def transform(self):
        body = self.soup.body
        table_list = body.find_all('table')
        map_list = body.find_all('map')

        level_soup = []
        for i in self.level_list:
            level_soup.append(self.new_soup.find('div', id=i))

        for i in range(len(map_list)):
            table_1 = table_list.pop(0)  # Phylum 제목 있는 테이블
            table_2 = table_list.pop(0)  # 'View Table (.txt)' 테이블
            table_3 = table_list.pop(0)  # Taxonomy 결과 테이블

            table_1.tr.extract()  # 붚필요한 테이블의 행 제거
            table_1.tr.extract()
            level_title = self._make_level_title(i)
            table_1.insert(3, level_title)

            level_soup[i].append(table_1)  # Phylum 제목 있는 테이블
            level_soup[i].append(map_list[i])  # <MAP> ... </MAP>
            level_soup[i].append(table_2)  # View Table (.txt)' 테이블
            level_soup[i].append(table_3)  # Taxonomy 결과 테이블

    def save(self):
        with open(self.new_html, "w") as output:
            output.write(str(self.new_soup))

    def del_html(self):
        cmd = f'rm {self.ori_html}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_exit=False,
        )

    def move_new_html(self):
        cmd = f'mv {self.new_html} {self.ori_html}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_exit=False,
        )


class TaxonomyAbundanceExcel:
    l_taxonomy_level_name = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    def __init__(self, kargs):
        """
        otu_table.L?.txt 파일을 읽어 Excel 파일로 변환한다.

        :param kargs:
                    db_type: [UNITE | Other]
                    level: [Genus | Species]
                    otu_table_dir_path: otu_table.L[2~7].txt
                    data_type: [count | ratio]
                    sort_type: [True | False | True_number | Custom]
                    custom_order: list()
        """
        self.db_type = kargs['db_type']
        self.level = kargs['level']
        self.otu_table_dir_path = kargs['otu_table_dir_path']
        self.data_type = kargs['data_type']
        self.sort_type = kargs['sort_type']
        self.l_custom_order = kargs['custom_order']
        if self.level == 'Genus':
            self.level_count = 7
        elif self.level == 'Species':
            self.level_count = 8
        else:
            raise ValueError(f'self.level: {self.level}. 지원하지 않는 Level입니다. 개발자에게 문의하세요.')
        self.l_otu_table_file = None
        self.l_level_count = list()
        self.l_read_file_data = list()

    def glob_otu_table_file(self):
        l_otu_table = glob_dir(self.otu_table_dir_path, 'otu_table.*.txt', p_mode='many')
        if not l_otu_table:
            raise RuntimeError(f'l_otu_table: {l_otu_table}. otu_table.*.txt 파일이 없습니다.')
        l_otu_table.sort()
        self.l_otu_table_file = l_otu_table

    @staticmethod
    def find_max(target1, target2):
        temp = []
        if not target1:
            return target2
        else:
            for index in range(len(target1)):
                if target1[index] > target2[index]:
                    temp.append(target1[index])
                elif target1[index] < target2[index]:
                    temp.append(target2[index])
                else:
                    temp.append(target1[index])
        return temp

    @staticmethod
    def sort_data(p_data, p_sample_start, p_sorting, p_custom_order=None):
        sorted_data = [[0 for i in p_data[0]] for j in p_data]
        taxon_data = p_data[0][:p_sample_start]
        sample_data = p_data[0][p_sample_start:]
        if p_sorting == "True_number":
            try:
                sample_data_numbered = [int(x) for x in sample_data]
            except ValueError:  # not int
                try:
                    sample_data_numbered = [float(x) for x in sample_data]
                except ValueError:
                    secho('Error: TypeError - The option of sorting is wrong!', fg='red')
                    exit(1)
            sample_data_numbered.sort()
            sample_data = [str(x) for x in sample_data_numbered]
        elif p_sorting == "True":
            sample_data.sort()
        elif p_sorting == 'Custom':
            if ','.join(sorted(sample_data)) == ','.join(sorted(p_custom_order)):
                sample_data = p_custom_order
            else:
                l_not_in_metadata = list()
                for ele in sample_data:
                    if ele in p_custom_order:
                        continue
                    else:
                        l_not_in_metadata.append(ele)
                else:
                    if l_not_in_metadata:
                        secho('Error: metadata.txt 파일에 기재되지 않은 시료들이 있습니다.', fg='red', blink=True)
                        secho(f'\t--> {l_not_in_metadata}', fg='magenta')
                        return False

        full_sample_data = taxon_data + sample_data
        for i in range(len(full_sample_data)):
            target = full_sample_data[i]
            sample_index = p_data[0].index(target)
            for j in range(len(p_data)):
                sorted_data[j][i] = p_data[j][sample_index]
        return sorted_data

    def read_otu_table(self):
        global l_taxonomy_level_name
        for level, table_file in enumerate(self.l_otu_table_file, 2):
            if f'L{level}' in table_file:
                pass
            else:
                raise RuntimeError('파일의 Level과 추측된 Level이 다릅니다. 개발자에게 문의하세요.')
            read_data_list = []
            self.l_level_count.append(level)
            with open(table_file, "r") as o_table:
                _ = o_table.readline()  # '# Constructed from biom file' 제거
                header = o_table.readline().rstrip().split("\t")
                final_header = self.l_taxonomy_level_name[:level] + header[1:]
                read_data_list.append(final_header)
                for element in o_table:
                    # ex) element -> 'Unassigned;Other\t626.0\t4101.0\t884.0\t2417.0\t861.0\t698.0\n'
                    if self.db_type == "UNITE":
                        element = element.replace("k__", "__")
                        element = element.replace("p__", "__")
                        element = element.replace("c__", "__")
                        element = element.replace("o__", "__")
                        element = element.replace("f__", "__")
                        element = element.replace("g__", "__")
                        element = element.replace("s__", "__")
                    else:
                        pass
                    # ex) element_split -> ['Unassigned', 'Other\t626.0\t4101.0\t884.0\t2417.0\t861.0\t698.0\n']
                    element_split = element.strip().split(";")
                    # ['Other', '626.0', '4101.0', '884.0', '2417.0', '861.0', '698.0']
                    element_split_2 = element_split[-1].split("\t")
                    final_element = element_split[:-1] + element_split_2
                    read_data_list.append(final_element)

                if (self.sort_type is True) or (self.sort_type == "True_number"):
                    self.l_read_file_data.append(self.sort_data(read_data_list, level, self.sort_type))
                elif self.sort_type == 'Custom':
                    self.l_read_file_data.append(self.sort_data(read_data_list, level,
                                                                self.sort_type, self.l_custom_order))
                elif self.sort_type is False:
                    self.l_read_file_data.append(read_data_list)
                else:
                    raise ValueError(f'self.sort_type: {self.sort_type}. 지원하지 않는 값입니다.')

    def save_excel(self):
        import xlsxwriter
        excel_file = os.path.join(self.otu_table_dir_path, f'Taxonomy_abundance_{self.data_type}.xlsx')
        workbook = xlsxwriter.Workbook(excel_file)
        worksheet_list = []

        if self.data_type == "count":
            data_format = workbook.add_format({'num_format': '#,##0'})
        elif self.data_type == "ratio":
            data_format = workbook.add_format({'num_format': '0.00%'})

        # Add worksheet
        for i in range(len(self.l_level_count)):
            # exclude Kingdom --> i+1
            worksheet_list.append(workbook.add_worksheet(self.l_taxonomy_level_name[i + 1]))

        for index in range(len(self.l_level_count)):
            column_size = len(self.l_read_file_data[index][0])
            target_data = self.l_read_file_data[index][1:]

            # Find long length by column
            column_large_length_list = []
            for element in target_data:
                length_element = list(map(len, element))
                column_large_length_list = self.find_max(column_large_length_list, length_element)
            # Set column widths
            for col_num, width in enumerate(column_large_length_list):
                worksheet_list[index].set_column(col_num, col_num, width)

            # make header dictionary
            header_list = []
            for i in self.l_read_file_data[index][0]:
                header_list.append({"header": i})

            # Add table
            worksheet_list[index].add_table(0, 0, len(target_data), column_size - 1,
                                            {'columns': header_list, 'autofilter': True,
                                             'style': 'Table Style Light 8'})

            # Write data
            for row_num, row_data in enumerate(self.l_read_file_data[index][1:], 1):
                for col_num, cell in enumerate(row_data):
                    try:
                        out_cell = float(cell)
                        worksheet_list[index].write_number(row_num, col_num, out_cell, data_format)
                    except ValueError:
                        out_cell = str(cell)
                        worksheet_list[index].write_string(row_num, col_num, out_cell)

        workbook.close()


class TaxonomySharedExcel:
    def __init__(self, kargs):
        """
        
        :param kargs:
                    db: [SILVA_111 | SILVA_123 | UNITE | NCBI_16S | NCBI_NT | RDP]
                    sort: [True | False | True_number | Custom]
                    custom_order: list()
                    method: [DENOVO | CLOSED]
                    shared_file:
                    tax_txt_file:
                    tax_log_file:
                    out_path:
        """
        if kargs['level'] == 'Genus':
            self.taxonomy_index = 7
        elif kargs['level'] == 'Species':
            self.taxonomy_index = 8

        self.db = kargs['db']
        self.sort = kargs['sort']
        self.custom_order = kargs['custom_order']
        self.otu_method = kargs['method']
        self.shared_file = kargs['shared_file']
        self.tax_txt_file = kargs['tax_txt_file']
        self.tax_log_file = kargs['tax_log_file']
        self.out_path = kargs['out_path']
        self.l_shared_data = None
        self.d_assignment_txt = None
        self.d_assignment_log = None

    def read_shared_file(self):
        with open(self.shared_file, 'r') as o_file:
            data = []
            for i in o_file:
                data.append(i.split())
            if self.sort is True:
                header = data[0]
                temp = data[1:]
                temp.sort(key=lambda x: x[1])
                data = [header] + temp
                del temp
                del header
            elif self.sort == 'True_number':
                header = data[0]
                temp = data[1:]
                for i in temp:
                    try:
                        i[1] = int(i[1])
                    except ValueError:
                        try:
                            i[1] = float(i[1])
                        except ValueError as err:
                            echo(err.message)
                            secho('Error: The option of sorting is wrong! Shared Excel 파일 생성 중단.',
                                  fg='red', blink=True)
                            secho('\t--> TaxonomySharedExcel', fg='red')
                            return False

                temp.sort(key=lambda x: x[1])
                for i in temp:
                    try:
                        i[1] = str(i[1])
                    except ValueError as err:
                        echo(err.message)
                        secho('Error: ValueError 발생. Shared Excel 파일 생성 중단.', fg='red', blink=True)
                        secho('\t--> TaxonomySharedExcel', fg='red')
                        return False
                data = [header] + temp
            elif self.sort == 'Custom':
                header = data[0]
                temp = data[1:]
                sorted_temp = sort_data_by_custom(temp, self.custom_order, p_key=1)
                data = [header] + sorted_temp
            elif self.sort is False:
                pass
        self.l_shared_data = data

    def read_tax_assignment_txt(self):
        with open(self.tax_txt_file, 'r') as o_tax_txt:
            assignment_dic = {}
            for i in o_tax_txt:
                i_split = i.rstrip().split('\t')
                taxonomy_data = i_split[:-2]
                score_infor = i_split[-2:]  # ['0.67', '3']
                taxon_split = [x.strip() for x in i_split[1].split(";")]
                # 'k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Atheliales; f__Atheliaceae'
                # --> ['k__Fungi', 'p__Basidiomycota', 'c__Agaricomycetes', 'o__Atheliales', 'f__Atheliaceae']
                taxonomy_data.extend(taxon_split)
                if len(taxonomy_data) > self.taxonomy_index + 1:
                    assignment_dic[i_split[0]] = taxonomy_data[:(self.taxonomy_index + 1 - len(taxonomy_data))] \
                                                 + score_infor
                elif len(taxonomy_data) == self.taxonomy_index + 1:
                    assignment_dic[i_split[0]] = taxonomy_data + score_infor
                else:
                    num = 1
                    add_data_number = self.taxonomy_index + 1 - len(taxonomy_data)
                    while num <= add_data_number:
                        taxonomy_data.append('-')
                        num += 1

                    assignment_dic[i_split[0]] = taxonomy_data + score_infor
        self.d_assignment_txt = assignment_dic

    def read_tax_assignment_log(self):
        with open(self.tax_log_file, 'r') as assignment_log_data:
            assignment_log_dic = {}
            echo('>>> otu_table_shared_with_tax_assignment.xlsx 데이터 정리 기준')
            if self.otu_method == 'DENOVO':
                secho('OTU Picking - DENOVO 방식에 대해서 데이터를 정리합니다.', fg='yellow')
            elif self.otu_method == 'CLOSED':
                secho('OTU Picking - CLOSED 방식에 대해서 데이터를 정리합니다.', fg='yellow')
            for i in assignment_log_data:
                if self.otu_method == 'DENOVO':
                    if 'denovo' in i:
                        i_split = i.split()
                        dic_value = assignment_log_dic.get(i_split[8])
                        if dic_value is None:
                            assignment_log_dic[i_split[8]] = [i_split[10]]
                        else:
                            assignment_log_dic[i_split[8]].append(i_split[10])
                    else:
                        continue
                if self.otu_method == 'CLOSED':
                    if i.startswith('H	'):
                        i_split = i.split()
                        dic_value = assignment_log_dic.get(i_split[8])
                        if dic_value is None:
                            assignment_log_dic[i_split[8]] = [i_split[10]]
                        else:
                            assignment_log_dic[i_split[8]].append(i_split[10])
                    else:
                        continue
        self.d_assignment_log = assignment_log_dic

    def save_excel(self):
        import xlsxwriter
        row = len(self.l_shared_data)
        column = len(self.l_shared_data[0])

        transpose_matrix = [[0 for x in range(row)] for y in range(column)]

        # Transpose
        for row_index in range(row):
            for column_index in range(column):
                transpose_matrix[column_index][row_index] = self.l_shared_data[row_index][column_index]

        l_taxonomy_level = ['Organism', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

        for i in range(self.taxonomy_index):
            transpose_matrix[1].insert(i + 1, l_taxonomy_level[i])
            transpose_matrix[2].insert(i + 1, l_taxonomy_level[i])

        # Insert taxonomy assignment information
        for index in range(len(transpose_matrix)):
            if index < 3:
                continue

            get_value = self.d_assignment_txt[transpose_matrix[index][0]]
            len_get_value = len(get_value[1:])
            if len_get_value != self.taxonomy_index:
                for j in range(self.taxonomy_index - len_get_value):
                    get_value.append(" ")
            transpose_matrix[index] = get_value + transpose_matrix[index][1:]

        # Save data to Excel file
        excel_file = os.path.join(self.out_path, 'otu_table_shared_with_tax_assignment.xlsx')
        workbook = xlsxwriter.Workbook(excel_file)
        worksheet = workbook.add_worksheet()

        # Set foramt
        # link_format = workbook.add_format({'color': 'red', 'underline': 1})

        # Find the number of column and row
        matrix_row = len(transpose_matrix[3:])
        # matrix_column = len(transpose_matrix[3])

        # Add table
        header_list = []
        for i in transpose_matrix[1]:
            header_list.append({'header': i})
            if i == l_taxonomy_level[self.taxonomy_index - 1]:
                header_list.append({'header': 'Score'})
                header_list.append({'header': 'Hit'})
                header_list.append({'header': 'Accession_Number_1'})
                header_list.append({'header': 'Accession_Number_2'})
                header_list.append({'header': 'Accession_Number_3'})

        for num, i in enumerate(transpose_matrix[1:], 1):
            if (num == 1) or (num == 2):
                i.insert(self.taxonomy_index+2, 'Hit')
                i.insert(self.taxonomy_index+1, 'Score')
                i.insert(self.taxonomy_index+3, 'Accession_Number_1')
                i.insert(self.taxonomy_index+4, 'Accession_Number_2')
                i.insert(self.taxonomy_index+5, 'Accession_Number_3')
            if num >= 3:
                dic_value = self.d_assignment_log.get(i[0])
                if len(dic_value) == 3:
                    i.insert(self.taxonomy_index+3, dic_value[0])
                    i.insert(self.taxonomy_index+4, dic_value[1])
                    i.insert(self.taxonomy_index+5, dic_value[2])
                if len(dic_value) == 2:
                    i.insert(self.taxonomy_index+3, dic_value[0])
                    i.insert(self.taxonomy_index+4, dic_value[1])
                    i.insert(self.taxonomy_index+5, 'NA')
                if len(dic_value) == 1:
                    if dic_value[0] == '*':
                        i.insert(self.taxonomy_index+3, 'NA')
                    else:
                        i.insert(self.taxonomy_index+3, dic_value[0])
                    i.insert(self.taxonomy_index+4, 'NA')
                    i.insert(self.taxonomy_index+5, 'NA')

        # Find column lengths
        def find_max(target1, target2):
            temp = list()
            if not target1:
                return target2
            else:
                for index in range(len(target1)):
                    if target1[index] > target2[index]:
                        temp.append(target1[index])
                    elif target1[index] < target2[index]:
                        temp.append(target2[index])
                    else:
                        temp.append(target1[index])
            return temp

        column_large_length_list = []
        for element in transpose_matrix[1:]:
            length_element = list(map(len, element))
            column_large_length_list = find_max(column_large_length_list, length_element)

        # Set column widths
        for col_num, width in enumerate(column_large_length_list):
            worksheet.set_column(col_num, col_num, width)

        worksheet.add_table(0, 0, matrix_row, len(header_list) - 1,
                            {'columns': header_list, 'style': "Table Style Medium 15"})

        # Write data
        for row_num, row_data in enumerate(transpose_matrix[3:], 1):
            for col_num, cell in enumerate(row_data):
                if (row_num > 0) and (col_num >= self.taxonomy_index+3) and (col_num <= self.taxonomy_index+3+2):
                    # taxonomy_index+3 --> index of Accession_Number_1,
                    # taxonomy_index+3+2 --> index of Accession_Number_3
                    if cell != 'NA':
                        if self.db == "SILVA_111":
                            url = "http://www.arb-silva.de/browser/ssu-104/" + cell
                            worksheet.write_url(row_num, col_num, url, string=cell)
                        elif self.db == "SILVA_123":
                            url = "http://www.arb-silva.de/browser/ssu-104/" + cell.split(".")[0]
                            worksheet.write_url(row_num, col_num, url, string=cell)
                        elif self.db == "UNITE":
                            url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=search&term=/' + cell
                            worksheet.write_url(row_num, col_num, url, string=cell)
                        else:
                            url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=search&term=/' + cell
                            worksheet.write_url(row_num, col_num, url, string=cell)
                    else:
                        worksheet.write_string(row_num, col_num, cell)
                else:
                    try:
                        out_cell = float(cell)
                        worksheet.write_number(row_num, col_num, out_cell)
                    except ValueError:
                        try:
                            out_cell = int(cell)
                            worksheet.write_number(row_num, col_num, out_cell)
                        except ValueError:
                            out_cell = str(cell)
                            worksheet.write_string(row_num, col_num, out_cell)
        workbook.close()


def run_make_shared(p_run_path, p_biom):
    """
    biom 파일을 복사하여, mothur make.shared() 명령어를 수행한다.
    shared 파일이 생성된다.

    * mothur 실행시 주의사항.
      - mothur를 실행하는 위치에 로그파일이 생성된다.
      - biom 파일 위치에 shared 파일이 생성된다.


    :type p_run_path: str
    :param p_run_path: 작업을 진행할 경로.
    :type p_biom: str
    :param p_biom: biom 파일
    :return:
    """
    config = parse_config()
    mothur = config['theCups']['mothur']

    copy_cmd = f'cp {p_biom} {p_run_path}'
    copy_run = run_cmd(copy_cmd)
    check_run_cmd(
        {
            'run': copy_run,
            'true_meg': None,
            'false_meg': 'copy biom',
        }, p_stdout=False, p_exit=False,
    )
    if copy_run.returncode == 0:
        pass
    else:
        return False

    biom_file = os.path.join(p_run_path, os.path.basename(p_biom))
    mothur_cmd = f'{mothur} "#make.shared(biom={biom_file})"'
    _ = run_move_and_cmd(p_run_path, mothur_cmd)
    l_globed_file = glob_dir(p_run_path, 'mothur*.logfile', p_mode='many')
    if len(l_globed_file) == 1:
        mothur_log = l_globed_file[0]
    else:
        l_globed_file.sort()
        mothur_log = l_globed_file[-1]
    with open(mothur_log, 'r') as o_log:
        for i in o_log:
            if i.startswith('[ERROR]'):
                echo(i)
                return False
        else:
            return True
