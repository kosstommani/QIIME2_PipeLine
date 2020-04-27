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
import xlsxwriter
from xlsxwriter.utility import xl_cell_to_rowcol
from pprint import pprint
from click import echo, secho, style
from theCups.report import Report as theCups_Report
from theCups.data_structure import ReportPath
from SpoON.util import check_file_type, parse_config

CONFIG = parse_config()
PROBIOTICS_CONFIG = CONFIG['theCups']['Probiotics']
SRC_PROBIOTICS_PATH = PROBIOTICS_CONFIG['src_probiotics_path']


class SampleTable:
    def __init__(self, p_sample_name):
        """
        시료명과 L7 결과를 클래스 변수로 가진다.
        
        :param p_sample_name: 시료명
        """
        self.sample_name = p_sample_name
        self.d_L7_table = dict()


class Data:
    def __init__(self, p_kargs):
        """
        고객이 섞은 유산균 종 정보를 입력받고, OTU분석결과를 파싱하여 보고서 생성에 필요한 데이터를 생성한다.

        :param p_kargs: 프로그램 옵션 값을 가지는 딕션너리
        """
        global PROBIOTICS_CONFIG, SRC_PROBIOTICS_PATH

        self.analysis_base_path = p_kargs['analysis_base_path']
        self.order_number = p_kargs['order_number']
        self.table_path_suffix = p_kargs['table_path_suffix']
        self.v1v2_dir = p_kargs['v1v2_dir']
        self.v1v2_num = p_kargs['v1v2_num']
        self.v3v4_dir = p_kargs['v3v4_dir']
        self.v3v4_num = p_kargs['v3v4_num']
        self.mixed_species_file = p_kargs['mixed_species']

        # Species List 파일
        self.species_file = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['species_list'])

        # 메서드 실행시 생성되는 클래스 변수
        self.info_file = None
        self.info = None
        self.mixed_species = None  # 고객이 첨가한 유산균
        self.d_mixed_species_from_file = None  # 고객이 첨가한 유산균 목록 (text file)
        self.species_list = None  # 순서가 매겨진 유산균 종 목록
        self.d_species_list = None
        self.d_species_list_index = None  # Key값 : species index. 식약처 19종 유산균 목록
        self.d_sample_product_name = None
        self.v1v2_table_file = None
        self.v3v4_table_file = None
        self.l_v1v2_sample_table_list = None
        self.l_v3v4_sample_table_list = None

        # Table_L7 경로 설정
        if (self.v1v2_dir is not None) and (self.v1v2_num is not None):
            self.v1v2_L7 = os.path.join(self.analysis_base_path, self.order_number,
                                        self.v1v2_dir, self.v1v2_num,
                                        'OTU_Results', self.table_path_suffix)
        else:
            self.v1v2_L7 = None
        if (self.v3v4_dir is not None) and (self.v3v4_num is not None):
            self.v3v4_L7 = os.path.join(self.analysis_base_path, self.order_number,
                                        self.v3v4_dir, self.v3v4_num,
                                        'OTU_Results', self.table_path_suffix)
        else:
            self.v3v4_L7 = None

    def check_info_file(self):
        for report_dir in [self.v1v2_dir, self.v3v4_dir]:
            analysis_dir = report_dir.rstrip('_Report')
            info_file = os.path.join(self.analysis_base_path, self.order_number, analysis_dir, 'info.txt')
            info_status, info_name = theCups_Report.exists(info_file)
            if info_status:
                self.info_file = info_file
                break

    def read_info_file(self):
        """
        info.txt 파일을 읽어 데이터를 리스트형태로 반환한다.

        -- info.txt 구조(yaml) --
        organization: 기관명
        name: 고객명
        order_number: 수주번호
        sample_count: 시료의 개수
        target_region: 목표영역
        e_mail: 본사 이메일
        
        # 이전 양식
        -- info.txt의 구조 --
        주영엔에스(주)         : 기관명
        이신영                : 고객명
        1706KMI-0081         : 수주번호
        1 (V1-V2)            : 시료의 개수
        ngskr@macrogen.com   : 본사 이메일

        값 할당 클래스 변수
        self.info : info.txt의 행별 문자열을 원소로 가지는 리스트

        :return:
        """
        from theCups.report import Report
        self.info = Report.parse_yaml(self.info_file)
        # 이전 형식
        # self.info = list()
        # with open(self.info_file, 'r') as o_input:
        #     for i in o_input:
        #         self.info.append(i.strip())
        secho('>>> info.txt 읽기 완료', fg='blue')
        echo(self.info_file)
        echo()

    def check_table_file(self):
        """
        역영별(V1V2, V3V4) table_L7 파일의 존재여부를 확인
        화면에 파일의 확인 여부 출력
        파일이 없는 경우 Error 메시지 출력 후, 프로그램 종료

        값 할당 클래스 변수
        self.v1v2_table_file
        self.v3v4_table_file

        :return:
        """
        if self.v1v2_L7:
            if check_file_type(self.v1v2_L7, 'exists'):
                secho('>>> V1V2 : table_L7 파일 확인', fg='blue')
                echo(self.v1v2_L7)
                echo()
            else:
                secho('Error : V1V2 table_L7 파일 확인 불가', fg='red', blink=True)
                echo(self.v1v2_L7)
        if self.v3v4_L7:
            if check_file_type(self.v3v4_L7, 'exists'):
                secho('>>> V3V4 : table_L7 파일 확인', fg='blue')
                echo(self.v3v4_L7)
                echo()
            else:
                secho('Error : V3V4 table_L7 파일 확인 불가', fg='red', blink=True)

    def _parse_table_data(self, p_table_file):
        """
        역영별(V1V2, V3V4) table_L7 파일을 읽은 후, 파싱한다.
        분석된(table_L7 파일 파싱) 종이 식약처에서 공시한 19종 유산균의 목록에 일치하지 않을 경우
        subspecies가 있는 아래의 유산균 종인지 확인한다.
            Lactobacillus delbrueckii ssp. bulgaricus
            Bifidobacterium animalis ssp. lactis
         probiotics DB의 경우 subspecies에 해당하는 서열들만 정리하였으므로, subspecies 일치 여부 확인 생략한다.

        :param p_table_file:
        :return: l_data_list : SampleTable 클래스의 인스턴스를 원소로 가지는 리스트
                    SampleTable.sample_name : 시료명
                    SampleTable.d_L7_table : 식약처 19종 유산균의 index를 Key로 가지는 딕션너리(str).
                                             value는 abundance ratio(float).
        """
        l_data_list = list()
        with open(p_table_file, 'r') as o_input:
            for i in o_input:
                if i.startswith('SampleID') or i.startswith('#OTU ID'):
                    for sample_name in i.strip().split('\t')[1:]:  # SampleID 제외
                        l_data_list.append(SampleTable(sample_name))
                    del sample_name
                elif i.startswith('# Constructed'):
                    continue
                else:
                    temp = i.strip().split('\t')
                    taxon = temp[0].split(';__')[-1]  # taxon의 species level
                    ratio_list = temp[1:]
                    species_index = self.d_species_list.get(taxon)

                    # Species 목록에 없는 종일 경우
                    # subspecies가 있는 경우  Lactobacillus delbrueckii ssp. bulgaricus,
                    #                       Bifidobacterium animalis ssp. lactis
                    # probiotics DB의 경우 subspecies에 해당하는 서열들만 정리하였음. 따라서 subspecies 일치 여부 확인 생략함.
                    if (species_index is None) and (taxon == 'Lactobacillus delbrueckii'):
                        species_index = self.d_species_list.get('Lactobacillus delbrueckii ssp. bulgaricus')
                    elif (species_index is None) and (taxon == 'Bifidobacterium animalis'):
                        species_index = self.d_species_list.get('Bifidobacterium animalis ssp. lactis')
                    elif species_index is None:
                        secho('Warning : {taxon} 식약처 종 목록에 없음!'.format(taxon=taxon), fg='red')
                        echo('--> {taxon_text}'.format(taxon_text=temp))
                        continue

                    for index, ratio in enumerate(ratio_list):
                        if float(ratio) > 0:
                            l_data_list[index].d_L7_table[species_index] = float(ratio)
                        elif float(ratio) == 0:
                            continue
                        else:
                            raise ValueError
                    del temp, taxon, ratio_list, index, ratio
        return l_data_list

    def read_table_file(self):  # check_mixed_species 또는 read_mixed_species 후 실행
        """
        self._parse_table_data() 메서드를 이용하여 table_L7파일을 파싱한 후, 클래스 변수에 파싱한 결과를 할당한다.
        파싱된 결과는 SampleTable 클래스의 인스턴스를 원소로 가지는 리스트 형태이다.

        값 할당 클래스 변수
        self.l_v1v2_sample_table_list
        self.l_v3v4_sample_table_list

        :return:
        """
        if self.d_species_list is None:
            raise SyntaxError(
                style('read_talbe_file 메서드 --> check_mixed_species 또는 read_mixed_species_list 메서드 실행 후 사용!',
                      fg='red', bg='yellow', blink=True))
        else:
            if self.v1v2_L7:
                secho(">>> V1V2 Table", fg="blue")
                self.l_v1v2_sample_table_list = self._parse_table_data(self.v1v2_L7)
            if self.v3v4_L7:
                secho("\n>>> V3V4 Table", fg="blue")
                self.l_v3v4_sample_table_list = self._parse_table_data(self.v3v4_L7)

    def input_product_name(self):
        echo()
        secho('유산균 제품명을 입력합니다.(입력완료 후, 수정 및 확인 가능)', fg='blue')
        self.d_sample_product_name = dict()
        if self.v1v2_L7:
            target_data_object_list = self.l_v1v2_sample_table_list
            target_region = 'v1v2'
        elif self.v3v4_L7:
            target_data_object_list = self.l_v3v4_sample_table_list
            target_region = 'v3v4'
        else:
            raise ValueError(style('v1v2 또는 v3v4 table list가 없습니다.', fg='red', blank=True))

        for o_table in target_data_object_list:
            sample_name = check_sample_name(o_table.sample_name, target_region)
            product_name = input('입력 > {sample} : '.format(sample=sample_name))
            if product_name == '':  # 제품명 미입력시
                product_name = '미기재'
            self.d_sample_product_name[sample_name] = product_name  #.decode('utf-8').encode('utf-8')

        answer = True
        while (answer != 'N') and (answer != 'n'):
            d_product_name_for_print = {}
            echo()
            for index, item in enumerate(self.d_sample_product_name.items(), 1):
                echo('{index}. {sample} : {product}'.format(index=index, sample=item[0], product=item[1]))
                d_product_name_for_print[index] = item
            answer = input('변경할 시료의 번호(선택완료: N or n) : ')
            if (answer != 'N') and (answer != 'n'):
                try:
                    selected_sample_name = d_product_name_for_print[int(answer)][0]
                except ValueError:  # invalid literal for int() with base 10 발생시 --> 한글입력
                    continue
                new_product_name = input(
                    '입력 > {index}. {sample} : '.format(index=answer, sample=selected_sample_name))
                if new_product_name == '':  # 제품명 미입력시
                    new_product_name = '미기재'
                self.d_sample_product_name[selected_sample_name] = new_product_name.decode('utf-8').encode('utf-8')

    def read_species_list(self):
        """
        19종 유산균 리스트는 SpeciesList.txt 파일을 파싱한다.

        값 할당 클래스 변수
        self.species_list         : 19종 유산균의 목록 리스트. ex) '1. Lactobacillus acidophilus'
        self.d_species_list       : 19종 유산균의 목록을 가지는 딕셔너리. key(str): 종명(str).
                                    ex) 'Lactobacillus acidophilus': '1'
        self.d_species_list_index : 19종 유산균의 목록 가지는 딕셔너리. key(str): index(str).
                                    ex) '1': 'Lactobacillus acidophilus'
                                    
        이전 경로
        /lustre/Tools/macrogen_analysis_toolbox2/OTU_customize/Probiotics/src/SpeciesList.txt

        :return:
        """
        self.species_list = list()
        self.d_species_list = dict()
        self.d_species_list_index = dict()
        with open(self.species_file, 'r') as o_input:
            for i in o_input:
                temp = i.strip().split('\t')
                self.species_list.append('%s. %s' % (temp[0], temp[1]))
                self.d_species_list[temp[1]] = temp[0]  # Lactobacillus acidophilus: 1
                self.d_species_list_index[temp[0]] = temp[1]  # 1: Lactobacillus acidophilus
        del temp, i

    def check_mixed_species(self):
        """
        19종 유산균 리스트를 화면에 출력하여, 고객이 섞은 유산균 종 정보를 키보드로 입력받는다.
        19종 유산균 리스트는 /lustre/Tools/macrogen_analysis_toolbox2/OTU_customize/Probiotics/src/SpeciesList.txt 파일을
        파싱하여 출력한다.

        값 할당 클래스 변수
        self.mixed_species : 숫자를 원소로 가지는 리스트. 숫자는 19종 유산균 목록의 index(int).

        Returns
        -------
        None
        """
        answer = True
        answer_list = list()
        while (answer != 'N') and (answer != 'n'):
            for index, element in enumerate(self.species_list, 1):
                if index in answer_list:
                    secho(element, fg='yellow')
                else:
                    echo(element)
            del index, element

            answer = input(style('투입된 종(species) 선택(선택완료: N or n, 제거: -) ex) 1 2 5-10 -2: ', fg='blue'))
            if (answer != 'N') and (answer != 'n'):
                if ',' in answer:
                    echo('"," 사용 금지!')
                    answer = 'N'
                    continue
                else:
                    temp = answer.split()
                    for i in temp:
                        if i.startswith('-'):  # 선택항목 제거일 경우
                            answer_list.remove(int(i.lstrip('-')))
                        elif '-' in i:  # 범위지정형일 경우
                            i_range = i.split('-')  # ex) '1-10'
                            answer_list.extend(range(int(i_range[0]), int(i_range[1]) + 1))  # 뒷범위는 포함X. +1 적용
                        else:  # 숫자형 문자일 경우
                            try:
                                i_num = int(i)
                                if i_num > len(self.species_list):
                                    continue
                                else:
                                    answer_list.append(int(i))
                            except ValueError:
                                continue
        self.mixed_species = list(set(answer_list))
        secho('>>> 선택되어진 종(Species)', fg='blue')
        echo(self.mixed_species)
        echo()

    def read_mixed_species_list(self):
        """
        고객이 섞은 종 목록을 정리한 파일을 파싱한다.
        목록파일에서 '#'문자로 시작하면 해당 열은 제외한다.
        제품명이 없을 경우 'None'으로 표기한다.
        시로명 입력시 v12, V1V2, v34, v3v4 등 목표 영역을 나타내는 접미사는 제거하고 입력한다.
        종 목록 파일의 구조는 다음과 같다.

        #시료명 제품명 투입종
        A1  유산균제품A  1   5   10
        B1  유산균제품B  2   5   7
        A2  None  1   5   7   9   15
        #B2  유산균제품D  4   6   12  15  19

        값 할당 클래스 변수
        self.d_mixed_species_from_file : 시료별로 유산균 목록을 가지는 딕션너리.
                                         Key: 시료명, Value: ('시료(제품)명', '1', '2', '3')

        Returns
        -------
        None
        """
        self.d_mixed_species_from_file = dict()
        self.d_sample_product_name = dict()
        secho('>>> 시료명 : 제품명', fg='blue')
        with open(self.mixed_species_file) as o_input:
            for i in o_input:
                if i.startswith('#'):
                    continue
                else:
                    temp = i.strip().split('\t')
                    sample_name = temp[0]
                    product_name = temp[1]
                    species_index = temp[2:]
                    if (product_name == 'None') or (product_name == ''):
                        self.d_sample_product_name[sample_name] = '미기재'
                    else:
                        self.d_sample_product_name[sample_name] = product_name
                    self.d_mixed_species_from_file[sample_name] = tuple([int(x) for x in species_index])
                    echo('{sample} : {product}'.format(sample=sample_name, product=product_name))
        echo()
        secho('>>> 투입된 종 목록 완료', fg='blue')
        echo(self.mixed_species_file)
        echo()


class Report:
    def __init__(self, kargs):
        """

        :param kargs: 다음을 Key로 가지는 딕션너리
                 o_v1v2_sample_table:
                 o_v3v4_sample_table:
                 product_name:
                 analysis_base_path:
                 order_number:
                 info:
                 d_species_list_index:
                 mixed_species:
        """
        global PROBIOTICS_CONFIG, SRC_PROBIOTICS_PATH

        self.analysis_base_path = kargs['analysis_base_path']
        self.order_number = kargs['order_number']
        self.v1v2_sample_table = kargs['o_v1v2_sample_table']
        self.v3v4_sample_table = kargs['o_v3v4_sample_table']
        self.product_name = kargs['product_name']
        if (self.v1v2_sample_table is not None) and (self.v3v4_sample_table is not None):
            self.region = 'both'
            self.sample_name = check_sample_name(self.v1v2_sample_table.sample_name, 'v1v2')
            if self.sample_name[-1] == '.':
                self.sample_name = self.sample_name[:-1]
        elif (self.v1v2_sample_table is not None) and (self.v3v4_sample_table is None):
            self.region = 'v1v2'
            self.sample_name = self.v1v2_sample_table.sample_name
        elif (self.v1v2_sample_table is None) and (self.v3v4_sample_table is not None):
            self.region = 'v3v4'
            self.sample_name = self.v3v4_sample_table.sample_name
        else:
            raise ValueError(
                style('Report 생성자의 v1v2 또는 v3v4의 sample_table 매개변수가 잘못되었습니다.', fg='red', blink=True)
                + '\n' + 'p_o_v1v2_sample_table : %s' % self.v1v2_sample_table
                + '\n' + 'p_o_v3v4_sample_table : %s' % self.v3v4_sample_table)

        self.report_dir_path = self.make_report_dir()
        self.workbook = xlsxwriter.Workbook(f'{self.report_dir_path}/{self.sample_name}.xlsx')
        self.workbook.set_properties({'title': 'Probiotics Product Analysis Report',
                                      'author': 'META분석팀',
                                      'manager': 'META분석팀',
                                      'company': '마크로젠(macrogen)'})

        self.info = kargs['info']
        self.d_species_list_index = kargs['d_species_list_index']
        self.mixed_species = kargs['mixed_species']

        # 그림 파일
        self.macrogen_logo_white = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['macrogen_logo_white'])
        self.process_jpg = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['process_jpg'])
        self.tree_png = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['tree_png'])
        self.cover_logo_jpg = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['cover_logo_jpg'])

        # 유산균 목록 파일
        self.reference_file = os.path.join(SRC_PROBIOTICS_PATH, PROBIOTICS_CONFIG['reference_list'])

        # 공통 서식 설정
        self.font_size = 10
        self.default_border = 1
        self.main_title_format = self.workbook.add_format(
            {
                'font_size': self.font_size + 14,
                'bold': True,
                'align': 'center',
                'valign': 'vcenter'
            })  # Info Sheet
        self.sheet_header_format = self.workbook.add_format(
            {
                'font_size': self.font_size + 4,
                'font_color': 'white',
                'bold': True,
                'align': 'right',
                'valign': 'vcenter'
            })  # Info Sheet
        self.title_format = self.workbook.add_format(
            {
                'font_size': self.font_size + 4,
                'border': self.default_border,
                'bg_color': '#f2f2f2',
                'bold': True,
                'align': 'center',
                'valign': 'vcenter'
            })  # Other Sheet

        # 메소드 실행시 생성되는 변수들

        self.report_data = None
        self.non_detection_counter = None
        self.outer_counter = None
        self.info_sheet = None
        self.probiotics_type_sheet = None
        self.taxonomy_sheet = None
        self.reference_sheet = None

    def make_report_dir(self):
        probiotics_name = 'Probiotics'
        i_report_path = ReportPath(self.analysis_base_path, self.order_number, probiotics_name)
        theCups_Report.make_dir(i_report_path.report_base_path)
        i_report_path.check_report_path_and_set(None)
        theCups_Report.make_dir(i_report_path.report_dir_path)
        return i_report_path.report_dir_path

    def _set_sheet_layout(self, p_o_sheet, p_last_col, p_format):
        """
        Excel Sheet에 대한 형태를 설정한다.
        설정 내용 : 한 페이지 맞춤, 페이지 가운데 맞춤, 가이드라인 제거, 행높이, 열넓이, 용지 설정

        Parameters
        ----------
        p_o_sheet: Excel Sheet Object
        p_last_col: 데이터가 입력되어지는 마지막 열명(문자)
        p_format: 첫행에 입력되어지는 보고서 제목을 위한 서식

        Returns
        -------
        None
        """
        p_o_sheet.set_paper(9)  # 9 : A4
        p_o_sheet.fit_to_pages(1, 1)  # 한 페이지 시트 맞춤
        p_o_sheet.center_horizontally()  # 페이지 가운데 맞춤
        p_o_sheet.hide_gridlines(2)  # 엑셀 가이드라인 제거
        p_o_sheet.set_row(0, 50)  # 1행 행높이 설정
        p_o_sheet.set_row(2, 43)  # 3행 행높이 설정
        p_o_sheet.merge_range('A1:{last}1'.format(last=p_last_col), 'NGS Analysis Report', p_format)
        # 상단 - 마크로젠 로고 삽입
        p_o_sheet.insert_image('A1', self.macrogen_logo_white)
        # p_o_sheet.set_column('{last_next}:XFD'.format(last_next=xl_col_to_name(last_column_index+2)),
        # options={'hidden': True})

    def _make_info_sheet(self):
        def return_region_text(p_region):
            v1v2_text = 'V1-V2 [340-370bp]'
            v3v4_text = 'V3-V4 [440-460bp]'
            both_text = 'V1-V2 [340-370bp]  /  V3-V4 [440-460bp]'
            if p_region == 'v1v2':
                return v1v2_text
            elif p_region == 'v3v4':
                return v3v4_text
            elif p_region == 'both':
                return both_text
            else:
                secho('Warning: 등록된 Region 정보가 아닙니다. 개발자 문의!', fg='red')
                secho('--> 등록되어진 정보로 대체합니다.')
                return both_text

        def return_taxon_text(p_key):
            if int(p_key) < 10:
                return '0{key}. {taxon}'.format(key=p_key, taxon=self.d_species_list_index[p_key])
            else:
                return '{key}. {taxon}'.format(key=p_key, taxon=self.d_species_list_index[p_key])

        # 서식 설정
        info_color = '#e7e1dc'
        border = 1
        indent = 2
        top_format = self.workbook.add_format({'top': border})
        bottom_format = self.workbook.add_format({'bottom': border})
        normal_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'valign': 'vcenter',
                'indent': indent
            })
        info_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'bg_color': info_color,
                'align': 'center',
                'valign': 'vcenter'
            })
        bottom_info_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'bg_color': info_color,
                'align': 'center',
                'valign': 'vcenter'
            })

        # 셀 설정
        last_column_name = 'E'
        last_column_index = xl_cell_to_rowcol(last_column_name + '1')[1]
        info_header_col_index = 2  # C열
        self.info_sheet = self.workbook.add_worksheet('Info')
        self._set_sheet_layout(self.info_sheet, last_column_name, self.sheet_header_format)
        self.info_sheet.set_column('A:C', 7)  # 열넓이 설정
        self.info_sheet.set_column('D:E', 36)  # 열넓이 설정

        # (서식, text1, text2)
        l_sheet_text = [(None, None, None),
                        ('title', '유산균 분석 보고서'),
                        ('bottom', None, None),
                        ('data', 'Company / Institution', self.info['organization']),
                        ('data', 'Customer', self.info['name']),
                        ('data', 'Order Number', self.info['order_number']),
                        ('data', 'Sample(Product) Name',
                         self.sample_name + "(" + self.product_name + ")"),
                        ('top', None, None),
                        ('bottom', None, None),
                        ('data', 'Type of Sequencer', 'MiSeq System'),
                        ('data', 'Sequencing Protocol', 'MiSeq System User Guide Part # 15027617 Rev. L'),
                        ('data', 'Application', 'Amplicon Metagenome'),
                        ('data', 'Target', return_region_text(self.region)),
                        ('data', 'Analysis Type', 'OTU analysis'),
                        ('data', 'Analysis Tool', 'OTU picking by CD-HIT-OTU'),
                        ('data', '', 'Taxonomic assignment by BLAST to "유산균 제품을 위한 유산균 19종 DB"'),
                        ('top', None, None),
                        ('bottom', None, None),
                        ('start', 'Probiotics Type'),
                        ('list', '[식약청고시 제 2011-68호]'),
                        ('list', ''),  # 3(species index)
                        ('list', ''),  # 4
                        ('list', ''),  # 5
                        ('list', ''),  # 6
                        ('list', ''),  # 7
                        ('list', ''),  # 8
                        ('list', ''),  # 9
                        ('list', ''),  # 10
                        ('list', ''),  # 11
                        ('top', None, None)
                        ]

        for row, element in enumerate(l_sheet_text, 1):
            if element[0] is None:
                pass
            elif element[0] == 'title':
                self.info_sheet.merge_range(row, 0,
                                            row, last_column_index,
                                            element[1], self.main_title_format)
            elif element[0] == 'bottom':
                self.info_sheet.merge_range(row, 0, row,
                                            info_header_col_index,
                                            element[1], bottom_format)
                self.info_sheet.merge_range(row, info_header_col_index + 1,
                                            row, info_header_col_index + 2,
                                            element[2], bottom_format)
            elif element[0] == 'top':
                self.info_sheet.merge_range(row, 0,
                                            row, info_header_col_index,
                                            element[1], top_format)
                self.info_sheet.merge_range(row, info_header_col_index + 1,
                                            row, info_header_col_index + 2,
                                            element[2], top_format)
            elif element[0] == 'data':
                self.info_sheet.merge_range(row, 0,
                                            row, info_header_col_index,
                                            element[1], info_format)
                self.info_sheet.merge_range(row, info_header_col_index + 1,
                                            row, info_header_col_index + 2,
                                            element[2], normal_format)
            elif element[0] == 'start':
                self.info_sheet.merge_range(row, 0,
                                            row, info_header_col_index,
                                            element[1], info_format)
                start_list_row = row
            elif element[0] == 'list':
                self.info_sheet.merge_range(row, 0,
                                            row, info_header_col_index,
                                            element[1], info_format)
            elif element[0] == 'last':
                self.info_sheet.merge_range(row, 0,
                                            row, info_header_col_index,
                                            element[1], bottom_format)

        # Probiotics Type 목록 출력
        for index in range(1, 12):  # 1 ~ 11번 유산균
            row_index = start_list_row + index - 1  # -1 : 위치 조정
            species_key_1 = str(index)  # species index
            species_key_2 = str(index + 11)
            self.info_sheet.write_string(row_index, info_header_col_index + 1,
                                         return_taxon_text(species_key_1), normal_format)
            if int(species_key_2) <= 19:
                self.info_sheet.write_string(row_index, info_header_col_index + 2,
                                             return_taxon_text(species_key_2), normal_format)
            else:
                self.info_sheet.write_string(row_index, info_header_col_index + 2,
                                             '', normal_format)
        end_list_row = row_index + 1

        # 하단 - 분석 단계 그림 삽입를 위한 영역 설정 및 그립 삽입
        for row_index in range(end_list_row + 1, 48):
            self.info_sheet.write(row_index, 0, '', normal_format)  # 0: 'A'
            self.info_sheet.write(row_index, last_column_index, '', normal_format)
        self.info_sheet.write(48, 0, '', normal_format)
        self.info_sheet.write_row(48, 1, ['' for x in range(1, last_column_index)], normal_format)
        self.info_sheet.write(48, last_column_index, '', normal_format)
        self.info_sheet.insert_image(end_list_row + 2, 0, self.process_jpg, {'x_scale': 1.1, 'y_scale': 1})  # A33

        # 프린트 영역 설정
        self.info_sheet.print_area(0, 0, row_index + 1, last_column_index)

    def _organize_data_for_report(self):
        """
        probiotics_type 시트와 taxonomy 시트를 구성하기 위한 데이터를 생성한다.
        데이터는 self.report_data 클래스 변수에 할당이 된다.
        self.report_data는 튜플을 원소로 가지는 리스트 구조이며, 튜플은 아래와 같은 형식으로 데이터를 구성한다.

        -*- self.report_data 데이터 구조 -*-
        튜플을 원소로 가지는 리스트 [(),(),...]
        index : type        : 설명
        0     : str         : 19종 유산균에 대한 인덱스(Key)
        1     : str         : 유산균 종명
        2     : str or None : 투입 여부 - O or None
        3     : str or None : 검출 여부 - O or O1 or O2 or None
        4     : float       : V1V2 영역의 abundance ratio
        5     : float       : V3V4 영역의 abundance ratio

        ex) ('1', 'Lactobacillus acidophilus', 'O', 'O', 0.801169395479, 0.908633915456)

        값 할당 클래스 변수
        self.non_detection_counter
        self.outer_counter
        self.report_data

        Returns
        -------
        None
        """
        temp = self.d_species_list_index.keys()
        species_indexes = [int(x) for x in temp]
        species_indexes.sort()
        del temp

        self.non_detection_counter = 0
        self.outer_counter = 0
        self.report_data = list()
        for species_index in species_indexes:
            temp_data = []
            mixed_species_indicator = None
            detection_indicator = None
            v1v2_indicator = None
            v3v4_indicator = None
            temp_data.append(str(species_index))  # 유산균 인덱스. ProbioticsType Sheet
            temp_data.append(self.d_species_list_index[str(species_index)])  # 유산균 종명. ProbioticsType Sheet
            # 투입종  확인
            if species_index in self.mixed_species:
                mixed_species_indicator = True
                temp_data.append('O')  # 투입여부(C열). ProbioticsType Sheet
            else:
                mixed_species_indicator = False
                temp_data.append(None)  # 투입여부(C열). ProbioticsType Sheet

            # V1V2, V3V4 영역에 대한 검출종 확인
            if self.v1v2_sample_table is not None:
                if str(species_index) in self.v1v2_sample_table.d_L7_table.keys():
                    v1v2_indicator = True
                else:
                    v1v2_indicator = False
            if self.v3v4_sample_table is not None:
                if str(species_index) in self.v3v4_sample_table.d_L7_table.keys():
                    v3v4_indicator = True
                else:
                    v3v4_indicator = False

            if v1v2_indicator and v3v4_indicator:  # 검출
                temp_data.append('O')  # both
                temp_data.append(self.v1v2_sample_table.d_L7_table[str(species_index)])
                temp_data.append(self.v3v4_sample_table.d_L7_table[str(species_index)])
                detection_indicator = True
            elif (v1v2_indicator is True) and (v3v4_indicator is False):  # 검출
                temp_data.append('O1')  # v1v2
                temp_data.append(self.v1v2_sample_table.d_L7_table[str(species_index)])
                temp_data.append(0)
                detection_indicator = True
            elif (v1v2_indicator is True) and (v3v4_indicator is None):  # 검출
                temp_data.append('O1')  # only v1v2
                temp_data.append(self.v1v2_sample_table.d_L7_table[str(species_index)])
                temp_data.append(None)
                detection_indicator = True
            elif (v1v2_indicator is False) and (v3v4_indicator is True):  # 검출
                temp_data.append('O2')  # v3v4
                temp_data.append(0)
                temp_data.append(self.v3v4_sample_table.d_L7_table[str(species_index)])
                detection_indicator = True
            elif (v1v2_indicator is None) and (v3v4_indicator is True):  # 검출
                temp_data.append('O2')  # only v3v4
                temp_data.append(None)
                temp_data.append(self.v3v4_sample_table.d_L7_table[str(species_index)])
                detection_indicator = True
            elif (v1v2_indicator is False) and (v3v4_indicator is False):  # 미검출
                temp_data.append('X')  # both
                temp_data.append(0)
                temp_data.append(0)
                detection_indicator = False
            elif (v1v2_indicator is False) and (v3v4_indicator is None):  # 미검출
                temp_data.append('X')  # only v1v2
                temp_data.append(0)
                temp_data.append(None)
                detection_indicator = False
            elif (v1v2_indicator is None) and (v3v4_indicator is False):  # 미검출
                temp_data.append('X')  # only v3v4
                temp_data.append(None)
                temp_data.append(0)
                detection_indicator = False
            else:
                secho('Error: 영역별 종의 검출을 확인하는 중 예상치 못한 조건이 실행되었습니다.',
                      fg='red', bg='yellow', blink=True)
                secho('알고리즘 검증이 필요합니다. 개발자에게 문의하세요.', fg='yellow')
                echo(f'v1v2_indicator : {v1v2_indicator}')
                echo(f'v3v4_indicator : {v3v4_indicator}')
                exit()

            # 검출종 개수 계산. ProbioticsType Sheet
            if (mixed_species_indicator is True) and (detection_indicator is False):
                self.non_detection_counter += 1
            elif (mixed_species_indicator is False) and (detection_indicator is True):
                self.outer_counter += 1

            # 영역별 종검출 데이터 보관
            self.report_data.append(tuple(temp_data))
        echo(f'>>> 유산균 검출 정보 : {self.sample_name}')
        pprint(self.report_data)
        echo(f'불검출 항목수: {self.non_detection_counter}')
        echo(f'투입 외 검출 항목수: {self.outer_counter}')
        echo()

    def _make_probiotics_type_sheet(self):
        # 서식 설정
        bg_color_no = '#d9d9d9'
        blue = '#0070c0'
        red = '#ff0000'
        indent = 2
        top_table_title_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'top': self.default_border + 1,
                'bg_color': '#e4dfec',
                'bold': True,
                'align': 'center',
                'valign': 'vcenter'
            })
        table_title_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'bg_color': '#e4dfec',
                'bold': True,
                'align': 'center',
                'valign': 'vcenter'
            })
        legend_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'bg_color': bg_color_no,
                'align': 'center',
                'valign': 'vcenter'
            })
        left_normal_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })
        thick_top_bottom_normal_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'top': self.default_border + 1,
                'bottom': self.default_border + 1,
                'left': self.default_border,
                'right': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })
        blue_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'font_color': blue,
                'align': 'center',
                'valign': 'vcenter'
            })
        blue_bold_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'font_color': blue,
                'bold': True,
                'align': 'center',
                'valign': 'vcenter'})
        blue_bold_no_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'font_color': blue,
                'bold': True,
                'bg_color': bg_color_no,
                'align': 'center',
                'valign': 'vcenter'
            })
        red_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'font_color': red,
                'align': 'center',
                'valign': 'vcenter'
            })
        self.no_mixed_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'bg_color': bg_color_no,
                'valign': 'vcenter'})
        self.species_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'italic': True,
                'indent': indent,
                'valign': 'vcenter'
            })

        # 셀 설정
        self.probiotics_type_sheet = self.workbook.add_worksheet('Probiotics Type')
        last_column_name = 'E'  # 데이터가 쓰여지는 마지막 열의 이름.
        self._set_sheet_layout(self.probiotics_type_sheet, last_column_name, self.sheet_header_format)
        self.probiotics_type_sheet.set_column('A:A', 11)  # B열 열넓이 설정 (A~E 열사이즈 합계 : 93)
        self.probiotics_type_sheet.set_column('B:C', 21)  # B열 열넓이 설정
        self.probiotics_type_sheet.set_column('D:E', 20)  # C,D열 열넓이 설정

        # 데이터 입력
        self.probiotics_type_sheet.merge_range('A3:E3', '유산균 종류 분석결과', self.title_format)
        self.probiotics_type_sheet.merge_range('D5:E5', '목표 영역 검출 여부', legend_format)
        contents = (('D6', '표기', legend_format), ('E6', '검출 영역', legend_format),
                    ('D7', 'O', blue_format), ('E7', 'V1V2 & V3V4', blue_format),
                    ('D8', 'O1', blue_format), ('E8', 'V1V2', blue_format),
                    ('D9', 'O2', blue_format), ('E9', 'V3V4', blue_format),
                    ('D10', 'X', red_format), ('E10', '미검출', red_format),
                    ('D13', '투입여부', table_title_format),
                    ('E13', '검출여부', table_title_format),
                    )
        for element in contents:
            self.probiotics_type_sheet.write(element[0], element[1], element[2])
        self.probiotics_type_sheet.merge_range('A12:C13', '유산균', top_table_title_format)
        self.probiotics_type_sheet.merge_range('D12:E12',
                                               self.sample_name + '(' + self.product_name + ')',
                                               top_table_title_format)  # 시료명 출력
        start = 13  # 14행
        temp = self.d_species_list_index.keys()
        species_indexes = [int(x) for x in temp]
        species_indexes.sort()
        del temp

        data_row_size = 24.8
        for row, cell_data in enumerate(self.report_data, start):
            # 행 높이 설정
            self.probiotics_type_sheet.set_row(row, data_row_size)
            # 데이터 입력
            self.probiotics_type_sheet.write_string(row, 0, cell_data[0], left_normal_format)  # 유산균 인덱스. A열.
            self.probiotics_type_sheet.merge_range(row, 1,
                                                   row, 2,
                                                   cell_data[1], self.species_format)  # 유산균 종명. B,C열
            cell_data[1]
            if cell_data[2] is not None:
                self.probiotics_type_sheet.write(row, 3, cell_data[2], left_normal_format)  # D열. 투입 여부
            else:
                self.probiotics_type_sheet.write_blank(row, 3, cell_data[2], self.no_mixed_format)  # D열. 투입 여부
            if cell_data[3] != 'X':
                self.probiotics_type_sheet.write_string(row, 4, cell_data[3], blue_bold_format)  # E열. 검출 여부
            else:
                self.probiotics_type_sheet.write_string(row, 4, cell_data[3], blue_bold_no_format)  # E열. 검출 여부
        self.probiotics_type_sheet.merge_range('A{row}:C{row}'.format(row=row + 2),
                                               '불검출 항목수',
                                               thick_top_bottom_normal_format)
        self.probiotics_type_sheet.merge_range('D{row}:{col}{row}'.format(col=last_column_name, row=row + 2),
                                               self.non_detection_counter,
                                               thick_top_bottom_normal_format)
        self.probiotics_type_sheet.merge_range('A{row}:C{row}'.format(row=row + 3),
                                               '투입 외 검출 항목수',
                                               thick_top_bottom_normal_format)
        self.probiotics_type_sheet.merge_range('D{row}:{col}{row}'.format(col=last_column_name, row=row + 3),
                                               self.outer_counter,
                                               thick_top_bottom_normal_format)
        # 행 높이 설정
        self.probiotics_type_sheet.set_row(row + 1, data_row_size)  # 행 - 불검출 항목수
        self.probiotics_type_sheet.set_row(row + 2, data_row_size)  # 행 - 투입 외 검출 항목수

        # 프린터 영역 설정
        self.probiotics_type_sheet.print_area('A1:{col}{row}'.format(col=last_column_name, row=row + 2 + 1))

    def _make_taxonomy_sheet(self):
        # 서식 설정
        table_title_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'bold': True,
                'bg_color': '#c5d9f1',
                'align': 'center',
                'valign': 'vcenter'
            })
        measure_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'align': 'right',
                'valign': 'vcenter'
            })
        left_normal_format = self.workbook.add_format(
            {
                'num_format': '0.00%',
                'font_size': self.font_size,
                'left': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })
        left_right_normal_format = self.workbook.add_format(
            {
                'num_format': '0.00%',
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })
        top_normal_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'top': self.default_border
            })

        # 셀 설정
        self.taxonomy_sheet = self.workbook.add_worksheet('Taxonomy')
        self._set_sheet_layout(self.taxonomy_sheet, 'E', self.sheet_header_format)
        self.taxonomy_sheet.hide_zero()  # 0 값 숨김 : 0인 데이터와 0.001이하의 데이터를 구분하기 위해서
        self.taxonomy_sheet.set_column('A:A', 11)  # (A~E 열사이즈 합계 : 93)
        self.taxonomy_sheet.set_column('B:B', 40)
        self.taxonomy_sheet.set_column('C:E', 14)

        # 데이터 입력
        self.taxonomy_sheet.merge_range('A3:E3', 'Taxonomy 분석결과', self.title_format)  # 1행
        self.taxonomy_sheet.write_string('E5', '단위 : %', measure_format)
        self.taxonomy_sheet.write_row('A6', ['No.', self.sample_name, 'V1-V2', 'V3-V4', '평균'],
                                      table_title_format)
        start = 6  # 7행 데이터 입력(제목 제외)
        for row, cell_data in enumerate(self.report_data, start):
            self.taxonomy_sheet.write_string(row, 0, cell_data[0], left_normal_format)  # 유산균 인덱스. A열.
            self.taxonomy_sheet.write_string(row, 1, cell_data[1], self.species_format)  # 유산균 종명. B열
            if cell_data[4] is not None:  # C열. V1V2 abundance ratio
                self.taxonomy_sheet.write_number(row, 2, cell_data[4], left_normal_format)
            else:
                self.taxonomy_sheet.write_blank(row, 2, cell_data[4], self.no_mixed_format)

            if cell_data[5] is not None:  # D열. V3V4 abundance ratio
                self.taxonomy_sheet.write_number(row, 3, cell_data[5], left_normal_format)
            else:
                self.taxonomy_sheet.write_blank(row, 3, cell_data[5], self.no_mixed_format)

            if (cell_data[4] is not None) and (cell_data[5] is not None):  # E열. 평균 abundance ratio
                self.taxonomy_sheet.write_number(row, 4, sum(cell_data[4:]) / 2, left_right_normal_format)
            else:
                self.taxonomy_sheet.write_blank(row, 4, None, self.no_mixed_format)
        self.taxonomy_sheet.write_row(row + 1, 0, [None for x in range(5)], top_normal_format)

        # 조건부 서식
        self.taxonomy_sheet.conditional_format(start, 2, row, 2,
                                               {'type': 'data_bar', 'bar_color': '#FFC000'})  # V1V2 : C열
        self.taxonomy_sheet.conditional_format(start, 3, row, 3,
                                               {'type': 'data_bar', 'bar_color': '#9ACD32'})  # V3V4 : D열
        self.taxonomy_sheet.conditional_format(start, 4, row, 4,
                                               {'type': 'data_bar', 'bar_color': '#DC1687'})  # 평균  : E열

        # Taxonomy Chart 생성
        taxonomy_chart = self.workbook.add_chart({'type': 'column', 'subtype': 'stacked'})
        for row in range(start, row + 1):
            if row == start:
                taxonomy_chart.add_series({
                    'categories': ['Taxonomy', start - 1, 2, start - 1, 4],
                    'values': ['Taxonomy', row, 2, row, 4],
                    'name': ['Taxonomy', row, 1]
                })
            else:
                taxonomy_chart.add_series({
                    'values': ['Taxonomy', row, 2, row, 4],
                    'name': ['Taxonomy', row, 1]
                })
        taxonomy_chart.set_y_axis({'num_format': '0%'})
        taxonomy_chart.set_size({'width': 677, 'height': 370})
        taxonomy_chart.set_style(10)
        taxonomy_chart.set_plotarea({'layout': {'width': 1, 'height': 0.7}})
        taxonomy_chart.set_legend(
            {
                'position': 'bottom',
                'font': {'size': 8},
                'layout': {'width': 1, 'height': 0.3}
            })
        self.taxonomy_sheet.insert_chart(row + 2, 0, taxonomy_chart)

        # 프린트 영역 설정 - 차트 삽입으로 인해 범위 지정
        self.taxonomy_sheet.print_area('A1:E45')

    def _make_reference_sheet(self):
        def read_reference_list(p_ref):
            """
            19종 유산균의 계통분류군을 정리한 파일을 읽는다.
            '#'으로 시작하는 행은 주석처리되어 읽지 않는다.
            첫행은 열제목으로 읽는다.

            이전 경로
            /lustre/Tools/macrogen_analysis_toolbox2/OTU_customize/Probiotics/src/ReferenceList.txt

            :param p_ref: 19종 유산균의 계통분류군을 정리한 파일
            :return data.header(list) : 열제목
                    data(nested list) : 19종 유산균의 계통분류군의 데이터
            """
            data = list()
            data_header = None
            with open(p_ref, 'r') as o_input:
                header = False
                for i in o_input:
                    if i.startswith('#'):
                        continue
                    elif header is False:
                        data_header = i.strip().split('\t')
                        header = True
                    elif header is True:
                        data.append(i.strip().split('\t'))
                    else:
                        raise ValueError(
                            style('예상치 못한 조건을 실행하였습니다. '
                                  'ReferenceList.txt의 파일 양식을 확인하거나 개발자에게 문의하세요.',
                                  fg='red', blink=True)
                        )
            return data_header, data

        reference_data_header, reference_data = read_reference_list(self.reference_file)

        def make_taxon_string(p_data):
            string = None
            for index, element in enumerate(p_data):
                if index == 0:
                    string = element
                else:
                    string = string + '__' + element
            return string

        # 서식 설정
        header_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'border': self.default_border,
                'bold': True,
                'bg_color': '#e7e1dc',
                'align': 'center',
                'valign': 'vcenter'
            })
        table_data_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'valign': 'vcenter'
            })
        center_table_data_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })
        last_row_table_data_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'bottom': self.default_border,
                'valign': 'vcenter'
            })
        last_row_center_table_data_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'left': self.default_border,
                'right': self.default_border,
                'bottom': self.default_border,
                'align': 'center',
                'valign': 'vcenter'
            })

        # 셀 설정
        self.reference_sheet = self.workbook.add_worksheet('Reference')
        self._set_sheet_layout(self.reference_sheet, 'E', self.sheet_header_format)
        self.reference_sheet.set_column('A:A', 3)
        self.reference_sheet.set_column('B:C', 26)
        self.reference_sheet.set_column('D:D', 30)
        self.reference_sheet.set_column('E:E', 7)

        # 데이터 입력
        self.reference_sheet.merge_range('A3:E3', '유산균 19종 Reference', self.title_format)
        self.reference_sheet.insert_image('C5', self.tree_png, {'x_scale': 0.4, 'y_scale': 0.3, 'x_offset': -30})
        start_row = 30
        self.reference_sheet.write_string(start_row, 0, reference_data_header[0], header_format)
        self.reference_sheet.merge_range(start_row, 1, start_row, 2, 'Phylum(문) ~ Family(과)',
                                         header_format)
        self.reference_sheet.write_string(start_row, 3, reference_data_header[-2], header_format)
        self.reference_sheet.write_string(start_row, 4, reference_data_header[-1], header_format)
        for row, data in enumerate(reference_data, start_row + 1):
            if row == len(reference_data) + start_row + 1 - 1:  # len() + stat_num - 1
                normal_cell_format = last_row_table_data_format
                center_cell_format = last_row_center_table_data_format
            else:
                normal_cell_format = table_data_format
                center_cell_format = center_table_data_format
            self.reference_sheet.write_string(row, 0, data[0], center_cell_format)
            self.reference_sheet.merge_range(row, 1, row, 2, make_taxon_string(data[2:-3]), normal_cell_format)
            self.reference_sheet.write_string(row, 3, data[-2], normal_cell_format)
            self.reference_sheet.write_number(row, 4, int(data[-1]), center_cell_format)

        # 프린트 영역 설정
        self.reference_sheet.print_area('A1:E{row}'.format(row=row + 2))

    def _make_contact_sheet(self):
        # 서식 설정
        header_format = self.workbook.add_format(
            {
                'font_size': self.font_size + 6,
                'bold': True,
                'valign': 'vcenter'
            })
        normal_format = self.workbook.add_format(
            {
                'font_size': self.font_size,
                'valign': 'vcenter'
            })

        # 셀 설정
        self.contact_sheet = self.workbook.add_worksheet('Contact')
        self._set_sheet_layout(self.contact_sheet, 'E', self.sheet_header_format)
        self.contact_sheet.set_column('A:B', 18)
        self.contact_sheet.set_column('C:C', 21)
        self.contact_sheet.set_column('D:E', 18)

        # 데이터 입력
        self.contact_sheet.insert_image('A10', self.cover_logo_jpg, {'x_offset': 25})
        self.contact_sheet.write_string('A43', 'Macrogen Korea', header_format)
        self.contact_sheet.write_string('D43', 'Contact', header_format)
        self.contact_sheet.write_string('A44', '10F, 254 Beotkkot-ro, Geumcheon-gu, Seoul, Rep. of Korea',
                                        normal_format)
        self.contact_sheet.write_string('D44', 'Web : www.macrogen.com', normal_format)
        self.contact_sheet.write_string('A45', 'Phone: +82-2113-7000', normal_format)
        self.contact_sheet.write_string('D45', 'Lims : http://dna.macrogen.com', normal_format)

        # 프린트 영역 설정
        self.contact_sheet.print_area('A1:E45')

    def _save_excel(self):
        self.info_sheet.activate()
        self.workbook.close()

    def make_report(self):
        self._make_info_sheet()
        self._organize_data_for_report()
        self._make_probiotics_type_sheet()
        self._make_taxonomy_sheet()
        self._make_reference_sheet()
        self._make_contact_sheet()
        self._save_excel()


def check_sample_name(p_name, p_region):
    if p_region == 'v1v2':
        a = 1
        b = 2
    elif p_region == 'v3v4':
        a = 3
        b = 4
    else:
        raise ValueError(
            style('check_sample_name 함수의 p_region 매개변수의 값이 잘못되었습니다.', fg='red', bg='yellow', blink=True) \
            + '\n' \
            + 'p_region : %s' % p_region)

    string = 'V{a}{b} v{a}{b} V{a}V{b} v{a}v{b} v{a}.{b} V{a}.{b} v{a}.v{b} V{a}.V{b}'.format(a=a, b=b)
    region_string = string.split(' ')
    for ele in region_string:
        stat = p_name.find(ele)
        if stat > -1:
            return p_name.replace(ele, '').rstrip('.')
    else:
        separator = '.'
        base, sep, suffix = p_name.rpartition(separator)
        if (base == '') and (sep == ''):
            return p_name
        elif (base == '') and (sep == separator):
            secho('Error: 시료명에 문제가 있습니다.', fg='red')
            secho(f'{p_name} --> ({base}, {sep}, {suffix})', fg='cyan')
            echo('base 값이 없습니다.')
            exit()
        else:
            return base
