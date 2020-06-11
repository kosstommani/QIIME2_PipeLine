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
from click import secho, echo, style
from beautifultable import BeautifulTable
from SpoON.util import parse_config, glob_dir, check_file_type


class Path:
    def __init__(self, p_base, p_order_number, p_analysis_number):
        """
        분석 디렉터리 경로를 반환한다.

        :param p_base: Amplicon Metagenome 분석 작업 디렉터리
        :param p_order_number: 수주번호.
        :param p_analysis_number: 분석 디렉터리. ex) Analysis_1
        """
        self.base = p_base
        self.order_number = p_order_number
        self.analysis_number = p_analysis_number
        self.analysis_dir_name = parse_config()['Analysis_Dir_Name']

    @property
    def order_number_path(self):
        return os.path.join(self.base, self.order_number)

    @property
    def analysis_number_path(self):
        return os.path.join(self.order_number_path, self.analysis_number)

    @property
    def read_assembly_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['read_assembly'])

    @property
    def alignment_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['alignment'])

    @property
    def cd_hit_otu_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['cd_hit_otu'])

    @property
    def closed_otu_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['closed'])

    @property
    def r_dada2_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['r_dada2'])

    @property
    def r_dada2_summary_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['r_dada2_summary'])

    @property
    def taxonomy_assignment_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['taxonomy_assignment'])

    @property
    def phylogeny_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['phylogeny'])

    @property
    def biom_path(self):
        return os.path.join(self.analysis_number_path, self.analysis_dir_name['biom'])


class ReportPath:
    def __init__(self, p_base, p_order_number, p_analysis_number):
        """
        보고서의 디렉터리 구조 및 경로를 생성한다.

        :param p_base:
        :param p_order_number:
        :param p_analysis_number:
        """
        self.base = p_base
        self.order_number = p_order_number
        self.analysis_number = p_analysis_number
        self.__report_dir_name = None

    @property
    def report_base_path(self):
        report_base_dir = f'{self.analysis_number}_Report'
        return os.path.join(self.base, self.order_number, report_base_dir)

    @property
    def report_dir_name(self):
        return self.__report_dir_name

    @report_dir_name.setter
    def report_dir_name(self, p_data):
        self.__report_dir_name = p_data

    @property
    def report_dir_path(self):
        return os.path.join(self.report_base_path, self.report_dir_name)

    @property
    def otu_results_path(self):
        return os.path.join(self.report_dir_path, 'OTU_Results')

    def get_report_dir_name(self):
        l_dir_name = glob_dir(self.report_base_path, 'Report_[0-9]*', 'many', False)
        if l_dir_name is None:
            return 'Report_1'
        else:
            new_num = max([int(x.split('Report_')[-1]) for x in l_dir_name]) + 1
            return f'Report_{new_num}'

    def check_report_path_and_set(self, p_report_dir):
        """
        Analysis_?_Report 디렉터리 내에 보고서를 생성할 보고서 번호(디렉터리)를 지정한다.

        :type p_report_dir: None or str
        :param p_report_dir: 보고서 번호(디렉터리)
                        None: Analysis_?_Report 디렉터리내에 보고서 번호(디렉터리)를 확인하고
                              있는 경우 그 다음 번호를 지정. 없을 경우 새번호 지정.
                        str: 전달된 보고서 번호로 설정.
        :return:
        """
        if p_report_dir is None:
            if check_file_type(self.report_base_path, 'exists'):
                self.report_dir_name = self.get_report_dir_name()
            else:
                self.report_dir_name = 'Report_1'
        else:
            self.report_dir_name = p_report_dir


class ReportPathV1(ReportPath):
    @property
    def pcoa_2d_path(self):
        return os.path.join(self.otu_results_path, '2D_PCoA')

    @property
    def check_list_for_pcoa_2d(self):
        item = ['weighted_unifrac_pc_2D_PCoA_plots.html', 'unweighted_unifrac_pc_2D_PCoA_plots.html']
        return item

    @property
    def alpha_rarefaction_path(self):
        return os.path.join(self.otu_results_path, 'Alpha_rarefaction')

    @property
    def check_list_for_alpha_rarefaction(self):
        item = ['alpha_div_collated', 'alpha_div_collated/PD_whole_tree.txt', 'alpha_div_collated/chao1.txt',
                'alpha_div_collated/observed_otus.txt', 'alpha_rarefaction_plots',
                'alpha_rarefaction_plots/average_plots', 'alpha_rarefaction_plots/rarefaction_plots.html']
        return item

    @property
    def beta_diversity_path(self):
        return os.path.join(self.otu_results_path, 'Beta_diversity')

    @property
    def check_list_for_beta_diversity(self):
        item = ['unweighted_unifrac_dm.txt', 'unweighted_unifrac_pc.txt', 'weighted_unifrac_dm.txt',
                'weighted_unifrac_pc.html', 'weighted_unifrac_pc.txt', 'unweighted_unifrac_emperor_pcoa_plot',
                'unweighted_unifrac_emperor_pcoa_plot/index.html', 'weighted_unifrac_emperor_pcoa_plot',
                'weighted_unifrac_emperor_pcoa_plot/index.html']
        return item

    @property
    def biom_path(self):
        return os.path.join(self.otu_results_path, 'biom')

    @property
    def diversity_index_path(self):
        return os.path.join(self.otu_results_path, 'Diversity_Index')

    @property
    def check_list_for_diversity_index(self):
        item = ['adiv.html', 'adiv.txt', 'adiv.xlsm']
        return item

    @property
    def src_path(self):
        return os.path.join(self.otu_results_path, 'src')

    @property
    def check_list_for_src(self):
        item = ['bootstrap.css', 'default.css', 'img']
        return item

    @property
    def taxonomy_assignment_path(self):
        return os.path.join(self.otu_results_path, f'Taxonomy_assignment')

    @property
    def check_list_for_taxonomy_assignment_uclust(self):
        item = ['Taxonomy_abundance_ratio.xlsx', 'otu_table_shared_with_tax_assignment.xlsx']
        return item

    @property
    def check_list_for_taxonomy_assignment_blast(self):
        item = ['OTU_BLAST.xlsx']
        return self.check_list_for_taxonomy_assignment_uclust + item

    @property
    def upgma_tree_path(self):
        return os.path.join(self.otu_results_path, 'UPGMAtree')

    @property
    def check_list_for_upgma_tree(self):
        item = ['weighted_unifrac_upgma_cluster.html', 'weighted_unifrac_upgma_cluster.tre',
                'unweighted_unifrac_upgma_cluster.html', 'unweighted_unifrac_upgma_cluster.tre']
        return item


class ReportPathV2(ReportPath):
    @property
    def data_path(self):
        return os.path.join(self.otu_results_path, 'Data')

    @property
    def biom_and_otus_rep_path(self):
        return os.path.join(self.otu_results_path, 'biom_and_otus_rep')

    @property
    def excel_path(self):
        return os.path.join(self.otu_results_path, 'Excel')

    @property
    def src_path(self):
        return os.path.join(self.data_path, 'src')

    @property
    def diversity_index_path(self):
        return os.path.join(self.data_path, 'Diversity_Index')

    @property
    def upgma_tree_path(self):
        return os.path.join(self.data_path, 'UPGMAtree')


class AnalysisDataTable:
    d_table_style = \
        {
            'default': BeautifulTable.STYLE_DEFAULT, 'none': BeautifulTable.STYLE_NONE,
            'dotted': BeautifulTable.STYLE_DOTTED, 'separated': BeautifulTable.STYLE_SEPARATED,
            'compact': BeautifulTable.STYLE_COMPACT, 'mysql': BeautifulTable.STYLE_MYSQL,
            'markdown': BeautifulTable.STYLE_MARKDOWN, 'rst': BeautifulTable.STYLE_RST,
            'box': BeautifulTable.STYLE_BOX, 'box_doubled': BeautifulTable.STYLE_BOX_DOUBLED,
            'box_rounded': BeautifulTable.STYLE_BOX_ROUNDED, 'grid': BeautifulTable.STYLE_GRID,
        }

    def __init__(self, kargs):
        """
        분석된 데이터(디렉터리, 파일)의 존재 유무를 딕셔너리에 저장하고, 표 형태로 출력한다.

        :param kargs:
                    'set_style':
                    'column_headers':
        """
        self.table = BeautifulTable()
        self.table_style = kargs['set_style']
        self.column_headers = kargs['column_headers']
        self.set_table_style(self.table_style)
        self.table.column_headers = self.column_headers
        self.__read_assembly = None
        self.__cd_hit_otu = None
        self.__closed = None
        self.__r_dada2 = None
        self.__r_dada2_summary = None
        self.__alignment = None
        self.__phylogeny = None
        self.__taxonomy_assignment = None
        self.__biom = None
        self.__info = None
        self.__metadata = None

    def set_table_style(self, p_style):
        self.table.set_style(self.d_table_style[p_style])

    def calculate_column_widths(self):
        self.table._calculate_column_widths()

    def get_table_width(self):
        return self.table.get_table_width()

    @staticmethod
    def get_style_string(p_l_name, p_l_state, p_true_color='cyan', p_false_color='red'):
        if len(p_l_name) == len(p_l_state):
            l_style_string = list()
            for i in range(len(p_l_name)):
                l_style_string.append(style(p_l_name[i], fg=p_true_color if p_l_state[i] else p_false_color))
            return '\n'.join(l_style_string)
        else:
            raise IndexError(p_l_name, p_l_state)

    def __set_data_only(self, l_dir, l_dir_state, l_file, l_file_state):
        """
        분석 디렉터리내에 파일만 있을 경우, 디렉터리 및 구성 파일을 예시의 딕셔너리 형태로 저장한다.

        ex)
        {
            'Read_Assembly':
                {
                    'dir': True,
                    'file':
                        {
                            'HN00108505.FLASH.xlsx': True,
                            'HN00108505.F_length.csv': True,
                            'HN00108505.pooled.fastq': True,
                            'STAT.txt': True
                        }
                }
        }

        :param l_dir:
        :param l_dir_state:
        :param l_file:
        :param l_file_state:
        :return:
        """
        if (len(l_dir) == 1) and (len(l_dir_state) == 1):
            data = {l_dir[0]: {'dir': l_dir_state[0]}}
            style_dir = self.get_style_string(l_dir, l_dir_state)
        else:
            raise IndexError(f'l_dir: {l_dir}\nl_dir_state: {l_dir_state}')

        if len(l_file) == len(l_file_state):
            d_files_state = {}
            for i in range(len(l_file)):
                d_files_state[l_file[i]] = l_file_state[i]
            data[l_dir[0]].update({'file': d_files_state})
            style_file = self.get_style_string(l_file, l_file_state, 'cyan')
        else:
            raise IndexError(f'l_file: {l_file}\nl_file_state: {l_file_state}')
        self.table.append_row([style_dir, style_file])
        return data

    def __set_data_many(self, l_dir, l_dir_state, l_file, l_file_state):
        """
        분석 디렉터리내에 디렉터리가 있을 경우, 디렉터리 및 구성 파일을 예시의 딕셔너리 형태로 저장한다.

        ex)
        {
            'Taxonomy_Assignment': {'dir': True},
            'blast_NCBI_16S':
                {
                    'dir': True,
                    'file':
                        {
                            'OTU_BLAST.xlsx': False,
                            'otus_rep_tax_assignments.log': True,
                            'otus_rep_tax_assignments.txt': True
                        }
                },
            'rdp_RDP':
                {
                    'dir': True,
                    'file':
                        {
                            'otus_rep_tax_assignments.log': False,
                            'otus_rep_tax_assignments.txt': True
                        }
                },
            'uclust_RDP':
                {
                    'dir': True,
                    'file':
                        {   'otus_rep_tax_assignments.log': True,
                            'otus_rep_tax_assignments.txt': True
                        }
                },
            'uclust_UNITE':
                {
                    'dir': True,
                    'file':
                        {
                            'otus_rep_tax_assignments.log': False,
                            'otus_rep_tax_assignments.txt': False
                        }
                }
        }

        :param l_dir:
        :param l_dir_state:
        :param l_file:
        :param l_file_state:
        :return:
        """
        data = dict()
        if len(l_dir) == len(l_dir_state):
            for i in range(len(l_dir)):
                data[l_dir[i]] = {'dir': l_dir_state[i]}
            style_dir = self.get_style_string(l_dir, l_dir_state)
        else:
            raise IndexError(f'l_dir: {l_dir}\nl_dir_state: {l_dir_state}')

        if len(l_file) == len(l_file_state):
            l_ele_state = [isinstance(x, list) for x in l_file]
            sub_table = BeautifulTable()
            sub_table.set_style(self.d_table_style.get('box'))
            if all(l_ele_state) is True:
                if len(l_dir) == len(l_file):
                    l_dir_list = l_dir
                else:
                    l_dir_list = l_dir[1:]
                for e_dir, files, states in zip(l_dir_list, l_file, l_file_state):
                    style_files = self.get_style_string(files, states, 'cyan')
                    sub_table.append_row([style_files])
                    d_file = dict()
                    for file, state in zip(files, states):
                        d_file[file] = state
                    data[e_dir].update({'file': d_file})
                self.table.append_row([style_dir, sub_table])
            elif all(l_ele_state) is False:
                d_file = dict()
                for i in range(len(l_file)):
                    d_file[l_file[i]] = l_file_state[i]
                data[l_dir[-1]].update({'file': d_file})
                style_file = self.get_style_string(l_file, l_file_state, 'cyan')
                self.table.append_row([style_dir, style_file])
        else:
            raise IndexError(f'l_file: {l_file}\nl_file_state: {l_file_state}')
        return data

    @property
    def read_assembly(self):
        return self.__read_assembly

    def set_read_assembly_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__read_assembly = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def cd_hit_otu(self):
        return self.__cd_hit_otu

    def set_cd_hit_otu_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__cd_hit_otu = self.__set_data_many(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def closed(self):
        return self.__closed

    def set_closed_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__closed = self.__set_data_many(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def r_dada2(self):
        return self.__r_dada2

    def set_r_dada2_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__r_dada2 = self.__set_data_many(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def r_dada2_summary(self):
        return self.__r_dada2_summary

    def set_r_dada2_summary_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__r_dada2_summary = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def alignment(self):
        return self.__alignment

    def set_alignment_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__alignment = self.__set_data_many(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)
        return

    @property
    def phylogeny(self):
        return self.__phylogeny

    def set_phylogeny(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__phylogeny = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)
        return

    @property
    def taxonomy_assignment(self):
        return self.__taxonomy_assignment

    def set_taxonomy_assignment_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__taxonomy_assignment = self.__set_data_many(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)
        return

    @property
    def biom(self):
        return self.__biom

    def set_biom_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__biom = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def info(self):
        return self.__info

    def set_info_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__info = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

    @property
    def metadata(self):
        return self.__metadata

    def set_metadata_data(self, p_l_dir, p_l_dir_state, p_l_file, p_l_file_state):
        self.__metadata = self.__set_data_only(p_l_dir, p_l_dir_state, p_l_file, p_l_file_state)

