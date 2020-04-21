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
__version__ = '1.0.3'

import os
from click import secho, echo, style
from SpoON.util import check_file_type, parse_html, read_metadata_for_sort, run_cmd, check_run_cmd, \
    read_data, glob_dir, read_all_data, parse_config
from theCups.report import Report
from theCups.data_structure import ReportPathV1
from theCups.read_assembly import StatExcel
from theCups.clustering_and_otu_picking import SummaryForCdHitOtu, SummaryForClosed
from theCups.diversity_index import run_alpha_diversity
from theCups.rarefaction import run_alpha_rarefaction
from theCups.beta_diversity import run_beta_diversity_through_plots
from theCups.pcoa import run_make_2d_pcoa
from theCups.upgma_tree import run_upgma_tree, TreeController
from theCups.taxonomy import run_summarize_taxa, run_plot_taxa_summary, TaxonHTML, \
    TaxonomyAbundanceExcel, TaxonomySharedExcel, run_make_shared
from humanize import intcomma
from collections import defaultdict

CONFIG = parse_config()
MiSeq_V1_TEMPLATE_PATH = CONFIG['theCups']['MiSeq_V1']['template_path']


class MiSeqReportV1(Report):
    def __init__(self, kargs):
        super().__init__(kargs, 'MiSeq')
        self.report_paths = ReportPathV1(self.base_path, self.order_number, self.analysis_number)

    @property
    def html_template_path(self):
        """
        HTML 보고서 작성에 필요한 HTML 템플릿의 경로를 지정합니다.

        :return:
        """
        global MiSeq_V1_TEMPLATE_PATH
        return MiSeq_V1_TEMPLATE_PATH

    def checklist_read_assembly_files(self, p_type):
        """
        Report.checklist_read_assembly_files : overriding

        :param p_type:
        :return:
        """
        l_files = [
            'STAT.txt',
            f'{self.order_number}.pooled.{p_type}',
            f'{self.order_number}.FLASH.xlsx',
            f'{self.order_number}.F_length.csv',
        ]
        return l_files

    def checklist_cd_hit_otu_files(self):
        """
        Report.checklist_cd_hit_otu_files : overriding

        :return:
        """
        l_files = ['trim.log', 'pick_otus.txt', 'otus_rep.fasta', 'chimaric.ids']
        return l_files

    def check_analysis_data(self):
        """
        HTML 보고서를 구성하기 위해 필요한 디렉터리 및 파일들이 존재하는지 확인합니다.
        확인된 디렉터리 및 파일들은 테이블 형태로 출력합니다.

        :return:
        """
        self.make_analysis_data_table(
            {
                'set_style': 'box',
                'column_headers': ['DIR', 'FILE'],
            })
        self.check_info_file()
        self.check_metadata_file()
        cd_hit_otu_dir_name, cd_hit_otu_dir_state = self.exists(self.analysis_paths.cd_hit_otu_path)
        closed_dir_name, closed_dir_state = self.exists(self.analysis_paths.closed_otu_path)

        # OTU Picking 방식 확인: Denovo(CD-HIT-OTU), CLOSED(UCLUST)
        if all([cd_hit_otu_dir_state, closed_dir_state]):
            secho(f'Error: 분석 디렉터리에 {cd_hit_otu_dir_name} 와 {closed_dir_name} 디렉터리가 모두 존재합니다.\n'
                  f'       설계 개념상 같이 존재하면 안 됩니다.', fg='red')
            secho('\t---> Analysis_? 디렉터리 생성하여, 두 디렉터리를 분리하세요.', fg='magenta')
            secho('\t---> 분석 데이터 검증을 위해 생성한 디렉터리일 경우 해당 디렉터리의 이름을 변경하세요.', fg='magenta')
            exit()
        elif cd_hit_otu_dir_state is True:
            self.check_read_assembly_dir(p_type='fastq')
            self.check_cd_hit_otu_dir()
            self.otu_picking_method = 'DENOVO'
        elif closed_dir_state is True:
            self.check_read_assembly_dir(p_type='fasta')
            self.check_closed_otu_dir()
            self.otu_picking_method = 'CLOSED'
        else:
            self.check_read_assembly_dir(p_type='fastq')
            self.check_read_assembly_dir(p_type='fasta')
            self.check_cd_hit_otu_dir()
            self.check_closed_otu_dir()
        self.check_alignment_dir()
        self.check_phylogeny_dir()
        self.check_taxonomy_assignment_dir()
        self.check_biom_dir()
        self.analysis_data.calculate_column_widths()
        echo('\n')
        secho(f'- {self.platform} 분석 데이터 확인 -'.center(self.analysis_data.get_table_width() - 8))
        echo(self.analysis_data.table)
        return

    def check_and_set_for_report_path(self):
        """
        Report Number를 지정합니다.

        :return:
        """
        if self.kargs['report_number'] is None:
            self.report_paths.check_report_path_and_set(None)
        else:
            self.report_paths.check_report_path_and_set(self.kargs['report_number'])

    def taste_cofi(self):
        i_data_table = self.get_analysis_data_table(
            {
                'set_style': 'dotted',
                'column_headers': ['주요 항목', '세부 항목', '비고']
            }
        )

        checked_data = list()
        metadata_status, metadata_file = self.confirm_metadata('taste')
        sample_count = 10
        group_count = 0
        if metadata_status:
            l_metadata = read_data(metadata_file, '\t')
            sample_count = len(l_metadata) - 1
            l_group_count = [len(x) for x in l_metadata]
            if len(set(l_group_count)) == 1:
                group_count = l_group_count[0] - 1
            else:
                secho('Error: metadata.txt 파일의 그룹의 개수가 다릅니다.', fg='red', blink=True)
                echo(f'{l_group_count}')
                secho('\t--> 행별 컬럼의 수를 확인하세요.', fg='magenta')
                secho('\t--> 시료의 개수와 그룹 정보에 기반한 결과 확인을 수행하지 않습니다.', fg='magenta')
        else:
            secho('Warning: metadata.txt 파일이 없습니다.', fg='yellow')
            secho('\t--> 시료의 개수와 그룹 정보에 기반한 결과 확인을 수행하지 않습니다.', fg='magenta')

        for target in ['OTUanalysis.html', 'OTU_Results']:
            target_path = os.path.join(self.report_paths.report_dir_path, target)
            returned_name, returned_status = self.exists(target_path)
            checked_data.append((returned_name, 'main', returned_status, ''))
            del returned_name, returned_status

        for target in ['STAT.html', 'otu_table_summary.html', 'STAT.xlsm', 'otus_rep.fasta']:
            target_path = os.path.join(self.report_paths.otu_results_path, target)
            returned_name, returned_status = self.exists(target_path)
            checked_data.append((returned_name, 'sub', returned_status, ''))
            del returned_name, returned_status

        l_dir = ['pcoa_2d_path', 'alpha_rarefaction_path', 'beta_diversity_path', 'diversity_index_path',
                 'src_path', 'upgma_tree_path']
        for dir_name in l_dir:
            dir_path = getattr(self.report_paths, dir_name)
            returned_name, returned_status = self.exists(dir_path)
            comment_text = ''
            if dir_name == 'pcoa_2d_path':
                if sample_count <= 2:
                    comment_text = style('2D PCoA 생략(<=2)', fg='yellow')
            elif dir_name == 'beta_diversity_path':
                if sample_count <= 3:
                    comment_text = style('3D PCoA 생략(<=3)', fg='yellow')
            elif dir_name == 'upgma_tree_path':
                if sample_count <= 2:
                    comment_text = style('UPGMA Tree 생략(<=2)', fg='yellow')
                else:
                    if group_count > 0:
                        header_temp = [x[1:] for x in l_metadata if x[0] == '#SampleID']
                        if len(header_temp) == 1:
                            l_header = header_temp[0]
                        else:
                            secho('Warning: metadata.txt 파일에 Header가 중복되어 있습니다.', fg='yellow')
                            secho('\t--> 첫 번째 Header가 적용됩니다.', fg='magenta')
                            l_header = header_temp[0]
            checked_data.append((returned_name, 'main', returned_status, comment_text))
            del returned_name, returned_status

            for target in getattr(self.report_paths, f'check_list_for_{dir_name.rstrip("_path")}'):
                target_path = os.path.join(dir_path, target)
                returned_name, returned_status = self.exists(target_path)
                checked_data.append((returned_name, 'sub', returned_status, ''))
                del returned_name, returned_status

            if dir_name == 'upgma_tree_path':
                if group_count > 0:
                    for group in l_header:
                        tre_file = os.path.join(dir_path, f'upgma_cluster.{group}.tre')
                        returned_name, returned_status = self.exists(tre_file)
                        checked_data.append((returned_name, 'sub', returned_status, ''))
                    del returned_name, returned_status

        # Taxonomy Assignment 결과 확인
        l_taxonomy_db = list()
        for key in self.analysis_data.taxonomy_assignment.keys():
            if key == self.analysis_paths.analysis_dir_name['taxonomy_assignment']:
                continue
            else:
                l_taxonomy_db.append(key)
        for tool_db in l_taxonomy_db:
            dir_path = f'{self.report_paths.taxonomy_assignment_path}_{tool_db}'
            returned_name, returned_status = self.exists(dir_path)
            if sample_count == 1:
                comment_text = 'Area Plot 생략(==1)'
            else:
                comment_text = ''
            checked_data.append((returned_name, 'main', returned_status, comment_text))

            tool, db = tool_db.split('_', 1)
            check_items1 = getattr(self.report_paths, f'check_list_for_taxonomy_assignment_{tool}')
            check_items2 = ['count_data', 'count_data/Taxonomy_abundance_count.xlsx']
            # f'otu_table.{tool_db}_L2.txt', f'otu_table.{tool_db}_L3.txt',
            # f'otu_table.{tool_db}_L4.txt', f'otu_table.{tool_db}_L6.txt',
            # f'otu_table.{tool_db}_L7.txt',
            # f'count_data/otu_table.{tool_db}_L2.txt', f'count_data/otu_table.{tool_db}_L3.txt',
            # f'count_data/otu_table.{tool_db}_L4.txt', f'count_data/otu_table.{tool_db}_L5.txt',
            # f'count_data/otu_table.{tool_db}_L6.txt', f'count_data/otu_table.{tool_db}_L7.txt']
            for item in check_items1 + check_items2:
                taxonomy_file = os.path.join(dir_path, item)
                returned_name, returned_status = self.exists(taxonomy_file)
                checked_data.append((returned_name, 'sub', returned_status, ''))
            l_globed_dir = glob_dir(dir_path, 'taxa_summary_plots_*', p_mode='many')
            if l_globed_dir is None:
                checked_data.append(('taxa_summary_plots_*', 'sub', False, '미검출'))
                continue
            for plot_dir in l_globed_dir:
                l_plot = ['area_charts.html', 'bar_charts.html']  # 'charts', 'css', 'js', 'raw_data']
                for ele in l_plot:
                    target_path = os.path.join(dir_path, plot_dir, ele)
                    returned_name, returned_status = self.exists(target_path)
                    checked_data.append((returned_name, 'sub', returned_status, ''))
                del plot_dir, ele
        del dir_path, returned_name, returned_status, tool_db, tool, db, item, comment_text, l_globed_dir

        # biom 파일 확인
        returned_name, returned_status = self.exists(self.report_paths.biom_path)
        checked_data.append((returned_name, 'main', returned_status, ''))
        biom_file = os.path.join(self.report_paths.biom_path, 'otu_table.biom')
        returned_name, returned_status = self.exists(biom_file)
        checked_data.append((returned_name, 'sub', returned_status, ''))
        for tool_db in l_taxonomy_db:
            biom_file = os.path.join(self.report_paths.biom_path, f'otu_table.{tool_db}.biom')
            returned_name, returned_status = self.exists(biom_file)
            checked_data.append((returned_name, 'sub', returned_status, ''))
        biom_metadata_file = os.path.join(self.report_paths.biom_path, 'metadata.txt')
        returned_name, returned_status = self.exists(biom_metadata_file)
        checked_data.append((returned_name, 'sub', returned_status, ''))
        del returned_name, returned_status, biom_file, tool_db

        # 확인된 결과 디렉터리 및 파일들 정리
        for item, item_type, status, comment in checked_data:
            if item_type == 'main':
                styled_item = i_data_table.get_style_string([item], [status])
                i_data_table.table.append_row([styled_item, item_type, comment])
            elif item_type == 'sub':
                styled_item = i_data_table.get_style_string([item], [status])
                i_data_table.table.append_row(['', styled_item, comment])
        i_data_table.calculate_column_widths()
        echo()  # 빈 줄
        secho(f'- {self.kargs["report_number"]}: 결과 확인 -'.center(i_data_table.get_table_width()),
              fg='white', bg='cyan', bold=True)
        echo(i_data_table.table)
        self.l_checked_data = checked_data

    def convert_checked_results_to_page_status(self):
        d_report_results = defaultdict(list)
        indicator = False
        for checked_data in self.l_checked_data:
            if checked_data[1] == 'main':
                if indicator is True:
                    d_report_results[main_result[0]] = l_sub_results
                else:
                    indicator = True
                l_sub_results = list()
                main_result = checked_data
            elif checked_data[1] == 'sub':
                l_sub_results.append(checked_data)
            else:
                raise ValueError(f'{checked_data}: main 또는 sub 문자가 원소에 있어야 됩니다.')
        del indicator, main_result, l_sub_results, checked_data

        for key in d_report_results.keys():
            l_sub_data = d_report_results.get(key)
            l_status = [x[2] for x in l_sub_data]
            if key == 'OTU_Results':
                for sub in l_sub_data:
                    if sub[0] == 'STAT.html':
                        stat_html = sub[2]
                    elif sub[0] == 'STAT.xlsm':
                        stat_xlsm = sub[2]
                    elif sub[0] == 'otu_table_summary.html':
                        if self.otu_picking_method == 'DENOVO':
                            if sub[2] is True:
                                self.page_status['summary.html'] = True
                            else:
                                self.page_status['summary.html'] = False
                        elif self.otu_picking_method == 'CLOSED':
                            if sub[2] is True:
                                self.page_status['summary_for_closed.html'] = True
                            else:
                                self.page_status['summary_for_closed.html'] = False
                else:
                    if all([stat_html, stat_xlsm]):
                        self.page_status['STAT.html'] = True
                    else:
                        self.page_status['STAT.html'] = False
            elif key == '2D_PCoA':
                if all(l_status):
                    self.page_status['weighted_unifrac_pc_2D_PCoA_plots.html'] = True
                else:
                    self.page_status['weighted_unifrac_pc_2D_PCoA_plots.html'] = False
            elif key == 'Alpha_rarefaction':
                if all(l_status):
                    self.page_status['rarefaction_plots.html'] = True
                else:
                    self.page_status['rarefaction_plots.html'] = False
            elif key == 'Beta_diversity':
                if all(l_status):
                    self.page_status['(DIR) Beta Diversity 3D PCoA HTML'] = True
                else:
                    self.page_status['(DIR) Beta Diversity 3D PCoA HTML'] = False
            elif key == 'Diversity_Index':
                if all(l_status):
                    self.page_status['adiv.html'] = True
                else:
                    self.page_status['adiv.html'] = False
            elif key == 'UPGMAtree':
                if all(l_status):
                    self.page_status['upgma_cluster.html'] = True
                else:
                    self.page_status['upgma_cluster.html'] = False
            elif 'Taxonomy_assignment' in key:
                tool_db = key.lstrip('Taxonomy_assignment')
                if all(l_status):
                    self.page_status[f'(DIR) Taxonomy Assignment HTML & Excel: {tool_db}'] = True
                else:
                    self.page_status[f'(DIR) Taxonomy Assignment HTML & Excel: {tool_db}'] = False
            else:
                continue

    def make_main_page(self):
        echo('>>> Main Page 생성 시작')
        d_info = self.parse_yaml(os.path.join(self.analysis_paths.analysis_number_path, 'info.txt'))
        # Length Filtering 조건 확인
        length_filtering_csv_file = os.path.join(self.analysis_paths.read_assembly_path,
                                                 f'{self.order_number}.F_length.csv')
        if check_file_type(length_filtering_csv_file, 'exists'):
            l_length_filtering = read_data(length_filtering_csv_file, ',')
            for stat in l_length_filtering:
                if stat[0] == 'lower_bp':
                    if len(set(stat[1:])) == 1:
                        lower_bp = stat[1]
                    else:
                        secho(f'Error: {stat}', fg='red')
                        raise RuntimeError(f'{length_filtering_csv_file} 파일의 길이 조건이 다른 시료가 있습니다.')
                elif stat[0] == 'upper_bp':
                    if len(set(stat[1:])) == 1:
                        upper_bp = stat[1]
                    else:
                        secho(f'Error: {stat}', fg='red')
                        raise RuntimeError(f'{length_filtering_csv_file} 파일의 길이 조건이 다른 시료가 있습니다.')
                else:
                    continue
        else:
            secho('Error: Length Filtering CSV 파일이 없습니다.', fg='red', blink=True)
            echo(length_filtering_csv_file)
            secho('\t--> None으로 설정합니다.', fg='magenta')
            lower_bp = None
            upper_bp = None
        length_filter = f'{lower_bp}bp <= Good Sequences <= {upper_bp}bp'

        # Taxonomy Assignment 확인
        l_taxonomy_dir = list()
        for key in self.page_status.keys():
            if 'Taxonomy Assignment HTML' in key:
                if self.page_status.get(key):
                    db_dir = key.split(':')[1]
                    l_taxonomy_dir.append(db_dir.strip())
        # DB Version
        l_db_version = list()
        for db_dir in l_taxonomy_dir:
            tool, db = db_dir.split('_', 1)
            db_version_file = os.path.join(self.analysis_paths.taxonomy_assignment_path, db_dir, 'db_version.txt')
            db_version = read_all_data(db_version_file)
            l_db_version.append(f'{db_version.strip()}({tool.upper()})')

        # OTU Picking Method & Tool
        if self.otu_picking_method == 'DENOVO':
            otu_picking_method = 'De novo(CD-HIT-OTU)'
        elif self.otu_picking_method == 'CLOSED':
            otu_picking_method = 'Closed-reference based(UCLUST)'
        else:
            raise RuntimeError(f'otu_picking_method: {self.otu_picking_method}')

        if self.template_env is None:
            raise RuntimeError('Report Template 설정이 필요합니다. set_report_templte() 실행')
        template = self.template_env.get_template('main.html')
        main_html = os.path.join(self.report_paths.report_dir_path, 'OTUanalysis.html')
        with open(main_html, 'w') as o_main_html:
            # DB 정보 데이터 추가
            o_main_html.write(
                template.render(
                    d_info,
                    length_filter=length_filter, db_version='<br>'.join(l_db_version),
                    page_status=self.page_status, skipped_page=self.skipped_page,
                    taxonomy_dir=l_taxonomy_dir, otu_method=otu_picking_method,
                )
            )
        secho('>>> Main Page 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def make_read_assembly_page(self):
        echo('>>> Read Assembly: STAT 생성 시작')
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        echo_html = 'STAT.html'
        stat_name = 'STAT.txt'

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            metadata_order = read_metadata_for_sort(metadata_file)
        else:
            return False

        stat_file_status = self.analysis_data.read_assembly[analysis_dir_name['read_assembly']]['file'][stat_name]
        if self.confirm_file(
            {
                'file_name': stat_name,
                'file_status': stat_file_status,
                'html_name': echo_html,
            }
        ):
            stat_file = os.path.join(self.analysis_paths.read_assembly_path, stat_name)
        else:
            return False

        # STAT.xlsm 생성
        try:
            i_stat_excel = StatExcel(
                {
                    'stat_file': stat_file,
                    'metadata_order': metadata_order,
                    'out_path': self.report_paths.otu_results_path,
                })
            i_stat_excel.read_assembly_stat()
            i_stat_excel.make_excel()
        except Exception as err:
            self.page_status[echo_html] = False
            secho('Error: STAT Excel 생성에 문제가 있습니다.', fg='red', blink=True)
            echo(f'메시지: {err.with_traceback(None)}')
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False
        
        # STAT.html 생성
        template = self.template_env.get_template('stat.html')
        stat_html = os.path.join(self.report_paths.otu_results_path, echo_html)
        with open(stat_html, 'w') as o_stat_html:
            o_stat_html.write(template.render(assembly_stat=i_stat_excel.l_stat))
        self.page_status[echo_html] = True
        secho('>>> Read Assembly: STAT 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def check_cd_hit_otu_data(self, p_echo_html, p_related_file=False):
        """
        CD-HIT-OTU 결과에 대한 결과 페이지를 생성하기 위한 결과파일들의 존재 유무를 확인한다.
        존재 유무가 확인된 파일들은 경로가 포함된 파일명 형태로 생성되며, 딕션너리 형태로 반환된다.

        :param p_echo_html:
        :param p_related_file: MiSeq_V1 에서 확인하는 파일이외에 다른 파일을 확인하고자 할 경우, 확인에 필요한 공통 데이터를 반환.
        :return: p_related_file의 상태값에 따라 반환되는 변수의 개수가 다름.
        """
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        cd_hit_otu_dir_status = self.analysis_data.cd_hit_otu[analysis_dir_name['cd_hit_otu']]['dir']
        if cd_hit_otu_dir_status is False:
            self.page_status[p_echo_html] = False
            secho(f'Error: {analysis_dir_name["cd_hit_otu"]} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {p_echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        if len(self.analysis_data.cd_hit_otu.keys()) > 2:
            raise RuntimeError('Clustering 디렉터리가 1개 이상입니다. 개발자에게 문의하세요.')

        # Clustering 디렉터리가 1개라고 가정.
        for keys in self.analysis_data.cd_hit_otu.keys():
            if keys == self.analysis_paths.analysis_dir_name['cd_hit_otu']:
                continue
            else:
                d_clusering_files = self.analysis_data.cd_hit_otu[keys]
                cluster_dir = keys

        # otus_rep.fasta 파일 확인
        rep_fasta_name = 'otus_rep.fasta'
        if self.confirm_file(
                {
                    'file_name': rep_fasta_name,
                    'file_status': d_clusering_files['file'][rep_fasta_name],
                    'html_name': p_echo_html,
                }
        ):
            otu = os.path.join(self.analysis_paths.cd_hit_otu_path, cluster_dir, rep_fasta_name)
            self.copy_otus_rep_fasta(otu)
        else:
            return False

        # chimaric.ids 파일 확인
        chimaric_name = 'chimaric.ids'
        if self.confirm_file(
                {
                    'file_name': chimaric_name,
                    'file_status': d_clusering_files['file'][chimaric_name],
                    'html_name': p_echo_html,
                }
        ):
            chimeric = os.path.join(self.analysis_paths.cd_hit_otu_path, cluster_dir, chimaric_name)
        else:
            return False

        # trim.log 파일 확인
        trim_log_name = 'trim.log'
        if self.confirm_file(
                {
                    'file_name': trim_log_name,
                    'file_status': d_clusering_files['file'][trim_log_name],
                    'html_name': p_echo_html,
                }
        ):
            trim_log = os.path.join(self.analysis_paths.cd_hit_otu_path, cluster_dir, trim_log_name)
        else:
            return False

        # STAT.txt 파일 확인(other count 계산시 필요)
        stat_name = 'STAT.txt'
        if self.confirm_file(
                {
                    'file_name': stat_name,
                    'file_status': self.analysis_data.read_assembly[analysis_dir_name['read_assembly']]['file'][
                        stat_name],
                    'html_name': p_echo_html,
                }
        ):
            stat = os.path.join(self.analysis_paths.read_assembly_path, stat_name)
        else:
            return False

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(p_echo_html)
        if metadata_status:
            pass
        else:
            return False

        # otu_table_summary.txt 파일 확인
        summary_name = 'otu_table_summary.txt'
        try:
            biom_summary_status = self.analysis_data.biom[analysis_dir_name['biom']]['file'][summary_name]
        except KeyError:
            biom_summary_status = False
        if self.confirm_file(
                {
                    'file_name': summary_name,
                    'file_status': biom_summary_status,
                    'html_name': p_echo_html,
                }
        ):
            biom_summary = os.path.join(self.analysis_paths.biom_path, summary_name)
        else:
            return False

        d_files = {
            'otu': otu,
            'chimeric': chimeric,
            'stat': stat,
            'trim_log': trim_log,
            'otu_table_summary': biom_summary,
            'metadata': metadata_file
        }
        cutoff = cluster_dir.split('_')[-1]
        try:
            int(cutoff)
        except ValueError:
            secho('Error: Clustering Cut-off 정보가 틀립니다.', fg='red', blink=True)
            echo(f'Clustering Dir: {cluster_dir}')
            echo(f'Cut-off: {cutoff}')
        if p_related_file is True:
            return d_files, cutoff, (d_clusering_files, cluster_dir)
        else:
            return d_files, cutoff

    def make_summary_page_for_cd_hit_otu(self):
        """
        CD-HIT-OTU 프로그램의 결과를 바탕으로 summary.html 페이지를 작성합니다.
        Read_Assembly 디렉터리와 CD-HIT-OTU 디렉터리(HN01XXXXX_9?)가 필요합니다.

        - 페이지 작성에 필요한 데이터 -
        Read_Assembly: STAT.txt
        CD-HIT-OTU: 수주번호_?? Clustering 디렉터리 - Cut-off Parsing
                    otus.rep.fasta
                    chimaric.ids
                    trim.log
        Analysis_?: metadata.txt
        BIOM: otu_table_summary.txt

        :return:
        """
        echo_html = 'summary.html'
        echo(f'>>> Summary for CD-HIT-OTU: {echo_html} 생성 시작')
        d_files, cutoff = self.check_cd_hit_otu_data(echo_html)

        summary = SummaryForCdHitOtu(d_files)
        summary.read_metadata()
        summary.read_stat()
        summary.read_trim_log()
        summary.read_otu_table_summary()
        summary.count_otus()
        summary.count_chimeric_id()

        # summary.html 생성
        template = self.template_env.get_template('summary.html')
        summary_html = os.path.join(self.report_paths.otu_results_path, 'otu_table_summary.html')
        l_sample_data = [[x[0], intcomma(x[1])] for x in summary.otu_table_summary.l_sample_data]
        if summary.gamma_diversity != summary.otu_table_summary.observation_count:
            secho('Warning: gamma_diversity 와 observation_count 가 다릅니다.', fg='yellow')
            echo(f'gamma_diversity({os.path.basename(d_files["otu"])}): {summary.gamma_diversity}')
            echo(f'observation_count({os.path.basename(d_files["otu_table_summary"])}): '
                 f'{summary.otu_table_summary.observation_count}')
            gamma_diversity = f'{intcomma(summary.otu_table_summary.observation_count)} of ' \
                              f'{intcomma(summary.gamma_diversity)}'
            secho(f'--> {gamma_diversity} 로 설정됩니다.', fg='magenta')
        else:
            gamma_diversity = intcomma(summary.gamma_diversity)

        with open(summary_html, 'w') as o_summary_html:
            o_summary_html.write(
                template.render(
                    clustering_tool='CD-HIT-OTU',
                    cutoff=cutoff,
                    sample_data_list=l_sample_data,
                    sample_count=intcomma(summary.otu_table_summary.sample_count),
                    read_count=intcomma(summary.otu_table_summary.total_count),
                    gamma_diversity=gamma_diversity,
                    min=intcomma(summary.otu_table_summary.min_count),
                    max=intcomma(summary.otu_table_summary.max_count),
                    median=intcomma(summary.otu_table_summary.median_count),
                    mean=intcomma(summary.otu_table_summary.mean_count),
                    ambiguous=intcomma(summary.ambiguous_base_calls),
                    wrong_prefix_or_primers=intcomma(summary.wrong_prefix_or_primers),
                    primer_seq=summary.primer_seq,
                    low_quality=intcomma(summary.low_quality_bases),
                    chimera=intcomma(summary.chimeric_count),
                    other=intcomma(summary.other_count)
                )
            )
        self.page_status['summary.html'] = True
        secho(f'>>> Summary for CD-HIT-OTU: {echo_html} 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def make_summary_page_for_closed(self):
        """
        pooled_failures.txt : 라인 수
        pooled_otus.log: Num failures 합계
        pooled_failures.txt = pooled_otus.log (동일)

        otu_table_summary.txt: total count + failures 합계 = Read_Assembly: F_length.csv 합계

        :return:
        """
        echo('>>> Summary for CLOSED 생성 시작')
        echo_html = 'summary_for_closed.html'
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        closed_dir_status = self.analysis_data.closed[analysis_dir_name['closed']]['dir']

        if closed_dir_status is False:
            self.page_status[echo_html] = False
            secho(f'Error: {analysis_dir_name["closed"]} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        # CLOSED - uclust_ref_picked_otus 디렉터리 및 파일 확인
        picked_otus_dir_name = 'uclust_ref_picked_otus'
        d_picked_otus_dir = self.analysis_data.closed[picked_otus_dir_name]
        if d_picked_otus_dir['dir'] is True:
            picked_txt_name = f'{self.order_number}.pooled_otus.txt'
            if self.confirm_file(
                {
                    'file_name': picked_txt_name,
                    'file_status': d_picked_otus_dir['file'][picked_txt_name],
                    'html_name': echo_html,
                }
            ):
                picked_txt_file = os.path.join(self.analysis_paths.closed_otu_path,
                                               picked_otus_dir_name, picked_txt_name)
            else:
                return False

            picked_log_name = f'{self.order_number}.pooled_otus.log'
            picked_failures_name = f'{self.order_number}.pooled_failures.txt'
            picked_log_status = d_picked_otus_dir['file'][picked_log_name]
            picked_failures_status = d_picked_otus_dir['file'][picked_failures_name]
            if all([picked_log_status, picked_failures_status]):
                picked_log_file = os.path.join(self.analysis_paths.closed_otu_path,
                                               picked_otus_dir_name, picked_log_name)
                picked_failures_file = os.path.join(self.analysis_paths.closed_otu_path,
                                                    picked_otus_dir_name, picked_failures_name)
            elif picked_log_status is True:
                picked_log_file = os.path.join(self.analysis_paths.closed_otu_path,
                                               picked_otus_dir_name, picked_log_name)
                picked_failures_file = None
                secho(f'Warning: {picked_failures_name} 파일이 없습니다.', fg='yellow')
                secho(f'\t--> OTU failures 의 개수는 {picked_log_name} 파일만 확인합니다.', fg='magenta')
            elif picked_failures_status is True:
                picked_log_file = None
                picked_failures_file = os.path.join(self.analysis_paths.clsoed_otu_path,
                                                    picked_otus_dir_name, picked_failures_name)
                secho(f'Warning: {picked_log_name} 파일이 없습니다.', fg='yellow')
                secho(f'\t--> OTU failures 의 개수는 {picked_failures_name} 파일만 확인합니다.', fg='magenta')
            else:
                self.page_status[echo_html] = False
                secho(f'Error: {picked_log_name}, {picked_failures_status}가 없습니다.', fg='red', blink=True)
                secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                echo()  # 빈 줄
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {picked_otus_dir_name} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        # CLOSED - rep_set 디렉터리 및 파일 확인
        rep_set_dir_name = 'rep_set'
        d_rep_set_dir = self.analysis_data.closed[rep_set_dir_name]
        if d_rep_set_dir['dir'] is True:
            rep_fasta_name = 'otus_rep.fasta'
            if self.confirm_file(
                {
                    'file_name': rep_fasta_name,
                    'file_status': d_rep_set_dir['file'][rep_fasta_name],
                    'html_name': echo_html,
                }
            ):
                rep_fasta_file = os.path.join(self.analysis_paths.closed_otu_path, rep_set_dir_name, rep_fasta_name)
                self.copy_otus_rep_fasta(rep_fasta_file)
            else:
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {rep_set_dir_name} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        filtered_stat_name = f'{self.order_number}.F_length.csv'
        filtered_stat_file_status = self.analysis_data.read_assembly[
            analysis_dir_name['read_assembly']]['file'][filtered_stat_name]
        if self.confirm_file(
            {
                'file_name': filtered_stat_name,
                'file_status': filtered_stat_file_status,
                'html_name': echo_html,
            }
        ):
            stat_file = os.path.join(self.analysis_paths.read_assembly_path, filtered_stat_name)
        else:
            return False

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            pass
        else:
            return False

        # otu_table_summary.txt 파일 확인
        summary_name = 'otu_table_summary.txt'
        biom_summary_status = self.analysis_data.biom[analysis_dir_name['biom']]['file'][summary_name]
        if self.confirm_file(
            {
                'file_name': summary_name,
                'file_status': biom_summary_status,
                'html_name': echo_html,
            }
        ):
            biom_summary_file = os.path.join(self.analysis_paths.biom_path, summary_name)
        else:
            return False

        # Data Parsing & 정리
        summary = SummaryForClosed(
            {
                'picked_log': picked_log_file,
                'picked_failures': picked_failures_file,
                'otu': rep_fasta_file,
                'stat': stat_file,
                'otu_table_summary': biom_summary_file,
                'metadata': metadata_file

            })
        summary.read_picked_log()
        summary.read_picked_failures()
        summary.read_metadata()
        summary.read_filtered_stat()
        summary.read_otu_table_summary()
        summary.count_otus()

        # summary.html 생성
        template = self.template_env.get_template('summary_for_closed.html')
        summary_html = os.path.join(self.report_paths.otu_results_path, 'otu_table_summary.html')
        l_sample_data = [[x[0], intcomma(x[1])] for x in summary.otu_table_summary.l_sample_data]
        if summary.gamma_diversity != summary.otu_table_summary.observation_count:
            secho('Warning: gamma_diversity 와 observation_count 가 다릅니다.', fg='yellow')
            echo(f'gamma_diversity({rep_fasta_name}): {summary.gamma_diversity}')
            echo(f'observation_count({summary_name}): {summary.otu_table_summary.observation_count}')
            gamma_diversity = f'{intcomma(summary.otu_table_summary.observation_count)} of ' \
                              f'{intcomma(summary.gamma_diversity)}'
            secho(f'\t--> {gamma_diversity} 로 설정됩니다.', fg='magenta')
        else:
            gamma_diversity = intcomma(summary.gamma_diversity)
        if summary.total_failures_from_log != summary.total_failures_from_failures:
            secho(f'Warning: {picked_failures_name} 와 {picked_log_name} 의 OTU failures 의 개수가 다릅니다.', fg='yellow')
            failures = f'{intcomma(summary.total_failures_from_failures)} of ' \
                       f'{intcomma(summary.total_failures_from_log)}'
            secho(f'\t--> {failures} 로 설정됩니다.', fg='magenta')
        else:
            failures = summary.total_failures_from_log

        with open(summary_html, 'w') as o_summary_html:
            o_summary_html.write(
                template.render(
                    clustering_tool='QIIME1 & UCLUST',
                    cutoff=summary.similarity,
                    sample_data_list=l_sample_data,
                    sample_count=intcomma(summary.otu_table_summary.sample_count),
                    read_count=intcomma(summary.otu_table_summary.total_count),
                    gamma_diversity=gamma_diversity,
                    min=intcomma(summary.otu_table_summary.min_count),
                    max=intcomma(summary.otu_table_summary.max_count),
                    median=intcomma(summary.otu_table_summary.median_count),
                    mean=intcomma(summary.otu_table_summary.mean_count),
                    failures=intcomma(failures)
                )
            )
        self.page_status[echo_html] = True
        secho('>>> Summary for CLOSED 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def make_alpha_diversity_page(self):
        """
        alpha_diversity.py 실행하여 adiv.txt 를 생성하고, adiv.xlsx 엑셀파일을 작성합니다.
        생성된 데이터를 바탕으로 adiv.html 페이지를 작성합니다.

        - 페이지 작성에 필요한 데이터 - 
        BIOM: otu_table.biom
        Phylogeny: rep_phylo.tre
        Analysis_?: metadata.txt

        :return:
        """
        echo('>>> Alpha Diversity 생성 시작')
        echo_html = 'adiv.html'
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        if self.analysis_data.biom[analysis_dir_name['biom']]['dir']:
            # otu_table.biom 파일 확인
            biom_status, biom_file = self.confirm_biom(echo_html)
            if biom_status:
                # BIOM(BASE) 파일 복사
                biom_name = os.path.split(biom_file)[-1]
                target_file = os.path.join(self.report_paths.biom_path, biom_name)
                if check_file_type(target_file, 'exists'):
                    pass
                else:
                    self.copy_biom(biom_file)
            else:
                return False

            # # rep_phylo.tre 파일 확인
            # tre_status, tre_file = self.confirm_tre(echo_html)
            # if tre_status:
            #     pass
            # else:
            #     return False
            tre_file = None  # PD_whole_tree 생성시 필요

            # metadata.txt 파일 확인 & 복사
            metadata_status, metadata_file = self.confirm_metadata(echo_html)
            if metadata_status:
                self.copy_metadata(metadata_file)
            else:
                return False

            out_path = self.report_paths.diversity_index_path
            self.make_dir(out_path)
            adiv = run_alpha_diversity(biom_file, out_path, tre_file, metadata_file)
            if adiv is False:
                self.page_status[echo_html] = False
                secho('\t --> adiv.html 생성 중단', fg='yellow')
                echo()  # 빈 줄
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {analysis_dir_name["biom"]} 디렉터리가 없습니다.', fg='red', blink=True)
            secho('\t --> adiv.html 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        # HTML 생성
        template = self.template_env.get_template(echo_html)
        adiv_html = os.path.join(self.report_paths.diversity_index_path, echo_html)
        with open(adiv_html, 'w') as o_adiv_html:
            o_adiv_html.write(template.render(adiv_data=adiv))
        self.page_status['adiv.html'] = True
        secho('>>> Alpha Diversity 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def make_alpha_rarefaction_page(self):
        echo('>>> Alpha Rarefaction 생성 시작')
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        echo_html = 'rarefaction_plots.html'
        if self.analysis_data.biom[analysis_dir_name['biom']]['dir']:
            # otu_table.biom 파일 확인
            biom_status, biom_file = self.confirm_biom(echo_html)
            if biom_status:
                # BIOM(BASE) 파일 복사
                biom_name = os.path.split(biom_file)[-1]
                target_file = os.path.join(self.report_paths.biom_path, biom_name)
                if check_file_type(target_file, 'exists'):
                    pass
                else:
                    self.copy_biom(biom_file)
            else:
                return False

            # rep_phylo.tre 파일 확인
            tre_status, tre_file = self.confirm_tre(echo_html)
            if tre_status:
                pass
            else:
                return False

            # metadata.txt 파일 확인 & 복사
            metadata_status, metadata_file = self.confirm_metadata(echo_html)
            if metadata_status:
                self.copy_metadata(metadata_file)
            else:
                return False

            out_path = self.report_paths.alpha_rarefaction_path
            self.make_dir(out_path)
            run_status = run_alpha_rarefaction(
                {
                    'biom': biom_file,
                    'out_dir': out_path,
                    'metadata': metadata_file,
                    'tre': tre_file,
                    'max_rare_depth': self.kargs['max_rare_depth']
                })
            if run_status:
                html_file = os.path.join(out_path, 'alpha_rarefaction_plots', 'rarefaction_plots.html')
                if check_file_type(html_file, 'exists'):
                    html_soup = parse_html(html_file)
                else:
                    self.page_status[echo_html] = False
                    secho(f'Error: {echo_html} 파일이 없습니다.', fg='red', blink=True)
                    echo(html_file)
                    echo()  # 빈 줄
                    return False

                # HTML 생성
                template = self.template_env.get_template('rarefaction_plots.html')
                with open(html_file, 'w') as o_rarefaction_html:
                    o_rarefaction_html.write(template.render(
                        style=html_soup.style,
                        script=html_soup.script,
                        form=html_soup.form,
                        table=html_soup.table,
                        curve_data=html_soup.div,
                    ))
                self.page_status[echo_html] = True
                secho('>>> Alpha Rarefaction 생성 완료', fg='cyan')
                echo()  # 빈 줄
                return True
            else:
                self.page_status[echo_html] = False
                echo()  # 빈 줄
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {analysis_dir_name["biom"]} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

    def make_beta_diversity_page(self, p_no_3d):
        echo('>>> Beta Diversity 3D PCoA 생성 시작')
        echo_html = '(DIR) Beta Diversity 3D PCoA HTML'
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        if self.analysis_data.biom[analysis_dir_name['biom']]['dir']:
            # otu_table.biom 파일 확인
            biom_status, biom_file = self.confirm_biom(echo_html)
            if biom_status:
                # BIOM(BASE) 파일 복사
                biom_name = os.path.split(biom_file)[-1]
                target_file = os.path.join(self.report_paths.biom_path, biom_name)
                if check_file_type(target_file, 'exists'):
                    pass
                else:
                    self.copy_biom(biom_file)
            else:
                return False

            # rep_phylo.tre 파일 확인
            tre_status, tre_file = self.confirm_tre(echo_html)
            if tre_status:
                pass
            else:
                return False

            # metadata.txt 파일 확인 & 복사
            metadata_status, metadata_file = self.confirm_metadata(echo_html)
            if metadata_status:
                self.copy_metadata(metadata_file)
            else:
                return False

            run_status = run_beta_diversity_through_plots(
                biom_file, metadata_file, self.report_paths.beta_diversity_path, tre_file, p_no_3d)
            if run_status:
                unifrac_pc_file = os.path.join(self.report_paths.beta_diversity_path, 'weighted_unifrac_pc.txt')
                l_pc_data = list()
                with open(unifrac_pc_file, 'r') as o_pc:
                    num = None
                    for i in o_pc:
                        if i.startswith('Site'):
                            num = int(i.strip().split()[1])
                            continue
                        if num is None:
                            continue
                        else:
                            if num == 0:
                                break
                            else:
                                l_pc_data.append(i.strip().split())
                                num -= 1
                # HTML 생성
                template = self.template_env.get_template('beta_diversity.html')
                pc_html = os.path.join(self.report_paths.beta_diversity_path, 'weighted_unifrac_pc.html')
                with open(pc_html, 'w') as o_pc_html:
                    o_pc_html.write(template.render(pc_list=l_pc_data))
                self.page_status[echo_html] = True
                secho('>>> Beta Diversity 3D PCoA 생성 완료', fg='cyan')
                echo()  # 빈 줄
                return True
            else:
                self.page_status[echo_html] = False
                echo()  # 빈 줄
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {analysis_dir_name["biom"]} 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

    def make_2d_pcoa_plot_page(self, p_unifrac='weighted'):
        """
        make_beta_diversity_page() 실행 후, 실행.

        :param p_unifrac:
        :return:
        """
        echo(f'>>> 2D PCoA Plot({p_unifrac}) 생성 시작')
        unifrac_pc = f'{p_unifrac}_unifrac_pc.txt'
        unifrac_html = f'{p_unifrac}_unifrac_pc_2D_PCoA_plots.html'
        unifrac_pc_file = os.path.join(self.report_paths.beta_diversity_path, unifrac_pc)
        if check_file_type(unifrac_pc_file, 'exists'):
            pass
        else:
            self.page_status[unifrac_html] = False
            secho(f'Error: {unifrac_pc} 파일이 없습니다.', fg='red', blink=True)
            secho(f'\t--> {unifrac_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(unifrac_html)
        if metadata_status:
            pass
        else:
            return False

        self.make_dir(self.report_paths.pcoa_2d_path)
        run_status = run_make_2d_pcoa(unifrac_pc_file, metadata_file, self.report_paths.pcoa_2d_path)

        # HTML 생성
        if run_status:
            html_soup = parse_html(os.path.join(self.report_paths.pcoa_2d_path, unifrac_html))
            head = html_soup.head.extract()
            html_soup.html.unwrap()
            html_soup.body.unwrap()
            for tag_table in html_soup.find_all('table'):
                tag_table['class'] = 'tb1'
            template = self.template_env.get_template('unifrac_pc_2D_PCoA_plots.html')
            pcoa_2d_html = os.path.join(self.report_paths.pcoa_2d_path, unifrac_html)
            with open(pcoa_2d_html, 'w') as o_pcoa_2d_html:
                o_pcoa_2d_html.write(template.render(
                    style=head.style,
                    script=head.script,
                    unifrac=p_unifrac,
                    pcoa_data=html_soup,
                ))
            self.page_status[unifrac_html] = True
            secho(f'>>> 2D PCoA Plot({p_unifrac}) 생성 완료', fg='cyan')
            echo()  # 빈 줄
            return True
        else:
            self.page_status[unifrac_html] = False
            echo()  # 빈 줄
            return False

    def make_beta_diversity_2d_3d_plot_page(self, p_sample_count):
        if p_sample_count <= 2:
            secho(f'Warning: 시료의 개수가 적어(<=2) Beta Diversity(2D & 3D PCoA)의 생성을 생략합니다.', fg='yellow')
            self.skipped_page.append('2D_PCoA')
            self.skipped_page.append('3D_PCoA')
        elif p_sample_count == 3:
            secho(f'Warning: 시료의 개수가 적어(==3) 3D PCoA의 생성을 생략합니다.', fg='yellow')
            self.skipped_page.append('3D_PCoA')
            self.make_beta_diversity_page(p_no_3d=True)
            self.make_2d_pcoa_plot_page('weighted')
        else:
            self.make_beta_diversity_page(p_no_3d=False)
            self.make_2d_pcoa_plot_page('weighted')

    def make_upgma_tree_page(self, p_unifrac='weighted'):
        echo(f'>>> UPGMA Tree({p_unifrac}) 생성 시작')
        unifrac_dm = f'{p_unifrac}_unifrac_dm.txt'
        echo_html = 'upgma_cluster.html'
        unifrac_dm_file = os.path.join(self.report_paths.beta_diversity_path, unifrac_dm)

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            pass
        else:
            return False

        if check_file_type(unifrac_dm_file, 'exists'):
            upgma_tre = os.path.join(self.report_paths.upgma_tree_path, 'upgma_cluster.tre')
            self.make_dir(self.report_paths.upgma_tree_path)
            run_status = run_upgma_tree(unifrac_dm_file, upgma_tre)
            if run_status:
                with open(metadata_file, 'r') as o_meta:
                    for i in o_meta:
                        if i.startswith('#'):
                            l_header = i.strip().split('\t')

                l_group = list()
                for header in l_header:
                    if 'sample' in header.lower():
                        continue
                    elif 'inputfile' in header.lower():
                        continue
                    else:
                        l_group.append(header)

                # UPGMA Tree - 그룹별 색지정
                l_tree = list()
                if len(l_group) == 0:
                    tree = TreeController(metadata_file, upgma_tre, 'SampleID')
                    tree.read_metadata_file_no_group()
                    tree.read_tree_file()
                    l_tree.append(tree)
                    del tree
                else:
                    for group in l_group:
                        tree = TreeController(metadata_file, upgma_tre, group)
                        tree.read_metadata_file()
                        tree.read_tree_file()
                        l_tree.append(tree)
                        del tree

                # HTML 생성
                template = self.template_env.get_template(echo_html)
                upgma_html_file = os.path.join(self.report_paths.upgma_tree_path, echo_html)
                with open(upgma_html_file, 'w') as o_upgma_html:
                    o_upgma_html.write(template.render(
                        l_tree=l_tree,
                        new_tre_file='',
                    ))
                self.page_status[echo_html] = True
                secho(f'>>> UPGMA Tree({p_unifrac}) 생성 완료', fg='cyan')
                echo()  # 빈 줄
                return True
            else:
                self.page_status[echo_html] = False
                echo()  # 빈 줄
                return False
        else:
            self.page_status[echo_html] = False
            secho(f'Error: {unifrac_dm} 파일이 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

    def make_taxonomy_assignment_page(self, p_chart):
        """

        :type p_chart: str
        :param p_chart: ['bar' | 'bar,area']
        :return:
        """
        echo('>>> Taxonomy Assignment Plot & Excel 생성 시작')
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        echo_html = '(DIR) Taxonomy Assignment HTML & Excel'

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            metadata_order = read_metadata_for_sort(metadata_file)
            self.copy_metadata(metadata_file)
        else:
            return False

        # Taxonomy Assignment 결과 확인
        l_taxonomy_data = list()
        for dir_name in self.analysis_data.taxonomy_assignment.keys():
            if dir_name == analysis_dir_name['taxonomy_assignment']:
                continue
            else:
                d_status = self.analysis_data.taxonomy_assignment.get(dir_name)
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

        # biom 파일 확인
        l_final_taxonomy_data = list()
        for dir_name in l_taxonomy_data:
            biom_name = f'otu_table.{dir_name}.biom'
            try:
                biom_status = self.analysis_data.biom[analysis_dir_name['biom']]['file'][biom_name]
            except KeyError:
                continue
            else:
                if biom_status is False:
                    secho(f'--- {dir_name} 제외', fg='yellow')
                    if biom_status is False:
                        secho(f'\t{biom_name} 가 없습니다.', fg='red', blink=True)
                else:
                    l_final_taxonomy_data.append(dir_name)
        try:
            del dir_name, biom_name, biom_status
        except UnboundLocalError:
            pass

        # Taxonomy Assignment Data 확인
        if len(l_final_taxonomy_data) == 0:
            self.page_status[echo_html] = False
            secho('Error: 사용할 수 있는 Taxonomy Assignment 디렉터리가 없습니다.', fg='red', blink=True)
            secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
            echo()  # 빈 줄

        for dir_name in l_final_taxonomy_data:
            secho(f'+++ {dir_name} 선택', fg='yellow')
            biom_name = f'otu_table.{dir_name}.biom'
            biom_file = os.path.join(self.analysis_paths.biom_path, biom_name)
            l_level = ['2', '3', '4', '5', '6', '7']
            l_level_name = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            out_dir_path = f'{self.report_paths.taxonomy_assignment_path}_{dir_name}'
            summarize_status = run_summarize_taxa(biom_file, out_dir_path, ','.join(l_level), p_absolute=False)
            summarize_status_count = run_summarize_taxa(biom_file, os.path.join(out_dir_path, 'count_data'),
                                                        ','.join(l_level), p_absolute=True)
            if all([summarize_status, summarize_status_count]) is False:
                self.page_status[echo_html] = False
                echo()  # 빈 줄
            l_table = [os.path.join(out_dir_path, f'otu_table.{dir_name}_L{level}.txt') for level in l_level]
            max_sample_name_length = max([len(x) for x in metadata_order])
            if max_sample_name_length < 6:
                l_y_size = [6]
            elif len(metadata_order) > 30:
                l_y_size = [6]
            else:
                l_y_size = [max_sample_name_length-1, max_sample_name_length]
            l_plot_run_status = list()
            for y_size in l_y_size:
                plot_dir = os.path.join(out_dir_path, f'taxa_summary_plots_{y_size}')
                run_status = run_plot_taxa_summary(','.join(l_table), ','.join(l_level_name),
                                                   plot_dir, y_size, p_chart)
                if run_status is True:
                    l_plot_html_file = list()
                    if 'bar' in p_chart:
                        l_plot_html_file.append(
                            (os.path.join(plot_dir, 'bar_charts.html'), os.path.join(plot_dir, 'new_bar_charts.html'))
                        )
                    if 'area' in p_chart:
                        l_plot_html_file.append(
                            (os.path.join(plot_dir, 'area_charts.html'), os.path.join(plot_dir, 'new_area_charts.html'))
                        )
                    for cur, new in l_plot_html_file:
                        plot_html = TaxonHTML(cur, new)
                        plot_html.transform()
                        plot_html.save()
                        plot_html.del_html()
                        plot_html.move_new_html()
                    secho(f'>>> Taxa Summary Plot {y_size}: HTML 변환 완료', fg='cyan')
                l_plot_run_status.append(run_status)

            # Excel File 생성
            if all(l_plot_run_status):
                if 'NCBI' in dir_name:
                    tool, _, db = dir_name.partition('_')
                elif 'UNITE' == dir_name:
                    tool, db = dir_name.split('_')
                elif 'UNITE_' in dir_name:
                    tool, _, db = dir_name.partition('_')
                elif 'GreenGenes' in dir_name:
                    tool, _, db = dir_name.partition('_')
                elif 'MIDORI' in dir_name:
                    tool, _, db = dir_name.partition('_')
                else:
                    tool, db = dir_name.split('_')

                # BLAST - Excel 복사
                if tool == 'blast':
                    blast_excel_file = os.path.join(self.analysis_paths.taxonomy_assignment_path, dir_name,
                                                    'OTU_BLAST.xlsx')
                    blast_copy_cmd = f'cp {blast_excel_file} {out_dir_path}'
                    run = run_cmd(blast_copy_cmd)
                    check_run_cmd(
                        {
                            'run': run,
                            'true_meg': None,
                            'false_meg': 'BLAST Excel 복사',
                        }, p_stdout=False, p_exit=False,
                    )
                    if run.returncode == 0:
                        pass
                    else:
                        secho('Error: BLAST Excel 파일 복사에 문제가 있습니다.', fg='red', blink=True)
                        secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                        echo()  # 빈 줄
                        return False

                # Taxonomy Assignment - Excel: count & ratio
                for data_type in ['count', 'ratio']:
                    if data_type == 'count':
                        out_table_dir_path = os.path.join(out_dir_path, 'count_data')
                    elif data_type == 'ratio':
                        out_table_dir_path = out_dir_path
                    try:
                        i_excel = TaxonomyAbundanceExcel(
                            {
                                'db_type': db,
                                'level': 'Species',
                                'otu_table_dir_path': out_table_dir_path,
                                'data_type': data_type,
                                'sort_type': 'Custom',
                                'custom_order': metadata_order,
                            }
                        )
                        i_excel.glob_otu_table_file()
                        i_excel.read_otu_table()
                        i_excel.save_excel()
                    except Exception as err:
                        self.page_status[f'{echo_html}: {dir_name}'] = False
                        secho('Error: Taxonomy Excel 생성에 문제가 있습니다.', fg='red', blink=True)
                        echo(f'메시지: {err}')
                        secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                        echo()  # 빈 줄
                        return False

                # OTU Table - Shared & Taxonomy Assignment Excel
                mothur_work_path = os.path.join(self.report_paths.report_dir_path, 'mothur_shared')
                self.make_dir(mothur_work_path)
                shared_file_status = run_make_shared(mothur_work_path, biom_file)
                if shared_file_status:
                    l_biom_file = os.path.splitext(os.path.split(biom_file)[1])  # 확장자 분리
                    shared_file = os.path.join(mothur_work_path, f'{l_biom_file[0]}.shared')
                else:
                    self.page_status[f'{echo_html}: {dir_name}'] = False
                    secho('Error: shared file 생성에 문제가 있습니다.', fg='red', blink=True)
                    secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                    echo()  # 빈 줄
                    return False

                for file_name, status in self.analysis_data.taxonomy_assignment.get(dir_name).get('file').items():
                    if status:
                        if file_name.endswith('.txt'):
                            tax_txt_file = os.path.join(self.analysis_paths.taxonomy_assignment_path,
                                                        dir_name, file_name)
                        elif file_name.endswith('.log'):
                            tax_log_file = os.path.join(self.analysis_paths.taxonomy_assignment_path,
                                                        dir_name, file_name)
                    else:
                        secho(f'Error: {file_name} 파일이 없습니다.', fg='red', blink=True)
                        secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                        echo()  # 빈 줄
                        return False
                else:
                    del file_name, status

                try:
                    i_shared_excel = TaxonomySharedExcel(
                        {
                            'level': 'Species',
                            'db': db,
                            'sort': 'Custom',
                            'custom_order': metadata_order,
                            'method': self.otu_picking_method,
                            'shared_file': shared_file,
                            'tax_txt_file': tax_txt_file,
                            'tax_log_file': tax_log_file,
                            'out_path': out_dir_path,
                        }
                    )
                    i_shared_excel.read_shared_file()
                    i_shared_excel.read_tax_assignment_txt()
                    i_shared_excel.read_tax_assignment_log()
                    i_shared_excel.save_excel()
                except Exception as err:
                    self.page_status[f'{echo_html}: {dir_name}'] = False
                    secho('Error: Shared Excel 생성에 문제가 있습니다.', fg='red', blink=True)
                    echo(f'메시지: {err.with_traceback(None)}')
                    secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                    echo()  # 빈 줄
                    return False
                else:
                    self.delete_dir(mothur_work_path)

                # taxonomy.biom 복사
                if self.copy_biom(biom_file) is False:
                    self.page_status[f'{echo_html}: {dir_name}'] = False
                    secho('Error: biom file 복사에 문제가 있습니다.', fg='red', blink=True)
                    echo(f'file: {biom_file}')
                    secho(f'\t--> {echo_html} 생성 중단', fg='yellow')
                    echo()  # 빈 줄
                    return False

                self.page_status[f'{echo_html}: {dir_name}'] = True
                secho('>>> Taxonomy Assignment Plot & Excel 생성 완료', fg='cyan')
                echo()  # 빈 줄 
            else:
                self.page_status[f'{echo_html}: {dir_name}'] = False
                return False

    def copy_metadata(self, p_metadata):
        """
        metadata.txt 파일을 biom 디렉터리에 복사한다.
        BIOM 디렉터리에 metadata.txt 파일이 존재할 경우 생략한다.

        :param p_metadata: meatadata.txt 파일
        :return:
        """
        if check_file_type(self.report_paths.biom_path, 'exists'):
            pass
        else:
            self.make_dir(self.report_paths.biom_path)

        biom_metadata_file = os.path.join(self.report_paths.biom_path, 'metadata.txt')
        if check_file_type(biom_metadata_file, 'exists'):
            secho('BIOM 디렉터리에 metadata.txt 파일이 존재합니다.', fg='yellow')
            secho('\t--> Metadata 복사를 생략합니다.', fg='magenta')
            return True
        else:
            cmd = f'cp {p_metadata} {self.report_paths.biom_path}'
            run = run_cmd(cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': 'Metadata 파일 복사 완료',
                    'false_meg': 'copy metadata',
                }, p_stdout=False, p_exit=False,
            )
            if run.returncode == 0:
                return True
            else:
                secho('Error: Metadata.txt 파일 복사 유무를 확인하세요,', fg='red', blink=True)
                secho(f'명령어: {cmd}', fg='magenta')
                return False

    def copy_biom(self, p_biom):
        """
        otu_table*.biom 파일을 분석 보고서의 지정된 위치(BIOM 디렉터리)에 복사한다.
        분석 보고서의 버전에 따라 위치가 달라진다.

        :param p_biom:
        :return:
        """
        if check_file_type(self.report_paths.biom_path, 'exists'):
            pass
        else:
            self.make_dir(self.report_paths.biom_path)
        cmd = f'cp {p_biom} {self.report_paths.biom_path}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'BIOM 파일 복사 완료',
                'false_meg': 'copy biom'
            }, p_stdout=False, p_exit=False,
        )
        if run.returncode == 0:
            return True
        else:
            secho('Error: BIOM 파일 복사 유무를 확인하세요.', fg='red', blink=True)
            secho(f'명령어: {cmd}', fg='magenta')
            return False

    def copy_src(self):
        """
        src 디렉터리를 분석 보고서의 지정된 위치에 복사한다.

        :return:
        """
        if check_file_type(self.report_paths.src_path, 'exists'):
            echo('>>> src 복사 생략(디렉터리 존재)')
        else:
            src_dir_path = os.path.join(self.html_template_path, 'src')
            cmd = f'cp -r {src_dir_path} {self.report_paths.otu_results_path}'
            run = run_cmd(cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': 'src dir 복사 완료',
                    'false_meg': 'copy src dir'
                }, p_stdout=False, p_exit=False,
            )
            if run.returncode == 0:
                return True
            else:
                secho('Error: src 디렉터리 복사 유무를 확인하세요.', fg='red', blink=True)
                secho(f'명령어: {cmd}', fg='magenta')
                return False

    def copy_otus_rep_fasta(self, p_file):
        """
        otus_rep.fasta 파일을 보고서 생성 디렉터리의 지정된 위치(OTU_Resluts)에 복사한다.
        분석 보고서의 버전에 따라 위치가 달라진다.

        :param p_file: otus_rep.fasta 파일
        :return:
        """
        fasta_name = os.path.split(p_file)[-1]
        fasta_file = os.path.join(self.report_paths.otu_results_path, fasta_name)
        if check_file_type(fasta_file, 'exists'):
            echo('>>> otus_rep.fasta 복사 생략(파일 존재)')
        else:
            cmd = f'cp {p_file} {fasta_file}'
            run = run_cmd(cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': 'otus_rep.fasta 복사 완료',
                    'false_meg': 'cp otus_rep.fasta',
                }, p_stdout=False, p_exit=False,
            )
            if run.returncode == 0:
                return True
            else:
                secho('Error: src 디렉터리 복사 유무를 확인하세요.', fg='red', blink=True)
                secho(f'명령어: {cmd}', fg='magenta')
                return False
