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
from click import secho, echo, style
from SpoON.util import glob_dir, check_file_type, parse_html, read_metadata_for_sort, read_data, \
    sort_data_by_custom, parse_config
from theCups.data_structure import ReportPathV2
from theCups.report_miseq_v1 import MiSeqReportV1
from theCups.read_assembly import transpose_data
from theCups.clustering_and_otu_picking import SummaryForCdHitOtu, SummaryForClosed
from theCups.diversity_index import run_alpha_diversity
from theCups.rarefaction import run_alpha_rarefaction
from theCups.beta_diversity import run_beta_diversity_through_plots
from theCups.pcoa import run_make_2d_pcoa
from theCups.upgma_tree import run_upgma_tree, TreeController
from theCups.taxonomy import run_summarize_taxa, run_plot_taxa_summary, TaxonHTML, \
    TaxonomyAbundanceExcel, TaxonomySharedExcel, run_make_shared
from jinja2 import FileSystemLoader
from jinja2.environment import Environment
from humanize import intcomma

CONFIG = parse_config()
MiSeq_V2_TEMPLATE_PATH = CONFIG['theCups']['MiSeq_V2']['template_path']


class MiSeqReportV2(MiSeqReportV1):
    def __init__(self, kargs):
        super().__init__(kargs)
        self.report_paths = ReportPathV2(self.base_path, self.order_number, self.analysis_number)

    @property
    def html_template_path(self):
        global MiSeq_V2_TEMPLATE_PATH
        return MiSeq_V2_TEMPLATE_PATH

    def checklist_read_assembly_files(self, p_type):
        """
        Report.checklist_read_assembly_files : overriding.

        :param p_type: 서열 파일 형식. fastq or fasta
        :return: 확인해야 할 파일명을 원소로 가지는 리스트
        """
        l_files = super().checklist_read_assembly_files(p_type)
        l_additional_files = [
            f'{self.order_number}.FLASH.csv',
            f'{self.order_number}.hist.csv',
        ]
        l_files.extend(l_additional_files)
        return l_files

    def checklist_cd_hit_otu_files(self):
        """
        Report.checklist_read_assembly_files : overriding

        :return: 확인해야 할 파일명을 원소로 가지는 리스트
        """
        l_files = super().checklist_cd_hit_otu_files()
        l_files.append('OTU.nr2nd.clstr.NAT.txt')
        return l_files

    def taste_cofi(self):
        # TODO ; overriding
        return

    def make_fastq_page(self):
        return

    def make_read_assembly_page(self):
        echo_html = 'Read_Assembly.html'
        echo(f'>>> Read Assembly: {echo_html} 생성 시작')
        analysis_dir_name = self.analysis_paths.analysis_dir_name

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            metadata_order = read_metadata_for_sort(metadata_file)
        else:
            return False

        # FLASH.csv 파일 확인
        # 파일 구조
        # ------------------
        # Sample
        # total_pairs
        # combined_pairs
        # uncombined_pairs
        # percent_combined
        # target_size
        # min_overlap
        # max_overlap
        # max_mismatch_densit
        # allow_outie_pairs
        # cap_mismatch_quals
        # ------------------
        flash_csv_name = f'{self.order_number}.FLASH.csv'
        flash_csv_file_status = self.analysis_data.read_assembly[
            analysis_dir_name['read_assembly']]['file'][flash_csv_name]
        if self.confirm_file(
            {
                'file_name': flash_csv_name,
                'file_status': flash_csv_file_status,
                'html_name': echo_html,
            }
        ):
            flash_csv_file = os.path.join(self.analysis_paths.read_assembly_path, flash_csv_name)
            l_flash_csv = read_data(flash_csv_file, p_sep=',')
            l_trans_flash_stat = transpose_data(l_flash_csv[:4])
            # Bar Chart Data - FLASH STAT
            l_sorted_flash_stat = [l_trans_flash_stat[0]]  # header
            l_sorted_flash_stat.extend(sort_data_by_custom(l_trans_flash_stat[1:], metadata_order))
            del l_trans_flash_stat
            # Bar Chart Data - FLASH Percentage
            l_percent_data = [x.replace('%', '') for x in l_flash_csv[4]]
            l_trans_percent = transpose_data([l_flash_csv[0], l_percent_data])
            l_sorted_percent = [l_trans_percent[0]]
            l_sorted_percent.extend(sort_data_by_custom(l_trans_percent[1:], metadata_order))
            del l_percent_data, l_trans_percent
            # Table Data - FLASH All
            l_trans_flash_all = transpose_data(l_flash_csv)
            l_sorted_flash_all = [l_trans_flash_all[0]]  # header
            l_sorted_flash_all.extend(sort_data_by_custom(l_trans_flash_all[1:], metadata_order))
        else:
            return False

        # hist.csv 파일 확인
        hist_csv_name = f'{self.order_number}.hist.csv'
        hist_csv_file_status = self.analysis_data.read_assembly[
            analysis_dir_name['read_assembly']]['file'][hist_csv_name]
        if self.confirm_file(
            {
                'file_name': hist_csv_name,
                'file_status': hist_csv_file_status,
                'html_name': echo_html,
            }
        ):
            # Line Chart Data - Length Distribution
            hist_csv_file = os.path.join(self.analysis_paths.read_assembly_path, hist_csv_name)
            l_hist_csv = read_data(hist_csv_file, p_sep=',')
            l_trans_hist_csv = transpose_data(l_hist_csv)
            l_sorted_trans_hist_csv = [l_trans_hist_csv[0]]  # header
            l_sorted_trans_hist_csv.extend(sort_data_by_custom(l_trans_hist_csv[1:], metadata_order))
            del l_hist_csv
            l_sorted_hist_csv = transpose_data(l_sorted_trans_hist_csv)
            l_sorted_sum = list()
            for l_length in l_sorted_trans_hist_csv[1:]:
                sum_data = sum([int(x) for x in l_length[1:]])
                l_sorted_sum.append([l_length[0], sum_data])
            l_sorted_trans_hist_percentage_data = list()
            l_sorted_trans_hist_percentage_data.append(l_sorted_trans_hist_csv[0])
            for index, l_sum in enumerate(l_sorted_sum, 1):
                if l_sum[0] == l_sorted_trans_hist_csv[index][0]:  # 시료명 확인
                    percentage = [(int(x)/l_sum[1])*100 for x in l_sorted_trans_hist_csv[index][1:]]
                    percentage.insert(0, l_sum[0])
                    l_sorted_trans_hist_percentage_data.append(percentage)
                else:
                    raise RuntimeError('개발자에게 문의하세요.')
            l_hist_header = [f"'{x}'" for x in l_sorted_hist_csv[0]]
            l_hist_data = l_sorted_hist_csv[1:]
            l_hist_percentage_data = transpose_data(l_sorted_trans_hist_percentage_data)[1:]
        else:
            return False

        # Read_Assembly.html 생성
        template = self.template_env.get_template('read_assembly.html')
        read_assembly_html = os.path.join(self.report_paths.otu_results_path, echo_html)
        with open(read_assembly_html, 'w') as o_read_assembly_html:
            o_read_assembly_html.write(
                template.render(
                    assembly_stat=l_sorted_flash_stat,
                    assembly_percent=l_sorted_percent,
                    hist_header=l_hist_header,
                    hist_count=l_hist_data,
                    hist_percentage=l_hist_percentage_data,
                    assembly_all=l_sorted_flash_all
                ))
        self.page_status[echo_html] = True
        secho(f'>>> Read Assembly: {echo_html} 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def make_length_filter_page(self):
        echo_html = 'Length_Filter.html'
        echo(f'>>> Length Filter: {echo_html} 생성 시작')
        analysis_dir_name = self.analysis_paths.analysis_dir_name

        # metadata.txt 파일 확인
        metadata_status, metadata_file = self.confirm_metadata(echo_html)
        if metadata_status:
            metadata_order = read_metadata_for_sort(metadata_file)
        else:
            return False

        # F_lenght.csv 파일 확인
        # 파일 구조
        # ---------------------
        # Name
        # Total_Read
        # Removed_Read
        # Removed_Read(Lower)
        # Removed_Read(Upper)
        # Final_Read
        # Removed_Rate
        # lower_bp
        # upper_bp
        # ---------------------
        length_csv_name = f'{self.order_number}.F_length.csv'
        length_csv_file_status = self.analysis_data.read_assembly[
            analysis_dir_name['read_assembly']]['file'][length_csv_name]
        if self.confirm_biom(
            {
                'file_name': length_csv_name,
                'file_status': length_csv_file_status,
                'html_name': echo_html,
            }
        ):
            length_csv_file = os.path.join(self.analysis_paths.read_assembly_path, length_csv_name)
            l_length_csv = read_data(length_csv_file, p_sep=',')
            l_length_csv[0] = [x.replace('.extendedFrags.fastq.gz', '') for x in l_length_csv[0]]
            l_trans_table = transpose_data(l_length_csv)
            # Table Data
            l_sorted_table = [l_trans_table[0]]  # header
            l_sorted_table.extend(sort_data_by_custom(l_trans_table[1:], metadata_order))
            # Bar Chart Data - STAT
            l_trans_stat = transpose_data(l_length_csv[0:3] + [l_length_csv[5]])
            l_sorted_stat = [l_trans_stat[0]]  # header
            l_sorted_stat.extend(sort_data_by_custom(l_trans_stat[1:], metadata_order))
            del l_trans_table, l_trans_stat
            # Bar Chart Data - percentage
            l_trans_percentage = transpose_data([l_length_csv[0]] + [l_length_csv[6]])
            l_sorted_percentage = [l_trans_percentage[0]]  # header
            l_sorted_percentage.extend(sort_data_by_custom(l_trans_percentage[1:], metadata_order))
            del l_length_csv, l_trans_percentage
        else:
            return False

        # Read_Assembly.html 생성
        template = self.template_env.get_template('length_filter.html')
        length_filter_html = os.path.join(self.report_paths.otu_results_path, echo_html)
        with open(length_filter_html, 'w') as o_length_filter_html:
            o_length_filter_html.write(
                template.render(
                    filter_stat=l_sorted_stat,
                    filter_percentage=l_sorted_percentage,
                    filter_table=l_sorted_table,
                ))
        self.page_status[echo_html] = True
        secho(f'>>> Length Filter: {echo_html} 생성 완료', fg='cyan')
        echo()
        return True

    def check_cd_hit_otu_data(self, p_echo_html):
        d_files, cutoff, (d_clustering_files, cluster_dir) = super().check_cd_hit_otu_data(p_echo_html,
                                                                                           p_related_file=True)

        # OTU.nr2nd.clstr.NAT.txt 확인 (MiSeq_V2)
        otu_nat_name = 'OTU.nr2nd.clstr.NAT.txt'
        if self.confirm_file(
                {
                    'file_name': otu_nat_name,
                    'file_status': d_clustering_files['file'][otu_nat_name],
                    'html_name': p_echo_html,
                }
        ):
            otu_nat = os.path.join(self.analysis_paths.cd_hit_otu_path, cluster_dir, otu_nat_name)
        else:
            return False
        d_files['otu_nat'] = otu_nat
        return d_files, cutoff

    def make_cd_hit_otu_page(self):
        echo_html = 'Clustering_OTU_Picking.html'
        echo(f'>>> Clustering & OTU Picking: {echo_html} 생성 시작')
        d_files, cutoff = self.check_cd_hit_otu_data(echo_html)
        clustering = SummaryForCdHitOtu(d_files)
        clustering.read_metadata()
        clustering.read_stat()
        clustering.read_trim_log()
        clustering.read_otu_table_summary()
        clustering.count_otus()
        clustering.count_chimeric_id()
        clustering.read_otu_nat()

        nat_header = [f'"{x}"' for x in clustering.otu_nat[0]]
        nat_data = clustering.otu_nat[1:]

        # Clustering_OTU_Picking.html 생성
        template = self.template_env.get_template('clustering_otu_picking.html')
        clustering_html = os.path.join(self.report_paths.otu_results_path, echo_html)
        # l_sample_data = [[x[0], intcomma(x[1])] for x in clustering.otu_table_summary.l_sample_data]

        if clustering.gamma_diversity != clustering.otu_table_summary.observation_count:
            secho('Warning: gamma_diversity 와 observation_count 가 다릅니다.', fg='yellow')
            echo(f'gamma_diversity({os.path.basename(d_files["otu"])}): {clustering.gamma_diversity}')
            echo(f'observation_count({os.path.basename(d_files["biom_summary"])}): '
                 f'{clustering.otu_table_summary.observation_count}')
            gamma_diversity = f'{intcomma(clustering.otu_table_summary.observation_count)} of ' \
                              f'{intcomma(clustering.gamma_diversity)}'
            secho(f'--> {gamma_diversity} 로 설정됩니다.', fg='magenta')
        else:
            gamma_diversity = intcomma(clustering.gamma_diversity)

        with open(clustering_html, 'w') as o_summary_html:
            o_summary_html.write(
                template.render(
                    clustering_tool='CD-HIT-OTU',
                    cutoff=cutoff,
                    sample_data_list=clustering.otu_table_summary.l_sample_data,
                    nat_header_list=nat_header,
                    nat_data_list=nat_data,
                    sample_count=intcomma(clustering.otu_table_summary.sample_count),
                    read_count=intcomma(clustering.otu_table_summary.total_count),
                    gamma_diversity=gamma_diversity,
                    min=intcomma(clustering.otu_table_summary.min_count),
                    max=intcomma(clustering.otu_table_summary.max_count),
                    median=intcomma(clustering.otu_table_summary.median_count),
                    mean=intcomma(clustering.otu_table_summary.mean_count),
                    ambiguous=intcomma(clustering.ambiguous_base_calls),
                    wrong_prefix_or_primers=intcomma(clustering.wrong_prefix_or_primers),
                    primer_seq=clustering.primer_seq,
                    low_quality=intcomma(clustering.low_quality_bases),
                    chimera=intcomma(clustering.chimeric_count),
                    other=intcomma(clustering.other_count)
                )
            )
        self.page_status[echo_html] = True
        secho(f'>>> Length Filter: {echo_html} 생성 완료', fg='cyan')
        echo()  # 빈 줄
        return True

    def summary_page(self):
        pass
