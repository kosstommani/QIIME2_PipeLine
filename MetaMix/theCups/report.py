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
__version__ = '1.1.0'

# ------------------------------------------
# Ver. 1.1.0 : 2020.06.04
# R_DADA2 추가

import os
import shutil
from click import secho, echo
from SpoON.util import glob_dir, check_file_type, parse_html, read_metadata_for_sort
from theCups.data_structure import Path, AnalysisDataTable
from jinja2 import FileSystemLoader
from jinja2.environment import Environment


class Report:
    def __init__(self, p_kargs, p_platform):
        """

        :param p_kargs: theCups.py 의 Argument & Option

        """
        self.platform = p_platform
        self.kargs = p_kargs
        self.base_path = p_kargs['analysis_base_path']
        self.order_number = p_kargs['order_number']
        self.analysis_number = p_kargs['analysis_number']
        self.analysis_paths = Path(self.base_path, self.order_number, self.analysis_number)
        self.analysis_data = None
        self.template_env = None
        self.page_status = dict()
        self.skipped_page = list()
        self.otu_picking_method = None

    def make_analysis_data_table(self, kargs):
        """

        :param kargs: 딕션너리
                    set_style:
                    column_headers: list()
        :return:
        """
        self.analysis_data = AnalysisDataTable(kargs)

    def get_analysis_data_table(self, kargs):
        return AnalysisDataTable(kargs)

    @staticmethod
    def exists(p_target):
        """
        디렉터리 또는 파일의 존재 유무를 확인하여 관련 값을 반환한다.

        :type p_target: str
        :param p_target: 존재 유무를 확인할 디렉터리
        :return Tuple - (디렉터리명 or 파일명, 존재유무)
        """
        state = os.path.exists(p_target)
        dir_name = os.path.basename(p_target)
        return dir_name, state

    def check_something(self, p_target):
        """
        특정 디렉터리 또는 파일의 존재유무를 확인하여,
        입력받은 리스트와 존재유무에 대한 bool 값을 원소로 가지는 리스트를 반환한다.

        :type p_target: list
        :param p_target: 존재유무를 확인할 디렉터리 경로 또는 파일의 경로를 원소로 가지는 리스트
            :rtype: list
        :return: 문자열을 원소로 가지는 리스트, bool 을 원소로 가지는 리스트
        """
        l_checked_something = list()
        l_states = list()
        for ele in p_target:
            name, state = self.exists(ele)
            l_checked_something.append(name)
            l_states.append(state)
        return l_checked_something, l_states

    def checklist_read_assembly_files(self, p_type):
        raise RuntimeError('Method checklist_read_assembly_files() 을 overriding 해야됩니다.')

    def check_read_assembly_dir(self, p_type='fastq'):
        l_dirs = [self.analysis_paths.read_assembly_path]
        l_files = self.checklist_read_assembly_files(p_type)
        l_path_files = map(lambda x: os.path.join(l_dirs[0], x), l_files)
        l_dirs_text, l_dirs_state = self.check_something(l_dirs)
        l_files_text, l_files_state = self.check_something(l_path_files)
        self.analysis_data.set_read_assembly_data(l_dirs_text, l_dirs_state, l_files_text, l_files_state)
        return

    def checklist_cd_hit_otu_files(self):
        raise RuntimeError('Method checklist_cd_hit_otu_files() 을 overriding 해야됩니다.')

    def check_cd_hit_otu_dir(self):
        # Clustering DIR 이 여러 개일 경우
        # l_dirs = list()
        # l_dirs.append(self.analysis_paths.cd_hit_otu_path)
        # l_cd_hit_list = glob_dir(self.analysis_paths.cd_hit_otu_path, f'{self.order_number}*',
        #                          p_mode='many', p_verbose=False)
        # clustering_dir, _ = check_file_type(l_cd_hit_list, 'isdir')
        # if len(clustering_dir) == 1:
        #     l_dirs.extend(clustering_dir)
        # elif len(clustering_dir) > 1:
        #     l_dirs.extend(clustering_dir)

        # l_files = ['trim.log']
        # if len(l_dirs) == 2:
        #     l_path_files = map(lambda x: os.path.join(l_dirs[1], x), l_files)
        # if len(l_dirs) > 2:
        #     l_path_files = list()
        #     for ele in l_dirs[1:]:
        #         l_path_files.extend(map(lambda x: os.path.join(ele, x), l_files))
        # p_o_table.append_row(self.__check_dir(l_dirs, l_path_files))

        # Clustering DIR 이 한 개일 경우
        cd_hit_dir_path = self.analysis_paths.cd_hit_otu_path
        clustering_dir_path = glob_dir(cd_hit_dir_path, f'{self.order_number}*[!.log]',
                                       p_mode='only', p_verbose=False)
        l_files = self.checklist_cd_hit_otu_files()
        l_path_files = map(lambda x: os.path.join(clustering_dir_path, x), l_files)
        if clustering_dir_path is None:
            clustering_dir_path = 'None'
        dirs_text, dirs_state = self.check_something([cd_hit_dir_path, clustering_dir_path])
        files_text, files_state = self.check_something(l_path_files)
        self.analysis_data.set_cd_hit_otu_data(dirs_text, dirs_state, [files_text], [files_state])
        return

    def check_closed_otu_dir(self):
        closed_dir_path = self.analysis_paths.closed_otu_path
        picked_otus_dir_path = os.path.join(closed_dir_path, 'uclust_ref_picked_otus')
        rep_set_dir_path = os.path.join(closed_dir_path, 'rep_set')
        # assigned_taxonomy_dir_path = os.path.join(closed_dir_path, 'uclust_assigned_taxonomy')
        picked_txt_file = f'{self.order_number}.pooled_otus.txt'
        picked_log_file = f'{self.order_number}.pooled_otus.log'
        picked_failures_file = f'{self.order_number}.pooled_failures.txt'
        rep_set_file = f'otus_rep.fasta'
        # assigned_file = f'{self.order_number}.pooled_rep_set_tax_assignments.txt'
        dirs_text, dirs_state = self.check_something([closed_dir_path, picked_otus_dir_path, rep_set_dir_path])
        l_target_files = [[os.path.join(picked_otus_dir_path, picked_txt_file),
                           os.path.join(picked_otus_dir_path, picked_log_file),
                           os.path.join(picked_otus_dir_path, picked_failures_file)],
                          [os.path.join(rep_set_dir_path, rep_set_file)]]
        l_files = list()
        l_files_state = list()
        for files in l_target_files:
            files_text, files_state = self.check_something(files)
            l_files.append(files_text)
            l_files_state.append(files_state)
        self.analysis_data.set_closed_data(dirs_text, dirs_state, l_files, l_files_state)
        return

    def check_r_dada2_dir(self):
        r_dada2_dir_path = self.analysis_paths.r_dada2_path
        fasta_file = 'all_ASVs.fasta'
        table_file = 'all_ASVs.tsv'
        rds_file = 'all_seqtab_noChimera.rds'
        dir_text, dir_state = self.check_something([r_dada2_dir_path])
        l_target_files = [[os.path.join(r_dada2_dir_path, fasta_file),
                           os.path.join(r_dada2_dir_path, table_file),
                           os.path.join(r_dada2_dir_path, rds_file)]]
        l_files = list()
        l_files_state = list()
        for file in l_target_files:
            file_text, file_state = self.check_something(file)
            l_files.append(file_text)
            l_files_state.append(file_state)
        self.analysis_data.set_r_dada2_data(dir_text, dir_state, l_files, l_files_state)
        return

    def check_r_dada2_summary_dir(self):
        r_dada2_summary_dir_path = self.analysis_paths.r_dada2_summary_path
        summary_file = os.path.join(r_dada2_summary_dir_path, 'all_dada2.summary')
        dir_text, dir_state = self.check_something([r_dada2_summary_dir_path])
        file_text, file_state = self.check_something([summary_file])
        self.analysis_data.set_r_dada2_summary_data(dir_text, dir_state, file_text, file_state)
        return 

    def check_alignment_dir(self):
        alignment_dir_path = self.analysis_paths.alignment_path
        filtered_dir_path = os.path.join(alignment_dir_path, 'filtered_alignment')
        filtered_fasta = os.path.join(filtered_dir_path, 'otus_rep_aligned_pfiltered.fasta')
        dirs_text, dirs_state = self.check_something([alignment_dir_path, filtered_dir_path])
        files_text, files_state = self.check_something([filtered_fasta])
        self.analysis_data.set_alignment_data(dirs_text, dirs_state, [files_text], [files_state])
        return

    def check_phylogeny_dir(self):
        phylogeny_dir_path = self.analysis_paths.phylogeny_path
        rep_tre_file = os.path.join(phylogeny_dir_path, 'rep_phylo.tre')
        dirs_text, dirs_state = self.check_something([phylogeny_dir_path])
        files_text, files_state = self.check_something([rep_tre_file])
        self.analysis_data.set_phylogeny(dirs_text, dirs_state, files_text, files_state)
        return

    def check_taxonomy_assignment_dir(self):
        taxa_assign_dir_path = self.analysis_paths.taxonomy_assignment_path
        l_assigned_work_list = glob_dir(taxa_assign_dir_path, '*', p_mode='many', p_verbose=False)
        if l_assigned_work_list is None:
            l_assigned_work_dir = ['None']
        else:
            l_assigned_work_dir, _ = check_file_type(l_assigned_work_list, 'isdir')

        l_files = list()
        l_files_state = list()
        for work_dir in l_assigned_work_dir:
            l_assigned_work_file = list()
            if os.path.basename(work_dir).startswith('blast'):
                l_assigned_work_file.append(os.path.join(work_dir, 'OTU_BLAST.xlsx'))
            l_assigned_work_file.append(os.path.join(work_dir, 'otus_rep_tax_assignments.txt'))
            l_assigned_work_file.append(os.path.join(work_dir, 'otus_rep_tax_assignments.log'))
            files_text, files_state = self.check_something(l_assigned_work_file)
            l_files.append(files_text)
            l_files_state.append(files_state)
        l_dir_path = list()
        l_dir_path.append(taxa_assign_dir_path)
        l_dir_path.extend(l_assigned_work_dir)
        dirs_text, dirs_state = self.check_something(l_dir_path)
        self.analysis_data.set_taxonomy_assignment_data(dirs_text, dirs_state, l_files, l_files_state)

    def check_biom_dir(self):
        """
        self.check_taxonomy_assignment_dir method 를 실행한 이후에 실행해야 한다.

        :return:
        """
        biom_dir_path = self.analysis_paths.biom_path
        dirs_text, dirs_state = self.check_something([biom_dir_path])
        if all(dirs_state):
            l_files = [os.path.join(biom_dir_path, 'otu_table.biom'),
                       os.path.join(biom_dir_path, 'otu_table_summary.txt')]
            try:
                taxonomy_assignment_keys = self.analysis_data.taxonomy_assignment.keys()
            except KeyError:
                secho('Error: check_taxonomy_assignment_dir 메소드를 먼저 실행해주세요.', fg='red', blink=True)
                exit()
            for dir_name in taxonomy_assignment_keys:
                if dir_name == self.analysis_paths.analysis_dir_name['taxonomy_assignment']:
                    continue
                else:
                    biom_name = f'otu_table.{dir_name}.biom'
                    summary_name = f'otu_table.{dir_name}_summary.txt'
                    l_glob_biom = glob_dir(biom_dir_path, biom_name, p_mode='many', p_verbose=False)
                    l_glob_summary = glob_dir(biom_dir_path, summary_name, p_mode='many', p_verbose=False)
                    if l_glob_biom is not None:
                        l_files.extend(l_glob_biom)
                    else:
                        l_files.append(os.path.join(biom_dir_path, biom_name))
                    if l_glob_summary is not None:
                        l_files.extend(l_glob_summary)
                    else:
                        l_files.append(os.path.join(biom_dir_path, summary_name))
            files_text, files_state = self.check_something(l_files)
            self.analysis_data.set_biom_data(dirs_text, dirs_state, files_text, files_state)
        else:
            self.analysis_data.set_biom_data(dirs_text, dirs_state, ['None'], [False])
        return

    def check_info_file(self):
        info_file = os.path.join(self.analysis_paths.analysis_number_path, 'info.txt')
        dirs_text, dirs_state = self.check_something([self.analysis_paths.analysis_number_path])
        files_text, files_state = self.check_something([info_file])
        self.analysis_data.set_info_data(dirs_text, dirs_state, files_text, files_state)

    def check_metadata_file(self):
        metadata_file = os.path.join(self.analysis_paths.analysis_number_path, 'metadata.txt')
        dirs_text, dirs_state = self.check_something([self.analysis_paths.analysis_number_path])
        files_text, files_state = self.check_something([metadata_file])
        self.analysis_data.set_metadata_data(dirs_text, dirs_state, files_text, files_state)

    @staticmethod
    def make_dir(p_dir_path):
        try:
            os.mkdir(p_dir_path)
        except FileExistsError as err:
            secho(f'Warning: {p_dir_path} 디렉터리가 존재합니다.', fg='yellow')
            echo(p_dir_path)
        else:
            echo(f'>>> {os.path.basename(p_dir_path)} 디렉터리 생성')
            echo(p_dir_path)

    @staticmethod
    def delete_dir(p_dir_path):
        try:
            shutil.rmtree(p_dir_path)
        except FileNotFoundError as err:
            echo('Warning 해당 디렉터리가 존재하지 않습니다.')
            echo(p_dir_path)
        except PermissionError as err:
            secho('PermissionError: 해당 디렉터리를 삭제할 수 없습니다.', fg='red', blink=True)
            echo(p_dir_path)
        else:
            echo(f'>>> {os.path.basename(p_dir_path)} 디렉터리 삭제')

    @staticmethod
    def parse_yaml(p_yaml):
        import yaml
        with open(p_yaml, 'r') as o_yaml:
            d_data = yaml.safe_load(o_yaml)
        return d_data

    @property
    def html_template_path(self):
        raise RuntimeError('Method html_template_path() 을 overriding 해야됩니다.')

    def set_report_template(self):
        self.template_env = Environment(lstrip_blocks=True, trim_blocks=True)
        self.template_env.loader = FileSystemLoader(self.html_template_path)

    def check_sample_count(self):
        """
        metadata.txt 파일과 STAT.txt 파일을 이용하여 시료의 개수를 확인하고, 확인된 시료의 개수를 반환한다.
        STAT.txt 파일이 없을 경우, metadata.txt 파일을 기준으로 시료의 개수를 산정한다.
        metadata.txt 파일이 없을 경우, 프로그램을 종료한다.

        :return:
        """
        echo('>>> 시료 개수 확인')
        metadata_file = os.path.join(self.analysis_paths.analysis_number_path, 'metadata.txt')
        metadata_name, metadata_status = self.exists(metadata_file)
        if metadata_status is True:
            l_sample_name = read_metadata_for_sort(metadata_file)
        else:
            secho(f'Error: {metadata_name} 파일이 없습니다.', fg='red', blink=True)
            echo(metadata_file)
        stat_file = os.path.join(self.analysis_paths.read_assembly_path, 'STAT.txt')
        stat_name, stat_status = self.exists(stat_file)
        sample_count = None
        if stat_status is True:
            from theCups.read_assembly import read_assembly_stat
            l_stat = read_assembly_stat(stat_file)
        else:
            secho(f'Error: {stat_name} 파일이 없습니다.', fg='red', blink=True)
            echo(stat_file)
        if all([metadata_status, stat_status]):
            if len(l_sample_name) == len(l_stat):
                sample_count = len(l_sample_name)
            else:
                secho(f'Error: {metadata_name}에 기재된 시료의 개수가 {stat_name}에 기재된 개수와 다릅니다.',
                      fg='red', blink=True)
                echo(f'{metadata_name}: {len(l_sample_name)}')
                echo(f'{stat_name}: {len(l_stat)}')
                exit()
        elif metadata_status is True:
            secho(f'Warning: {stat_name} 파일이 없습니다.', fg='yellow')
            secho(f'\t---> {metadata_name} 파일을 기준으로 시료의 개수가 산정됩니다.', fg='yellow')
            secho(f'{metadata_name}: {len(l_sample_name)}', fg='yellow')
            sample_count = len(l_sample_name)
        else:
            exit()
        secho(f'>>> 시료 개수: {sample_count}개 - 확인 완료', fg='cyan')
        return sample_count

    def confirm_file(self, kargs):
        """
        HTML Page 생성을 위해 필요한 파일의 존재 유무를 확인하고, 존재하지 않을 경우 안내 메시지를 출력한다.

        :param kargs:
                    file_name:
                    file_status:
                    html_name:
        :return:
        """
        if kargs['file_status']:
            return True
        else:
            self.page_status[kargs['html_name']] = False
            secho(f'Error: {kargs["file_name"]} 파일이 없습니다.', fg='red', blink=True)
            secho(f'\t --> {kargs["html_name"]} 생성 중단', fg='yellow')
            echo()  # 빈 줄
            return False

    def confirm_metadata(self, p_html):
        """
        metadata.txt 파일의 존재 유무를 확인하고, metadata.txt 파일의 위치를 반환한다.
        metadata.txt 파일이 없을 경우 안내 메시지를 출력한다.
        
        :param p_html: HTML 이름 또는 작업명.
        :return: (bool, [str | None])
        """
        metadata_name = 'metadata.txt'
        status = self.confirm_file(
            {
                'file_name': metadata_name,
                'file_status': self.analysis_data.metadata[self.analysis_number]['file'][metadata_name],
                'html_name': p_html,
            })
        if status:
            metadata_file = os.path.join(self.analysis_paths.analysis_number_path, metadata_name)
            return True, metadata_file
        else:
            return False, None

    def confirm_biom(self, p_html):
        biom_name = 'otu_table.biom'
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        biom_status = self.analysis_data.biom[analysis_dir_name['biom']]['file']['otu_table.biom']
        status = self.confirm_file(
            {
                'file_name': biom_name,
                'file_status': biom_status,
                'html_name': p_html,
            })
        if status:
            biom_file = os.path.join(self.analysis_paths.biom_path, biom_name)
            return True, biom_file
        else:
            return False, None

    def confirm_tre(self, p_html):
        tre_name = 'rep_phylo.tre'
        analysis_dir_name = self.analysis_paths.analysis_dir_name
        tre_status = self.analysis_data.phylogeny[analysis_dir_name['phylogeny']]['file'][tre_name]
        status = self.confirm_file(
            {
                'file_name': tre_name,
                'file_status': tre_status,
                'html_name': p_html,
            })
        if status:
            tre_file = os.path.join(self.analysis_paths.phylogeny_path, tre_name)
            return True, tre_file
        else:
            return False, None
