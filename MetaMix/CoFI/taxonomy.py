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
__version__ = '1.0.4'
# ----------------------------
# 2020.04.17
# Taxonomy Paser 변경 - XML
# ----------------------------
# 2020.05.20
# Queue Job 확인 간격 변경 - NT
# ----------------------------


import os
import sys
from SpoON.util import glob_dir, check_file_type, run_cmd, check_run_cmd, \
    parse_config, run_move_and_cmd, read_all_data, check_queue_node
from SpoON.run_time import compute_run_time
from click import secho, echo, style

CONFIG = parse_config()
QIIME2 = CONFIG['CoFI_QIIME2']
BLAST_PATH = CONFIG['CoFI_TAXONOMY']['blast_path']
BLAST_OPTION = CONFIG['CoFI_TAXONOMY']['blast_option']
BLAST_DB = CONFIG['CoFI_TAXONOMY']['blast_db']
UCLUST_DB = CONFIG['CoFI_TAXONOMY']['uclust_db']
BLAST_PARSER_PATH = CONFIG['CoFI_TAXONOMY']['blast_parser_path']
BLAST_PARSER_EXCEL_PATH = CONFIG['CoFI_TAXONOMY']['blast_parser_excel_path']


# NCBI 16S blast
# qiime tools import --input-path taxonomy.tsv --type FeatureData[Taxonomy] --source-format TSVTaxonomyFormat --output-path blast_taxonomy
# qiime taxa barplot --i-table ../table.qza --i-taxonomy blast_taxonomy.qza --m-metadata-file ../mappingfile.txt --o-visualization taxa-bar-plots.NCBI_16S.qzv


# 2018.11
# """
# Usage: qiime feature-classifier [OPTIONS] COMMAND [ARGS]...
#
#   Description: This QIIME 2 plugin supports taxonomic classification
#   of features using a variety of methods, including Naive Bayes,
#   vsearch, and BLAST+.
#
#   Plugin website: https://github.com/qiime2/q2-feature-classifier
#
#   Getting user support: Please post to the QIIME 2 forum for help
#   with this plugin: https://forum.qiime2.org
#
# Options:
#   --version    Show the version and exit.
#   --citations  Show citations and exit.
#   --help       Show this message and exit.
#
# Commands:
#   classify-consensus-blast    BLAST+ consensus taxonomy classifier
#   classify-consensus-vsearch  VSEARCH consensus taxonomy classifier
#   classify-sklearn            Pre-fitted sklearn-based taxonomy
#                               classifier
#   extract-reads               Extract reads from reference
#   fit-classifier-naive-bayes  Train the naive_bayes classifier
#   fit-classifier-sklearn      Train an almost arbitrary scikit-learn
#                               classifier

# """
# Usage: qiime feature-classifier classify-consensus-blast
#            [OPTIONS]
#
#   Assign taxonomy to query sequences using BLAST+. Performs BLAST+
#   local alignment between query and reference_reads, then assigns
#   consensus taxonomy to each query sequence from among maxaccepts
#   hits, min_consensus of which share that taxonomic assignment. Note
#   that maxaccepts selects the first N hits with > perc_identity
#   similarity to query, not the top N matches. For top N hits, use
#   classify-consensus-vsearch.
#
# Options:
#   --i-query ARTIFACT PATH FeatureData[Sequence]
#                                   Sequences to classify taxonomically.
#                                   [required]
#   --i-reference-reads ARTIFACT PATH FeatureData[Sequence]
#                                   reference sequences.  [required]
#   --i-reference-taxonomy ARTIFACT PATH FeatureData[Taxonomy]
#                                   reference taxonomy labels.
#                                   [required]
#   --p-maxaccepts INTEGER RANGE    Maximum number of hits to keep for
#                                   each query. Must be in range [1,
#                                   infinity]. BLAST will choose the
#                                   first N hits in the reference
#                                   database that exceed perc_identity
#                                   similarity to query.  [default: 10]
#   --p-perc-identity FLOAT         Reject match if percent identity to
#                                   query is lower. Must be in range
#                                   [0.0, 1.0].  [default: 0.8]
#   --p-strand [plus|minus|both]    Align against reference sequences in
#                                   forward ("plus"), reverse ("minus"),
#                                   or both directions ("both").
#                                   [default: both]
#   --p-evalue FLOAT                BLAST expectation value (E)
#                                   threshold for saving hits.
#                                   [default: 0.001]
#   --p-min-consensus FLOAT         Minimum fraction of assignments must
#                                   match top hit to be accepted as
#                                   consensus assignment. Must be in
#                                   range (0.5, 1.0].  [default: 0.51]
#   --p-unassignable-label TEXT     [default: Unassigned]
#   --o-classification ARTIFACT PATH FeatureData[Taxonomy]
#                                   Taxonomy classifications of query
#                                   sequences.  [required if not passing
#                                   --output-dir]
#   --output-dir DIRECTORY          Output unspecified results to a
#                                   directory
#   --cmd-config FILE               Use config file for command options
#   --verbose                       Display verbose output to stdout
#                                   and/or stderr during execution of
#                                   this action.  [default: False]
#   --quiet                         Silence output if execution is
#                                   successful (silence is golden).
#                                   [default: False]
#   --citations                     Show citations and exit.
#   --help                          Show this message and exit.
# """

def get_q2_blast_cmd(p_args):
    cmd = '{qiime2} feature-classifier classify-consensus-blast ' \
        '--i-query {query_qza}' \
        '--i-reference-reads {ref_seq_qza} ' \
        '--i-reference-taxonomy {ref_tax_qza} ' \
        '--p-maxaccepts 10 ' \
        '--p-perc-identity 0.85' \
        '--p-strand both ' \
        '--p-evalue 0.001 ' \
        '--p-min-consensus 0.51 ' \
        '--p-unassignable-label Unassigned ' \
        '--output-dir {output_dir}'.format(
            qiime2=QIIME2,
            query_qza='',
            ref_seq_qza='',
            ref_tax_qza='',
            output_dir=''
        )
    return cmd


# """
# Usage: qiime feature-classifier classify-consensus-vsearch
#            [OPTIONS]
#
#   Assign taxonomy to query sequences using VSEARCH. Performs VSEARCH
#   global alignment between query and reference_reads, then assigns
#   consensus taxonomy to each query sequence from among maxaccepts
#   top hits, min_consensus of which share that taxonomic assignment.
#   Unlike classify-consensus-blast, this method searches the entire
#   reference database before choosing the top N hits, not the first N
#   hits.
#
# Options:
#   --i-query ARTIFACT PATH FeatureData[Sequence]
#                                   Sequences to classify taxonomically.
#                                   [required]
#   --i-reference-reads ARTIFACT PATH FeatureData[Sequence]
#                                   reference sequences.  [required]
#   --i-reference-taxonomy ARTIFACT PATH FeatureData[Taxonomy]
#                                   reference taxonomy labels.
#                                   [required]
#   --p-maxaccepts INTEGER RANGE    Maximum number of hits to keep for
#                                   each query. Set to 0 to keep all
#                                   hits > perc_identity similarity.
#                                   Must be in range [0, infinity].
#                                   [default: 10]
#   --p-perc-identity FLOAT         Reject match if percent identity to
#                                   query is lower. Must be in range
#                                   [0.0, 1.0].  [default: 0.8]
#   --p-strand [plus|both]          Align against reference sequences in
#                                   forward ("plus") or both directions
#                                   ("both").  [default: both]
#   --p-min-consensus FLOAT         Minimum fraction of assignments must
#                                   match top hit to be accepted as
#                                   consensus assignment. Must be in
#                                   range (0.5, 1.0].  [default: 0.51]
#   --p-unassignable-label TEXT     [default: Unassigned]
#   --p-threads INTEGER             [default: 1]
#   --o-classification ARTIFACT PATH FeatureData[Taxonomy]
#                                   The resulting taxonomy
#                                   classifications.  [required if not
#                                   passing --output-dir]
#   --output-dir DIRECTORY          Output unspecified results to a
#                                   directory
#   --cmd-config FILE               Use config file for command options
#   --verbose                       Display verbose output to stdout
#                                   and/or stderr during execution of
#                                   this action.  [default: False]
#   --quiet                         Silence output if execution is
#                                   successful (silence is golden).
#                                   [default: False]
#   --citations                     Show citations and exit.
#   --help                          Show this message and exit.
# """
def get_q2_vsearch_cmd(p_args):
    cmd = '{qiime2} feature-classifier classify-consensus-vsearch ' \
        '--i-query ' \
        '--i-reference-reads ' \
        '--i-reference-taxonomy ' \
        '--p-maxaccepts ' \
        '--p-perc-identity ' \
        '--p-strand both ' \
        '--p-min-consensus 0.51 ' \
        '--p-unassignable-label Unassigned ' \
        '--p-threads 0 ' \
        '--output-dir '.format(
            qiime2=QIIME2
        )
    return cmd


# """
# Usage: qiime taxa [OPTIONS] COMMAND [ARGS]...
#
#   Description: This QIIME 2 plugin provides functionality for
#   working with and visualizing taxonomic annotations of features.
#
#   Plugin website: https://github.com/qiime2/q2-taxa
#
#   Getting user support: Please post to the QIIME 2 forum for help
#   with this plugin: https://forum.qiime2.org
#
# Options:
#   --version    Show the version and exit.
#   --citations  Show citations and exit.
#   --help       Show this message and exit.
#
# Commands:
#   barplot       Visualize taxonomy with an interactive bar plot
#   collapse      Collapse features by their taxonomy at the specified
#                 level
#   filter-seqs   Taxonomy-based feature sequence filter.
#   filter-table  Taxonomy-based feature table filter.
# """

# """
# Usage: qiime taxa barplot [OPTIONS]
#
#   This visualizer produces an interactive barplot visualization of
#   taxonomies. Interactive features include multi-level sorting, plot
#   recoloring, sample relabeling, and SVG figure export.
#
# Options:
#   --i-table ARTIFACT PATH FeatureTable[Frequency]
#                                   Feature table to visualize at
#                                   various taxonomic levels.
#                                   [required]
#   --i-taxonomy ARTIFACT PATH FeatureData[Taxonomy]
#                                   Taxonomic annotations for features
#                                   in the provided feature table. All
#                                   features in the feature table must
#                                   have a corresponding taxonomic
#                                   annotation. Taxonomic annotations
#                                   that are not present in the feature
#                                   table will be ignored.  [required]
#   --m-metadata-file MULTIPLE FILE
#                                   Metadata file or artifact viewable
#                                   as metadata. This option may be
#                                   supplied multiple times to merge
#                                   metadata. The sample metadata.
#                                   [required]
#   --o-visualization VISUALIZATION PATH
#                                   [required if not passing --output-
#                                   dir]
#   --output-dir DIRECTORY          Output unspecified results to a
#                                   directory
#   --cmd-config FILE               Use config file for command options
#   --verbose                       Display verbose output to stdout
#                                   and/or stderr during execution of
#                                   this action.  [default: False]
#   --quiet                         Silence output if execution is
#                                   successful (silence is golden).
#                                   [default: False]
#   --citations                     Show citations and exit.
#   --help                          Show this message and exit.
# """
def get_taxa_barplot(p_args):
    cmd = '{qiime2} taxa barplot ' \
        '--i-table ' \
        '--i-taxonomy ' \
        '--m-metadata-file '.format(qiime2=QIIME2)
    return cmd


def get_taxonomy_qiime1_cmd(p_fasta, p_out_dir, p_tool, p_db_fasta, p_db_taxa):
    cmd = 'assign_taxonomy.py ' \
        '-i {fasta} ' \
        '-m {tool} ' \
        '-r {db_fasta} ' \
        '-t {db_taxa} ' \
        '-o {out_dir} '.format(
            fasta=p_fasta,
            tool=p_tool,
            db_fasta=p_db_fasta,
            db_taxa=p_db_taxa,
            out_dir=p_out_dir)
    return cmd


def run_taxonomy_qiime1(p_fasta, p_db_tool, p_db, p_path):
    """
    QIIME1 프로그램을 이용하여 Taxonomy Assignment 을 진행한다.
    db_tool을 blast 로 선택시 stand

    :param p_fasta:
    :param p_db_tool:
    :param p_db:
    :param p_path: TAXONOMY 디렉터리 경로
    :return:
    """
    global UCLUST_DB, BLAST_DB
    if p_db_tool == 'blast':
        d_db = BLAST_DB[p_db]
        path = os.path.join(p_path,  f'blast_{p_db}')
    elif p_db_tool == 'uclust':
        d_db = UCLUST_DB[p_db]
        path = os.path.join(p_path, f'uclust_{p_db}')
    else:
        raise ValueError('사용할 DB가 설정되지 않았습니다.')
    cmd = get_taxonomy_qiime1_cmd(p_fasta, path, p_db_tool, d_db['fasta'], d_db['taxon'])
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': f'QIIME1 Taxonomy Assignment 완료: {p_db_tool}, {p_db}\n'
                        f'DB Version: {d_db["version"]}',
            'false_meg': 'QIIME1 assign_taxonomy'
        }
    )
    with open(os.path.join(path, 'db_version.txt'), 'w') as o_db_version:
        o_db_version.write(d_db['version'])
        o_db_version.write('\n')
    with open(os.path.join(path, f'TAXONOMY_CMD.{p_db_tool}_{p_db}.recipe'), 'w') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')


class BLAST:
    def __init__(self, p_kargs):
        """

         blast_option
            query_coverage: 85
            identity_percentage: 85
            queue: bi3m.q
            e_value: 1.0E-3
            num_alignment: 5

        :param p_kargs:
                tool: ex) blastn
                db: NCBI_16S, NCBI_NT 둥
                fasta: otus_rep.fasta
                out_path:
                remove_taxon:
                keep_taxon:
                query_coverage:
                identity_percentage:
                read_per_job:
                nt_max_job:
                queue:
                tool:
                db:
        """
        global BLAST_PATH, BLAST_OPTION, BLAST_PARSER_PATH, BLAST_PARSER_EXCEL_PATH
        self.node = None
        self.queue_state = None
        self.blast_path = BLAST_PATH
        self.blast_option = BLAST_OPTION
        self.blast_parser_path = BLAST_PARSER_PATH
        self.blast_parser_excel_path = BLAST_PARSER_EXCEL_PATH
        self.read_per_job = p_kargs['read_per_job']
        self.nt_max_job = p_kargs['nt_max_job']
        self.queue = p_kargs['queue']
        self.tool = p_kargs['tool']
        self.db = p_kargs['db']
        self.d_db = None
        self.query_file = p_kargs['fasta']
        self.otus_count_for_query = None
        self.out_path = p_kargs['out_path']
        self.remove_taxon = p_kargs['remove_taxon']
        self.keep_taxon = p_kargs['keep_taxon']
        self.query_coverage = p_kargs['query_coverage']
        self.identity_percentage = p_kargs['identity_percentage']
        self.out_tmp = 'blast_tmp'
        self.out_archive = 'blast_archive'
        self.out_xml = 'blast_xml'
        self.out_result = 'blast_result'
        self.db_text = "db_version.txt"
        self.l_split_query = list()

    def check_otus_count(self):
        from humanize import intcomma
        cmd = f'grep -c ">" {self.query_file}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': f'OTU 개수 확인: {intcomma(run.stdout.strip())} 개',
                'false_meg': 'check_otus_count',
            }
        )
        self.otus_count_for_query = int(run.stdout.strip())

    def set_db(self):
        import math
        self.d_db = BLAST_DB[self.db]
        if self.db == 'NCBI_16S':
            self.blast_option['qstat_check_time'] = 15
        elif self.db == 'NCBI_NT':
            self.blast_option['qstat_check_time'] = 60 * 1
            self.read_per_job = math.ceil(self.otus_count_for_query / self.nt_max_job)
            secho(f'>>> read_per_job 설정 : {self.read_per_job}', fg='cyan')
        else:
            self.blast_option['qstat_check_time'] = 15

    def make_job_dir(self, p_version):
        """
        BLAST 실행에 필요한 디렉터리(blast_archive, blast_tmp, blast_xml)들을 생성한다.
        :return:
        """
        if p_version == 'V2':
            l_dirs = [self.out_tmp, self.out_archive, self.out_result]
        elif p_version == 'V3':
            l_dirs = [self.out_tmp, self.out_archive, self.out_xml]
        elif p_version == 'JSON':
            l_dirs = [self.out_tmp, self.out_archive, 'blast_json']
        else:
            l_dirs = [self.out_tmp, self.out_archive, self.out_result]

        for dir_name in l_dirs:
            try:
                dir_path = os.path.join(self.out_path, dir_name)
                os.mkdir(dir_path)
            except FileExistsError as err:
                secho(f'Error: {dir_name} 디렉터리가 이미 있습니다.', fg='red', err=True, blink=True)
                echo(err, err=True)
                exit(1)
            else:
                echo(f'>>> {dir_name} 디렉터리 생성')
                echo(dir_path)
        echo()  # 빈 줄

    def split_otus_fasta(self):
        """
        otus_rep.fasta 파일에 기재된 서열을 특정 개수로 나누어서 저장한다.

        :return:
        """
        with open(self.query_file, 'r') as o_query:
            file_count = 1
            split_data = list()
            line_num = 0
            text = o_query.readline()
            while text:
                if text.startswith('>'):
                    line_num += 1
                if line_num <= self.read_per_job * file_count:
                    split_data.append(text)
                else:
                    split_file = os.path.join(self.out_path, self.out_tmp, f'query.fasta_{file_count}')
                    self.l_split_query.append(split_file)
                    with open(split_file, 'w') as o_outfile:
                        o_outfile.writelines(split_data)
                    file_count += 1
                    split_data = list()
                    split_data.append(text)
                text = o_query.readline()
            else:
                split_file = os.path.join(self.out_path, self.out_tmp, f'query.fasta_{file_count}')
                self.l_split_query.append(split_file)
                with open(split_file, 'w') as o_outfile:
                    o_outfile.writelines(split_data)
            echo(f'>>> 분할된 대표서열 파일 생성: {file_count}개')

    def get_blastn_code(self, p_query_file):
        """

        :param p_query_file:
        :return:
        """
        out_name = os.path.basename(p_query_file)
        out_file = os.path.join(self.out_path, self.out_archive, f'Job_{out_name}.archive')
        cmd = f'{self.blast_path}/{self.tool} ' \
              f'-db {self.d_db["db_name"]} ' \
              f'-query {p_query_file} ' \
              f'-out {out_file} ' \
              f'-evalue {self.blast_option["e_value"]} ' \
              f'-num_alignments {self.blast_option["num_alignments"]} ' \
              '-outfmt 11'
        return cmd

    def get_blast_formatter_code(self, p_query_file, p_type):
        """

        :param p_query_file: blast_tmp 디렉터리에 생성된 query_?.fasta
        :param p_type: blast outfmt
                    0: text
                    5: XML
                    12, 13 ,15 : JSON
        :return:
        """
        if p_type == 0:
            extension = '.blast'
            out_dir = self.out_result
        elif p_type == 5:
            extension = '.xml'
            out_dir = self.out_xml
        elif p_type in [12, 13, 15]:
            extension = '.json'
            out_dir = 'blast_json'
        else:
            extension = '.txt'
            out_dir = self.out_result

        out_name = os.path.basename(p_query_file)
        archive_file = os.path.join(self.out_path, self.out_archive, f'Job_{out_name}.archive')
        result_file = os.path.join(self.out_path, out_dir, f'Job_{out_name}{extension}')
        cmd = f'{self.blast_path}/blast_formatter ' \
              f'-archive {archive_file} ' \
              f'-out {result_file} ' \
              f'-outfmt {p_type}'
        return cmd

    def save_db_version_file(self):
        version_file = os.path.join(self.out_path, self.db_text)
        with open(version_file, 'w') as o_version:
            o_version.write(self.d_db['version'])
            o_version.write('\n')
        echo(f'>>> {self.db_text} 생성')
        echo(version_file)

    def save_blast_bash(self, p_query_file, p_type):
        """

        :param p_query_file: blast_tmp 디렉터리에 생성된 query_?.fasta
        :param p_type: blast outfmt
                    0: text
                    5: XML
                    12, 13 ,15 : JSON
        :return:
        """
        query_name = os.path.basename(p_query_file)
        blast_bash_file = os.path.join(self.out_path, self.out_tmp, f'{query_name}.sh')
        with open(blast_bash_file, 'w') as o_bash:
            o_bash.write(self.get_blastn_code(p_query_file))
            o_bash.write('\n')
            o_bash.write('\n')
            o_bash.write(self.get_blast_formatter_code(p_query_file, p_type))
            o_bash.write('\n')
            o_bash.write('\n')
            o_bash.write(f'date > {blast_bash_file}.done')
            o_bash.write('\n')
        return blast_bash_file

    def run_blast_qsub(self, p_blast_bash):
        """
        BLAST 실행을 위한 BASH 파일을 qsub 명령어를 이용하여 Queue 에 등록한다.

        :param p_blast_bash: query_?.bash
        :return:
        """
        cmd = 'qsub -cwd -j y ' \
              f'-q {self.queue} ' \
              '-S /bin/sh ' \
              f'-o {self.out_path}/qsub.err ' \
              f'{p_blast_bash} >> {self.out_path}/qsub.log'

        self.queue_state, self.node = check_queue_node()
        if self.queue_state is True:
            run = run_move_and_cmd(self.out_path, cmd)
            check_run_cmd(
                {
                    'run': run,
                    'true_meg': f'{os.path.basename(p_blast_bash)} 작업 등록 완료',
                    'false_meg': 'run_blast_qsub'
                }
            )
        else:
            secho('Error: Queue 사용할 수 없는 서버입니다.', fg='red', blink=True, err=True)
            echo(f'서버: {self.node}', err=True)
            secho('\t--> qsub_blast.recipe만 생성합니다.', fg='magenta')
        recipe_file = os.path.join(self.out_path, 'qsub_blast.recipe')
        with open(recipe_file, 'a') as o_recipe:
            o_recipe.write(cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')

    def check_jobs(self):
        """
        Queue 에 등록된 작업의 완료 여부를 확인한다.
        완료여부는 blast_tmp 디렉터리의 *.sh.done 파일의 생성 여부로 확인한다.

        :return:
        """
        from time import time, sleep
        start_time = time()
        echo(f'>>> BLAST 진행 중 (Queue) : 확인주기 {self.blast_option["qstat_check_time"]}초')
        echo('!! Queue 에 등록된 작업은 CoFI 프로그램과 연결이 끊기더라도 동작합니다.')
        target_path = os.path.join(self.out_path, 'blast_tmp')
        l_glob_files = glob_dir(target_path, '*.sh', 'many', False)
        while l_glob_files is None:
            sleep(5)
            l_glob_files = glob_dir(target_path, '*.sh', 'many', False)
        l_target_files = [f'{x}.done' for x in l_glob_files]
        job_status = [check_file_type(x, 'exists') for x in l_target_files]
        timer_up = 0
        status_check = 1
        while all(job_status) is False:
            check_time = time()
            run_time = compute_run_time(start_time, check_time)
            secho(f'\r{job_status.count(True)} / {len(job_status)} 완료 '
                  f'- {run_time} 경과, 상태확인: {status_check}회', nl=False, fg='yellow')
            sleep(1)
            timer_up += 1
            if timer_up % self.blast_option['qstat_check_time'] == 0:
                job_status = [check_file_type(x, 'exists') for x in l_target_files]
                status_check += 1
        else:
            echo()  # 줄바꿈
            secho('>>> blast qsub 완료', fg='cyan')

    def check_qsub_error(self):
        """
        qsub.err 파일을 확인하여 오류 여부가 있는지 확인한다.

        :return:
        """
        error_log_file = os.path.join(self.out_path, 'qsub.err')
        if check_file_type(error_log_file, 'exists'):
            error_log = read_all_data(error_log_file, 'r')
            if error_log == '':
                pass
            else:
                secho('Warning: qsub.err 내용을 확인하세요.', fg='yellow')
                echo('---- qsub.err ----')
                echo(error_log)
                echo()
        else:
            secho('Warning: qsub.err 파일이 없습니다.', fg='yellow')
            echo(error_log_file)
            secho('\t--> qsub 명령어가 정상적으로 완료되었는지 확인하세요.', fg='magenta')

    def copy_blast_parser_scripts(self):
        cmd = f'cp -t {self.out_path} {os.path.join(self.blast_parser_path, "BLUEARA.jar")} ' \
              f'{os.path.join(self.blast_parser_path, "BlastPlusParser2*")} ' \
              f'{os.path.join(self.blast_parser_path, "MakeBlastNParsing.jar")}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'BLAST Parser Scripts 복사 완료'
            }
        )

    def get_blastplus_parser2_cmd(self):
        cmd = 'java -cp BLUEARA.jar:. BlastPlusParser2 ' \
              f'-d {os.path.join(self.out_path, self.out_result)} ' \
              f'-o {os.path.join(self.out_path, "blast_results.txt")} ' \
              '-l F ' \
              f'-jn Job_query.fasta'
        return cmd

    def get_make_blastn_parsing_cmd(self):
        cmd = 'java -jar MakeBlastNParsing.jar blast_results.txt'
        return cmd

    def get_add_taxonomy_cmd(self):
        cmd = f'{os.path.join(self.blast_parser_path, "AddTaxonomy.py")} ' \
              f'{self.db_text} ' \
              f'{self.d_db["taxon"]} ' \
              f'blast_results.txt.final.txt'
        return cmd

    def get_blast_xml_parser_cmd(self):
        cmd = f'{os.path.join(self.blast_parser_excel_path, "BLAST_XML_Parser.py")} ' \
              f'{os.path.join(self.out_path, self.out_xml)} ' \
              f'--output {os.path.join(self.out_path, "blast_results.txt")}'
        return cmd

    def get_add_taxonomy2_cmd(self):
        cmd = f'{os.path.join(self.blast_parser_path, "AddTaxonomy_V2.py")} ' \
              f'{os.path.join(self.out_path,self.db_text)} ' \
              f'{self.d_db["taxon"]} ' \
              f'{os.path.join(self.out_path, "blast_results.txt")}'
        return cmd

    def run_blast_parse(self):
        echo('==== BLAST Parser 실행 ====')
        blast_parser_cmd = self.get_blastplus_parser2_cmd()
        blastn_parsing_cmd = self.get_make_blastn_parsing_cmd()
        add_taxonomy_cmd = self.get_add_taxonomy_cmd()
        with open(os.path.join(self.out_path, 'blast_parser.recipe'), 'w') as o_recipe:
            o_recipe.write(blast_parser_cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')
            o_recipe.write(blastn_parsing_cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')
            o_recipe.write(add_taxonomy_cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')

        blast_parser_run = run_move_and_cmd(self.out_path, blast_parser_cmd)
        check_run_cmd(
            {
                'run': blast_parser_run,
                'true_meg': 'BLAST+ Parser2 완료',
                'false_meg': 'blastplus_parser2',
            }
        )
        blastn_parsing_run = run_move_and_cmd(self.out_path, blastn_parsing_cmd)
        check_run_cmd(
            {
                'run': blastn_parsing_run,
                'true_meg': 'BLASTN Paring 완료',
                'false_meg': 'make_blastn_parsing',
            }
        )
        add_taxonomy_run = run_move_and_cmd(self.out_path, add_taxonomy_cmd)
        check_run_cmd(
            {
                'run': add_taxonomy_run,
                'true_meg': 'Add_Taxonomy 완료',
                'false_meg': 'add_taxonomy',
            }
        )

    def run_blast_parse2(self):
        echo('==== BLAST Parser 실행(XML) ====')
        blast_xml_parser_cmd = self.get_blast_xml_parser_cmd()
        add_taxonomy_cmd = self.get_add_taxonomy2_cmd()
        with open(os.path.join(self.out_path, 'blast_parser.recipe'), 'w') as o_recipe:
            o_recipe.write(blast_xml_parser_cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')
            o_recipe.write(add_taxonomy_cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')

        blast_parser_run = run_cmd(blast_xml_parser_cmd)
        check_run_cmd(
            {
                'run': blast_parser_run,
                'true_meg': 'BLAST XML Parser 완료',
                'false_meg': 'BLAST_XML_Parser',
            }
        )
        add_taxonomy_run = run_cmd(add_taxonomy_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': add_taxonomy_run,
                'true_meg': 'Add_Taxonomy2 완료',
                'false_meg': 'add_taxonomy2',
            },  p_stdout=False,
        )

    def make_excel_and_taxonomy_file(self, p_version):
        """

        :param p_version: [V2, V3]
        :return:
        """
        if p_version == 'V2':
            parser = 'BLAST_Taxonomy_parser_V2.py'
            input_name = 'blast_results.txt.final.txt_taxonomy.txt'
        elif p_version == 'V3':
            parser = 'BLAST_Taxonomy_parser_V3.py'
            input_name = 'blast_results_taxonomy.txt'
        else:
            raise ValueError(f'p_version: {p_version}')
        cmd = f'{os.path.join(self.blast_parser_excel_path, parser)} ' \
              f'-i {os.path.join(self.out_path, input_name)} ' \
              f'--out_path {self.out_path} ' \
              f'{"" if self.remove_taxon is None else self.remove_taxon} ' \
              f'{"" if self.keep_taxon is None else self.keep_taxon} ' \
              f'-qc {self.query_coverage} ' \
              f'-ip {self.identity_percentage}'
        with open(os.path.join(self.out_path, 'blast_parser.recipe'), 'a') as o_recipe:
            o_recipe.write(cmd)
            o_recipe.write('\n')
            o_recipe.write('\n')
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'BLAST Excel & QIIME Taxonomy Files 생성 완료',
                'false_meg': 'make_excel_and_taxonomy_file'
            }
        )

    def confirm_otus_count_from_blast_results(self, p_version):
        if p_version == 'V2':
            blast_taxonomy_file = os.path.join(self.out_path, 'blast_results.txt.final.txt_taxonomy.txt')
        elif p_version == 'V3':
            blast_taxonomy_file = os.path.join(self.out_path, 'otus_rep_tax_assignments.txt')
        else:
            raise ValueError(f'p_version: {p_version}')
        # cmd = f'grep -c "denovo" {blast_taxonomy_file}'
        cmd = f'wc -l {blast_taxonomy_file}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'OTU 개수 확인: blast__taxonomy.txt',
                'false_meg': 'confirm_otus_count_from_blast_results'
            }, p_exit=False,
        )
        if p_version == 'V2':
            otus_count_for_run = int(run.stdout.split()[0]) - 2
        elif p_version == 'V3':
            otus_count_for_run = int(run.stdout.split()[0])
        else:
            raise ValueError(f'p_version: {p_version}')

        if self.otus_count_for_query == otus_count_for_run:
            secho('>>>OTU 개수 확인: query == blast', fg='cyan')
        else:
            secho('Error: query != blast', fg='red', blink=True, err=True)
            echo(f'Query: {self.otus_count_for_query}', err=True)
            echo(f'\t--> {self.query_file}', err=True)
            echo(f'BLAST: {run.stdout}', err=True)
            echo(f'\t--> {blast_taxonomy_file}', err=True)
            exit()


def run_standalone_blast(p_kargs):
    """

    :param p_kargs: 아래의 항목을 키로 가지는 딕션너리.
            fasta: otus_rep.fasta
            db: 분석에 사용할 DB. config.yaml 파일 참조.
            out_path:
            remove_taxon: taxon 제거 조건.
                - no_hit, filtered, uncultured, environmental, bacterium, marine, unidentified, eukaryote, all
                - CoFI Option 과 BLAST_Taxonomy_parser_V2.py 참조
            keeep_taoxon: taxon 유지 조건. 나머지 모두 제거.
            query_coverage:
            identity_percentage:
            read_per_job: Queue Job 당 read 의 개수. otus_rep.fasta 를 나눌 때 파일당 read의 개수.
            nt_max_job: NT DB를 사용할 경우 Queue 에 넣을 수 있는 Job의 최대 개수.
            queue: Queue 이름.
            mode: [blast, parse, all]
    :return:
    """
    global BLAST_DB
    i_blast = BLAST(
        {
            'tool': p_kargs['tool'],
            'fasta': p_kargs['fasta'],
            'db': p_kargs['db'],
            'out_path': p_kargs['out_path'],
            'remove_taxon': p_kargs['remove_taxon'],
            'keep_taxon': p_kargs['keep_taxon'],
            'query_coverage': p_kargs['query_coverage'],
            'identity_percentage': p_kargs['identity_percentage'],
            'read_per_job': p_kargs['read_per_job'],
            'nt_max_job': p_kargs['nt_max_job'],
            'queue': p_kargs['queue'],
        }
    )
    if p_kargs['mode'] == 'parse':
        i_blast.check_otus_count()
        i_blast.set_db()
        i_blast.copy_blast_parser_scripts()
        i_blast.run_blast_parse()
        i_blast.make_excel_and_taxonomy_file('V2')
        i_blast.confirm_otus_count_from_blast_results('V2')
    else:  # blast, all
        i_blast.check_otus_count()
        i_blast.set_db()
        i_blast.make_job_dir('V2')
        i_blast.split_otus_fasta()
        for query in i_blast.l_split_query:
            i_blast.save_blast_bash(query, 0)
            i_blast.run_blast_qsub(f'{query}.sh')
        if i_blast.queue_state is True:
            i_blast.save_db_version_file()
            i_blast.check_jobs()
            i_blast.check_qsub_error()
            if p_kargs['mode'] == 'all':
                i_blast.copy_blast_parser_scripts()
                i_blast.run_blast_parse()
                i_blast.make_excel_and_taxonomy_file('V2')
                i_blast.confirm_otus_count_from_blast_results('V2')
        else:
            secho('Warning: qsub는 현재 denovo06, denovo10 서버에서만 가능합니다.', fg='yellow', bold=True)
            secho('\t---> denovo6, denovo10 서버에서 qsub_blast.recipe 를 실행하세요.', fg='magenta')
            echo(f'{os.path.join(i_blast.out_path, "qsub_blast.recipe")}')


def run_standalone_blast2(p_kargs):
    """

    :param p_kargs: 아래의 항목을 키로 가지는 딕션너리.
            fasta: otus_rep.fasta
            db: 분석에 사용할 DB. config.yaml 파일 참조.
            out_path:
            remove_taxon: taxon 제거 조건.
                - no_hit, filtered, uncultured, environmental, bacterium, marine, unidentified, eukaryote, all
                - CoFI Option 과 BLAST_Taxonomy_parser_V2.py 참조
            keeep_taoxon: taxon 유지 조건. 나머지 모두 제거.
            query_coverage:
            identity_percentage:
            read_per_job: Queue Job 당 read 의 개수. otus_rep.fasta 를 나눌 때 파일당 read의 개수.
            nt_max_job: NT DB를 사용할 경우 Queue 에 넣을 수 있는 Job의 최대 개수.
            queue: Queue 이름.
            mode: [blast, parse, all]
    :return:
    """
    global BLAST_DB
    i_blast = BLAST(
        {
            'tool': p_kargs['tool'],
            'fasta': p_kargs['fasta'],
            'db': p_kargs['db'],
            'out_path': p_kargs['out_path'],
            'remove_taxon': p_kargs['remove_taxon'],
            'keep_taxon': p_kargs['keep_taxon'],
            'query_coverage': p_kargs['query_coverage'],
            'identity_percentage': p_kargs['identity_percentage'],
            'read_per_job': p_kargs['read_per_job'],
            'nt_max_job': p_kargs['nt_max_job'],
            'queue': p_kargs['queue'],
        }
    )
    if p_kargs['mode'] == 'parse':
        i_blast.check_otus_count()
        i_blast.set_db()
        i_blast.run_blast_parse2()
        i_blast.make_excel_and_taxonomy_file('V3')
        i_blast.confirm_otus_count_from_blast_results('V3')
    else:  # blast, all
        i_blast.check_otus_count()
        i_blast.set_db()
        i_blast.make_job_dir('V3')
        i_blast.split_otus_fasta()
        for query in i_blast.l_split_query:
            i_blast.save_blast_bash(query, 5)
            i_blast.run_blast_qsub(f'{query}.sh')
        if i_blast.queue_state is True:
            i_blast.save_db_version_file()
            i_blast.check_jobs()
            i_blast.check_qsub_error()
            if p_kargs['mode'] == 'all':
                i_blast.run_blast_parse2()
                i_blast.make_excel_and_taxonomy_file('V3')
                i_blast.confirm_otus_count_from_blast_results('V3')
        else:
            secho('Warning: qsub는 현재 denovo06, denovo10 서버에서만 가능합니다.', fg='yellow', bold=True)
            secho('\t---> denovo6, denovo10 서버에서 qsub_blast.recipe 를 실행하세요.', fg='magenta')
            echo(f'{os.path.join(i_blast.out_path, "qsub_blast.recipe")}')