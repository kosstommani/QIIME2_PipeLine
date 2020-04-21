#!/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
# ----------------------------------------------------------------------------------------------------------------------
#                        888b     d888          888             888b     d888 d8b
#                        8888b   d8888          888             8888b   d8888 Y8P
#                        88888b.d88888          888             88888b.d88888
# 88888b.d88b.  888  888 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b 888  888 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 Y88b 888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888  "Y88888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#                    888
#               Y8b d88P
#                "Y88P"
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '0.2.1'

import click
import sys
import os
import time
from pprint import pprint
from SpoON.util import (check_order_number_system, run_cmd, check_run_cmd,
                        parse_config, glob_dir, launcher_cmd, check_file_type)
from SpoON.fastq_handler import exist_something, exist_assembled_fastq
from myBiomeStory.my_util import (print_run_time, check_sample_list, set_sample_list, get_sample_info,
                                  read_csv_to_dict, read_csv, check_db_version)


# 기본값 설정
CONFIG = parse_config()
ANALYSIS_BASE_PATH = '/garnet/Analysis/BI/AmpliconMetaGenome_MAM'
R1_SUFFIX = CONFIG['CoFI']['R1_suffix']
R2_SUFFIX = CONFIG['CoFI']['R2_suffix']
ASSEMBLY_SUFFIX = '.extendedFrags.fasta'
ADAPTER_TRIM_TOOL = CONFIG['CoFI']['adapter_trim_tool']
TARGET_DIR_SUFFIX = CONFIG['CoFI']['target_dir_suffix']
POOL_WORKER = 10

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.argument('order_number')
@click.option('--copy_rawdata', '-c',
              is_flag=True,
              help='RawData 복사.')
@click.option('--no_prema',
              is_flag=True,
              help='PreMA 실행 안함.')
@click.option('--full_course',
              is_flag=True,
              help='풀코스 요리를 먹을 수 있음(PreMA 별도).')
@click.option('--flash',
              is_flag=True,
              help='Read Assembly 실행.')
@click.option('--closed',
              is_flag=True,
              help='closed-reference OTU Picking 실행.')
@click.option('--diversity',
              is_flag=True,
              help='Alpha Diversity 계산(Shannon, Simpson)')
@click.option('--summarize_taxa',
              is_flag=True,
              help='biom 파일로부터 L7_Table 생성.')
@click.option('--make_story',
              is_flag=True,
              help='보고서에 사용할 데이터 생성.')
@click.option('--insert',
              is_flag=True,
              help='DB에 데이터 삽입.')
@click.option('--include',
              help='분석에 포함할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용.'
                   'ex) --include "A B C"')
@click.option('--exclude',
              help='분석에 제외할 시료명. --closed, --diversity, --summarize_taxa, --make_story 에만 적용'
                   'ex) --exclude "D E F"')
@click.option('--no_check',
              is_flag=True,
              help='Analysis_? 디렉터리의 개수가 여러 개일 경우 종료하지 않음. 사용할 디렉터리 선택.')
def my_main(**kargs):
    """
    \b
                           888b     d888          888             888b     d888 d8b
                           8888b   d8888          888             8888b   d8888 Y8P
                           88888b.d88888          888             88888b.d88888
    88888b.d88b.  888  888 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
    888 "888 "88b 888  888 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
    888  888  888 888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
    888  888  888 Y88b 888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
    888  888  888  "Y88888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
                       888
                  Y8b d88P
                   "Y88P"

    \b
    개인 유전체 - 장내 미생물 모니터링 서비스.
    MY BIOMESTORY

    !! 301bp 기준으로 설정됨.

    !! Closed-reference OTU Picking 진행시 QIIME config file 설정을 해야 됩니다.
    """
    pprint(kargs)
    global ANALYSIS_BASE_PATH, R1_SUFFIX, R2_SUFFIX, ADAPTER_TRIM_TOOL, ASSEMBLY_SUFFIX
    start_time = time.time()
    kargs.update({
        'analysis_base_path': ANALYSIS_BASE_PATH,
        'r1_suffix': R1_SUFFIX,
        'r2_suffix': R2_SUFFIX,
        'adapter_trim_tool': ADAPTER_TRIM_TOOL,
        'target_dir_suffix': TARGET_DIR_SUFFIX,
    })

    if (kargs['include'] is not None) and (kargs['exclude'] is not None):
        click.secho('--include 와 --exclude 를 동시에 사용할 수 없습니다.', fg='red')
        exit(1)
    else:
        if kargs['include'] is not None:
            l_include = kargs['include'].strip().split()
        else:
            l_include = None
        if kargs['exclude'] is not None:
            l_exclude = kargs['exclude'].strip().split()
        else:
            l_exclude = None

    # 수주번호 형식 확인
    if check_order_number_system(kargs['order_number']) == 'LIMS2' or 'LIMS3':
        pass

    # PreMA - Data Copy & QC
    if kargs['no_prema'] is False:
        prema_cmd = \
            'python3 /garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/PreMA.py ' \
            'CORE ' \
            '{copy} ' \
            '--analysis_base_path {analysis_base_path} ' \
            '--adapter_trim_tool {adapter_trim_tool} ' \
            '--my_story ' \
            '--no_multiqc ' \
            '{order_number} '.format(
                copy='-c' if kargs['copy_rawdata'] else '',
                analysis_base_path=kargs['analysis_base_path'],
                adapter_trim_tool=kargs['adapter_trim_tool'],
                order_number=kargs['order_number'])
        run = run_cmd(prema_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'PreMA CORE 완료',
                'false_meg': 'PreMA CORE',
            }, p_exit=True)
        del prema_cmd, run
    prema_end_time = time.time()

    # Read Assembly
    if kargs['flash'] or kargs['full_course']:
        cofi_cmd = \
            'python3 /garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/CoFI.py ' \
            'FLASH ' \
            '--analysis_base_path {analysis_base_path} ' \
            '--read_length 251 ' \
            '--target_size 465 ' \
            '--length_min 430 ' \
            '--length_max 475 ' \
            '--pooled_sample False ' \
            '--stat False ' \
            '--no_compress ' \
            '--out_seq_type FASTA ' \
            '{order_number}'.format(
                analysis_base_path=ANALYSIS_BASE_PATH,
                order_number=kargs['order_number'])
        run = run_cmd(cofi_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'CoFI FLASH 완료',
                'false_meg': 'CoFI FLASH',
            }, p_exit=True)
        del cofi_cmd, run
    cofi_end_time = time.time()

    # Closed-reference OTU Picking
    if kargs['closed'] or kargs['full_course']:
        # 시료 목록 작성 및 Assemled FASTQ 파일 존재 유무 확인
        l_samples, analysis_path = check_sample_list(kargs['analysis_base_path'], kargs['order_number'],
                                                     kargs['target_dir_suffix'], kargs['no_check'])
        # 시료 설정
        if l_include or l_exclude:
            l_samples = set_sample_list(l_samples, l_include, l_exclude)
        assembly_path = glob_dir(analysis_path, 'Read_Assembly', 'only')
        exist_assembled_fastq(assembly_path,
                              ['F_{}'.format(sample) for sample in l_samples],
                              ASSEMBLY_SUFFIX)

        from myBiomeStory.qiime1_cmd import get_closed_reference_otu_cmd
        l_closed_otu_cmd = list()
        for sample in l_samples:
            assembled_fasta = os.path.join(*(assembly_path, 'F_' + sample + ASSEMBLY_SUFFIX))
            out_dir = os.path.join(*(analysis_path, sample))
            cmd = get_closed_reference_otu_cmd(assembled_fasta, out_dir)
            l_closed_otu_cmd.append(cmd)

        global POOL_WORKER
        launcher_cmd(l_closed_otu_cmd, 'Closed-reference OTU Picking', POOL_WORKER, False)
    closed_end_time = time.time()

    if kargs['diversity'] or kargs['summarize_taxa'] or kargs['full_course']:
        # 시료 목록화, 디렉터리 및 파일 존재 유무 확인.
        l_samples, analysis_path = check_sample_list(kargs['analysis_base_path'], kargs['order_number'],
                                                     kargs['target_dir_suffix'], kargs['no_check'])
        # 시료 설정
        if l_include or l_exclude:
            l_samples = set_sample_list(l_samples, l_include, l_exclude)
        exist_something(analysis_path, l_samples,
                        '', '>>> Closed-Reference OTU 디렉터리 확인')
        exist_something(analysis_path, l_samples,
                        '/otu_table.biom', '>>> biom 파일 확인')
        exist_something(analysis_path, l_samples,
                        '/rep_set', '>>> rep_set 디렉터리 확인')
        exist_something(analysis_path, l_samples,
                        '/uclust_assigned_taxonomy', '>>> uclust_assigned_taxonomy 디렉터리 확인')
        exist_something(analysis_path, l_samples,
                        '/uclust_ref_picked_otus', '>>> uclust_ref_picked_otus 디렉터리 확인')

    # alpha-diversity : Simpson, Shannon
    if kargs['diversity'] or kargs['full_course']:
        from myBiomeStory.qiime1_cmd import get_cmd_alpha_diversity_qiime1
        l_alpha_diversity_cmd = list()
        for sample in l_samples:
            biom = os.path.join(analysis_path, sample + '/otu_table.biom')
            output = os.path.join(analysis_path, sample + '/Alpha_Diversity.txt')
            cmd = get_cmd_alpha_diversity_qiime1(biom, output)
            l_alpha_diversity_cmd.append(cmd)
            o_recipe = open(os.path.join(analysis_path, sample + '/ALPHA_DIVERSITY.recipe'), 'w')
            click.echo(cmd, file=o_recipe)
            o_recipe.close()
        launcher_cmd(l_alpha_diversity_cmd, 'Alpha-Diversity', 30, False)
        del l_alpha_diversity_cmd, biom, output, cmd, o_recipe
    diversity_end_time = time.time()

    if kargs['summarize_taxa'] or kargs['full_course']:
        from myBiomeStory.qiime1_cmd import get_cmd_summarize_taxa_qiime1
        l_summarize_taxa_cmd = list()
        for sample in l_samples:
            biom = os.path.join(analysis_path, sample + '/otu_table.biom')
            output = os.path.join(analysis_path, sample + '/Taxonomy_Assignment')
            cmd = get_cmd_summarize_taxa_qiime1(biom, output)
            l_summarize_taxa_cmd.append(cmd)
            o_recipe = open(os.path.join(analysis_path, sample + '/SUMMARIZE_TAXA.recipe'), 'w')
            click.echo(cmd, file=o_recipe)
            o_recipe.close()
        launcher_cmd(l_summarize_taxa_cmd, 'Summarize_Taxa', 30, False)
    summarize_end_time = time.time()

    if kargs['make_story'] or kargs['full_course']:
        # 디렉터리 구조
        # analysis_base_path/Order_Number/시료명/Alpha_Diversity.txt
        #                                      /Taxonomy_Assignment/otu_table_L7.txt
        alpha_txt = 'Alpha_Diversity.txt'
        table_dir = 'Taxonomy_Assignment'
        table_txt = 'otu_table_L7.txt'
        l_samples, analysis_path = check_sample_list(kargs['analysis_base_path'], kargs['order_number'],
                                                     kargs['target_dir_suffix'], kargs['no_check'])
        # 시료 설정
        if l_include or l_exclude:
            l_samples = set_sample_list(l_samples, l_include, l_exclude)

        exist_something(analysis_path, l_samples,
                        f'/{alpha_txt}', '>>> Alpha_Diversity.txt 확인')
        exist_something(analysis_path, l_samples,
                        f'/{table_dir}/{table_txt}', '>>> otu_table_L7.txt 확인')

        from myBiomeStory.data_struct import MyBiomeStory, StoryIndex
        MyBiomeStory.read_bacteria_list()
        MyBiomeStory.read_mybiome_db()
        StoryIndex.read_tax_id()
        with click.progressbar(l_samples, label='Making StoryIndex',
                               fill_char=click.style('#', fg='green')) as bar_l_samples:
            for sample in bar_l_samples:
                alpha_file = os.path.join(analysis_path, sample, alpha_txt)
                table_file = os.path.join(analysis_path, sample, table_dir, table_txt)
                mybiome = MyBiomeStory(sample, table_file, alpha_file)
                mybiome.build_data_structure()
                mybiome.compute_data()
                mybiome.story_index.transform_index()
                mybiome.story_index.arrange_sample_microbiome(mybiome.d_bact_list)
                story_index_path = os.path.join(analysis_path, sample, 'Story_Index')
                os.mkdir(story_index_path)
                mybiome.story_index.save_points_to_csv(os.path.join(story_index_path, f'{sample}_points.csv'))
                mybiome.story_index.save_biome_to_csv(os.path.join(story_index_path, f'{sample}_biome.csv'))
    make_story_end_time = time.time()

    if kargs['insert'] or kargs['full_course']:
        l_samples, analysis_path = check_sample_list(kargs['analysis_base_path'], kargs['order_number'],
                                                     kargs['target_dir_suffix'], kargs['no_check'])
        # 시료 설정
        if l_include or l_exclude:
            l_samples = set_sample_list(l_samples, l_include, l_exclude)
        # Input File 확인
        exist_something(analysis_path, [f'{sample}/Story_Index/{sample}' for sample in l_samples],
                        '_points.csv', '>>> points.csv 확인')
        exist_something(analysis_path, [f'{sample}/Story_Index/{sample}' for sample in l_samples],
                        '_biome.csv', '>>> biome.csv 확인')
        sample_info_file = os.path.join(os.path.dirname(analysis_path), f'{kargs["order_number"]}_SampleInfo.csv')
        if check_file_type(sample_info_file, 'exists'):
            click.secho('>>> SampleInfo.csv 확인')
            click.echo(sample_info_file)
            l_sample_info = read_csv_to_dict(sample_info_file)
        else:
            click.secho('Warning: SampleInfo.csv 파일이 없습니다.', fg='yellow')
            if get_sample_info(kargs['order_number'], sample_info_file):
                click.secho('>>> SampleInfo.csv 생성 완료', fg='cyan')
                click.echo(sample_info_file)
                l_sample_info = read_csv_to_dict(sample_info_file)
            else:
                exit(1)

        # 성별 데이터 변환
        d_filtered_sample_info = dict()
        for info in l_sample_info:
            if info['KitId'] in l_samples:
                if info['Sex'] == '남성':
                    info['Sex'] = 'male'
                elif info['Sex'] == '여성':
                    info['Sex'] = 'female'
                else:
                    click.secho('SampleInfo.csv 파일의 성별 항목에 사용할 수 없는 값이 존재합니다.',
                                fg='red', err=True)
                    click.echo(f'Sex: {info["Sex"]}')
                d_filtered_sample_info[info['KitId']] = info

        from DBcontrol.dbMyBiomeStory import DBmyBiomeStory
        db_story = DBmyBiomeStory({'ip': '172.19.87.50',
                                   'user': 'mybiome',
                                   'password': 'Qnstjr3qn!',
                                   'db': 'my-biomestory'
                                   })
        date_fmt = '%Y-%m-%d %H:%M:%S'
        analysis_time = time.strftime(date_fmt, time.localtime())
        db_story.start_transaction()
        with click.progressbar(l_samples, label='Inserting StoryIndex',
                               fill_char=click.style('#', fg='green')) as bar_l_samples:
            for sample in bar_l_samples:
                points_csv = os.path.join(analysis_path, sample, 'Story_Index', f'{sample}_points.csv')
                biome_csv = os.path.join(analysis_path, sample, 'Story_Index', f'{sample}_biome.csv')
                try:
                    l_csv = read_csv(points_csv, p_verbose=False)
                    db_version = check_db_version(l_csv)
                    if db_version is False:
                        click.secho('Error: points.csv 파일에 여러 개의 DB 버전이 존재합니다.', fg='red')
                        click.echo(points_csv)
                        exit(1)
                except Exception as err:
                    click.secho('Error: DB Version 확인 중 Exception 발생', fg='red', blink=True)
                    click.echo(err)
                    db_story.rollback()
                try:
                    db_story.insert_analysis_version((0, sample, db_version, analysis_time, 0))
                    db_story.insert_sample_info(d_filtered_sample_info[sample], kargs['order_number'], time.time())
                    db_story.insert_point(points_csv)
                    db_story.insert_biome(biome_csv)
                except KeyError as err:
                    click.secho('Error: 데이터베이스에 데이터 입력 중 KeyError 발생', fg='red', blink=True)
                    click.echo(err)
                    db_story.rollback()
                except Exception as err:
                    click.secho('Error: 데이터베이스에 데이터 입력 중 Exception 발생', fg='red', blink=True)
                    click.echo(err)
                    db_story.rollback()
        db_story.commit()
        db_story.close()
    insert_end_time = time.time()
    print_run_time(start_time, prema_end_time, cofi_end_time, closed_end_time,
                   diversity_end_time, summarize_end_time, make_story_end_time, insert_end_time)


if __name__ == '__main__':
    my_main()
