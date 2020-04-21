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
from SpoON.util import run_cmd, check_run_cmd
from click import secho

RECIPE = 'BIOM_CMD.recipe'


def get_make_otu_table_cmd(p_otu_map, p_biom, p_taxonomy=None, p_metadata=None):
    option_cmd = ''
    if p_taxonomy is not None:
        option_cmd += f'--taxonomy {p_taxonomy} '
    if p_metadata is not None:
        option_cmd += f'--mapping_fp {p_metadata}'

    cmd = 'make_otu_table.py ' \
        '--otu_map_fp {otu_map} ' \
        '--output_biom_fp {biom} ' \
        '{option_cmd}'.format(
            otu_map=p_otu_map,
            biom=p_biom,
            option_cmd=option_cmd)
    return cmd


def get_sort_otu_table_cmd(p_biom, p_sorted_biom, p_order):
    cmd = 'sort_otu_table.py ' \
        '-i {biom} ' \
        '-o {sorted_biom} ' \
        '-l {order_list}'.format(
            biom=p_biom,
            sorted_biom=p_sorted_biom,
            order_list=p_order,
        )
    return cmd


def get_biom_convert_cmd(p_in, p_out):
    cmd = 'biom convert ' \
        '-i {hdf5} ' \
        '-o {json} ' \
        '--to-json'.format(
            hdf5=p_in,
            json=p_out
        )
    return cmd


def get_biom_summarize_cmd(p_in, p_out):
    cmd = 'biom summarize-table ' \
        '-i {biom} ' \
        '-o {summary}'.format(
            biom=p_in,
            summary=p_out
    )
    return cmd


def make_sample_id_list_for_sorting(p_metadata):
    path = os.path.split(p_metadata)[0]
    sample_id_file = os.path.join(path, 'sample_list.txt')
    cmd = f'cut -f 1 {p_metadata} > {sample_id_file} '
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': None,
            'false_meg': 'Sample List 파일(고객맞춤정렬) 생성'
        }, p_stdout=False, p_exit=False,
    )
    if run.returncode == 0:
        return sample_id_file
    else:
        return False


def run_make_otu_table(p_path, p_otu_map, p_biom_base, p_taxonomy, p_metadata):
    global RECIPE
    hdf5_biom = f'{p_biom_base}.HDF5.biom'
    hdf5_biom_fp = os.path.join(p_path, hdf5_biom)
    hdf5_biom_cmd = get_make_otu_table_cmd(p_otu_map, hdf5_biom_fp, p_taxonomy, p_metadata)
    run = run_cmd(hdf5_biom_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME1 Make OTU Table(BIOM) 완료',
            'false_meg': 'QIIME1 make_otu_table',
        }
    )

    sample_id_file = make_sample_id_list_for_sorting(p_metadata)
    if sample_id_file is False:
        secho('Error: 고객맞춤정렬을 위한 Sample List 파일이 생성되지 않았습니다.', fg='red', blink=True)
        exit()
    sorted_biom = f'{p_biom_base}.sorted.biom'
    sorted_biom_fp = os.path.join(p_path, sorted_biom)
    sorted_biom_cmd = get_sort_otu_table_cmd(hdf5_biom_fp, sorted_biom_fp, p_metadata)
    run = run_cmd(sorted_biom_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME1 Sort OTU Table(BIOM) 완료',
            'false_meg': 'QIIME1 sort_otu_table',
        }
    )

    json_biom = f'{p_biom_base}.biom'
    json_biom_fp = os.path.join(p_path, json_biom)
    json_biom_cmd = get_biom_convert_cmd(sorted_biom_fp, json_biom_fp)
    run = run_cmd(json_biom_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'biom convert(JSON) 완료',
            'false_meg': 'biom convert',
        }
    )
    summary_txt = f'{p_biom_base}_summary.txt'
    summary_txt_fp = os.path.join(p_path, summary_txt)
    summary_txt_cmd = get_biom_summarize_cmd(sorted_biom_fp, summary_txt_fp)
    run = run_cmd(summary_txt_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'biom summarize-table 완료',
            'false_meg': 'biom summarize-table'
        }
    )

    with open(os.path.join(p_path, RECIPE), 'a') as o_recipe:
        o_recipe.write(f'#------------------------------ {p_biom_base} ------------------------------')
        o_recipe.write('\n')
        o_recipe.write(hdf5_biom_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
        o_recipe.write(sorted_biom_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
        o_recipe.write(json_biom_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
        o_recipe.write(summary_txt_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')


def get_collapse_samples_cmd(p_kargs):
    """

    ex) collapse_samples.py --collapse_fields Test,OK
    #SampleID	Test	OK
    200	1	A
    305	1	B
    Free	1	A
    MLE	2	A
    Methanol	2	B
    SMM.IAR	2	A

    #SampleID	original-sample-ids
    1.A	(200, Free)
    1.B	305
    2.A	(MLE, SMM.IAR)
    2.B	Methanol

    :type p_kargs: dict
    :param p_kargs: 다음을 키로 가지는 딕션너리
            biom
            metadtaa
            field
            mode
            out_biome
            out_metadata
    """
    cmd = 'collapse_samples.py ' \
        '-b {biom} ' \
        '-m {metadata} ' \
        '--collapse_fields {field} ' \
        '--collapse_mode {mode} ' \
        '--output_biom_fp {out_biom} ' \
        '--output_mapping_fp {out_metadata}'.format(
            biom=p_kargs['biom'],
            metadata=p_kargs['metadata'],
            field=p_kargs['field'],
            mode=p_kargs['mode'],
            out_biom=p_kargs['out_biom'],
            out_metadata=p_kargs['out_metadata'],
        )
    return cmd


def run_collapse_samples(p_kargs):
    """
    각 그룹의 해당되는 시료들을 그룹명으로 통합한다.
    통합하는 방법에는 mean|sum|random|median|first 가 있다.
    collapsed biom 파일을 만들어지는 절차는 다음과 같다.
    biom(BIOM 디렉터리) --> collapsed HDF5 biom --> collapsed JSON biom
                                               --> summary.txt

    :param p_kargs:
                biom: otu_table.biom
                metadata: Metadata.txt 파일
                field: Metadata 파일의 그룹 열명
                mode: 통합 방법. mean|sum|random|median|first
                out_path:
                out_biom_base: collapsed biom 파일명의 기준이름.
                out_metadata
    :return:
    """
    global RECIPE
    out_path = p_kargs.pop('out_path')
    out_biom_base = p_kargs.pop('out_biom_base')
    collapsed_biom = f'{out_biom_base}.HDF5.biom'
    collapsed_biom_fp = os.path.join(out_path, collapsed_biom)
    p_kargs.update({'out_biom': collapsed_biom_fp})
    collapsed_biom_cmd = get_collapse_samples_cmd(p_kargs)
    run = run_cmd(collapsed_biom_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME1 Collapse Sample(BIOM) 완료',
            'false_meg': 'QIIME1 collapse_samples'
        }
    )

    json_biom = f'{out_biom_base}.biom'
    json_biom_fp = os.path.join(out_path, json_biom)
    json_biom_cmd = get_biom_convert_cmd(collapsed_biom_fp, json_biom_fp)
    run = run_cmd(json_biom_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'biom convert(JSON) 완료',
            'false_meg': 'biom conver',
        }
    )

    summary_txt = f'{out_biom_base}_summary.txt'
    summary_txt_fp = os.path.join(out_path, summary_txt)
    summary_txt_cmd = get_biom_summarize_cmd(json_biom_fp, summary_txt_fp)
    run = run_cmd(summary_txt_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'biom summarize-table 완료',
            'false_meg': 'biom summarize-table',
        }
    )

    out_path = os.path.dirname(p_kargs['out_biom'])
    out_biom_name = os.path.basename(p_kargs['out_biom'])
    with open(os.path.join(out_path, RECIPE), 'a') as o_recipe:
        o_recipe.write(f'#-------------------- {out_biom_name} --------------------')
        o_recipe.write('\n')
        o_recipe.write(collapsed_biom_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
        o_recipe.write(json_biom_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
        o_recipe.write(summary_txt_cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
