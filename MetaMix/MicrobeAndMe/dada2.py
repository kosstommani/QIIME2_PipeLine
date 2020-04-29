# ----------------------------------------------------------------------------------------------------------------------
#                                      888b     d888          888             888b     d888 d8b
#                                      8888b   d8888          888             8888b   d8888 Y8P
#                                      88888b.d88888          888             88888b.d88888
# 88888b.d88b.   8888b.  88888b.d88b.  888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b     "88b 888 "888 "88b 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 .d888888 888  888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 888  888 888  888  888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888 "Y888888 888  888  888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.0'

import os
import sys
from SpoON.util import run_cmd, check_run_cmd, parse_config, \
    check_queue_node, run_cmd_qsub, check_jobs_in_queue, parse_job_id
from click import secho, echo
from time import sleep

CONFIG = parse_config()
R_SCRIPT = CONFIG['mamMetaMix']['Rscript']


def get_dada2_r_cmd(p_kargs):
    """

    :param p_kargs:
            order_number:
            analysis_number:
            run_dir:
            analysis_base_path:
            multithread:
            sample_name:
    :return:
    """
    global R_SCRIPT
    dada2_r = '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/dada2.R'
    cmd = f'{R_SCRIPT} {dada2_r} ' \
          f'{p_kargs["order_number"]} ' \
          f'{p_kargs["analysis_number"]} ' \
          f'{p_kargs["run_dir"]} ' \
          f'--analysis_base_path {p_kargs["analysis_base_path"]} ' \
          f'--multithread {p_kargs["multithread"]} ' \
          f'--sample_name {p_kargs["sample_name"]} ' \
          f'--trimLeftF 17 ' \
          f'--trimLeftR 21 ' \
          f'--truncLenF 245 ' \
          f'--truncLenR 245'
    return cmd


def run_cmd_qsub_parallel(p_kargs):
    """

    :param p_kargs:
            queue: [bi3m.q|meta.q]
            queue_mode: [auto|fix]
            sample_name
            run_file:
            log_out_path:
            slots:
    :return: 
    """
    if p_kargs['queue'] == 'bi3m.q':
        if p_kargs['slots'] > 56:
            slots = 56
        else:
            slots = p_kargs['slots']
    elif p_kargs['queue'] == 'meta.q':
        if p_kargs['slots'] > 48:
            slots = 48
        else:
            slots = p_kargs['slots']
    qsub_run = os.path.join(p_kargs['log_out_path'], f'{p_kargs["sample_name"]}_qsub.run')
    qsub_log = os.path.join(p_kargs['log_out_path'], f'{p_kargs["sample_name"]}_qsub.log')
    cmd = 'qsub ' \
          '-j y ' \
          '-V ' \
          '-S /bin/sh ' \
          f'-N Job_{p_kargs["sample_name"]} ' \
          f'-pe peBWA {slots} ' \
          f'-q {p_kargs["queue"]} ' \
          f'-o {qsub_run} ' \
          f'{p_kargs["run_file"]} >> {qsub_log}'
    qsub_recipe_file = os.path.join(p_kargs['log_out_path'], f'{p_kargs["sample_name"]}_qsub.recipe')
    with open(qsub_recipe_file, 'w') as o_qsub:
        o_qsub.write(cmd)

    queue_state, node = check_queue_node()
    if queue_state is True:
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_stderr=False, p_exit=False,
        )
        return run.returncode
    else:
        secho('Error: Queue 사용할 수 없는 서버입니다.', fg='red', blink=True, err=True)
        echo(f'서버: {node}', err=True)
        exit(1)


def run_dada2_r(p_kargs):
    """

    :param p_kargs:
            order_number:
            analysis_number:
            run_dir:
            analysis_base_path:
            sample_name:
            queue: [bi3m.q|meta.q]
            no_queue:
            slots:
    :return:
    """
    dada2_r_cmd = get_dada2_r_cmd(
        {
            'order_number': p_kargs['order_number'],
            'analysis_number': p_kargs['analysis_number'],
            'run_dir': p_kargs['run_dir'],
            'analysis_base_path': p_kargs['analysis_base_path'],
            'multithread': p_kargs['slots'],
            'sample_name': p_kargs['sample_name']
        }
    )
    analysis_number_path = os.path.join(p_kargs['analysis_base_path'], p_kargs['order_number'],
                                        p_kargs['analysis_number'])
    dada2_r_file = os.path.join(analysis_number_path, f'{p_kargs["sample_name"]}_dada2.recipe')
    with open(dada2_r_file, 'w') as o_dada2:
        o_dada2.write(dada2_r_cmd)
        o_dada2.write('\n')

    if p_kargs['no_queue'] is True:
        echo('\n')
        secho('------ R - DADA2 ------ ', fg='yellow')
        run = run_cmd(dada2_r_cmd, p_stdout=sys.stdout, p_stderr=sys.stderr)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_stderr=False, p_exit=False,
        )
        return run.returncode
    else:
        return_code = run_cmd_qsub_parallel(
            {
                'queue': p_kargs['queue'],
                'sample_name': p_kargs['sample_name'],
                'run_file': dada2_r_file,
                'log_out_path': analysis_number_path,
                'slots': p_kargs['slots'],
            }
        )
        return return_code


def get_merge_asvs_cmd(p_rds_dir: str, p_rds_name: str) -> str:
    global R_SCRIPT
    merge_asvs_r = '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/merge_ASVs.R'
    cmd = f'{R_SCRIPT} {merge_asvs_r} ' \
          f'{p_rds_dir} ' \
          f'--rds_name {p_rds_name}'
    return cmd


def run_merge_asvs(p_rds_dir: str, p_rds_name: str, p_queue: str, p_no_queue: bool):
    cmd = get_merge_asvs_cmd(p_rds_dir, p_rds_name)
    out_path = os.path.dirname(p_rds_dir)
    cmd_file = os.path.join(out_path, 'merge_ASVs.recipe')
    with open(cmd_file, 'w') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
    if p_no_queue is False:
        run_cmd_qsub(
            {
                'stream': 'y',
                'queue': p_queue,
                'interpreter': '/bin/bash',
                'out_path': out_path,
                'cmd': cmd_file,
            }
        )
        # sleep(5)
        job_id = parse_job_id(out_path, 'qsub.log', 'only')
        check_jobs_in_queue([job_id])
        echo()  # 빈 줄
    else:
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'ASVs 통합 완료',
                'false_meg': 'run_merge_asvs',
            }
        )
