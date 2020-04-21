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
__version__ = '1.0.0'

import os
from time import sleep
from click import echo
from SpoON.util import run_cmd, check_run_cmd, check_file_type, \
    run_cmd_qsub, check_jobs_in_queue, parse_job_id

# QIIME2 = '/crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/envs/qiime2-2018.11/bin/qiime'
QIIME2 = 'qiime'
QIIME1_ALIGN = 'align_seqs.py'
QIIME1_FILTER_ALIGN = 'filter_alignment.py'

# qiime alignment
"""
Usage: qiime alignment [OPTIONS] COMMAND [ARGS]...

  Description: This QIIME 2 plugin provides support for generating
  and manipulating sequence alignments.

  Plugin website: https://github.com/qiime2/q2-alignment

  Getting user support: Please post to the QIIME 2 forum for help
  with this plugin: https://forum.qiime2.org

Options:
  --version    Show the version and exit.
  --citations  Show citations and exit.
  --help       Show this message and exit.

Commands:
  mafft  De novo multiple sequence alignment with MAFFT
  mask   Positional conservation and gap filtering.
"""


def get_alignment_mafft_cmd(p_args):
    # """
    # Usage: qiime alignment mafft [OPTIONS]
    #
    #   Perform de novo multiple sequence alignment using MAFFT.
    #
    # Options:
    #   --i-sequences ARTIFACT PATH FeatureData[Sequence]
    #                                   The sequences to be aligned.
    #                                   [required]
    #   --p-n-threads INTEGER RANGE     The number of threads. (Use 0 to
    #                                   automatically use all available
    #                                   cores)  [default: 1]
    #   --p-parttree / --p-no-parttree  This flag is required if the number
    #                                   of sequences being aligned are
    #                                   larger than 1000000. Disabled by
    #                                   default  [default: False]
    #   --o-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   The aligned sequences.  [required if
    #                                   not passing --output-dir]
    #   --output-dir DIRECTORY          Output unspecified results to a
    #                                   directory
    #   --cmd-config PATH               Use config file for command options
    #   --verbose                       Display verbose output to stdout
    #                                   and/or stderr during execution of
    #                                   this action.  [default: False]
    #   --quiet                         Silence output if execution is
    #                                   successful (silence is golden).
    #                                   [default: False]
    #   --citations                     Show citations and exit.
    #   --help                          Show this message and exit.
    # """
    global QIIME2
    cmd = '{qiime2} alignment mafft ' \
        '--i-sequences {qza} ' \
        '--output-dir {output_dir} ' \
        '--p-n-threads 0'.format(
            qiime2=QIIME2,
            qza=p_args['qza'],
            qzv=os.path.join(*(p_args['output_path'], 'alignment_MAFFT')),
        )
    return cmd


def run_alignment_mafft(p_args):
    cmd = get_alignment_mafft_cmd(p_args)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2 alignment - MAFFT 완료',
            'false_meg': 'QIIME2 alignment mafft'
        }
    )


def get_alignment_mask_cmd(p_args):
    # """
    # Usage: qiime alignment mask [OPTIONS]
    #
    #   Mask (i.e., filter) unconserved and highly gapped columns from an
    #   alignment. Default min_conservation was chosen to reproduce the
    #   mask presented in Lane (1991).
    #
    # Options:
    #   --i-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   The alignment to be masked.
    #                                   [required]
    #   --p-max-gap-frequency FLOAT     The maximum relative frequency of
    #                                   gap characters in a column for the
    #                                   column to be retained. This relative
    #                                   frequency must be a number between
    #                                   0.0 and 1.0 (inclusive), where 0.0
    #                                   retains only those columns without
    #                                   gap characters, and 1.0 retains all
    #                                   columns regardless of gap character
    #                                   frequency.  [default: 1.0]
    #   --p-min-conservation FLOAT      The minimum relative frequency of at
    #                                   least one non-gap character in a
    #                                   column for that column to be
    #                                   retained. This relative frequency
    #                                   must be a number between 0.0 and 1.0
    #                                   (inclusive). For example, if a value
    #                                   of 0.4 is provided, a column will
    #                                   only be retained if it contains at
    #                                   least one character that is present
    #                                   in at least 40% of the sequences.
    #                                   [default: 0.4]
    #   --o-masked-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   The masked alignment.  [required if
    #                                   not passing --output-dir]
    #   --output-dir DIRECTORY          Output unspecified results to a
    #                                   directory
    #   --cmd-config PATH               Use config file for command options
    #   --verbose                       Display verbose output to stdout
    #                                   and/or stderr during execution of
    #                                   this action.  [default: False]
    #   --quiet                         Silence output if execution is
    #                                   successful (silence is golden).
    #                                   [default: False]
    #   --citations                     Show citations and exit.
    #   --help                          Show this message and exit.
    # """
    global QIIME2
    cmd = '{qiime2} alignment mask ' \
        '--i-alignment {qza} ' \
        '--o-masked-alignment {qzv}'.format(
            qiime2=QIIME2,
            qza=p_args['qza'],
            qzv=p_args['qzv'],
        )
    return cmd


def run_alignment_mask(p_args):
    cmd = get_alignment_mask_cmd(p_args)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2 alignment - mask 완료',
            'false_meg': 'QIIME2 alignment mask'
        }
    )


def get_alignment_muscle_cmd(p_otus_rep, p_out_path):
    """
    QIIME1 프로그램의 Multiple Alignment(method: muscle) 의 실행 명령어를 반환한다.

    :type p_otus_rep: str
    :param p_otus_rep: OTU 대표 서열의 파일명. ex) otus_rep.fasta
    :type p_out_path: str
    :param p_out_path: 결과를 생성할 경로
    :rtype: str
    :return: 실행 명려어
    """
    # align_seqs.py -i otus_rep.fasta -m muscle -o muscle_alignment
    global QIIME1_ALIGN
    cmd = '{qiime1} ' \
        '-i {otus_rep_fasta} ' \
        '-m muscle ' \
        '-o {out_dir}'.format(
            qiime1=QIIME1_ALIGN,
            otus_rep_fasta=p_otus_rep,
            out_dir=p_out_path,
        )
    return cmd


def get_filter_alignment_cmd(p_aligned_fasta, p_out_dir):
    # filter_alignment.py -i muscle_alignment/otus_rep_aligned.fasta
    # -o filtered_alignment/ -e 0.10 --suppress_lane_mask_filter
    global QIIME1_FILTER_ALIGN
    cmd = '{qiime1} ' \
        '-i {aligned_fasta} ' \
        '-o {out_dir} ' \
        '-e 0.10 ' \
        '--suppress_lane_mask_filter'.format(
            qiime1=QIIME1_FILTER_ALIGN,
            aligned_fasta=p_aligned_fasta,
            out_dir=p_out_dir,
        )
    return cmd


def run_alignment_muscle(p_otus_rep: str, p_out_path: str, p_queue: str, p_no_queue: bool) -> str:
    """
    OTU 대표 서열들에 대해서 Multiple Alignment을 진행한다.
    QIIME1 프로그램의 align_seqs.py 실행한다.


    :type p_otus_rep: str
    :param p_otus_rep: OTU 대표 서열의 파일명. ex) otus_rep.fasta
    :type p_out_path: str
    :param p_out_path: 결과를 생성할 경로
    :param p_queue: Queue 이름
    :param p_no_queue: Queue 사용 여부
    :rtype: str
    :return: Multiple Alignment가 생성된 경로
    """
    out_dir = os.path.join(p_out_path, 'muscle_alignment')
    cmd = get_alignment_muscle_cmd(p_otus_rep, out_dir)
    recipe_file = os.path.join(p_out_path, 'ALIGNMENT_CMD.recipe')
    with open(recipe_file, 'w') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')

    if p_no_queue is True:
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'QIIME1 Multiple Alignment(muscle) 완료',
                'false_meg': 'QIIME1 align_seqs'
            }
        )
    else:
        run_cmd_qsub(
            {
                'stream': 'y',
                'queue': p_queue,
                'interpreter': '/bin/bash',
                'out_path': p_out_path,
                'cmd': recipe_file,
            }
        )
        # sleep(5)
        job_id = parse_job_id(p_out_path, 'qsub.log', 'only')
        check_jobs_in_queue([job_id])
        echo('\n')
    return out_dir


def run_filter_alignment(p_aligned, p_out_path):
    """

    :type p_aligned: str
    :param p_aligned:
    :type p_out_path: str
    :param p_out_path: 결과를 생성할 경로
    :return:
    """
    out_dir = os.path.join(p_out_path, 'filtered_alignment')
    cmd = get_filter_alignment_cmd(p_aligned, out_dir)
    with open(os.path.join(p_out_path, 'ALIGNMENT_CMD.recipe'), 'a') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')
        o_recipe.write('\n')
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME1 Filtered Multiple Alignment(muscle) 완료',
            'false_meg': 'QIIME1 filtered_alignment'
        }
    )

