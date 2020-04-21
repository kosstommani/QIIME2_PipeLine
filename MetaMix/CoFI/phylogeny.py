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
from SpoON.util import run_cmd, check_run_cmd
from click import secho, echo, style
import pdb



# QIIME2 = '/crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/envs/qiime2-2018.08/bin/qiime'
QIIME2 = 'qiime'
QIIME1 = 'make_phylogeny.py'

# qiime phylogeny
"""
Usage: qiime phylogeny [OPTIONS] COMMAND [ARGS]...

  Description: This QIIME 2 plugin supports generating and
  manipulating phylogenetic trees.

  Plugin website: https://github.com/qiime2/q2-phylogeny

  Getting user support: Please post to the QIIME 2 forum for help
  with this plugin: https://forum.qiime2.org

Options:
  --version    Show the version and exit.
  --citations  Show citations and exit.
  --help       Show this message and exit.

Commands:
  align-to-tree-mafft-fasttree  Build a phylogenetic tree using
                                fasttree and mafft alignment
  fasttree                      Construct a phylogenetic tree with
                                FastTree.
  filter-table                  Remove features from table if they're
                                not present in tree.
  iqtree                        Construct a phylogenetic tree with IQ-
                                TREE.
  iqtree-ultrafast-bootstrap    Construct a phylogenetic tree with IQ-
                                TREE with bootstrap supports.
  midpoint-root                 Midpoint root an unrooted phylogenetic
                                tree.
  raxml                         Construct a phylogenetic tree with
                                RAxML.
  raxml-rapid-bootstrap         Construct a phylogenetic tree with
                                bootstrap supports using RAxML.
"""


def get_align_to_tree_mafft_fasttree_cmd(kargs):
    """
    Alignment(MAFFT) & Phylogeny(FastTree)

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                i_sequences: str - FeatureData[Sequence]
                         The sequences to be used for creating a fasttree based rooted phylogenetic tree.
                out_dir: None or str
    :return: cmd
    """
    # """
    # Usage: qiime phylogeny align-to-tree-mafft-fasttree [OPTIONS]
    #
    #   This pipeline will start by creating a sequence alignment using
    #   MAFFT, after which any alignment columns that are phylogenetically
    #   uninformative or ambiguously aligned will be removed (masked). The
    #   resulting masked alignment will be used to infer a phylogenetic
    #   tree and then subsequently rooted at its midpoint. Output files
    #   from each step of the pipeline will be saved. This includes both
    #   the unmasked and masked MAFFT alignment from q2-alignment methods,
    #   and both the rooted and unrooted phylogenies from q2-phylogeny
    #   methods.
    #
    # Options:
    #   --i-sequences ARTIFACT PATH FeatureData[Sequence]
    #                                   The sequences to be used for
    #                                   creating a fasttree based rooted
    #                                   phylogenetic tree.  [required]
    #   --p-n-threads INTEGER RANGE     The number of threads. (Use 0 to
    #                                   automatically use all available
    #                                   cores) This value is used when
    #                                   aligning the sequences and creating
    #                                   the tree with fasttree.  [default:
    #                                   1]
    #   --p-mask-max-gap-frequency FLOAT
    #                                   The maximum relative frequency of
    #                                   gap characters in a column for the
    #                                   column to be retained. This relative
    #                                   frequency must be a number between
    #                                   0.0 and 1.0 (inclusive), where 0.0
    #                                   retains only those columns without
    #                                   gap characters, and 1.0 retains all
    #                                   columns  regardless of gap character
    #                                   frequency. This value is used when
    #                                   masking the aligned sequences.
    #                                   [default: 1.0]
    #   --p-mask-min-conservation FLOAT
    #                                   The minimum relative frequency of at
    #                                   least one non-gap character in a
    #                                   column for that column to be
    #                                   retained. This relative frequency
    #                                   must be a number between 0.0 and 1.0
    #                                   (inclusive). For example, if a value
    #                                   of  0.4 is provided, a column will
    #                                   only be retained  if it contains at
    #                                   least one character that is present
    #                                   in at least 40% of the sequences.
    #                                   This value is used when masking the
    #                                   aligned sequences.  [default: 0.4]
    #   --o-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   The aligned sequences.  [required if
    #                                   not passing --output-dir]
    #   --o-masked-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   The masked alignment.  [required if
    #                                   not passing --output-dir]
    #   --o-tree ARTIFACT PATH Phylogeny[Unrooted]
    #                                   The unrooted phylogenetic tree.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-rooted-tree ARTIFACT PATH Phylogeny[Rooted]
    #                                   The rooted phylogenetic tree.
    #                                   [required if not passing --output-
    #                                   dir]
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
    #
    # """
    cmd = '{qiime2} phylogeny align-to-tree-mafft-fasttree ' \
        '--i-sequences {qza} ' \
        '--output-dir {output_dir} ' \
        '--p-n-threads 0'.format(
            qiime2=QIIME2,
            qza=kargs['i_sequences'],
            output_dir=kargs['out_dir'],
        )
    return cmd


def run_align_to_tree_mafft_fasttree(kargs):
    """
    This pipeline will start by creating a sequence alignment using MAFFT, after which any alignment columns
    that are phylogenetically uninformative or ambiguously aligned will be removed (masked).
    The resulting masked alignment will be used to infer a phylogenetic tree and then subsequently rooted at
    its midpoint. Output files from each step of the pipeline will be saved. This includes both the unmasked
    and masked MAFFT alignment from q2-alignment methods, and both the rooted and unrooted phylogenies from
    q2-phylogeny methods.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                i_sequences: str - FeatureData[Sequence]
                         The sequences to be used for creating a fasttree based rooted phylogenetic tree.
                out_dir: None or str
    :return:
    """
    global QIIME2
    echo('>>> Alignment(MAFFT) & Phylogeny(FastTree) 시작')
    cmd = get_align_to_tree_mafft_fasttree_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: MAFFT & FastTree 완료',
            'false_meg': 'QIIME2 - phylogeny align-to-tree-mafft-fasttree',
        }
    )


def get_fasttree_cmd(p_args):
    """

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime phylogeny fasttree [OPTIONS]
    #
    #   Construct a phylogenetic tree with FastTree.
    #
    # Options:
    #   --i-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   Aligned sequences to be used for
    #                                   phylogenetic reconstruction.
    #                                   [required]
    #   --p-n-threads INTEGER RANGE     The number of threads. Using more
    #                                   than one thread runs the non-
    #                                   deterministic variant of `FastTree`
    #                                   (`FastTreeMP`), and may result in a
    #                                   different tree than single-
    #                                   threading. See http://www.microbeson
    #                                   line.org/fasttree/#OpenMP for
    #                                   details. (Use 0 to automatically use
    #                                   all available cores)  [default: 1]
    #   --o-tree ARTIFACT PATH Phylogeny[Unrooted]
    #                                   The resulting phylogenetic tree.
    #                                   [required if not passing --output-
    #                                   dir]
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
    cmd = '{qiime2} phylogeny fasttree ' \
        '--i-alignment {qza} ' \
        '--output-dir {output_dir} ' \
        '--p-n-threads 0'.format(
            qiime2=QIIME2,
            qza=p_args['qza'],
            output_dir=os.path.join(*(p_args['output_path'], 'phylogeny')),
        )
    return cmd


def get_midpoint_root_cmd(p_args):
    """

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime phylogeny midpoint-root [OPTIONS]
    #
    #   Midpoint root an unrooted phylogenetic tree.
    #
    # Options:
    #   --i-tree ARTIFACT PATH Phylogeny[Unrooted]
    #                                   The phylogenetic tree to be rooted.
    #                                   [required]
    #   --o-rooted-tree ARTIFACT PATH Phylogeny[Rooted]
    #                                   The rooted phylogenetic tree.
    #                                   [required if not passing --output-
    #                                   dir]
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
    cmd = '{qiime2} phylogeny midpoint-root ' \
        '--i-tree {qza} ' \
        '--o-rooted-tree {qzv}'.format(
            qiime2=QIIME2,
            qza=p_args['qza'],
            qzv=os.path.join(*(p_args['output_path'], 'rooted_tree.qza')),
        )
    return cmd


def get_iqtree_cmd(p_args):
    """

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime phylogeny iqtree [OPTIONS]
    #
    #   Construct a phylogenetic tree using IQ-TREE
    #   (http://www.iqtree.org/) with automatic model selection.
    #
    # Options:
    #   --i-alignment ARTIFACT PATH FeatureData[AlignedSequence]
    #                                   Aligned sequences to be used for
    #                                   phylogenetic reconstruction.
    #                                   [required]
    #   --p-seed INTEGER                Random number seed. If not set,
    #                                   program defaults will be used. See
    #                                   IQ-TREE manual for details.
    #                                   [optional]
    #   --p-n-cores INTEGER RANGE       The number of cores to use for
    #                                   parallel processing. Use '0' to let
    #                                   IQ-TREE automatically determine the
    #                                   optimal number of cores to use.
    #                                   [default: 1]
    #   --p-substitution-model [TN+I|TVM+R4|HKY+I+G|K81u+I+G|TPM3+R2|K80+R6|TIM+I+G|GTR+R7|JC+I|TVM+R5|
    #                           HKY+R2|F81+R7|TPM3+R9|TVM+R6|TPM2u+R7|JC+R4|TIM+R9|GTR+R4|F81+R2|TIM2+R2|
    #                           SYM+R7|TPM2|K81u+R7|K80+R9|TPM3+I|GTR+G|K81+R10|K80|HKY+R7|TVM+G|
    #                           TVM+R7|GTR+R5|TIM2e+R6|HKY+R8|TIM2|K81+G|TNe+R3|TIM3+R9|TIMe+R8|TVMe+R10|
    #                           TPM2u+I+G|TIM3+I+G|F81+R10|K81u+R2|TPM2u+R8|TPM2+R3|TIM+R7|TIM+R2|TIM+R10|TIM2e+R2|
    #                           TIM2+G|K80+R10|K80+R4|K81u+R9|K81+R6|TIM2+I+G|K81u+R10|TPM2u+R3|TIM|TPM2+I+G|
    #                           TIM2e+R7|TVMe+R8|TVMe+R6|GTR+I|TPM3+I+G|TIMe+R4|TNe|TPM2+R4|F81+R3|TNe+R2|
    #                           TN+R4|TPM3u|K80+G|TIM3e+R4|K81u+G|K81+R7|TVMe+G|TVM+R2|SYM+I|TIMe+R2|
    #                           TIM+R3|K81+I+G|K81+R3|TIM2+R7|TIM3e+R8|K81u+R8|F81+R9|TIM3e|TIM3e+R9|F81+R8|
    #                           K81+R5|K80+I|TPM3+R10|TPM3u+I|TIM3e+I|HKY|GTR+R10|SYM+R9|K80+R8|TN+R10|
    #                           TPM2u+R9|TIM2e+R10|TIM3+R2|TPM2u+G|TPM2u+R4|TIM3e+R10|HKY+R9|SYM+R5|TN+R9|TIM2+R3|
    #                           K80+R2|TIM2+R4|SYM+G|TPM3u+R3|TN+R8|TPM2+I|TIM3|TVM+I|JC+I+G|TPM3+R6|
    #                           TIM2e+R5|TVM+R3|TIM3e+R6|TIM+I|TIMe+I|TPM3+R4|TIM+R5|TIM2+R10|K81u+R6|TPM2+R10|
    #                           TVMe+R3|SYM+R2|TPM3u+R10|JC+R9|TIM2e+R8|TVM+R10|F81+G|JC|HKY+R3|HKY+G|
    #                           TIM2e+G|TIM3e+R3|TIM3e+R7|TPM2+G|K80+R5|HKY+R4|TVMe|TIM+G|TN+R5|K81u|
    #                           TPM3u+R9|TIM2e+R3|TPM2u+R5|SYM+R3|F81+R6|TN|TN+I+G|TVM|F81+I+G|TPM2+R9|
    #                           K81u+R5|TIM3+G|TIM3e+I+G|TPM2u+R2|GTR+R8|K81u+I|JC+R3|GTR+R3|JC+R6|TIM3+R5|
    #                           JC+R5|TIMe+R3|GTR+R2|TIM3+R8|TNe+I|TPM3u+R7|TIM2e|TNe+R9|TPM2u+R6|GTR+R6|
    #                           TVM+R8|K80+I+G|TIM+R6|GTR+R9|HKY+R10|K81+I|K81u+R4|F81+R5|TIMe+G|TNe+R8|
    #                           TIM2+R8|SYM+I+G|GTR+I+G|JC+G|HKY+I|K81+R4|K81+R9|TNe+R4|TPM3u+R5|TPM3u+R2|
    #                           TIMe|TIM+R4|TPM3+R7|TIM2e+I|TIM2e+R4| HKY+R5|TVMe+R4|TPM2u|MFP|TPM3+R5|
    #                           TPM2u+I|TPM2u+R10|TIM3+I|TPM3u+R8|TIM2+R5|TVMe+R5|TPM3+R8|TIM3e+R5|SYM+R8|TPM3u+I+G|
    #                           SYM+R10|TPM3u+G|TPM2+R2|TVM+R9|SYM+R4|TIMe+R5|TIM2e+R9|TEST|TPM3+R3|TIMe+I+G|
    #                           TIM2e+I+G|TN+R3|TN+R7|TVM+I+G|F81|TIMe+R6|HKY+R6|TIM3+R6|TIM3+R10|K81+R2|
    #                           GTR|TN+R6|TVMe+I|TIM3+R4|TPM2+R8|TIM3e+R2|SYM+R6|JC+R7|F81+I|TPM3u+R4|
    #                           TN+G|TIM2+I|SYM|TIM3e+G|F81+R4|JC+R2|JC+R8|K80+R7|TPM2+R7|K81|
    #                           K81+R8|TVMe+I+G|K81u+R3|TIMe+R9|TNe+R6|TVMe+R9|TPM3|TPM2+R6|TIM3+R7|TVMe+R2|
    #                           TIM2+R9|JC+R10|TNe+I+G|K80+R3|TN+R2|TPM2+R5|TPM3+G|TPM3u+R6|TIMe+R10|TIM3+R3|
    #                           TIMe+R7|TVMe+R7|TIM+R8|TIM2+R6|TNe+R10|TNe+G|TNe+R7|TNe+R5]
    #                                   Model of Nucleotide Substitution. If
    #                                   not provided, IQ-TREE will determine
    #                                   the best fit substitution model
    #                                   automatically.  [default: MFP]
    #   --p-n-init-pars-trees INTEGER RANGE
    #                                   Number of initial parsimony trees.
    #                                   If not set, program defaults will be
    #                                   used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-n-top-init-trees INTEGER RANGE
    #                                   Number of top initial trees. If not
    #                                   set, program defaults will be used.
    #                                   See IQ-TREE manual for details.
    #                                   [optional]
    #   --p-n-best-retain-trees INTEGER RANGE
    #                                   Number of best trees retained during
    #                                   search. If not set, program defaults
    #                                   will be used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-n-iter INTEGER RANGE        Fix number of iterations to stop. If
    #                                   not set, program defaults will be
    #                                   used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-stop-iter INTEGER RANGE     Number of unsuccessful iterations to
    #                                   stop. If not set, program defaults
    #                                   will be used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-perturb-nni-strength FLOAT  Perturbation strength for randomized
    #                                   NNI. If not set, program defaults
    #                                   will be used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-spr-radius INTEGER RANGE    Radius for parsimony SPR search. If
    #                                   not set, program defaults will be
    #                                   used. See IQ-TREE manual for
    #                                   details.  [optional]
    #   --p-allnni / --p-no-allnni      Perform more thorough NNI search.
    #                                   [default: False]
    #   --p-safe / --p-no-safe          Safe likelihood kernel to avoid
    #                                   numerical underflow  [default:
    #                                   False]
    #   --o-tree ARTIFACT PATH Phylogeny[Unrooted]
    #                                   The resulting phylogenetic tree.
    #                                   [required if not passing --output-
    #                                   dir]
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
    #
    # """
    cmd = '{qiime2} phylogeny iqtree ' \
        '--i-alignment {qza} ' \
        '--output-dir {output_dir} ' \
        '--p-n-cores 0'.format(
            qiime2=QIIME2,
            qza=p_args['qza'],
            output_dir=os.path.join(*(p_args['output_path'], 'phylogeny')),
        )
    return cmd


def get_phylogeny_qiime1_cmd(p_fasta, p_out_file):
    # make_phylogeny.py -i filtered_alignment/otus_rep_aligned_pfiltered.fasta -o phylogeny/rep_phylo.tre
    cmd = '{qiime1} ' \
        '-i {filtered_align} ' \
        '-o {out_file}'.format(
            qiime1=QIIME1,
            filtered_align=p_fasta,
            out_file=p_out_file,
        )
    return cmd


def run_phylogeny_qiime1(p_fasta, p_out_path):
    """
    OTU 대표 서열들간의 Multiple Alignment 결과를 바탕으로 Phylogeny Tree 를 생성한다.
    Tree 는 FastTree method를 이용한다(Price, Dehal, & Arkin, 2009).

    :type p_fasta: str
    :param p_fasta: filtered_alignment 디렉터리의 otus_rep_aligned_pfiltered.fasta
    :type p_out_path: str
    :param p_out_path: 결과를 생성할 경로
    :return:
    """
    out_file = os.path.join(p_out_path, 'rep_phylo.tre')
    cmd = get_phylogeny_qiime1_cmd(p_fasta, out_file)
    with open(os.path.join(p_out_path, 'PHYLOGENY_CMD.recipe'), 'w') as o_recipe:
        o_recipe.write(cmd)
        o_recipe.write('\n')
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME1 Phylogeny 완료',
            'false_meg': 'QIIME1 make_phylogeny',
        }
    )

