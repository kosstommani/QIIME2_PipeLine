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
from SpoON.util import run_cmd, check_run_cmd, parse_config
from click import secho, echo, style
import pdb

CONFIG = parse_config()
QIIME2 = CONFIG['CoFI_QIIME2']


# qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file mappingfile.txt --o-visualization faith_pd_group_significance.qzv
# qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file mappingfile.txt --m-metadata-column Kit --o-visualization core-metrics-results/unweighted_unifrac_Kit-significance.qzv --p-pairwise
# qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza --m-metadata-file mappingfile.txt --m-metadata-column Kit --o-visualization core-metrics-results/weighted_unifrac_Kit-significance.qzv --p-pairwise
# qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny phylogeny/rooted-tree.qza --p-max-depth 892 --m-metadata-file mappingfile.txt --o-visualization alpha-rarefaction.qzv


# """
# Usage: qiime diversity [OPTIONS] COMMAND [ARGS]...
#
#   Description: This QIIME 2 plugin supports metrics for calculating
#   and exploring community alpha and beta diversity through
#   statistics and visualizations in the context of sample metadata.
#
#   Plugin website: https://github.com/qiime2/q2-diversity
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
#   alpha                      Alpha diversity
#   alpha-correlation          Alpha diversity correlation
#   alpha-group-significance   Alpha diversity comparisons
#   alpha-phylogenetic         Alpha diversity (phylogenetic)
#   alpha-rarefaction          Alpha rarefaction curves
#   beta                       Beta diversity
#   beta-correlation           Beta diversity correlation
#   beta-group-significance    Beta diversity group significance
#   beta-phylogenetic          Beta diversity (phylogenetic)
#   beta-rarefaction           Beta diversity rarefaction
#   bioenv                     bioenv
#   core-metrics               Core diversity metrics (non-phylogenetic)
#   core-metrics-phylogenetic  Core diversity metrics (phylogenetic and
#                              non-phylogenetic)
#   filter-distance-matrix     Filter samples from a distance matrix.
#   mantel                     Apply the Mantel test to two distance
#                              matrices
#   pcoa                       Principal Coordinate Analysis
#   pcoa-biplot                Principal Coordinate Analysis Biplot
#   procrustes-analysis        Procrustes Analysis
# """

def get_core_metrics_phylogenetic_cmd(p_args):
    """

    ex) qiime diversity core-metrics-phylogenetic
            --i-phylogeny phylogeny/rooted-tree.qza
            --i-table table.qza
            --p-sampling-depth 892
            --m-metadata-file mappingfile.txt
            --output-dir core-metrics-results

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime diversity core-metrics-phylogenetic [OPTIONS]
    #
    #   Applies a collection of diversity metrics (both phylogenetic and
    #   non-phylogenetic) to a feature table.
    #
    # Options:
    #   --i-table ARTIFACT PATH FeatureTable[Frequency]
    #                                   The feature table containing the
    #                                   samples over which diversity metrics
    #                                   should be computed.  [required]
    #   --i-phylogeny ARTIFACT PATH Phylogeny[Rooted]
    #                                   Phylogenetic tree containing tip
    #                                   identifiers that correspond to the
    #                                   feature identifiers in the table.
    #                                   This tree can contain tip ids that
    #                                   are not present in the table, but
    #                                   all feature ids in the table must be
    #                                   present in this tree.  [required]
    #   --p-sampling-depth INTEGER RANGE
    #                                   The total frequency that each sample
    #                                   should be rarefied to prior to
    #                                   computing diversity metrics.
    #                                   [required]
    #   --m-metadata-file MULTIPLE FILE
    #                                   Metadata file or artifact viewable
    #                                   as metadata. This option may be
    #                                   supplied multiple times to merge
    #                                   metadata. The sample metadata to use
    #                                   in the emperor plots.  [required]
    #   --p-n-jobs INTEGER RANGE        [beta/beta-phylogenetic methods
    #                                   only, excluding weighted_unifrac] -
    #                                   The number of jobs to use for the
    #                                   computation. This works by breaking
    #                                   down the pairwise matrix into n_jobs
    #                                   even slices and computing them in
    #                                   parallel. If -1 all CPUs are used.
    #                                   If 1 is given, no parallel computing
    #                                   code is used at all, which is useful
    #                                   for debugging. For n_jobs below -1,
    #                                   (n_cpus + 1 + n_jobs) are used. Thus
    #                                   for n_jobs = -2, all CPUs but one
    #                                   are used. (Description from
    #                                   sklearn.metrics.pairwise_distances)
    #                                   [default: 1]
    #   --o-rarefied-table ARTIFACT PATH FeatureTable[Frequency]
    #                                   The resulting rarefied feature
    #                                   table.  [required if not passing
    #                                   --output-dir]
    #   --o-faith-pd-vector ARTIFACT PATH SampleData[AlphaDiversity]
    #                                   Vector of Faith PD values by sample.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-observed-otus-vector ARTIFACT PATH SampleData[AlphaDiversity]
    #                                   Vector of Observed OTUs values by
    #                                   sample.  [required if not passing
    #                                   --output-dir]
    #   --o-shannon-vector ARTIFACT PATH SampleData[AlphaDiversity]
    #                                   Vector of Shannon diversity values
    #                                   by sample.  [required if not passing
    #                                   --output-dir]
    #   --o-evenness-vector ARTIFACT PATH SampleData[AlphaDiversity]
    #                                   Vector of Pielou's evenness values
    #                                   by sample.  [required if not passing
    #                                   --output-dir]
    #   --o-unweighted-unifrac-distance-matrix ARTIFACT PATH DistanceMatrix
    #                                   Matrix of unweighted UniFrac
    #                                   distances between pairs of samples.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-weighted-unifrac-distance-matrix ARTIFACT PATH DistanceMatrix
    #                                   Matrix of weighted UniFrac distances
    #                                   between pairs of samples.  [required
    #                                   if not passing --output-dir]
    #   --o-jaccard-distance-matrix ARTIFACT PATH DistanceMatrix
    #                                   Matrix of Jaccard distances between
    #                                   pairs of samples.  [required if not
    #                                   passing --output-dir]
    #   --o-bray-curtis-distance-matrix ARTIFACT PATH DistanceMatrix
    #                                   Matrix of Bray-Curtis distances
    #                                   between pairs of samples.  [required
    #                                   if not passing --output-dir]
    #   --o-unweighted-unifrac-pcoa-results ARTIFACT PATH PCoAResults
    #                                   PCoA matrix computed from unweighted
    #                                   UniFrac distances between samples.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-weighted-unifrac-pcoa-results ARTIFACT PATH PCoAResults
    #                                   PCoA matrix computed from weighted
    #                                   UniFrac distances between samples.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-jaccard-pcoa-results ARTIFACT PATH PCoAResults
    #                                   PCoA matrix computed from Jaccard
    #                                   distances between samples.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-bray-curtis-pcoa-results ARTIFACT PATH PCoAResults
    #                                   PCoA matrix computed from Bray-
    #                                   Curtis distances between samples.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-unweighted-unifrac-emperor VISUALIZATION PATH
    #                                   Emperor plot of the PCoA matrix
    #                                   computed from unweighted UniFrac.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-weighted-unifrac-emperor VISUALIZATION PATH
    #                                   Emperor plot of the PCoA matrix
    #                                   computed from weighted UniFrac.
    #                                   [required if not passing --output-
    #                                   dir]
    #   --o-jaccard-emperor VISUALIZATION PATH
    #                                   Emperor plot of the PCoA matrix
    #                                   computed from Jaccard.  [required if
    #                                   not passing --output-dir]
    #   --o-bray-curtis-emperor VISUALIZATION PATH
    #                                   Emperor plot of the PCoA matrix
    #                                   computed from Bray-Curtis.
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
    global QIIME2
    cmd = '{qiime2} diversity core-metrics-phylogenetic ' \
        '--i-table ' \
        '--i-phylogeny ' \
        '--p-sampling-depth ' \
        '--m-metadata-file ' \
        '--p-n-jobs -1' \
        '--output-dir '.format(
            qiime2=QIIME2
        )
    return cmd


def run_core_metrics_phylogenetic(kargs):
    echo('>>> Diveristy(core-metrics-phylogenetic) 시작')
    cmd = get_core_metrics_phylogenetic_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd({
        'run': run,
        'true_meg': 'QIIME2: core-metrics-phylogenetic 완료',
        'false_meg': 'QIIME2 - diversity core-metrics-phylogenetic',
    })
    pass


def get_emperor_pcoa_cmd(p_args):
    """

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime emperor plot [OPTIONS]
    #
    #   Generates an interactive ordination plot where the user can
    #   visually integrate sample metadata.
    #
    # Options:
    #   --i-pcoa ARTIFACT PATH PCoAResults
    #                                   The principal coordinates matrix to
    #                                   be plotted.  [required]
    #   --m-metadata-file MULTIPLE FILE
    #                                   Metadata file or artifact viewable
    #                                   as metadata. This option may be
    #                                   supplied multiple times to merge
    #                                   metadata. The sample metadata.
    #                                   [required]
    #   --p-custom-axes MULTIPLE TEXT   Numeric sample metadata columns that
    #                                   should be included as axes in the
    #                                   Emperor plot.  [optional]
    #   --p-ignore-missing-samples / --p-no-ignore-missing-samples
    #                                   This will suppress the error raised
    #                                   when the coordinates matrix contains
    #                                   samples that are not present in the
    #                                   metadata. Samples without metadata
    #                                   are included by setting all metadata
    #                                   values to: "This sample has no
    #                                   metadata". This flag is only applied
    #                                   if at least one sample is present in
    #                                   both the coordinates matrix and the
    #                                   metadata.  [default: False]
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
    global QIIME2
    cmd = '{qiime2} emperor plot ' \
        '--i-pcoa ' \
        '--m-metadata-file ' \
        '--output-dir '
    return cmd



def get_emperor_bioplot_cmd(p_args):
    """

    :param p_args:
    :return:
    """
    # """
    # Usage: qiime emperor biplot [OPTIONS]
    #
    #   Generates an interactive ordination biplot where the user can
    #   visually integrate sample and feature metadata.
    #
    # Options:
    #   --i-biplot ARTIFACT PATH PCoAResults % Properties(['biplot'])
    #                                   The principal coordinates matrix to
    #                                   be plotted.  [required]
    #   --m-sample-metadata-file MULTIPLE FILE
    #                                   Metadata file or artifact viewable
    #                                   as metadata. This option may be
    #                                   supplied multiple times to merge
    #                                   metadata. The sample metadata
    #                                   [required]
    #   --m-feature-metadata-file MULTIPLE FILE
    #                                   Metadata file or artifact viewable
    #                                   as metadata. This option may be
    #                                   supplied multiple times to merge
    #                                   metadata. The feature metadata
    #                                   (useful to manipulate the arrows in
    #                                   the plot).  [optional]
    #   --p-ignore-missing-samples / --p-no-ignore-missing-samples
    #                                   This will suppress the error raised
    #                                   when the coordinates matrix contains
    #                                   samples that are not present in the
    #                                   metadata. Samples without metadata
    #                                   are included by setting all metadata
    #                                   values to: "This sample has no
    #                                   metadata". This flag is only applied
    #                                   if at least one sample is present in
    #                                   both the coordinates matrix and the
    #                                   metadata.  [default: False]
    #   --p-number-of-features INTEGER RANGE
    #                                   The number of most important
    #                                   features (arrows) to display in the
    #                                   ordination.  [default: 5]
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
    global QIIME2
    cmd = '{qiime2} emperor biplot ' \
        '--i-biplot ' \
        '--m-sample-metadata-file ' \
        '--p-number-of-features ' \
        '--output-dir '.format(
            qiime2=QIIME2
        )
    return cmd

