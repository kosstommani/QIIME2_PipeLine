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
from SpoON.util import check_file_type, run_cmd, check_run_cmd, parse_config
from click import secho, echo, style
import pdb


CONFIG = parse_config()
QIIME2 = CONFIG['CoFI_QIIME2']

# QIIME2 Version 2018.08
# qiime --help
"""

Usage: qiime [OPTIONS] COMMAND [ARGS]...

  QIIME 2 command-line interface (q2cli)
  --------------------------------------

  To get help with QIIME 2, visit https://qiime2.org.

  To enable tab completion in Bash, run the following command or add it to
  your .bashrc/.bash_profile:

      source tab-qiime

  To enable tab completion in ZSH, run the following commands or add them to
  your .zshrc:

      autoload bashcompinit && bashcompinit && source tab-qiime

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  info                Display information about current deployment.
  tools               Tools for working with QIIME 2 files.
  dev                 Utilities for developers and advanced users.
  alignment           Plugin for generating and manipulating alignments.
  composition         Plugin for compositional data analysis.
  cutadapt            Plugin for removing adapter sequences, primers, and
                      other unwanted sequence from sequence data.
  dada2.R               Plugin for sequence quality control with DADA2.
  deblur              Plugin for sequence quality control with Deblur.
  demux               Plugin for demultiplexing & viewing sequence quality.
  diversity           Plugin for exploring community diversity.
  emperor             Plugin for ordination plotting with Emperor.
  feature-classifier  Plugin for taxonomic classification.
  feature-table       Plugin for working with sample by feature tables.
  gneiss              Plugin for building compositional models.
  longitudinal        Plugin for paired sample and time series analyses.
  metadata            Plugin for working with Metadata.
  phylogeny           Plugin for generating and manipulating phylogenies.
  quality-control     Plugin for quality control of feature and sequence data.
  quality-filter      Plugin for PHRED-based filtering and trimming.
  sample-classifier   Plugin for machine learning prediction of sample
                      metadata.
  taxa                Plugin for working with feature taxonomy annotations.
  vsearch             Plugin for clustering and dereplicating with vsearch.
"""

# QIIME tools import
"""
--show-importable-types
    DeblurStats
    DistanceMatrix
    EMPPairedEndSequences
    EMPSingleEndSequences
    FeatureData[AlignedSequence]
    FeatureData[Importance]
    FeatureData[PairedEndSequence]
    FeatureData[Sequence]
    FeatureData[Taxonomy]
    FeatureTable[Balance]
    FeatureTable[Composition]
    FeatureTable[Frequency]
    FeatureTable[PercentileNormalized]
    FeatureTable[PresenceAbsence]
    FeatureTable[RelativeFrequency]
    Hierarchy
    MultiplexedPairedEndBarcodeInSequence
    MultiplexedSingleEndBarcodeInSequence
    PCoAResults
    Phylogeny[Rooted]
    Phylogeny[Unrooted]
    QualityFilterStats
    RawSequences
    SampleData[AlphaDiversity]
    SampleData[BooleanSeries]
    SampleData[ClassifierPredictions]
    SampleData[DADA2Stats]
    SampleData[FirstDifferences]
    SampleData[JoinedSequencesWithQuality]
    SampleData[PairedEndSequencesWithQuality]
    SampleData[RegressorPredictions]
    SampleData[SequencesWithQuality]
    SampleData[Sequences]
    SampleEstimator[Classifier]
    SampleEstimator[Regressor]
    TaxonomicClassifier
    UchimeStats
    
--show-importable-formats    
    AlignedDNAFASTAFormat
    AlignedDNASequencesDirectoryFormat
    AlphaDiversityDirectoryFormat
    AlphaDiversityFormat
    BIOMV100DirFmt
    BIOMV100Format
    BIOMV210DirFmt
    BIOMV210Format
    BooleanSeriesDirectoryFormat
    BooleanSeriesFormat
    CasavaOneEightLanelessPerSampleDirFmt
    CasavaOneEightSingleLanePerSampleDirFmt
    DADA2StatsDirFmt
    DADA2StatsFormat
    DNAFASTAFormat
    DNASequencesDirectoryFormat
    DeblurStatsDirFmt
    DeblurStatsFmt
    DistanceMatrixDirectoryFormat
    EMPPairedEndCasavaDirFmt
    EMPPairedEndDirFmt
    EMPSingleEndCasavaDirFmt
    EMPSingleEndDirFmt
    FastqGzFormat
    FirstDifferencesDirectoryFormat
    FirstDifferencesFormat
    HeaderlessTSVTaxonomyDirectoryFormat
    HeaderlessTSVTaxonomyFormat
    ImportanceDirectoryFormat
    ImportanceFormat
    LSMatFormat
    MultiplexedPairedEndBarcodeInSequenceDirFmt
    MultiplexedSingleEndBarcodeInSequenceDirFmt
    NewickDirectoryFormat
    NewickFormat
    OrdinationDirectoryFormat
    OrdinationFormat
    PairedDNASequencesDirectoryFormat
    PairedEndFastqManifestPhred33
    PairedEndFastqManifestPhred64
    PredictionsDirectoryFormat
    PredictionsFormat
    QIIME1DemuxDirFmt
    QIIME1DemuxFormat
    QualityFilterStatsDirFmt
    QualityFilterStatsFmt
    SampleEstimatorDirFmt
    SingleEndFastqManifestPhred33
    SingleEndFastqManifestPhred64
    SingleLanePerSamplePairedEndFastqDirFmt
    SingleLanePerSampleSingleEndFastqDirFmt
    TSVTaxonomyDirectoryFormat
    TSVTaxonomyFormat
    TaxonomicClassiferTemporaryPickleDirFmt
    UchimeStatsDirFmt
    UchimeStatsFmt
"""

# qiime metadata --help
"""
Usage: qiime metadata [OPTIONS] COMMAND [ARGS]...

  Description: This QIIME 2 plugin provides functionality for working with
  and visualizing Metadata.

  Plugin website: https://github.com/qiime2/q2-metadata

  Getting user support: Please post to the QIIME 2 forum for help with this
  plugin: https://forum.qiime2.org

Options:
  --version    Show the version and exit.
  --citations  Show citations and exit.
  --help       Show this message and exit.

Commands:
  distance-matrix  Create a distance matrix from a numeric Metadata column
  tabulate         Interactively explore Metadata in an HTML table
"""

# qiime feature-table
"""
Usage: qiime feature-table [OPTIONS] COMMAND [ARGS]...

  Description: This is a QIIME 2 plugin supporting operations on
  sample by feature tables, such as filtering, merging, and
  transforming tables.

  Plugin website: https://github.com/qiime2/q2-feature-table

  Getting user support: Please post to the QIIME 2 forum for help
  with this plugin: https://forum.qiime2.org

Options:
  --version    Show the version and exit.
  --citations  Show citations and exit.
  --help       Show this message and exit.

Commands:
  core-features       Identify core features in table
  filter-features     Filter features from table
  filter-samples      Filter samples from table
  filter-seqs         Filter features from sequences
  group               Group samples or features by a metadata column
  heatmap             Generate a heatmap representation of a feature
                      table
  merge               Combine multiple tables
  merge-seqs          Combine collections of feature sequences
  merge-taxa          Combine collections of feature taxonomies
  presence-absence    Convert to presence/absence
  rarefy              Rarefy table
  relative-frequency  Convert to relative frequencies
  subsample           Subsample table
  summarize           Summarize table
  tabulate-seqs       View sequence associated with each feature
"""


def make_manifest_file(kargs):
    """
    FASTQ 파일들을 Artifact로 변환하기 위해서 필요한 manifest file을 생성한다.

    ex) manifest.csv
    sample-id,absolute-filepath,direction
    Mock.01,$PWD/RawData/Mock.01_1.fastq.gz,forward
    Mock.01,$PWD/RawData/Mock.01_2.fastq.gz,reverse
    Mock.02,$PWD/RawData/Mock.02_1.fastq.gz,forward
    Mock.02,$PWD/RawData/Mock.02_2.fastq.gz,reverse

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            run_list: list
            R1_suffix: str
            R2_suffix: str
            output_path: str

    :rtype: str
    :return: csv_file
    """
    l_header = ['sample-id', 'absolute-filepath', 'direction']
    csv_file = os.path.join(*(kargs['output_path'], 'manifest_file.csv'))
    with open(csv_file, 'w') as manifest_file:
        manifest_file.write(','.join(l_header))
        manifest_file.write('\n')
        for run in kargs['run_list']:
            for sample in run.samples:
                r1_text = [sample.name, os.path.join(*(sample.path, sample.name + kargs['R1_suffix'])), 'forward']
                r2_text = [sample.name, os.path.join(*(sample.path, sample.name + kargs['R2_suffix'])), 'reverse']
                manifest_file.write(','.join(r1_text))
                manifest_file.write('\n')
                manifest_file.write(','.join(r2_text))
                manifest_file.write('\n')
    echo('>>> manifest file 생성')
    return csv_file


def get_import_fastq_cmd(p_manifest, p_output):
    """
    FASTQ 파일들을 QIIME2 Artifact 로 변환하는 QIIME2 명령어를 생성한다.

    :type p_manifest: str
    :param p_manifest: manifest file

    :type p_output: str
    :param p_output: Artifact 생성 경로

    :return: cmd
    """
    # """
    # Usage: qiime tools import [OPTIONS]
    #
    #   Import data to create a new QIIME 2 Artifact. See https://docs.qiime2.org/
    #   for usage examples and details on the file types and associated semantic
    #   types that can be imported.
    #
    # Options:
    #   --type TEXT                The semantic type of the artifact that will be
    #                              created upon importing. Use --show-importable-
    #                              types to see what importable semantic types are
    #                              available in the current deployment.  [required]
    #   --input-path PATH          Path to file or directory that should be
    #                              imported.  [required]
    #   --output-path PATH         Path where output artifact should be written.
    #                              [required]
    #   --input-format TEXT        The format of the data to be imported. If not
    #                              provided, data must be in the format expected by
    #                              the semantic type provided via --type.
    #   --show-importable-types    Show the semantic types that can be supplied to
    #                              --type to import data into an artifact.
    #   --show-importable-formats  Show formats that can be supplied to --input-
    #                              format to import data into an artifact.
    #   --help                     Show this message and exit.
    # """
    global QIIME2
    cmd = '{qiime2} tools import ' \
        '--input-path {manifest} ' \
        '--type SampleData[PairedEndSequencesWithQuality] ' \
        '--input-format PairedEndFastqManifestPhred33 ' \
        '--output-path {output}'.format(
            qiime2=QIIME2,
            manifest=p_manifest,
            output=p_output,
        )
    return cmd


def run_import_fastq(kargs):
    """
    FASTQ 파일들을 QIIME2 Artifact 로 변환한다.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                run_list: list - named tuple인 RunData 를 원소로 가지는 리스트
                                 RunData.path
                                        .run_info
                                        .samples: named tuple인 SampleList을 원소로 가지는 리스트
                                            SampleList.path: 경로
                                                      .name: 시료명
                R1_suffix: str
                R2_suffix: str
                output_path: str
                order_number: str
    :return: artifact(qza) --> SampleData[PairedEndSequencesWithQuality]
    """
    echo('>>> Import FASTQ 시작')
    path = os.path.join(kargs['output_path'], 'FASTQ_artifact')
    if check_file_type(kargs['output_path'], 'exists'):
        os.mkdir(path)
    else:
        secho('Error: 디렉터리을 찾을 수 없습니다.', fg='red', blink=True)
        echo(kargs['output_path'])
        exit()

    manifest_file = make_manifest_file(
        {
            'run_list': kargs['run_list'],
            'R1_suffix': kargs['R1_suffix'],
            'R2_suffix': kargs['R2_suffix'],
            'output_path': path,
        }
    )
    fastq_artifact = os.path.join(path, kargs['order_number'] + '.qza')
    import_cmd = get_import_fastq_cmd(manifest_file, fastq_artifact)
    run = run_cmd(import_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: import FASTQ 완료',
            'false_meg': 'QIIME2 - import FASTQ'
        }
    )
    return fastq_artifact


def get_demux_summarize_cmd(p_artifact, p_output, p_random_read):
    """
    QIIME2의 demux -> summarize 의 실행 명령어를 반환한다.
    시료별 Read count, Quality Score Plot 을 생성한다.
    --p-n 값에 의해 Quality Score Plot 에 반영되는 Read 의 개수가 정해진다.

    :param p_artifact: str - qza
    :param p_output: str
    :param p_random_read: int
    :return: cmd
    """
    # """
    # Usage: qiime demux summarize [OPTIONS]
    #
    #   Summarize counts per sample for all samples, and generate interactive
    #   positional quality plots based on `n` randomly selected sequences.
    #
    # Options:
    #   --i-data ARTIFACT PATH SampleData
    #           [
    #                          JoinedSequencesWithQuality | PairedEndSequencesWithQuality | SequencesWithQuality
    #           ]
    #                                   The demultiplexed sequences to be
    #                                   summarized.  [required]
    #   --p-n INTEGER                   The number of sequences that should be
    #                                   selected at random for quality score plots.
    #                                   The quality plots will present the average
    #                                   positional qualities across all of the
    #                                   sequences selected. If input sequences are
    #                                   paired end, plots will be generated for both
    #                                   forward and reverse reads for the same `n`
    #                                   sequences.  [default: 10000]
    #   --o-visualization VISUALIZATION PATH
    #                                   [required if not passing --output-dir]
    #   --output-dir DIRECTORY          Output unspecified results to a directory
    #   --cmd-config PATH               Use config file for command options
    #   --verbose                       Display verbose output to stdout and/or
    #                                   stderr during execution of this action.
    #                                   [default: False]
    #   --quiet                         Silence output if execution is successful
    #                                   (silence is golden).  [default: False]
    #   --citations                     Show citations and exit.
    #   --help                          Show this message and exit.
    # """
    global QIIME2
    cmd = '{qiime2} demux summarize ' \
        '--i-data {qza} ' \
        '--o-visualization {qzv} ' \
        '--p-n {random_read}'.format(
            qiime2=QIIME2,
            qza=p_artifact,
            qzv=p_output,
            random_read=p_random_read,
        )
    return cmd


# matplotlib 관련 error : 설정 필요함.
def run_demux_summarize(p_artifact, p_output, p_random_read):
    """
    QIIME2의 demux -> summarize 의 실행한다.
    시료별 Read count, Quality Score Plot 을 생성한다.

    :type p_artifact: str
    :param p_artifact: qza

    :type p_output: str
    :param p_output:
    
    :type p_random_read: int
    :param p_random_read: Quality Score Plot에 사용할 Read의 개수.

    :return: None
    """
    echo('>>> FASTQ Summary 시작')
    cmd = get_demux_summarize_cmd(p_artifact, p_output, p_random_read)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: FASTQ Summary 완료',
            'false_meg': 'QIIME2 - demux summarize',
        }
    )


# TODO: DADA2 옵션 설정
# TODO: DADA2 매뉴얼 학습 필요
# DADA2 1.9 --> PacBio 데이터. 논문 revision
# --p-chimera-method[pooled|consensus|none] default: consensus --> 결과 차이?
# --p-n-reads-learn default: 1000000 --> 얼만큼?
#  run_dada_paired.R : /crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/envs/qiime2-2018.08/bin 에 있음.
def get_dada2_cmd(kargs):
    """
    DADA2의 실행명령어를 반환한다.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                qza: str
                output_path: str
                p_trunc_len_f: int
                p_trunc_len_r: int
                p_trim_left_f: int
                p_trim_left_r: int
    :return: cmd
    """
    # """
    # Usage: qiime dada2.R denoise-paired [OPTIONS]
    #
    #   This method denoises paired-end sequences, dereplicates them, and filters
    #   chimeras.
    #
    # Options:
    #   --i-demultiplexed-seqs ARTIFACT PATH SampleData[PairedEndSequencesWithQuality]
    #                                   The paired-end demultiplexed sequences to be
    #                                   denoised.  [required]
    #   --p-trunc-len-f INTEGER         Position at which forward read sequences
    #                                   should be truncated due to decrease in
    #                                   quality. This truncates the 3' end of the of
    #                                   the input sequences, which will be the bases
    #                                   that were sequenced in the last cycles.
    #                                   Reads that are shorter than this value will
    #                                   be discarded. After this parameter is
    #                                   applied there must still be at least a 20
    #                                   nucleotide overlap between the forward and
    #                                   reverse reads. If 0 is provided, no
    #                                   truncation or length filtering will be
    #                                   performed  [required]
    #   --p-trunc-len-r INTEGER         Position at which reverse read sequences
    #                                   should be truncated due to decrease in
    #                                   quality. This truncates the 3' end of the of
    #                                   the input sequences, which will be the bases
    #                                   that were sequenced in the last cycles.
    #                                   Reads that are shorter than this value will
    #                                   be discarded. After this parameter is
    #                                   applied there must still be at least a 20
    #                                   nucleotide overlap between the forward and
    #                                   reverse reads. If 0 is provided, no
    #                                   truncation or length filtering will be
    #                                   performed  [required]
    #   --p-trim-left-f INTEGER         Position at which forward read sequences
    #                                   should be trimmed due to low quality. This
    #                                   trims the 5' end of the input sequences,
    #                                   which will be the bases that were sequenced
    #                                   in the first cycles.  [default: 0]
    #   --p-trim-left-r INTEGER         Position at which reverse read sequences
    #                                   should be trimmed due to low quality. This
    #                                   trims the 5' end of the input sequences,
    #                                   which will be the bases that were sequenced
    #                                   in the first cycles.  [default: 0]
    #   --p-max-ee FLOAT                Reads with number of expected errors higher
    #                                   than this value will be discarded.
    #                                   [default: 2.0]
    #   --p-trunc-q INTEGER             Reads are truncated at the first instance of
    #                                   a quality score less than or equal to this
    #                                   value. If the resulting read is then shorter
    #                                   than `trunc_len_f` or `trunc_len_r`
    #                                   (depending on the direction of the read) it
    #                                   is discarded.  [default: 2]
    #   --p-chimera-method [pooled|consensus|none]
    #                                   The method used to remove chimeras. "none":
    #                                   No chimera removal is performed. "pooled":
    #                                   All reads are pooled prior to chimera
    #                                   detection. "consensus": Chimeras are
    #                                   detected in samples individually, and
    #                                   sequences found chimeric in a sufficient
    #                                   fraction of samples are removed.  [default:
    #                                   consensus]
    #   --p-min-fold-parent-over-abundance FLOAT
    #                                   The minimum abundance of potential parents
    #                                   of a sequence being tested as chimeric,
    #                                   expressed as a fold-change versus the
    #                                   abundance of the sequence being tested.
    #                                   Values should be greater than or equal to 1
    #                                   (i.e. parents should be more abundant than
    #                                   the sequence being tested). This parameter
    #                                   has no effect if chimera_method is "none".
    #                                   [default: 1.0]
    #   --p-n-threads INTEGER           The number of threads to use for
    #                                   multithreaded processing. If 0 is provided,
    #                                   all available cores will be used.  [default:
    #                                   1]
    #   --p-n-reads-learn INTEGER       The number of reads to use when training the
    #                                   error model. Smaller numbers will result in
    #                                   a shorter run time but a less reliable error
    #                                   model.  [default: 1000000]
    #   --p-hashed-feature-ids / --p-no-hashed-feature-ids
    #                                   If true, the feature ids in the resulting
    #                                   table will be presented as hashes of the
    #                                   sequences defining each feature. The hash
    #                                   will always be the same for the same
    #                                   sequence so this allows feature tables to be
    #                                   merged across runs of this method. You
    #                                   should only merge tables if the exact same
    #                                   parameters are used for each run.  [default:
    #                                   True]
    #   --o-table ARTIFACT PATH FeatureTable[Frequency]
    #                                   The resulting feature table.  [required if
    #                                   not passing --output-dir]
    #   --o-representative-sequences ARTIFACT PATH FeatureData[Sequence]
    #                                   The resulting feature sequences. Each
    #                                   feature in the feature table will be
    #                                   represented by exactly one sequence, and
    #                                   these sequences will be the joined paired-
    #                                   end sequences.  [required if not passing
    #                                   --output-dir]
    #   --o-denoising-stats ARTIFACT PATH SampleData[DADA2Stats]
    #                                   [required if not passing --output-dir]
    #   --output-dir DIRECTORY          Output unspecified results to a directory
    #   --cmd-config PATH               Use config file for command options
    #   --verbose                       Display verbose output to stdout and/or
    #                                   stderr during execution of this action.
    #                                   [default: False]
    #   --quiet                         Silence output if execution is successful
    #                                   (silence is golden).  [default: False]
    #   --citations                     Show citations and exit.
    #   --help                          Show this message and exit.
    #
    # """
    global QIIME2
    cmd = '{qiime2} dada2.R denoise-paired ' \
        '--i-demultiplexed-seqs {qza} ' \
        '--output-dir {output_dir} ' \
        '--p-trunc-len-f {trunc_len_f} ' \
        '--p-trunc-len-r {trunc_len_r} ' \
        '--p-trim-left-f {trim_left_f} ' \
        '--p-trim-left-r {trim_left_r} ' \
        '--p-chimera-method consensus ' \
        '--p-n-threads 0'.format(
            qiime2=QIIME2,
            qza=kargs['qza'],
            output_dir=os.path.join(*(kargs['output_path'], 'DADA2')),
            trunc_len_f=kargs['p_trunc_len_f'],
            trunc_len_r=kargs['p_trunc_len_r'],
            trim_left_f=kargs['p_trim_left_f'],
            trim_left_r=kargs['p_trim_left_r'],
        )
    return cmd


def run_dada2(kargs):
    """
    QIIME2의 DADA2를 실행한다.

    :param kargs:
            qza: str
            output_path: str
            p_trunc_len_f: int
            p_trunc_len_r: int
            p_trim_left_f: int
            p_trim_left_r: int

    :return:
    """
    echo('>>> DADA2 시작')
    cmd = get_dada2_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: DADA2 완료',
            'false_meg': 'QIIME2 - dada2.R denoise-paired'
        }
    )


def get_deblur_cmd():
    # TODO: deblur cmd
    cmd = '{qiime2} deblur'
    return cmd


def get_metadata_tabulate_cmd(kargs):
    """
    QIIME2 metadata tabulate 실행 명령어를 반환한다.
    DADA2의 STAT 정보를 시각화하는 명령어를 반환한다.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            qza: str - denoising_stats.qza 파일
            qzv: str - denoising_stats.qzv 파일
    :return: cmd
    """
    # """
    # Usage: qiime metadata tabulate [OPTIONS]
    #
    #   Generate a tabular view of Metadata. The output visualization supports
    #   interactive filtering, sorting, and exporting to common file formats.
    #
    # Options:
    #   --m-input-file MULTIPLE PATH    Metadata file or artifact viewable as
    #                                   metadata. This option may be supplied
    #                                   multiple times to merge metadata. The
    #                                   metadata to tabulate.  [required]
    #   --p-page-size INTEGER           The maximum number of Metadata records to
    #                                   display per page  [default: 100]
    #   --o-visualization VISUALIZATION PATH
    #                                   [required if not passing --output-dir]
    #   --output-dir DIRECTORY          Output unspecified results to a directory
    #   --cmd-config PATH               Use config file for command options
    #   --verbose                       Display verbose output to stdout and/or
    #                                   stderr during execution of this action.
    #                                   [default: False]
    #   --quiet                         Silence output if execution is successful
    #                                   (silence is golden).  [default: False]
    #   --citations                     Show citations and exit.
    #   --help                          Show this message and exit.
    # """
    global QIIME2
    cmd = '{qiime2} metadata tabulate ' \
        '--m-input-file {qza} ' \
        '--o-visualization {qzv} '.format(
            qiime2=QIIME2,
            qza=kargs['qza'],
            qzv=kargs['qzv'],
        )
    return cmd


def run_dada2_stats_qzv(kargs):
    """
    DADA2의 STAT 정보를 시각화한다.
    qzv 파일 생성.

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            qza: str - denoising_stats.qza 파일
            qzv: str - denoising_stats.qzv 파일
    :return:
    """
    echo('>>> DADA2 STATs : Tabulate 시작')
    cmd = get_metadata_tabulate_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: DADA2 Denoising STATs - Visualization 완료',
            'false_meg': 'QIIME2 - metadata tabulate'
        }
    )


def get_feature_table_summarize_cmd(kargs):
    """
    DADA2 실행 후 생성된 feature table(table.qza)를 시각화하는 명령어를 반환한다.
    table.qza --> table.qzv

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                qza: str - FeatureTable[Frequency | PresenceAbsence | RelativeFrequency]
                qzv: str -
                metadata: str -
    :return: cmd
    """
    #  """
    #  Usage: qiime feature-table summarize [OPTIONS]
    #
    #   Generate visual and tabular summaries of a feature table.
    #
    # Options:
    #   --i-table ARTIFACT PATH FeatureTable[Frequency | PresenceAbsence | RelativeFrequency]
    #                                   The feature table to be summarized.
    #                                   [required]
    #   --m-sample-metadata-file MULTIPLE PATH
    #                                   Metadata file or artifact viewable
    #                                   as metadata. This option may be
    #                                   supplied multiple times to merge
    #                                   metadata. The sample metadata.
    #                                   [optional]
    #   --o-visualization VISUALIZATION PATH
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
    #  """
    global QIIME2
    cmd = '{qiime2} feature-table summarize ' \
        '--i-table {qza} ' \
        '--o-visualization {qzv} ' \
        '--m-sample-metadata-file {metadata}'.format(
            qiime2=QIIME2,
            qza=kargs['qza'],
            qzv=kargs['qzv'],
            metadata=kargs['metadata']
        )
    return cmd


def get_feature_table_tabulate_seqs_cmd(kargs):
    """
    DADA2 실행 후 생성된 feature sequences(representative_sequences.qza)를 시각화하는 명령어를 반환한다.
    representative_sequences.qza --> representative_sequences.qzv

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                qza: str - FeatureData[Sequence]
                qzv: str -
    :return:
    """
    # """
    # Usage: qiime feature-table tabulate-seqs [OPTIONS]
    #
    #   Generate tabular view of feature identifier to sequence mapping,
    #   including links to BLAST each sequence against the NCBI nt
    #   database.
    #
    # Options:
    #   --i-data ARTIFACT PATH FeatureData[Sequence]
    #                                   The feature sequences to be
    #                                   tabulated.  [required]
    #   --o-visualization VISUALIZATION PATH
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
    global QIIME2
    cmd = '{qiime2} feature-table tabulate-seqs ' \
        '--i-data {qza} ' \
        '--o-visualization {qzv} '.format(
            qiime2=QIIME2,
            qza=kargs['qza'],
            qzv=kargs['qzv'],
        )
    return cmd


def run_feature_table_summarize(kargs):
    """
    DADA2 실행 후 생성된 feature table(table.qza)를 시각화한다.
    table.qza --> table.qzv

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                qza: str -
                qzv: str -
                metadata: str -
    :return:
    """
    echo('>>> Feature Table Summarize 시작')
    cmd = get_feature_table_summarize_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: feature-table summarize 완료',
            'false_meg': 'QIIME2 - feature-table summarize'
        }
    )


def run_feature_table_tabulate_seqs(kargs):
    """
    DADA2 실행 후 생성된 feature sequences(representative_sequences.qza)를 시각화한다.
    representative_sequences.qza --> representative_sequences.qzv

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
                qza: str -
                qzv: str -
    :return:
    """
    echo('>>> Feature Table Tabulate Seqs 시작')
    cmd = get_feature_table_tabulate_seqs_cmd(kargs)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'QIIME2: feature-table tabulate-seqs 완료',
            'false_meg': 'QIIME2 - feature-table tabulate-seqs'
        }
    )
