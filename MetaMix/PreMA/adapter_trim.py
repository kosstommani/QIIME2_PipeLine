# ----------------------------------------------------------------------------------------------------------------------
# 888b     d888          888             888b     d888 d8b
# 8888b   d8888          888             8888b   d8888 Y8P
# 88888b.d88888          888             88888b.d88888
# 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#  _____         _____ _____
# |  _  |___ ___|     |  _  |
# |   __|  _| -_| | | |     |
# |__|  |_| |___|_|_|_|__|__|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------
__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.1'

from multiprocessing import Pool
from click import secho, echo, style, progressbar
import os
from SpoON.util import typeWrite, run_cmd, check_run_cmd, parse_config, start_process
# from tqdm import tqdm
# import progressbar
import pdb


CONFIG = parse_config()
FASTP = CONFIG['PreMA_Adapter_Trim']['fastp']
SEQPURGE = CONFIG['PreMA_Adapter_Trim']['SeqPurge']
SCYTHE = CONFIG['PreMA_Adapter_Trim']['Scythe']
POOL_WORKER = CONFIG['PreMA_PoolWorker']
FASTP_THREAD = CONFIG['PreMA_Adapter_Trim']['fastp_thread']


# FASTP
# """
#     options:
#   -i, --in1                            read1 input file name (string)
#   -o, --out1                           read1 output file name (string [=])
#   -I, --in2                            read2 input file name (string [=])
#   -O, --out2                           read2 output file name (string [=])
#   -6, --phred64                        indicate the input is using phred64 scoring (it'll be converted to phred33,
#                                        so the output will still be phred33)
#   -z, --compression                    compression level for gzip output (1 ~ 9).
#                                        1 is fastest, 9 is smallest, default is 4. (int [=4])
#       --stdout                         stream passing-filters reads to STDOUT. This option will result in interleaved
#                                        FASTQ output for paired-end input.Disabled by defaut.
#       --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and
#                                        read2. Disabled by defaut.
#       --reads_to_process               specify how many reads/pairs to be processed. Default 0 means process
#                                        all reads. (int [=0])
#       --dont_overwrite                 don't overwrite existing files. Overwritting is allowed by default.
#   -V, --verbose                        output verbose log information (i.e. when every 1M reads are processed).
#   -A, --disable_adapter_trimming       adapter trimming is enabled by default. If this option is specified,
#                                        adapter trimming is disabled
#   -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be
#                                        auto-detected. For PE data, this is used if R1/R2 are found not overlapped.
#                                        (string [=auto])
#       --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found
#                                        not overlapped. If not specified, it will be the same as <adapter_sequence>
#                                        (string [=])
#   -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
#   -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
#   -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified,
#                                        it will follow read1's settings (int [=0])
#   -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified,
#                                        it will follow read1's settings (int [=0])
#   -g, --trim_poly_g                    force polyG tail trimming, by default trimming is automatically enabled for
#                                        Illumina NextSeq/NovaSeq data
#       --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
#   -G, --disable_trim_poly_g            disable polyG tail trimming, by default trimming is automatically enabled
#                                        for Illumina NextSeq/NovaSeq data
#   -x, --trim_poly_x                    enable polyX trimming in 3' ends.
#       --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
#   -5, --cut_by_quality5                enable per read cutting by quality in front (5'), default is disabled
#                                        (WARNING: this will interfere deduplication for both PE/SE data)
#   -3, --cut_by_quality3                enable per read cutting by quality in tail (3'), default is disabled
#                                        (WARNING: this will interfere deduplication for SE data)
#   -W, --cut_window_size                the size of the sliding window for sliding window trimming,
#                                        default is 4 (int [=4])
#   -M, --cut_mean_quality               the bases in the sliding window with mean quality below cutting_quality
#                                        will be cut, default is Q20 (int [=20])
#   -Q, --disable_quality_filtering      quality filtering is enabled by default. If this option is specified,
#                                        quality filtering is disabled
#   -q, --qualified_quality_phred        the quality value that a base is qualified.
#                                        Default 15 means phred quality >=Q15 is qualified. (int [=15])
#   -u, --unqualified_percent_limit      how many percents of bases are allowed to be unqualified (0~100).
#                                        Default 40 means 40% (int [=40])
#   -n, --n_base_limit                   if one read's number of N base is >n_base_limit,
#                                        then this read/pair is discarded. Default is 5 (int [=5])
#   -L, --disable_length_filtering       length filtering is enabled by default. If this option is specified,
#                                        length filtering is disabled
#   -l, --length_required                reads shorter than length_required will be discarded, default is 15.
#                                        (int [=15])
#       --length_limit                   reads longer than length_limit will be discarded,
#                                        default 0 means no limitation. (int [=0])
#   -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of
#                                        base that is different from its nextbase (base[i] != base[i+1]).
#   -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30,
#                                        which means 30% complexity is required. (int [=30])
#       --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out,
#                                        one barcode per line (string [=])
#       --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out,
#                                        one barcode per line (string [=])
#       --filter_by_index_threshold      the allowed difference of index barcode for index filtering,
#                                        default 0 means completely identical. (int [=0])
#   -c, --correction                     enable base correction in overlapped regions (only for PE data),
#                                        default is disabled
#       --overlap_len_require            the minimum length of the overlapped region for overlap analysis based
#                                        adapter trimming and correction. 30 by default. (int [=30])
#       --overlap_diff_limit             the maximum difference of the overlapped region for overlap analysis based
#                                        adapter trimming and correction. 5 by default. (int [=5])
#   -U, --umi                            enable unique molecular identifer (UMI) preprocessing
#       --umi_loc                        specify the location of UMI, can be
#                                        (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
#       --umi_len                        if the UMI is in read1/read2, its length should be provided (int [=0])
#       --umi_prefix                     if specified, an underline will be used to connect prefix and UMI
#                                        (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG).
#                                        No prefix by default (string [=])
#       --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI,
#                                        default is 0 (int [=0])
#   -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
#   -P, --overrepresentation_sampling    one in (--overrepresentation_sampling) reads will be computed for
#                                        overrepresentation analysis (1~10000), smaller is slower,
#                                        default is 20. (int [=20])
#   -j, --json                           the json format report file name (string [=fastp.json])
#   -h, --html                           the html format report file name (string [=fastp.html])
#   -R, --report_title                   should be quoted with ' or ",
#                                        default is "fastp report" (string [=fastp report])
#   -w, --thread                         worker thread number, default is 2 (int [=2])
#   -s, --split                          split output by limiting total split file number with this option (2~999),
#                                        a sequential number prefix will be added to output name
#                                        ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
#   -S, --split_by_lines                 split output by limiting lines of each file with this option(>=1000),
#                                        a sequential number prefix will be added tooutput name
#                                        ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
#   -d, --split_prefix_digits            the digits for the sequential number padding (1~10), default is 4,
#                                        so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
#   -?, --help                           print this message
# """
def get_fastp_cmd(kargs):
    """
    fastp 실행 명령어를 반환한다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary.
                - in1(str)
                - in2(str)
                - out1(str)
                - out2(str)
                - adap_seq1(str)
                - adap_seq2(str)
                - trim_tail(int)
                - html(str)
                - json(str)
                - title(str)
                - log(str)

    :return: fastp_cmd
    """
    global FASTP, FASTP_THREAD
    if kargs['trim_tail'] is None:
        trim_tail_text = ''
    else:
        trim_tail_text = '--trim_tail1 {base} --trim_tail2 {base} '.format(base=kargs['trim_tail'])
    fastp_cmd = \
        '{fastp} ' \
        '-i {input1} ' \
        '-I {input2} ' \
        '-o {output1} ' \
        '-O {output2} ' \
        '--adapter_sequence {adapter_seq1} ' \
        '--adapter_sequence_r2 {adapter_seq2} ' \
        '{trim_tail}' \
        '--html {html} ' \
        '--json {json} ' \
        '--report_title {title} ' \
        '--thread {thread} ' \
        '--disable_trim_poly_g ' \
        '--disable_length_filtering ' \
        '--correction |& tee {log}'.format(
            fastp=FASTP,
            input1=kargs['in1'],
            input2=kargs['in2'],
            output1=kargs['out1'],
            output2=kargs['out2'],
            adapter_seq1=kargs['adap_seq1'],
            adapter_seq2=kargs['adap_seq2'],
            trim_tail=trim_tail_text,
            html=kargs['html'],
            json=kargs['json'],
            title=kargs['title'],
            thread=FASTP_THREAD,
            log=kargs['log'],
        )
    return fastp_cmd


def get_fastp_my_cmd(kargs):
    """
    개인 유전체 - 장내모니터링 서비스인 my Biomestory 에 적용되는 fastp 실행 명령어를 반환한다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary.
                - in1(str)
                - in2(str)
                - out1(str)
                - out2(str)
                - adap_seq1(str)
                - adap_seq2(str)
                - html(str)
                - json(str)
                - title(str)
                - log(str)

    :return: fastp_cmd
    """
    global FASTP
    fastp_cmd = \
        '{fastp} ' \
        '-i {input1} ' \
        '-I {input2} ' \
        '-o {output1} ' \
        '-O {output2} ' \
        '--adapter_sequence {adapter_seq1} ' \
        '--adapter_sequence_r2 {adapter_seq2} ' \
        '--html {html} ' \
        '--json {json} ' \
        '--report_title {title} ' \
        '--thread 12 ' \
        '--trim_tail1 50 ' \
        '--trim_tail2 50 ' \
        '--trim_front1 0 ' \
        '--trim_front2 0 ' \
        '--disable_trim_poly_g ' \
        '--length_required 200 ' \
        '--correction |& tee {log}'.format(
            fastp=FASTP,
            input1=kargs['in1'],
            input2=kargs['in2'],
            output1=kargs['out1'],
            output2=kargs['out2'],
            adapter_seq1=kargs['adap_seq1'],
            adapter_seq2=kargs['adap_seq2'],
            html=kargs['html'],
            json=kargs['json'],
            title=kargs['title'],
            log=kargs['log'],
        )
    return fastp_cmd


def get_fastp_mam_cmd(kargs):
    """
    개인 유전체 - 장내모니터링 서비스인 MicrobeAndMe 에 적용되는 fastp 실행 명령어를 반환한다.
    
    적용 옵션
    '--trim_tail1 50 ' 
    '--trim_tail2 70 ' 
    '--trim_front1 17 '
    '--trim_front2 21 '

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary.
                - in1(str)
                - in2(str)
                - out1(str)
                - out2(str)
                - adap_seq1(str)
                - adap_seq2(str)
                - html(str)
                - json(str)
                - title(str)
                - log(str)

    :return: fastp_cmd
    """
    global FASTP
    fastp_cmd = \
        '{fastp} ' \
        '-i {input1} ' \
        '-I {input2} ' \
        '-o {output1} ' \
        '-O {output2} ' \
        '--adapter_sequence {adapter_seq1} ' \
        '--adapter_sequence_r2 {adapter_seq2} ' \
        '--html {html} ' \
        '--json {json} ' \
        '--report_title {title} ' \
        '--thread 12 ' \
        '--trim_tail1 0 ' \
        '--trim_tail2 0 ' \
        '--trim_front1 17 ' \
        '--trim_front2 21 ' \
        '--disable_trim_poly_g ' \
        '--correction |& tee {log}'.format(
            fastp=FASTP,
            input1=kargs['in1'],
            input2=kargs['in2'],
            output1=kargs['out1'],
            output2=kargs['out2'],
            adapter_seq1=kargs['adap_seq1'],
            adapter_seq2=kargs['adap_seq2'],
            html=kargs['html'],
            json=kargs['json'],
            title=kargs['title'],
            log=kargs['log'],
        )
    return fastp_cmd


# SeqPurge
# """
# SeqPurge (2018_06)
#
# Removes adapter sequences from paired-end sequencing data.
#
# Mandatory parameters:
#   -in1 <filelist>     Forward input gzipped FASTQ file(s).
#   -in2 <filelist>     Reverse input gzipped FASTQ file(s).
#   -out1 <file>        Forward output gzipped FASTQ file.
#   -out2 <file>        Reverse output gzipped FASTQ file.
#
# Optional parameters:
#   -a1 <string>        Forward adapter sequence (at least 15 bases).
#                       Default value: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
#   -a2 <string>        Reverse adapter sequence (at least 15 bases).
#                       Default value: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
#   -match_perc <float> Minimum percentage of matching bases for sequence/adapter matches.
#                       Default value: '80'
#   -mep <float>        Maximum error probability of insert and adapter matches.
#                       Default value: '9.9999999999999995e-07'
#   -qcut <int>         Quality trimming cutoff for trimming from the end of reads using a sliding window approach.
#                       Set to 0 to disable.
#                       Default value: '15'
#   -qwin <int>         Quality trimming window size.
#                       Default value: '5'
#   -qoff <int>         Quality trimming FASTQ score offset.
#                       Default value: '33'
#   -ncut <int>         Number of subsequent Ns to trimmed using a sliding window approach from the front of reads.
#                       Set to 0 to disable.
#                       Default value: '7'
#   -min_len <int>      Minimum read length after adapter trimming. Shorter reads are discarded.
#                       Default value: '30'
#   -threads <int>      The number of threads used for trimming (an additional thread is used for reading data).
#                       Default value: '1'
#   -out3 <file>        Name prefix of singleton read output files (if only one read of a pair is discarded).
#                       Default value: ''
#   -summary <file>     Write summary/progress to this file instead of STDOUT.
#                       Default value: ''
#   -qc <file>          If set, a read QC file in qcML format is created (just like ReadQC).
#                       Default value: ''
#   -prefetch <int>     Maximum number of reads that may be pre-fetched to speed up trimming
#                       Default value: '1000'
#   -ec                 Enable error-correction of adapter-trimmed reads (only those with insert match).
#                       Default value: 'false'
#   -debug              Enables debug output (use only with one thread).
#                       Default value: 'false'
#   -progress           Enables progress output.
#                       Default value: 'false'
#
# Special parameters:
#   --help              Shows this help and exits.
#   --version           Prints version and exits.
#   --changelog         Prints changeloge and exits.
#   --tdx               Writes a Tool Definition Xml file. The file name is the application name with the suffix '.tdx'.
#
# """
def get_seqpurge_cmd(kargs):
    """
    SeqPurge 실행 명령어를 반환한다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                - in1(str)
                - in2(str)
                - out1(str)
                - out2(str)
                - adap_seq1(str)
                - adap_seq2(str)
                - summary(str)

    :return: seqpurge_cmd
    """
    global SEQPURGE
    seqpurge_cmd = \
        '{seqpurge} ' \
        '-in1 {input1} ' \
        '-in2 {input2} ' \
        '-out1 {output1} ' \
        '-out2 {output2} ' \
        '-a1 {adapter_seq1} ' \
        '-a2 {adapter_seq2} ' \
        '-qcut 0 ' \
        '-ncut 0 ' \
        '-threads 20 ' \
        '-summary {summary_file} ' \
        '-ec'.format(
            seqpurge=SEQPURGE,
            input1=kargs['in1'],
            input2=kargs['in2'],
            output1=kargs['out1'],
            output2=kargs['out2'],
            adapter_seq1=kargs['adap_seq1'],
            adapter_seq2=kargs['adap_seq2'],
            summary_file=kargs['summary']
        )
    return seqpurge_cmd


# scythe
# """
# Usage: scythe -a adapter_file.fasta sequence_file.fastq
# Trim 3'-end adapter contaminants off sequence files. If no output file
# is specified, scythe will use stdout.
#
# Options:
#   -p, --prior		    prior (default: 0.300)
#   -q, --quality-type	quality type, either illumina, solexa, or sanger (default: sanger)
#   -m, --matches-file	matches file (default: no output)
#   -o, --output-file 	output trimmed sequences file (default: stdout)
#   -t, --tag	        	add a tag to the header indicating Scythe cut a sequence (default: off)
#   -n, --min-match   	smallest contaminant to consider (default: 5)
#   -M, --min-keep    	filter sequnces less than or equal to this length (default: 35)
#   --quiet	           	don't output statistics about trimming to stdout (default: off)
#   --help		        display this help and exit
#   --version		        output version information and exit
#
#   Information on quality schemes:
#   phred			PHRED quality scores (e.g. from Roche 454). ASCII with no offset, range: [4, 60].
#   sanger		Sanger are PHRED ASCII qualities with an offset of 33, range: [0, 93]. From
#                 NCBI SRA, or Illumina pipeline 1.8+.
#   solexa		Solexa (also very early Illumina - pipeline < 1.3). ASCII offset of
#                 64, range: [-5, 62]. Uses a different quality-to-probabilities conversion than other
#                 schemes.
#   illumina		Illumina output from pipeline versions between 1.3 and 1.7. ASCII offset of 64,
#                 range: [0, 62]
# """
def get_scythe_cmd(kargs):
    """
    Scythe 실행 명령어를 반환한다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                - adap_seq(str)
                - out(str)
                - in(str)
                - log(str)

    :return: scythe_cmd
    """
    global SCYTHE
    scythe_cmd = \
        '{scythe} ' \
        '-a {adapter_seq} ' \
        '-o {output} {input} -p 0.01 |& tee {log}'.format(
            scythe=SCYTHE,
            adapter_seq=kargs['adap_seq'],
            output=kargs['out'],
            input=kargs['in'],
            log=kargs['log'],
        )
    return scythe_cmd


def get_adapter_seq(p_method, p_kit):
    """
    Adapter Trimming Tool과 index kit 정보에 따라 적절한 Adataper Sequences을 반환한다.

    :type p_method: str
    :param p_method: Adapter Trimming Tool('fastp', 'SeqPurge', 'Scythe')

    :type p_kit: str
    :param p_kit: index kit 정보('Nextera', 'TruSeq')

    :return: seq1,seq2
    """
    nextera_seq1 = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
    nextera_seq2 = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
    truseq_seq1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    truseq_seq2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
    scythe_nextera_seq = '/garnet/Tools/Amplicon_MetaGenome/Scythe/nextera_enrichment.fa'
    # scythe_nextera_seq = '/lustre/Tools/macrogen_analysis_toolbox/scythe/vsbuffalo-scythe-20d3cff/nextera_enrichment.fa'
    scythe_truseq_seq = '/garnet/Tools/Amplicon_MetaGenome/Scythe/truseq.fa'
    # scythe_truseq_seq = '/lustre/Tools/macrogen_analysis_toolbox/scythe/vsbuffalo-scythe-20d3cff/Scythe/truseq.fa'

    if p_method == 'Scythe':
        if p_kit.lower() == 'nextera':
            return scythe_nextera_seq
        elif p_kit.lower() == 'truseq':
            return scythe_truseq_seq
    else:
        if p_kit.lower() == 'nextera':
            return nextera_seq1, nextera_seq2
        elif p_kit.lower() == 'truseq':
            return truseq_seq1, truseq_seq2


def check_fastp_run_error(p_run):
    if 'Filtering result' in p_run.stdout:
        return True
    else:
        return False


def check_scythe_run_error(p_run):
    if 'Adapter Trimming Complete' in p_run.stdout:
        return True
    else:
        return False


def run_adapter_trim(kargs, p_tool):
    """
    fastq파일에 대해서 Adapter Trimming 을 진행한다.
    index kit에 따라 Adapter Sequences 가 달리 적용된다.
    Tool: fastp, SeqPurge, Scythe
    index_kit: Nextera, TruSeq


    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                sample_list(list): namedtuple인 SampleList 를 원소로 가지는 리스트.
                             SampleList.rawdata_path_in_analysis
                                       .cur_name
                                       .new_name
                R1_suffix(str): Read1에 대한 fastq 파일명의 접미사.
                R2_suffix(str): Read2에 대한 fastq 파일명의 접미사.
                index_kit(str): index kit 정보
                my_story(bool): 개인 유전체 - 장내모니터링 서비스에 대한 FASTQ QC 기준 적용(fastp).
                trim_tail(int): fastp 전용. Read의 끝에서부터 제거할 base의 개수.

    :type p_tool: str
    :param p_tool: Adapter Trim Tool

    :rtype: list
    :return: l_summary_files.
             Adapter Trimming Tool 에서 생성하는 Summary 또는 Log 파일의 경로를 포함한 파일명을 원소로 가지는 리스트.
             Scythe의 경우 tuple을 원소로 가지는 리스트. tuple은 2개의 log파일명을 원소로 가진다. [tuple(log1,log2), ...]
    """
    global POOL_WORKER
    if p_tool == 'fastp':
        adapter_seqs = get_adapter_seq(p_tool, kargs['index_kit'])
    elif p_tool == 'SeqPurge':
        adapter_seqs = get_adapter_seq(p_tool, kargs['index_kit'])
    elif p_tool == 'Scythe':
        adapter_seqs = get_adapter_seq(p_tool, kargs['index_kit'])
    else:
        raise ValueError(style(p_tool, fg='red', blink=True))

    l_cmd = list()
    l_summary_files = list()
    for nt_sample in kargs['sample_list']:
        src_path = os.path.join(*(nt_sample.path, nt_sample.new_name))
        r1_file = '{path}/{name}{suffix}'.format(
            path=src_path,
            name=nt_sample.new_name,
            suffix=kargs['R1_suffix'])
        r2_file = '{path}/{name}{suffix}'.format(
            path=src_path,
            name=nt_sample.new_name,
            suffix=kargs['R2_suffix'])
        out1_file = '{path}/{name}{suffix}'.format(path=src_path,
                                                   name=nt_sample.new_name,
                                                   suffix='_1.{method}.fastq'.format(method=p_tool))
        out2_file = '{path}/{name}{suffix}'.format(path=src_path,
                                                   name=nt_sample.new_name,
                                                   suffix='_2.{method}.fastq'.format(method=p_tool))

        if p_tool == 'fastp':
            html = '{path}/{name}{suffix}'.format(path=src_path,
                                                  name=nt_sample.new_name,
                                                  suffix='.fastp.html')
            json = '{path}/{name}{suffix}'.format(path=src_path,
                                                  name=nt_sample.new_name,
                                                  suffix='.fastp.json')
            log = '{path}/{name}{suffix}'.format(path=src_path,
                                                 name=nt_sample.new_name,
                                                 suffix='.fastp.log')
            title = '"fastp report : {name}"'.format(name=nt_sample.new_name)
            if all([kargs['my_story'], kargs['microbe_and_me']]) is True:
                secho('Option Error: --my_stroy 와 --microbe_and_me 옵션 모두 적용 불가', fg='red')
                exit()
            elif any([kargs['my_story'], kargs['microbe_and_me']]) is False:
                cmd = get_fastp_cmd(
                    {
                        'in1': r1_file,
                        'in2': r2_file,
                        'out1': out1_file+'.gz',
                        'out2': out2_file+'.gz',
                        'adap_seq1': adapter_seqs[0],
                        'adap_seq2': adapter_seqs[1],
                        'trim_tail': kargs['trim_tail'],
                        'html': html,
                        'json': json,
                        'title': title,
                        'log': log,
                    }
                )
            elif kargs['my_story'] is True:
                cmd = get_fastp_my_cmd(
                    {
                        'in1': r1_file,
                        'in2': r2_file,
                        'out1': out1_file + '.gz',
                        'out2': out2_file + '.gz',
                        'adap_seq1': adapter_seqs[0],
                        'adap_seq2': adapter_seqs[1],
                        'html': html,
                        'json': json,
                        'title': title,
                        'log': log,
                    }
                )
            elif kargs['microbe_and_me'] is True:
                cmd = get_fastp_mam_cmd(
                    {
                        'in1': r1_file,
                        'in2': r2_file,
                        'out1': out1_file + '.gz',
                        'out2': out2_file + '.gz',
                        'adap_seq1': adapter_seqs[0],
                        'adap_seq2': adapter_seqs[1],
                        'html': html,
                        'json': json,
                        'title': title,
                        'log': log,
                    }
                )
            else:
                raise ValueError('설정값이 잘못되었습니다.\n'
                                 f'--my_story : {kargs["my_story"]}\n'
                                 f'--microbe_and_me : {kargs["microbe_and_me"]}\n')
            l_cmd.append(cmd)
        elif p_tool == 'SeqPurge':
            summary = '{path}/{name}{suffix}'.format(path=src_path,
                                                     name=nt_sample.new_name,
                                                     suffix='.SeqPurge.summary')
            cmd = get_seqpurge_cmd(
                {
                    'in1': r1_file,
                    'in2': r2_file,
                    'out1': out1_file+'.gz',
                    'out2': out2_file+'.gz',
                    'adap_seq1': adapter_seqs[0],
                    'adap_seq2': adapter_seqs[1],
                    'summary': summary,
                }
            )
            l_cmd.append(cmd)
            l_summary_files.append(summary)
        elif p_tool == 'Scythe':
            log1 = '{path}/{name}{suffix}'.format(path=src_path,
                                                  name=nt_sample.new_name,
                                                  suffix='_1.scythe.log'
                                                  )
            log2 = '{path}/{name}{suffix}'.format(path=src_path,
                                                  name=nt_sample.new_name,
                                                  suffix='_2.scythe.log'
                                                  )
            cmd1 = get_scythe_cmd(
                {
                    'in': r1_file,
                    'out': out1_file,
                    'adap_seq': adapter_seqs,
                    'log': log1
                }
            )
            cmd2 = get_scythe_cmd(
                {
                    'in': r2_file,
                    'out': out2_file,
                    'adap_seq': adapter_seqs,
                    'log': log2
                }
            )
            l_cmd.append(cmd1)
            l_cmd.append(cmd2)
            l_summary_files.append((log1, log2))
        else:
            raise ValueError(style(p_tool, fg='red', blink=True))

    # Adapter Trimming 명령어 저장
    o_trim_cmd = open(os.path.join(*(nt_sample.path, 'ADAPTER_TRIM.recipe')), 'w')
    _ = [echo(x, file=o_trim_cmd) for x in l_cmd]

    # TODO multiple progress bar 구현(어려움)
    # TODO multiprocessing or concurrent.futures 학습 필요
    # progressbar2, tqdm 라이브러리

    process = POOL_WORKER
    if len(l_cmd) < POOL_WORKER:
        process = len(l_cmd)
    secho('>>> Adapter Trimming 시작 : {tool}, {kit}'.format(tool=p_tool, kit=kargs['index_kit']), fg='cyan')
    with Pool(process, initializer=start_process) as pool:
        pool_outputs = pool.map(run_cmd, l_cmd)
        pool.close()
        pool.join()

    done_count = 0
    error_count = 0
    for run in pool_outputs:
        if p_tool == 'fastp':
            status = check_fastp_run_error(run)
            if status:
                done_count += 1
            else:
                error_count += 1
        elif p_tool == 'Scythe':
            status = check_scythe_run_error(run)
            if status:
                done_count += 1
            else:
                error_count += 1
        else:
            if (run.returncode == 0) and (run.stderr.decode('utf-8').strip() is ''):
                done_count += 1
            else:
                error_count += 1

    done_text = style('완료: {}'.format(done_count), fg='cyan')
    error_text = style('에러: {}'.format(error_count),
                       fg='red' if (error_count != 0) else 'cyan',
                       blink=True if (error_count != 0) else False)
    echo('{msg} {done}, {error}'.format(
        msg=style('>>> Adapter Trimming', fg='cyan'),
        done=done_text,
        error=error_text,
        )
    )
    if error_count != 0:
        # TODO error cmd 처리를 어떻게 할 것인지. cmd 를 파일로 저장
        for run in pool_outputs:
            if (p_tool == 'fastp') or (p_tool == 'Scythe'):
                run.returncode = 1
                run.stderr = run.stdout
            check_run_cmd({
                'run': run,
                'true_meg': None,
                'false_meg': 'Adapter Trimming 실행'
            }, p_exit=False)
        exit(1)
    return l_summary_files


def make_excel_for_scythe_log(p_files, p_path, p_order_number):
    """
    Scythe Log 파일을 엑셀로 변환

    :type p_files: list
    :param p_files: tuple(str, str)을 원소로 가지는 리스트. Read1 & Read2 에 대한 Log 파일.

    :type p_path: str
    :param p_path: 엑셀파일 생성 경로.

    :type p_order_number: str
    :param p_order_number: 수주 번호

    :type: None
    :return: None
    """

    def read_log(p_file):
        """
        Scythe Log 파일을 읽는다.

        :type p_file: str
        :param p_file: 경로를 포함한 Scythe Log 파일명

        :return: result.
                 Log 파일의 내용을 포함하는 중첩 리스트
        """
        with open(p_file, 'r') as input_data:
            result = []
            for i in input_data:
                if i.startswith('Adapter Trimming Complete'):
                    continue
                elif i.startswith('contaminated'):
                    i_split = i.strip().split(',')
                    i_data = []
                    for j in i_split:
                        j_split = [x.strip() for x in j.strip().split(':')]
                        i_data.append(j_split)
                    result.extend(i_data)
                elif i.strip() == '':
                    continue
                else:
                    i_data = [x.strip() for x in i.strip().split(':')]
                    result.append(i_data)
        return result

    log_data_r1_list = []
    log_data_r2_list = []
    with progressbar(p_files, label='Scythe Log 파일 읽는 중') as bar__p_list:
        for log_file in bar__p_list:
            log_data_r1_list.append(read_log(log_file[0]))
            log_data_r2_list.append(read_log(log_file[1]))

    import xlsxwriter
    workbook = xlsxwriter.Workbook(
        os.path.join(*(p_path, p_order_number + '.scythe.log.xlsx')))
    worksheet = workbook.add_worksheet()

    # 출력 형식
    number_format = workbook.add_format({'num_format': '#,##0'})

    # 항목 이름 출력
    for row, cell in enumerate(log_data_r1_list[0], 1):
        worksheet.write_string(row, 0, cell[0])  # Read1
        worksheet.write_string(len(log_data_r1_list[0]) + row + 2, 0, cell[0])  # Read2

    # Read1 항목별 값 출력
    for col1, data in enumerate(log_data_r1_list, 1):
        for row1, cell in enumerate(data, 1):
            typeWrite(worksheet, row1, col1, cell[1], number_format)

    # 표 : 시료명 출력
    header_list = [{'header': 'Read1'}]
    for file in p_files:
        header_list.append({'header': os.path.split(file[0])[1].replace('_1.scythe.log', '')})
    # 표 지정 : Read1
    worksheet.add_table(0, 0, row1, col1,
                        {'columns': header_list, 'style': 'Table Style Light 8', 'autofilter': False})

    # Read2 항목별 값 출력
    for col2, data in enumerate(log_data_r2_list, 1):
        for row2, cell in enumerate(data, row1 + 3):
            typeWrite(worksheet, row2, col2, cell[1], number_format)
    # 표 지정 : Read2
    header_list[0] = {'header': 'Read2'}
    worksheet.add_table(row1 + 2, 0, row2, col2,
                        {'columns': header_list, 'style': 'Table Style Light 8', 'autofilter': False})

    # 그래프 : adapter contamination rate
    contamination_rate_chart = workbook.add_chart({'type': 'column'})
    contamination_rate_chart.add_series({
        'values': ['Sheet1', row1, 1, row1, col1],
        'categories': ['Sheet1', 0, 1, 0, col1],
        'name': ['Sheet1', 0, 0],
    })
    contamination_rate_chart.add_series({
        'values': ['Sheet1', row2, 1, row2, col2],
        'categories': ['Sheet1', row1 + 2, 1, row1 + 2, col1],
        'name': ['Sheet1', row1 + 2, 0],
    })
    contamination_rate_chart.set_x_axis({'name': 'Sample Name'})
    contamination_rate_chart.set_title({'name': 'Scythe - Rate of Adapter Contamination'})
    contamination_rate_chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
    worksheet.insert_chart(row2 + 3, 0, contamination_rate_chart)

    workbook.close()
    secho('>>> Scythe.log.xlsx 생성', fg='cyan')


def make_excel_for_seqpurge_summary(p_files, p_path, p_order_number):
    """
    SeqPurge.summary.txt 파일들을 엑셀파일로 변환

    :type p_files: list
    :param p_files: 경로를 포함하는 SeqPurge Summary 파일명을 원소로 가지는 리스트.

    :type p_path: str
    :param p_path: 엑셀파일 생성 경로.

    :type p_order_number: str
    :param p_order_number: 수주번호.

    :rtype: None
    :return: None
    """
    data_1 = dict()
    data_2 = dict()  # length_distribution
    data_3 = dict()  # read1_error
    data_4 = dict()  # read2_error
    data_5 = dict()  # error_distribution
    for summary_file in p_files:
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        temp5 = []
        sample_name = os.path.split(summary_file)[1].replace('.SeqPurge.summary', '')
        with open(summary_file, 'r') as o_summary:
            indicator = 1
            for i in o_summary:
                i_split = i.strip().split(':')
                if i == '\n':
                    continue
                elif 'Read length distribution after trimming' in i:
                    indicator = 2
                    temp2.append(i_split)
                elif 'Read error per cycle (read 1)' in i:
                    indicator = 3
                    temp3.append(i_split)
                elif 'Read error per cycle (read 2)' in i:
                    indicator = 4
                    temp4.append(i_split)
                elif 'Read error count distribution' in i:
                    indicator = 5
                    temp5.append(i_split)
                else:
                    if indicator == 1:
                        temp1.append(i_split)
                    elif indicator == 2:
                        temp2.append(i_split)
                    elif indicator == 3:
                        temp3.append(i_split)
                    elif indicator == 4:
                        temp4.append(i_split)
                    elif indicator == 5:
                        temp5.append(i_split)
                    else:
                        raise RuntimeError(style('파싱 알고리즘 오류!', fg='red', blink=True))
        data_1[sample_name] = temp1
        data_2[sample_name] = temp2
        data_3[sample_name] = temp3
        data_4[sample_name] = temp4
        data_5[sample_name] = temp5

    def find_max_item_name(p_dict):
        max = (None, 0)
        for name in p_dict.keys():
            temp = len(p_dict.get(name))
            if temp > max[1]:
                max = (name, temp)
        return max[0]

    data_2_max_name = find_max_item_name(data_2)
    data_3_max_name = find_max_item_name(data_3)
    data_4_max_name = find_max_item_name(data_4)
    data_5_max_name = find_max_item_name(data_5)

    import xlsxwriter
    workbook = xlsxwriter.Workbook(os.path.join(*(p_path, p_order_number + '.SeqPurge.summary.xlsx')))
    worksheet = workbook.add_worksheet('Summary')
    worksheet2 = workbook.add_worksheet('Trimmed Reads Chart')

    # 출력 형식
    number_format = workbook.add_format({'num_format': '#,##0'})

    sample_name_list = list(data_1.keys())
    sample_name_list.sort()

    # 항목별 이름 출력
    for clause1, cell in enumerate(data_1.get(sample_name_list[0]), 1):
        worksheet.write_string(clause1, 0, cell[0])
    for clause2, cell in enumerate(data_2.get(data_2_max_name), clause1 + 1):
        worksheet.write_string(clause2, 0, cell[0])
    for clause3, cell in enumerate(data_3.get(data_3_max_name), clause2 + 1):
        worksheet.write_string(clause3, 0, cell[0])
    for clause4, cell in enumerate(data_4.get(data_4_max_name), clause3 + 1):
        worksheet.write_string(clause4, 0, cell[0])
    for clause5, cell in enumerate(data_5.get(data_5_max_name), clause4 + 1):
        worksheet.write_string(clause5, 0, cell[0])

    # 데이터 출력
    for col, name in enumerate(sample_name_list, 1):
        for row1, cell in enumerate(data_1.get(name), 1):
            typeWrite(worksheet, row1, col, cell[1], number_format)
        for row2, cell in enumerate(data_2.get(name), clause1 + 1):
            typeWrite(worksheet, row2, col, cell[1], number_format)
        for row3, cell in enumerate(data_3.get(name), clause2 + 1):
            typeWrite(worksheet, row3, col, cell[1], number_format)
        for row4, cell in enumerate(data_4.get(name), clause3 + 1):
            typeWrite(worksheet, row4, col, cell[1], number_format)
        for row5, cell in enumerate(data_5.get(name), clause4 + 1):
            typeWrite(worksheet, row5, col, cell[1], number_format)

        # Trimmed Reads Chart Sheet
        worksheet2.write_string(1, 0, 'Reads trimmed by insert match')
        worksheet2.write_string(2, 0, 'Reads trimmed by adapter match')
        worksheet2.write_string(3, 0, 'Trimmed reads')
        worksheet2.write_string(4, 0, 'Removed reads')
        worksheet2.write_string(5, 0, 'Removed bases')
        for index, data in enumerate(data_1.get(name)):
            if data[0] == 'Reads trimmed by insert match':
                typeWrite(worksheet2, 1, col, data[1], number_format)
            elif data[0] == 'Reads trimmed by adapter match':
                typeWrite(worksheet2, 2, col, data[1], number_format)
            elif data[0] == 'Trimmed reads':
                cell = data[1].strip().split('(')[1].split('%')[0]  # ex)  162 of 408154 (0.04%)
                typeWrite(worksheet2, 3, col, cell, number_format)
            elif data[0] == 'Removed reads':
                cell = data[1].strip().split('(')[1].split('%')[0]  # ex)  26 of 408154 (0.01%)
                typeWrite(worksheet2, 4, col, cell, number_format)
            elif data[0] == 'Removed bases':
                cell = data[1].strip().split('%')[0]  # ex)  0.01%
                typeWrite(worksheet2, 5, col, cell, number_format)  # ex)  0.01%
            else:
                continue

    # 표 지정
    header_list = [{'header': 'Summary'}]
    for name in sample_name_list:
        header_list.append({'header': name})
    worksheet.add_table(0, 0, clause1, col,
                        {'columns': header_list, 'style': 'Table Style Light 8', 'autofilter': False})
    worksheet2.add_table(0, 0, 5, col,
                         {'columns': header_list, 'style': 'Table Style Light 8', 'autofilter': False})
    header_list[0] = {'header': data_2.get(data_2_max_name)[0][0]}
    worksheet.add_table(clause1 + 1, 0, clause2, col, {'columns': header_list, 'style': 'Table Style Light 8'})
    header_list[0] = {'header': data_3.get(data_3_max_name)[0][0]}
    worksheet.add_table(clause2 + 1, 0, clause3, col, {'columns': header_list, 'style': 'Table Style Light 8'})
    header_list[0] = {'header': data_4.get(data_4_max_name)[0][0]}
    worksheet.add_table(clause3 + 1, 0, clause4, col, {'columns': header_list, 'style': 'Table Style Light 8'})
    header_list[0] = {'header': data_5.get(data_5_max_name)[0][0]}
    worksheet.add_table(clause4 + 1, 0, clause5, col, {'columns': header_list, 'style': 'Table Style Light 8'})

    # 차트 추가
    length_distribution_chartsheet = workbook.add_chartsheet('Length Chart')
    error_distribution_chartsheet = workbook.add_chartsheet('Error Chart')

    trimmed_reads_bar_chart = workbook.add_chart({'type': 'column'})
    trimmed_reads_line_chart = workbook.add_chart({'type': 'line'})
    length_distribution_chart = workbook.add_chart({'type': 'line'})
    error_distribution_chart = workbook.add_chart({'type': 'line'})

    for col in range(1, len(sample_name_list) + 1):
        length_distribution_chart.add_series(
            {
                'values': ['Summary', clause1 + 2, col, clause2, col],
                'categories': ['Summary', clause1 + 2, 0, clause2, 0],
                'name': ['Summary', clause1 + 1, col],
            }
        )
        error_distribution_chart.add_series(
            {
                'values': ['Summary', clause4 + 2, col, clause5, col],
                'categories': ['Summary', clause4 + 2, 0, clause5, 0],
                'name': ['Summary', clause4 + 1, col],
            }
        )
    for row in range(1, 3):
        trimmed_reads_bar_chart.add_series(
            {
                'values': ['Trimmed Reads Chart', row, 1, row, len(sample_name_list)],
                'categories': ['Trimmed Reads Chart', 0, 1, 0, len(sample_name_list)],
                'name': ['Trimmed Reads Chart', row, 0],
            }
        )
    for row in range(3, 6):
        trimmed_reads_line_chart.add_series(
            {
                'values': ['Trimmed Reads Chart', row, 1, row, len(sample_name_list)],
                'categories': ['Trimmed Reads Chart', 0, 1, 0, len(sample_name_list)],
                'name': ['Trimmed Reads Chart', row, 0],
            }
        )

    trimmed_reads_bar_chart.set_x_axis({'name': 'sample name'})
    trimmed_reads_bar_chart.set_y_axis({'name': 'read count'})
    trimmed_reads_bar_chart.set_title({'name': 'SeqPurge - The number of trimmed reads by methods'})
    trimmed_reads_bar_chart.set_legend({'position': 'top'})
    trimmed_reads_line_chart.set_x_axis({'name': 'sample name'})
    trimmed_reads_line_chart.set_y_axis({'name': 'percentage(%)'})
    trimmed_reads_line_chart.set_title({'name': 'SeqPurge - Percentage of trimmed data'})
    trimmed_reads_line_chart.set_legend({'position': 'top'})
    length_distribution_chart.set_x_axis({'name': 'base position'})
    length_distribution_chart.set_y_axis({'name': 'read count'})
    length_distribution_chart.set_title({'name': 'SeqPurge - Read length distribution after trimming'})
    error_distribution_chart.set_x_axis({'name': 'count'})
    error_distribution_chart.set_y_axis({'name': 'the number of error'})
    error_distribution_chart.set_title({'name': 'SeqPurge - Read error count distribution'})
    trimmed_reads_bar_chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
    trimmed_reads_line_chart.set_size({'x_scale': 3.2, 'y_scale': 1.1})
    worksheet2.insert_chart('A8', trimmed_reads_bar_chart)
    worksheet2.insert_chart('A25', trimmed_reads_line_chart)
    length_distribution_chartsheet.set_chart(length_distribution_chart)
    error_distribution_chartsheet.set_chart(error_distribution_chart)
    workbook.close()
    secho('>>> 엑셀 파일 생성 : SeqPurge.Summary', )


def make_excel(kargs, p_tool):
    """
    Adapter Trimming 에 대한 Summary 또는 Log 파일을 정리하여 엑셀 파일로 만든다.
    fastp tool은 HTML 로 제공되므로 엑셀 파일이 생성되지 않는다.

    :type kargs: dict
    :param kargs: 다음을 key로 가지는 Dictionary
                files: 경로를 포함하는 SeqPurge Summary 파일명을 원소로 가지는 리스트.
                output_path: 엑셀파일 생성 경로.
                order_number: 수주번호.

    :type p_tool: str
    :param p_tool: Adapter Trimming Tool('fastp', 'SeqPurge', 'Scythe')

    :rtype : None
    :return: None
    """
    if p_tool == 'fastp':
        pass
    elif p_tool == 'SeqPurge':
        make_excel_for_seqpurge_summary(kargs['files'], kargs['output_path'], kargs['order_number'])
    elif p_tool == 'scythe':
        make_excel_for_scythe_log(kargs['files'], kargs['output_path'], kargs['order_number'])

