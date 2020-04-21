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
__version__ = '1.0.2'

from click import secho, echo, style
import os
from SpoON.util import check_file_type, parse_config, run_cmd, check_run_cmd, get_file_copy_cmd, \
                       launcher_cmd, glob_dir

# 기본값 설정
CONFIG = parse_config()
CD_HIT_OTU_SINGLE = os.path.join(
    CONFIG['CoFI_OTU']['cd-hit-otu']['path'],
    CONFIG['CoFI_OTU']['cd-hit-otu']['single'])
CD_HIT_OTU_SINGLE_SINGLETON = os.path.join(
    CONFIG['CoFI_OTU']['cd-hit-otu']['path'],
    CONFIG['CoFI_OTU']['cd-hit-otu']['single_singleton'])
CD_HIT_OTU_SINGLE_DOUBLETON = os.path.join(
    CONFIG['CoFI_OTU']['cd-hit-otu']['path'],
    CONFIG['CoFI_OTU']['cd-hit-otu']['single_doubleton'])
CD_HIT_OTU_CLSTR_SAMPLE_COUNT = os.path.join(
    CONFIG['CoFI_OTU']['cd-hit-otu']['path'],
    CONFIG['CoFI_OTU']['cd-hit-otu']['clstr_sample_count'])
CD_HIT_OTU_CLSTR_SAMPLE_COUNT_MATRIX = os.path.join(
    CONFIG['CoFI_OTU']['cd-hit-otu']['path'],
    CONFIG['CoFI_OTU']['cd-hit-otu']['clstr_sample_count_matrix'])
CLOSED_OTU_POOLWORKER = CONFIG['CoFI_OTU_Closed_PoolWorker']
CLOSED_PICKING_DB = CONFIG['CoFI_OTU']['closed_otu']['database']
ANALYSIS_DIR_NAME = CONFIG['Analysis_Dir_Name']


def get_primer_seq_file(p_seq, p_path):
    """
    Primer Sequence 을 CD-HIT-OTU 의 Primer File 형식으로 변환한다.

    :type p_seq: str
    :param p_seq: Primer Sequence

    :type p_path: str
    :param p_path: 파일 생셩 경로

    :rtype: str
    :return: temp_file.name. 파일 이름
    """
    from tempfile import NamedTemporaryFile
    temp_file = NamedTemporaryFile('wt', dir=p_path, prefix='primer.', suffix='.seq', encoding='utf-8', delete=False)
    os.chmod(temp_file.name, 0o0660)
    for index, symbol in enumerate(p_seq):
        if symbol.upper() == 'W':    # weak
            base = '[AT]'
        elif symbol.upper() == 'S':  # strong
            base = '[CG]'
        elif symbol.upper() == 'M':  # amino
            base = '[AC]'
        elif symbol.upper() == 'K':  # keto
            base = '[GT]'
        elif symbol.upper() == 'R':  # purine
            base = '[AG]'
        elif symbol.upper() == 'Y':  # pyrimidine
            base = '[CT]'
        elif symbol.upper() == 'B':  # not A (B comes after A)
            base = '[CGT]'
        elif symbol.upper() == 'D':  # not C (D comes after C)
            base = '[AGT]'
        elif symbol.upper() == 'H':  # not G (H comes after G)
            base = '[ACT]'
        elif symbol.upper() == 'V':  # not T (V comes after T and U)
            base = '[ACG]'
        elif symbol.upper() == 'N':  # any base (not a gap)
            base = '[ACGT]'
        elif symbol.upper() in ['A', 'C', 'G', 'T']:
            base = symbol.upper()
        else:
            secho('>>> Primer Sequences Error: "{}" 문자는 염기에 해당하지 않습니다.'.format(symbol), fg='red')
            echo('{front}{error}{back}'.format(
                front=p_seq[:index],
                error=style(p_seq[index], fg='red'),
                back=p_seq[index+1:])
            )
            exit()
        temp_file.write(base)
    temp_file.write('\n')
    temp_file.close()
    echo('>>> Primer Sequences File 생성')
    echo(temp_file.name)
    return temp_file.name


def get_cd_hit_otu_single_cmd(kargs):
    """
    CD-HIT-OTU 프로그램의 single 스크립트 실행 명령어를 반환한다.

    ex) 실행 코드 예시
    cd-hit-otu-all-single.pl -i $1 -o $1"_97" -t $2 -c 0.97 -p 0 -m true -e 0.01
    cd-hit-otu-all-single.cutoff_1.pl : Singleton
    cd-hit-otu-all-single.cutoff_2.pl : Doubleton

    :type kargs: dict
    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            pooled_sample: str
            prefix_length_or_primer: int or str
            out_dir: str
            trim_cutoff: int or float
            otu_cutoff: float
            prefix_primer: int or float
            check_chimera: bool
            pcr_error: float

    :return: cmd
    """
    # """
    # cd-hit-otu-illumina-0.0.1/cd-hit-otu-all-single.pl
    # Usage:
    #
    # /lustre/Tools/CD-HIT-OTU/cd-hit-otu-illumina-0.0.1/cd-hit-otu-all-single.pl -i input_fasta_file
    # -o output_dir -f length-cutoff-lower-fraction  -p prefix-length/primers_file -e sequence_error_rate -c OTU_cutoff
    #
    # Parameters:
    # -i input fastq file
    # -o output dir
    # -t trim_cutoff, default 1.0 (means no trimming)
    #         if cutoff is a integer number > 1 (like 200), the program will trim reads to this length
    #         if cutoff is a fraction (like 0.8), the program will trim reads to its length times this fraction
    # -p prefix-length/primers_file default 6
    #         if a primers_file is provided,
    #           read primers from this file, remove the reads don't match the primers
    #         if a prefix-length (a digit number) is provided
    #           get the consensus of prefix of the all reads
    #           remove the reads without this consensus
    # -c OTU cutoff default 0.97
    # -m whether to perform chimera checking (true/false), default true
    # -e per-base PCR error, default 0.01
    #
    # format of primer_file, each line is a primer sequence, where [AT] mean either A or T
    # AGTGCGTAGTG[ACTG]CAGC[AC]GCCGCGGTAA
    # """
    global CD_HIT_OTU_SINGLE, CD_HIT_OTU_SINGLE_SINGLETON, CD_HIT_OTU_SINGLE_DOUBLETON
    # Prefix_Length or Primer Sequences or Primer File 확인.
    if kargs['prefix_length_or_primer'].isdigit():
        prefix_primer = int(kargs['prefix_length_or_primer'])
    elif kargs['prefix_length_or_primer'].isalpha():
        path = os.path.split(kargs['out_dir'])[0]
        prefix_primer = get_primer_seq_file(kargs['prefix_length_or_primer'], path)
    elif check_file_type(kargs['prefix_length_or_primer'], 'isfile'):
        prefix_primer = kargs['prefix_length_or_primer']
    else:
        raise ValueError('{} : 지원하지 않는 값입니다.'.format(style(kargs['prefix_length_or_primer'], fg='red')))

    # Singleton 또는 Doubleton 포함 여부 확인
    if kargs['small_cluster'] is None:
        cd_hit_otu = CD_HIT_OTU_SINGLE
    elif kargs['small_cluster'] is 'singleton':
        cd_hit_otu = CD_HIT_OTU_SINGLE_SINGLETON
    elif kargs['small_cluster'] is 'doubleton':
        cd_hit_otu = CD_HIT_OTU_SINGLE_DOUBLETON
    else:
        raise ValueError('{} : 지원하지 않는 값입니다.'.format(style(kargs['small_cluster'], fg='red')))

    cmd = '{cd_hit_otu} ' \
        '-i {pooled_sample} ' \
        '-o {out_dir} ' \
        '-t {trim_cutoff} ' \
        '-c {otu_cutoff} ' \
        '-p {prefix_primer} ' \
        '-m {chimera} ' \
        '-e {pcr_error} '.format(
            cd_hit_otu=cd_hit_otu,
            pooled_sample=kargs['pooled_sample'],
            out_dir=kargs['out_dir'],
            trim_cutoff=kargs['trim_cutoff'],
            otu_cutoff=kargs['otu_cutoff'],
            prefix_primer=prefix_primer,
            chimera='true',  # TODO : True(Python) --> true(Perl)
            pcr_error=kargs['pcr_error']
            )
    return cmd


def get_clstr_sample_count_cmd(p_path):
    """
    CD-HIT-OTU 프로그램의 clstr_sample_count.pl 스크립트를 실행하는 명령어를 반환한다.
    OTU Cluster별 시료의 개수 및 Read의 개수를 산출하여 정리한다.

    생성파일
     - OTU.nr2nd.clstr.core.txt
     - OTU.nr2nd.clstr.smp

    ex) /lustre/Tools/CD-HIT-OTU/cd-hit-otu-illumina-0.0.1/clstr_sample_count.pl _ $1"_97/OTU.nr2nd.clstr"
    "_" is the delimiter between sample name and sequence name

    :type p_path: str
    :param p_path: OTU.nr2nd.clstr 파일
    :return: cmd
    """
    global CD_HIT_OTU_CLSTR_SAMPLE_COUNT
    cmd = '{cd_hit_otu} ' \
        '_ ' \
        '{file} '.format(
            cd_hit_otu=CD_HIT_OTU_CLSTR_SAMPLE_COUNT,
            file=os.path.join(p_path, 'OTU.nr2nd.clstr'))
    return cmd


def get_clstr_sample_count_matrix(p_path):
    """
    CD-HIT-OTU 프로그램의 clstr_sample_count_matrix.pl 스크립트를 실행하는 명령어를 반환한다.
    OTU Cluster에 해당되는 각 시료의 Read 개수를 매트릭스 형태로 정리한다.

    생성파일
     - OTU.nr2nd.clstr.NAT.txt
     - OTU.nr2nd.clstr.otu.txt
     - OTU.nr2nd.clstr.sample.txt

    eX) /lustre/Tools/CD-HIT-OTU/cd-hit-otu-illumina-0.0.1/clstr_sample_count_matrix.pl _ $1"_97/OTU.nr2nd.clstr"
    "_" is the delimiter between sample name and sequence name

    :param p_path:
    :return
    """
    global CD_HIT_OTU_CLSTR_SAMPLE_COUNT_MATRIX
    cmd = '{cd_hit_otu} ' \
        '_ ' \
        '{file} '.format(
            cd_hit_otu=CD_HIT_OTU_CLSTR_SAMPLE_COUNT_MATRIX,
            file=os.path.join(p_path, 'OTU.nr2nd.clstr'))
    return cmd


def run_cd_hit_otu(kargs):
    """
    CD-HIT-OTU 프로그램의 single 스크립트를 이용하여 denovo 방식의 OTU Picking을 진행한다.
    CD-HIT-OTU 작업 디렉터리 생성
    CD-HIT-OTU 실행
    STAT.txt 파일 복사

    :param kargs: 다음의 항목을 키로 가지는 딕션너리
            pooled_sample: str
            out_dir: None or str
            otu_cutoff: float
            prefix_length_or_primer: int or str
            trim_cutoff: int or float
            small_cluster: str
            check_chimera: bool
            pcr_error: float
            analysis_number: str

    :rtype: str
    :return: out_dir : CD-HIT-OTU 생성 디렉터리
    """
    global ANALYSIS_DIR_NAME
    if (kargs['order_number'] != 'None') and (kargs['analysis_number'] != 'None'):
        analysis_number_path = os.path.join(kargs['analysis_base_path'], kargs['order_number'],
                                            kargs['analysis_number'])
        read_assembly_path = os.path.join(analysis_number_path, ANALYSIS_DIR_NAME['read_assembly'])
        pooled_file = glob_dir(read_assembly_path, f"*{kargs['order_number']}*.pooled.*fastq*", 'only')
        if pooled_file is not None:
            pooled_file_name = os.path.basename(pooled_file)
            kargs['pooled_sample'] = pooled_file
        cd_hit_otu_dir_path = os.path.join(analysis_number_path, ANALYSIS_DIR_NAME['cd_hit_otu'])
    elif (kargs['pooled_sample'] != 'None') or (kargs['pooled_sample'] is not None):
        # pooled_sample 파일로부터 작업 경로 파싱
        if os.path.isabs(kargs['pooled_sample']):
            pooled_sample_file = kargs['pooled_sample']
        else:
            pooled_sample_file = os.path.abspath(kargs['pooled_sample'])
        read_assembly_path, pooled_file_name = os.path.split(pooled_sample_file)
        if read_assembly_path is '':
            cwd = os.getcwd()
            cd_hit_otu_dir_path = os.path.join(os.path.split(cwd)[0], ANALYSIS_DIR_NAME['cd_hit_otu'])
        elif read_assembly_path.lower().startswith(ANALYSIS_DIR_NAME['read_assembly'].lower()):
            cwd = os.getcwd()
            cd_hit_otu_dir_path = os.path.join(cwd, ANALYSIS_DIR_NAME['cd_hit_otu'])
        else:
            cd_hit_otu_dir_path = os.path.join(os.path.split(read_assembly_path)[0], ANALYSIS_DIR_NAME['cd_hit_otu'])

        dir_count = cd_hit_otu_dir_path.upper().count(ANALYSIS_DIR_NAME['cd_hit_otu'])
        if dir_count > 1:
            secho('Error: 경로에 CD-HIT-OTU 디렉터리명이 중복({})으로 존재합니다.'.format(dir_count), fg='red')
            echo('개발자에게 문의하세요. 또는 -ps 옵션의 값을 절대경로로 입력하세요.', fg='yellow')
    else:
        raise RuntimeError('예기치 않은 오류. Argument & Option 조합 오류.')

    if kargs['out_dir'] is None:
        out_dir_name = '{file}_{otu_cutoff}'.format(file=pooled_file_name, otu_cutoff=int(kargs['otu_cutoff']*100))
    else:
        out_dir_name = kargs['out_dir']

    # CD-HIT-OTU 디렉터리 생성
    if check_file_type(cd_hit_otu_dir_path, 'exists'):
        out_dir = os.path.join(cd_hit_otu_dir_path, out_dir_name)
        if check_file_type(out_dir, 'exists'):
            secho('>>> Error: 같은 이름의 디렉터리가 존재합니다.', fg='red', blink=True)
            secho('Warning: {}'.format(cd_hit_otu_dir_path), fg='yellow')
            secho('Error: {}'.format(out_dir), fg='red')
            exit()
        else:
            secho('>>> Warning: 같은 이름의 디렉터리가 존재합니다.', fg='yellow')
            echo(cd_hit_otu_dir_path)
    else:
        os.mkdir(cd_hit_otu_dir_path)
        echo('>>> CD-HIT-OTU 디렉터리 생성')
        echo(cd_hit_otu_dir_path)
        out_dir = os.path.join(cd_hit_otu_dir_path, out_dir_name)

    kargs['out_dir'] = out_dir
    l_cmd = list()
    l_cmd.append(get_cd_hit_otu_single_cmd(kargs))
    l_cmd.append(get_clstr_sample_count_cmd(kargs['out_dir']))
    l_cmd.append(get_clstr_sample_count_matrix(kargs['out_dir']))
    log_file = os.path.join(*(cd_hit_otu_dir_path, '{}.log'.format(out_dir)))
    with open(log_file, 'w') as log:
        echo('>>> CD-HIT-OTU 시작')
        for count, cmd in enumerate(l_cmd, 1):
            run = run_cmd(cmd)
            small_cluster = '' if kargs['small_cluster'] is None else '({})'.format(kargs['small_cluster'])
            check_run_cmd({
                'run': run,
                'true_meg': 'CD-HIT-OTU-Single{small} {count}완료'.format(small=small_cluster, count=count),
                'false_meg': 'CD-HIT-OTU-Single{small}'.format(small=small_cluster),
            }, p_exit=True, p_stdout=False, p_stderr=True)
            if run.returncode == 0:
                log.write('='*30 + ' CMD ' + '='*30)
                log.write('\n')
                log.write(run.args)
                log.write('\n')
                log.write(run.stdout)
                log.write('*' * 30 + ' ERROR Message ' + '*' * 30)
                log.write('\n')
                log.write(run.stderr)
            elif run.returncode == 1:
                log.write('='*30 + ' CMD ' + '='*30)
                log.write('\n')
                log.write(run.args)
                log.write('\n')
                log.write(run.stdout)
                log.write('\n')
                log.write('*'*30 + ' ERROR Message ' + '*'*30)
                log.write('\n')
                log.write(run.stderr)
    del run

    # STAT.txt 복사
    copy_cmd = get_file_copy_cmd(os.path.join(read_assembly_path, 'STAT.txt'), cd_hit_otu_dir_path)
    run = run_cmd(copy_cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'STAT.txt 복사 완료',
            'false_meg': 'STAT.txt 복사',
        }, p_exit=False)
    return out_dir


def make_pick_otus_rep_file(p_out_dir):
    """
    QIIME 프로그램에서 BIOM 파일을 생성하기 위해 필요한 pick_otus.txt 파일을 생성한다.
    OTU 대표 서열에 대한 fasta 파일(otus_rep.fasta)을 생성한다.
    OTU.nr2nd.clstr --> pick_outs.txt 생성
    OTU             --> otus_rep.fasta 생성.

    :type p_out_dir: str
    :param p_out_dir: 대상 파일이 있는 경로
    :return: None
    """
    echo('>>> pick_otus.txt & otus_rep.fasta 파일 생성 시작')
    otu_clstr_file = os.path.join(p_out_dir, 'OTU.nr2nd.clstr')
    echo('>>> OTU.nr2nd.clstr 파일 연결')
    echo(otu_clstr_file)
    pick_otus_file = os.path.join(p_out_dir, 'pick_otus.txt')
    otus_rep_file = os.path.join(p_out_dir, 'otus_rep.fasta')
    d_otu_id = dict()
    with open(pick_otus_file, 'w') as o_pick_otus:
        with open(otu_clstr_file, 'r') as o_clstr:
            for row in o_clstr:
                if row.startswith('>'):
                    _, num = row.strip().lstrip('>').split()
                    if int(num) > 0:
                        o_pick_otus.write('\n')
                    o_pick_otus.write('denovo{}'.format(num))
                else:
                    # 데이터 예시
                    # 0 443nt, >control1.72h.1_14160... *
                    # 10	443nt, >control1.72h.1_82095... at 1:443:1:443/+/100.00%
                    sample_read_id = row.strip().split('>')[1].split('...')[0]
                    if row.strip().endswith('*'):
                        d_otu_id[sample_read_id] = 'denovo{}'.format(num)
                    o_pick_otus.write('\t{}'.format(sample_read_id))
    secho('>>> pick_otus.txt 생성 완료', fg='cyan')
    echo(pick_otus_file)
    del sample_read_id, num, row

    otu_file = os.path.join(p_out_dir, 'OTU')
    echo('>>> OTU 파일 연결')
    echo(otu_file)
    with open(otus_rep_file, 'w') as o_otus_rep:
        with open(otu_file, 'r') as o_otu:
            for row in o_otu:
                if row.startswith('>'):
                    sample_read_id = row.strip().lstrip('>')
                    try:
                        o_otus_rep.write('>{otu_id}\t{read_id}\n'.format(
                            otu_id=d_otu_id[sample_read_id],
                            read_id=sample_read_id))
                    except KeyError as err:
                        secho('OTU 파일의 sample_read_id인 {}가 '
                              'OTU.nr2nd.clstr 파일에서 생성한 딕션너리에 없습니다.'.format(err), fg='red')
                        secho('Error: 개발자와 상의하세요!', fg='red', blink=True)
                        exit()
                else:
                    o_otus_rep.write(row)
    secho('>>> otus_rep.fasta 생성 완료', fg='cyan')
    echo(otus_rep_file)
    return


def get_closed_reference_otu_cmd(kargs):
    """
    Closed-Reference OTU Picking 에 관련된 Pipeline Script를 실행한다.
    다음의 단계를 자동으로 실행한다.
    parallel_pick_otus_uclust_ref.py
    pick_rep_set.py
    parallel_assign_taxonomy_uclust.py
    make_otu_table.py

    :param kargs:
            fasta: str. Assembled Read(fasta)의 파일명.
            out_dir: str. 결과 생성 디렉터리명.
            parameter: str. parameter 파일명.
            assign: bool. Taxonomy Assignment 실행 여부.
            jobs: int. 0 ~ . 작업의 개수. 0일 경우 병렬처리 안함.
    :return:
    """
    cmd = \
        'pick_closed_reference_otus.py ' \
        '-i {fasta} ' \
        '-o {out_dir} ' \
        '--parameter_fp {para} ' \
        '{assign} ' \
        '{parallel}'.format(
            fasta=kargs['fasta'],
            out_dir=kargs['out_dir'],
            para=kargs['parameter'],
            assign='--assign_taxonomy' if kargs['assign'] else '',
            parallel=f'--parallel --jobs_to_start {kargs["jobs"]}' if kargs['jobs'] else '')
    return cmd


def get_parallel_pick_otus_uclust_ref_cmd(kargs):
    """

    :param kargs:
                pooled_fasta: full path to input_fasta_fp
                out_dir: path to store output files
                ref_fasta: full path to reference collection
                similarity: Sequence similarity threshold
                jobs: Number of jobs to start
    :return:
    """
    cmd = \
        'parallel_pick_otus_uclust_ref.py ' \
        '--input_fasta_fp {pooled_fasta} ' \
        '--output_dir {out_dir} ' \
        '--refseqs_fp {ref_fasta} ' \
        '--similarity {similarity} ' \
        '--poll_directly ' \
        '--jobs_to_start {jobs}'.format(
            pooled_fasta=kargs['pooled_fasta'],
            out_dir=kargs['out_dir'],
            ref_fasta=kargs['ref_fasta'],
            similarity=kargs['similarity'],
            jobs=kargs['jobs'],
        )
    return cmd


def get_pick_rep_set_cmd(kargs):
    """
    
    :param kargs:
                pick_otus: Path to input otu mapping file
                pooled_fasta: Path to input fasta file
                log: Path to store log file
                out_file:  Path to store result file
    :return:
    """
    cmd = \
        'pick_rep_set.py ' \
        '--input_file {pick_otus} ' \
        '--fasta_file {pooled_fasta} ' \
        '--log_fp {log} ' \
        '--result_fp {out_file} ' \
        '--rep_set_picking_method most_abundant'.format(
            pick_otus=kargs['pick_otus'],
            pooled_fasta=kargs['pooled_fasta'],
            log=kargs['log'],
            out_file=kargs['out_file'],
        )
    return cmd


def get_parallel_assign_taxonomy_uclust(kargs):
    """

    :param kargs:
                rep_fasta: full path to fasta file containing query sequences
                out_dir: path to store output files
                ref_fasta:
                ref_taxon:
                jobs: Number of jobs to start
                max_accepts: Number of database hits to consider when making an assignment
    :return:
    """
    cmd = \
        'parallel_assign_taxonomy_uclust.py ' \
        '--input_fasta_fp {rep_fasta} ' \
        '--output_dir {out_dir} ' \
        '--reference_seqs_fp {ref_fasta} ' \
        '--id_to_taxonomy_fp {ref_taxon} ' \
        '--poll_directly ' \
        '--jobs_to_start {jobs} ' \
        '--uclust_max_accepts {max_accepts}'.format(
            rep_fasta=kargs['rep_fasta'],
            out_dir=kargs['out_dir'],
            ref_fasta=kargs['ref_fasta'],
            ref_taxon=kargs['ref_taxon'],
            jobs=kargs['jobs'],
            max_accepts=kargs['max_accepts'],
        )
    return cmd


def make_closed_parameter_file(kargs):
    """

    :param kargs:
                out_file:
                db_fasta:
                db_taxon:
                max_accepts:
                similarity:
    :return:
    """
    prameter_text = f"""pick_otus:refseqs_fp\t{kargs['db_fasta']}
pick_otus:similarity\t{kargs['similarity']}
pick_rep_set:rep_set_picking_method\tmost_abundant
assign_taxonomy:uclust_max_accepts\t{kargs['max_accepts']}
assign_taxonomy:id_to_taxonomy_fp\t{kargs['db_taxon']}
assign_taxonomy:reference_seqs_fp\t{kargs['db_fasta']}
"""
    with open(kargs['out_file'], 'w') as o_out_file:
        o_out_file.write(prameter_text)
    secho('>>> Closed-Ref. Parameter 파일 생성 완료.', fg='cyan')


def run_closed_reference_otu(kargs, p_mode):
    """
    QIIME1.9 의 Closed-Reference OTU Picking을 실행한다.
    입력되는 서열은 FASTA 형식이다.
    p_mode가 pool일 경우, Pick OTUs, Pick Representative Set, Assign Taxonomy 를 각각 진행한다.
    p_mode가 sample 일 경우 pick_closed_reference_otus.py 를 시료별로 실행이 하며, 시료명으로 디렉터리가 생성된다.

    :type kargs: dict
    :param kargs: 아래와 같이 구성된 딕션너리.
            p_mode가 pool일 경우.
                pooled_sample: Pooled Assembled Read 파일(fasta).
            p_mode가 sample일 경우.
                sample_list: named tuple인 SampleList을 원소로 가지는 리스트
                             SampleList.path: 경로
                                       .name: 시료명
                # analysis_path: Analysis_? 디렉터리 경로. myMetaMix용(삭제)
                closed_path: CLOSED 디렉터리 경로.
                assembly_path: Read_Assembly 디렉터리 경로.
                assembly_suffix: Assembled Read 파일의 접미사(fasta).
            공통항목.
                no_parallel: 병렬 처리 안함.
                no_assign: Taxonomy Assignment 실행 안함.
                jobs: 병렬 처리에 사용할 작업의 개수.
                database: OTU Picking 과 Taxonomy Assignment 에 사용할 DB
                similarity: Sequence similarity threshold

    :type p_mode: str
    :param p_mode: pool or sample
    :return:
    """
    global ANALYSIS_DIR_NAME
    echo('>>> Closed-reference OTU Picking 시작')
    if kargs['no_parallel']:
        jobs = 0
    else:
        jobs = kargs['jobs']

    if kargs['no_assign']:
        assign = False
    else:
        assign = True

    db_fasta = CLOSED_PICKING_DB[kargs['database']]['fasta']
    db_taxon = CLOSED_PICKING_DB[kargs['database']]['taxon']
    db_version = CLOSED_PICKING_DB[kargs['database']]['version']
    secho(f'>>> Database: {db_version} 설정', fg='yellow')
    secho(db_fasta)
    secho(db_taxon)
    if kargs['database'] == 'NCBI_16S':
        max_accepts = 1
    elif kargs['database'] == 'NCBI_NT':
        max_accepts = 1
    elif kargs['database'] == 'MIDORI_LONGEST':
        max_accepts = 1
    elif kargs['database'] == 'NCBI_Probiotics':
        max_accepts = 1
    elif kargs['database'].startswith('UNITE'):
        max_accepts = 1
    elif kargs['database'].startswith('GreenGenes'):
        max_accepts = 3
    elif kargs['database'] == 'RDP':
        max_accepts = 3
    elif kargs['database'] == 'MIDORI_UNIQUE':
        max_accepts = 3
    else:
        raise ValueError(f'kargs["database"]: {kargs["database"]}\nmax_accepts 설정 필요!')

    secho(f'>>> UCLUST max_accepts: {max_accepts} 설정', fg='yellow')

    if p_mode == 'pool':
        analysis_path = os.path.dirname(os.path.dirname(kargs['pooled_sample']))
        closed_dir_path = os.path.join(analysis_path, ANALYSIS_DIR_NAME['closed'])
        with open(os.path.join(closed_dir_path, 'db_version.txt'), 'w') as o_version:
            o_version.write(db_version)
            o_version.write('\n')

        # Pick OTUs command
        picked_otu_dir_path = os.path.join(closed_dir_path, 'uclust_ref_picked_otus')
        pick_otus_cmd = get_parallel_pick_otus_uclust_ref_cmd(
            {
                'pooled_fasta': kargs['pooled_sample'],
                'out_dir': picked_otu_dir_path,
                'ref_fasta': db_fasta,
                'similarity': kargs['similarity'],
                'jobs': jobs,
            }
        )
        pick_otus_run = run_cmd(pick_otus_cmd)
        check_run_cmd(
            {
                'run': pick_otus_run,
                'true_meg': 'Pick OTUs 완료',
                'false_meg': 'parallel_pick_otus_uclust_ref',
            }, p_exit=True,
        )

        # Pick representative set command
        file_base_name = os.path.splitext(os.path.basename(kargs['pooled_sample']))[0]
        pick_otus_name = f'{file_base_name}_otus.txt'
        pick_otus_file = os.path.join(closed_dir_path, 'uclust_ref_picked_otus', pick_otus_name)
        rep_set_dir_path = os.path.join(closed_dir_path, 'rep_set')
        log_file_name = f'{file_base_name}_rep_set.log'
        log_file = os.path.join(rep_set_dir_path, log_file_name)
        rep_set_fasta_name = 'otus_rep.fasta'
        rep_set_fasta_file = os.path.join(rep_set_dir_path, rep_set_fasta_name)
        os.makedirs(rep_set_dir_path)
        pick_rep_set_cmd = get_pick_rep_set_cmd(
            {
                'pick_otus': pick_otus_file,
                'pooled_fasta': kargs['pooled_sample'],
                'log': log_file,
                'out_file': rep_set_fasta_file,
            }
        )
        pick_rep_set_run = run_cmd(pick_rep_set_cmd)
        check_run_cmd(
            {
                'run': pick_rep_set_run,
                'true_meg': 'Pick representative set 완료',
                'false_meg': 'pick_rep_set',
            }, p_exit=True,
        )

        # Assign taxonomy command
        if assign is True:
            taxonomy_dir_path = os.path.join(analysis_path, ANALYSIS_DIR_NAME['taxonomy_assignment'])
            taxonomy_results_dir_name = f'uclust_{kargs["database"]}'
            assigned_taxonomy_dir_path = os.path.join(taxonomy_dir_path, taxonomy_results_dir_name)
            os.makedirs(assigned_taxonomy_dir_path)
            assign_taxonomy_cmd = get_parallel_assign_taxonomy_uclust(
                {
                    'rep_fasta': rep_set_fasta_file,
                    'out_dir': assigned_taxonomy_dir_path,
                    'ref_fasta': db_fasta,
                    'ref_taxon': db_taxon,
                    'jobs': kargs['jobs'],
                    'max_accepts': max_accepts,
                }
            )
            assign_taxonomy_run = run_cmd(assign_taxonomy_cmd)
            check_run_cmd(
                {
                    'run': assign_taxonomy_run,
                    'true_meg': 'Assign Taxonomy 완료',
                    'false_meg': 'parallel_assign_taxonomy_uclust',
                }, p_exit=True,
            )
        with open(os.path.join(closed_dir_path, 'CLOSED_CMD.recipe'), 'w') as o_recipe:
            o_recipe.write('# Pick OTUs command\n')
            o_recipe.write(pick_otus_cmd)
            o_recipe.write('\n\n')
            o_recipe.write('# Pick representative set command\n')
            o_recipe.write(pick_rep_set_cmd)
            o_recipe.write('\n\n')

        if assign is True:
            recipe_name = f'CLOSED_TAXONOMY_CMD.{kargs["database"]}.recipe'
            taxonomy_results_dir_path = os.path.join(taxonomy_dir_path, taxonomy_results_dir_name)
            with open(os.path.join(taxonomy_results_dir_path, recipe_name), 'w') as o_recipe:
                o_recipe.write('# Assign taxonomy command\n')
                o_recipe.write(assign_taxonomy_cmd)
                o_recipe.write('\n\n')
            with open(os.path.join(taxonomy_results_dir_path, 'db_version.txt'), 'w') as o_db_version:
                o_db_version.write(db_version)
                o_db_version.write('\n')
        secho('>>> Closed-Reference OTU Picking 완료(QIIME 1.9)', fg='cyan')
    elif p_mode == 'sample':
        global CLOSED_OTU_POOLWORKER
        l_closed_otu_cmd = list()
        closed_parameter_file = os.path.join(kargs['closed_path'], f'closed_parameter.{kargs["database"]}.txt')
        make_closed_parameter_file(
            {
                'out_file': closed_parameter_file,
                'db_fasta': db_fasta,
                'db_taxon': db_taxon,
                'max_accepts': max_accepts,
                'similarity': kargs['similarity'],
            }
        )
        with open(os.path.join(kargs['closed_path'], 'db_version.txt'), 'w') as o_version:
            o_version.write(db_version)
            o_version.write('\n')

        for sample in kargs['sample_list']:
            assembled_fasta = os.path.join(kargs['assembly_path'], 'F_' + sample + kargs['assembly_suffix'])
            out_dir = os.path.join(kargs['closed_path'], sample)
            cmd = get_closed_reference_otu_cmd(
                {
                    'fasta': assembled_fasta,
                    'out_dir': out_dir,
                    'assign': assign,
                    'jobs': jobs,
                    'parameter': closed_parameter_file,
                }
            )
            l_closed_otu_cmd.append(cmd)
        launcher_cmd(l_closed_otu_cmd, 'Closed-reference OTU Picking', CLOSED_OTU_POOLWORKER, False)
    else:
        raise ValueError(f'p_mode: {p_mode}. pool or sample')