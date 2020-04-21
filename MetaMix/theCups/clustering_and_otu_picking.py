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

from SpoON.util import run_cmd, check_run_cmd, cast_num, read_metadata_for_sort, sort_data_by_custom
from theCups.read_assembly import transpose_data
from collections import namedtuple
from click import secho, echo
import os

OTU_Table_Summary = namedtuple('otu_table_summary',
                               ['sample_count', 'observation_count', 'total_count', 'min_count',
                                'max_count', 'median_count', 'mean_count', 'l_sample_data'])


class Summary:
    def __init__(self, p_kargs):
        """
        de novo 방식(CD-HIT-OTU) 과 Closed-Reference 방식(UCLUST)에 대해서
        summary.html 페이지 작성에 필요한 Base Class 입니다.

        :param p_kargs:
                otu: OTU or otus.rep.fasta
                stat: STAT.txt
                otu_table_summary: otu_table_summary.txt
                metadata: metadata.txt
        """
        self.kargs = p_kargs
        self.__metadata_order = None
        self.__total_sample_read_count = None
        self.__gamma_diversity = None
        self.__otu_table_summary = None

    @property
    def metadata_order(self):
        return self.__metadata_order

    def read_metadata(self):
        metadata_order = read_metadata_for_sort(self.kargs['metadata'])
        self.__metadata_order = metadata_order

    @property
    def total_sample_read_count(self):
        return self.__total_sample_read_count

    def read_stat(self):
        with open(self.kargs['stat'], 'r') as o_stat:
            read_count_sum = 0
            for i in o_stat:
                if i.startswith('==>') or i.startswith('\n'):
                    continue
                else:
                    i_split = i.split()
                    read_count_sum += int(i_split[2])
        self.__total_sample_read_count = read_count_sum

    @property
    def otu_table_summary(self):
        return self.__otu_table_summary

    def read_otu_table_summary(self):
        with open(self.kargs['otu_table_summary'], 'r') as o_otu_summary:
            sample_data_indicator = False
            for i in o_otu_summary:
                if "Num samples:" in i:
                    _, sample_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Num observations:" in i:
                    _, observation_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Total count:" in i:
                    _, total_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Min:" in i:
                    _, min_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Max:" in i:
                    _, max_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Median:" in i:
                    _, median_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Mean:" in i:
                    _, mean_count = cast_num(i.strip().split(":")[1].lstrip().replace(',', ''))
                elif "Counts/sample detail:" in i:
                    sample_data_list = []
                    sample_data_indicator = True
                elif sample_data_indicator is True:
                    temp_data = [x.strip() for x in i.strip().split(":")]
                    if len(temp_data) == 2:
                        _, float_count = cast_num(temp_data[1].replace(',', ''))
                        temp_data[1] = int(float(float_count))
                    else:
                        raise IndexError('otu_table_summary 파일의 파싱 조건이 다릅니다.')
                    try:
                        sample_data_list.append(temp_data)
                    except NameError:
                        continue
                else:
                    continue

        sorted_sample_data_list = sort_data_by_custom(sample_data_list, self.metadata_order)
        self.__otu_table_summary = OTU_Table_Summary(
            sample_count=sample_count,
            observation_count=observation_count,
            total_count=total_count,
            min_count=min_count,
            max_count=max_count,
            median_count=median_count,
            mean_count=mean_count,
            l_sample_data=sorted_sample_data_list)

    @property
    def gamma_diversity(self):
        return self.__gamma_diversity

    def count_otus(self):
        cmd = f'grep -c ">" {self.kargs["otu"]}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_exit=False)
        self.__gamma_diversity = int(run.stdout)


class SummaryForCdHitOtu(Summary):
    def __init__(self, p_kargs):
        """
        CD-HIT-OTU 프로그램의 결과를 바탕으로 summary.html 페이지 작성에 필요한 데이터를 확인합니다.

        :param p_kargs: 아래의 이름을 Key 로 가지는 딕션너리.
                otu: OTU or otus.rep.fasta
                chimeric: chimaric.ids
                stat: STAT.txt
                trim_log: trim.log
                otu_table_summary: otu_table_summary.txt
                metadata: metadata.txt
                otu_nat: OTU.nr2nd.clstr.NAT.txt (MiSeq_V2일 경우)
        """
        super().__init__(p_kargs)
        self.__ambiguous_base_calls = None
        self.__wrong_prefix_or_primers = None
        self.__primer_seq = None
        self.__low_quality_bases = None
        self.__chimeric_count = None
        self.__otu_nat = None

    @property
    def ambiguous_base_calls(self):
        return self.__ambiguous_base_calls

    @property
    def wrong_prefix_or_primers(self):
        return self.__wrong_prefix_or_primers

    @property
    def primer_seq(self):
        return self.__primer_seq

    @property
    def low_quality_bases(self):
        return self.__low_quality_bases

    def read_trim_log(self):
        with open(self.kargs['trim_log'], 'r') as o_trim:
            for i in o_trim:
                if "filtered seqs with ambiguous base calls" in i:
                    text = i.strip().split(':')[1].lstrip()
                    _, ambiguous_base_calls = cast_num(text)
                elif "filtered seqs with wrong prefix or primers" in i:
                    text = i.strip().split(':')[1].lstrip()
                    _, wrong_prefix_or_primers = cast_num(text)
                elif "      primers" in i:
                    text = i.strip().split(':')[1].lstrip()
                    primer_seq = text
                elif "filtered seqs with low-quality bases" in i:
                    text = i.strip().split(':')[1].lstrip()
                    _, low_quality_bases = cast_num(text)
        self.__ambiguous_base_calls = int(ambiguous_base_calls)
        self.__wrong_prefix_or_primers = int(wrong_prefix_or_primers)
        self.__primer_seq = primer_seq
        self.__low_quality_bases = int(low_quality_bases)

    @property
    def chimeric_count(self):
        return self.__chimeric_count

    def count_chimeric_id(self):
        cmd = f'wc -l {self.kargs["chimeric"]}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_exit=False)
        self.__chimeric_count = int(run.stdout.split()[0])  # '48629 /crystal/~~/chimaric.ids'

    @property
    def other_count(self):
        other = self.total_sample_read_count \
            - self.otu_table_summary.total_count \
            - self.ambiguous_base_calls \
            - self.wrong_prefix_or_primers \
            - self.low_quality_bases \
            - self.chimeric_count
        return other

    @property
    def otu_nat(self):
        return self.__otu_nat

    def read_otu_nat(self):
        with open(self.kargs['otu_nat'], 'r') as o_nat:
            l_data = list()
            for i in o_nat:
                l_data.append(i.strip().split())
        l_data[0].insert(0, 'NAT')
        l_trans_data = transpose_data(l_data)
        l_trans_header = l_trans_data[0]
        l_sorted_trans_data = sort_data_by_custom(l_trans_data[1:], self.metadata_order)
        l_sorted_trans_data.insert(0, l_trans_header)
        l_sorted_data = transpose_data(l_sorted_trans_data)
        self.__otu_nat = l_sorted_data


class SummaryForClosed(Summary):
    def __init__(self, p_kargs):
        """
        QIIME1 Closed-Ref. Method(UCLUST) 프로그램의 결과를 바탕으로 summary.html 페이지 작성에 필요한 데이터를 확인합니다.

        :param p_kargs:
                    picked_log:
                    picked_failures:
                    otu: XXX_rep_set.fasta
                    stat: XXX.F_length.csv
                    otu_table_summary: otu_table_summary.txt
                    metadata: metadata.txt or mappingfile.txt
        """
        super().__init__(p_kargs)
        self.__total_failures_from_log = None
        self.__total_failures_from_failures = None
        self.__similarity = None

    @property
    def total_failures_from_log(self):
        return self.__total_failures_from_log

    @property
    def similarity(self):
        return self.__similarity

    def read_picked_log(self):
        if self.kargs['picked_log'] is None:
            return
        with open(self.kargs['picked_log'], 'r') as o_log:
            total_failures = 0
            s_similarity = set()
            for i in o_log:
                if i.startswith('Num failures'):
                    total_failures += int(i.strip().split(':')[1])
                elif i.startswith('Similarity'):
                    s_similarity.update([i.strip().split(':')[1]])

        if len(s_similarity) == 1:
            similarity = list(s_similarity)[0]
        else:
            secho(f'Error: {os.path.basename(self.kargs["picked_log"])} 파일에 기재된 '
                  f'Similarity(Clustering Cut-off)의 값이 모두 동일하지 않습니다.', fg='red', blink=True)
            secho('\t -->  발생할 수 없는 사항입니다. 개발자와 상의하세요.', fg='magenta')
            echo(self.kargs['picked_log'])
        self.__total_failures_from_log = total_failures
        self.__similarity = int(float(similarity) * 100)

    @property
    def total_failures_from_failures(self):
        return self.__total_failures_from_failures

    def read_picked_failures(self):
        if self.kargs['picked_failures'] is None:
            return
        cmd = f'wc -l {self.kargs["picked_failures"]}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': None,
            }, p_stdout=False, p_exit=False
        )
        # wc -l 결과 예시:  '555542 HN00112054.pooled_failures.txt'
        self.__total_failures_from_failures = int(run.stdout.strip().split()[0])

    def read_filtered_stat(self):
        with open(self.kargs['stat'], 'r') as o_stat:
            for i in o_stat:
                if i.startswith('Final_Read'):
                    text = i.strip().split(',')
                    read_count_sum = sum([int(x) for x in text[1:]])
        self.__total_sample_read_count = read_count_sum

