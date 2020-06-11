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
__version__ = '1.0.2'

from collections import namedtuple, defaultdict
from click import secho, echo
import os
from SpoON.util import run_cmd, check_run_cmd


RawDataPath = namedtuple('RawDataPath', 'rundata_path sample_list stat_file nread_file')
SampleList = namedtuple('SampleList', 'path cur_name new_name')


class PreMA_Pack:
    def __init__(self, kargs, p_mode='core'):
        """

        :param kargs:
        :type p_mode: str
        :param p_mode: core or mini
        """
        # No Setter
        self.__order_number = kargs['order_number']
        self.__rawdata_base_path = kargs['rawdata_base_path']
        self.__target_dir_suffix = kargs['target_dir_suffix']
        self.__r1_suffix = kargs['r1_suffix']
        self.__r2_suffix = kargs['r2_suffix']
        self.__dic_custom_sample_name = None
        self.__dic_metadata_group = None
        self.__sample_name_mode = kargs['sample_name_mode']
        if p_mode == 'core':
            self.__analysis_base_path = kargs['analysis_base_path']
            self.__copy_fastqc = False if kargs['no_copy_fastqc'] else True
            self.__index_kit = kargs['index_kit']
            self.__sample_read = kargs['sample_read']
            self.__q30 = kargs['q30']
            self.__n_base = kargs['n_base']
            self.__n_read = kargs['n_read']
            self.__adapter_trim_tool = kargs['adapter_trim_tool']
            self.__trim_tail = kargs['trim_tail']
            self.__copy_rawdata = kargs['copy_rawdata']
            self.__my_story = kargs['my_story']
            self.__microbe_and_me = kargs['microbe_and_me']

        # find_rawdata 함수 반환
        self.__nt_rundata_path = None
        # self.__dic_flowcell_id = None

        # check_rawdata 함수 반환
        self.__index_kit_info_form_rawdata = None
        self.__l_nt_sample_list = None

        # copy_rawdata 함수 반환
        self.__l_nt_sample_list_in_analysis = None
        self.__l_nt_sample_list_in_analysis_for_fastqc = None

        # run_adapter_trim 함수 반환
        self.__l_summary_files = None

    def parse_custom_sample_name(self, p_file):
        """
        입력받은 파일을 읽어 시료명을 변경하고, 그룹 정보를 설정한다.
        그룹정보는 metadata.txt 파일을 작성할 때 이용된다.
        그룹정보를 입력할 경우, Header을 입력해야 한다.
        Header는 #SampleID 로 시작해야 되며, 어느 위치에 있어도 된다.
        다만, 중복되는 경우 마지막에 위치한 Header를 사용한다.
        파일의 열 구분자는 공백 또는 탭으로 한다.

        - 파일 양식 -
            #SampleID   New_name   Group1   Group2
            Sample1     Test1      A        A_1
            Sample2     Sample2    B        B_1     : 시료명 변경 안함.
            Sample3     Test3      A        A_2
            #Sample4    Test4      A        A_3     : 제외
            Sample5     Test5      A        A_2

        - Header
            #SampleID : 이름 변경 불가
            New_name: 이름 변경 가능. 2번째 열.
            Group: 이름 변경 가능. 3번째 열과 그 이후에 존재하는 열들은 모두 Group 정보로 간주.

        - 그룹정보 없고, 시료명 변경하는 경우: 변경할 시료명만 입력
            Sample1     Test1
            Sample5     Test5

        - 그룹정보 있고, 시료명 변경하지 않는 경우: 변경할 시료명에 동일한 이름을 입력.
            #SampleID   New_name   Group1   Group2
            Sample1     Sample1       A        A_1
            Sample2     Sample2       B        B_1
            Sample3     Sample3       A        A_2
            #Sample4    Sample4       A        A_3     : 제외
            Sample5     Sample5       A        A_2

        :param p_file: txt 파일
        :return:
        """
        self.__dic_custom_sample_name = defaultdict(str)
        self.__dic_metadata_group = defaultdict(str)
        l_name = list()
        l_new_name = list()
        l_group = list()
        header = None
        dup_indicator1 = False
        dup_indicator2 = False
        no_new_name_indicator = False
        no_group_indicator = False
        with open(p_file, 'r') as file:
            for line in file:
                if line.startswith('#SampleID'):
                    header = line.strip().split()
                    continue
                elif line.startswith('#'):
                    continue
                elif line.strip() is '':
                    continue
                l_cleaned_data = line.strip().split()
                if len(l_cleaned_data) == 2:
                    name, new_name = l_cleaned_data
                    group = [None]
                elif len(l_cleaned_data) > 2:
                    name = l_cleaned_data[0]
                    new_name = l_cleaned_data[1]
                    group = l_cleaned_data[2:]
                l_name.append(name)
                l_new_name.append(new_name)
                l_group.append(group)
            del name, new_name, group
        
        # 시료명 중복 및 시료명 미변경 확인
        set_name = set(l_name)
        set_new_name = set(l_new_name)
        temp_group = list()
        for group in l_group:  # l_group: [[], []]
            temp_group += group
        set_group = set(temp_group)
        if len(set_name) != len(l_name):
            dup_indicator1 = True
        if len(set_new_name) != len(l_new_name):
            if (len(set_new_name) == 1) and (list(set_new_name)[0] is None):
                no_new_name_indicator = True
            else:
                dup_indicator2 = True
        if (len(set_group) == 1) and (list(set_group)[0] is None):
            no_group_indicator = True

        # 시료명 중복일 경우
        if dup_indicator1 or dup_indicator2:
            l_error_name = l_name.copy()
            l_error_new_name = l_new_name.copy()  # 얕은 복사
            try:
                for name in set_name:
                    l_error_name.remove(name)
                for new_name in set_new_name:
                    l_error_new_name.remove(new_name)
            except TypeError:
                pass
            secho('>>> Error: 변경할 시료명 중복!({})'.format(p_file), fg='red', blink=True)
            echo('name: {}'.format(l_error_name))
            echo('new name: {}'.format(l_error_new_name))
            exit()

        # 시로명을 변경하지 않는 경우.
        if no_new_name_indicator is True:
            if no_group_indicator is True:
                pass
            else:
                for name, new_name, group in zip(l_name, l_new_name, l_group):
                    if self.__dic_metadata_group.get(name) is None:
                        self.__dic_metadata_group[name] = group
                else:
                    del header[1]  # New_name 삭제
                    self.__dic_metadata_group['header'] = header
        # 시료명 변경 & 시료명 중복이 없을 경우
        else:
            if no_group_indicator is True:
                for name, new_name in zip(l_name, l_new_name):
                    if self.__dic_custom_sample_name.get(name) is None:
                        self.__dic_custom_sample_name[name] = new_name
                    else:
                        raise RuntimeError('알고리즘 에러. 개발자에게 문의하세요.')
            else:
                for name, new_name, group in zip(l_name, l_new_name, l_group):
                    if self.__dic_custom_sample_name.get(name) is None:
                        self.__dic_custom_sample_name[name] = new_name
                    else:
                        raise RuntimeError('알고리즘 에러(name). 개발자에게 문의하세요.')
                    if self.__dic_metadata_group.get(name) is None:
                        self.__dic_metadata_group[name] = group
                    else:
                        raise RuntimeError('알고리즘 에러(group). 개발자에게 문의하세요.')
                else:
                    if header is None:
                        self.__dic_metadata_group['header'] = ['SampleID']
                    else:
                        del header[1]
                        self.__dic_metadata_group['header'] = header

    def save_metadata_group_file(self, p_out_path):
        if not self.__dic_metadata_group:
            return
        else:
            group_file = os.path.join(p_out_path, 'metadata_group.txt')
            with open(group_file, 'w') as o_group:
                o_group.write('\t'.join(self.__dic_metadata_group['header']))
                o_group.write('\n')
                for nt_sample in self.l_nt_sample_list:
                    group_text = self.__dic_metadata_group[nt_sample.cur_name]
                    if group_text == '':
                        secho('Error: 시료명 & 그룹정보 입력 파일에 기재된 시료명과 실제 생산된 시료명이 다릅니다.',
                              fg='red', blink=True, err=True)
                        echo(f'name: {nt_sample.cur_name}')
                        exit()
                    if group_text[0] is None:
                        o_group.write(nt_sample.new_name)
                        o_group.write('\n')
                    else:
                        group_text = '\t'.join(group_text)
                        o_group.write(f"{nt_sample.new_name}\t{group_text}")
                        o_group.write('\n')
            secho('>>> metadata_group 파일 생성', fg='cyan')
            echo(group_file)

    def copy_custom_sample_name_file(self, p_file):
        initial_metadata_name = f'{self.order_number}_name_group.txt'
        initial_meatadata_file = os.path.join(self.analysis_base_path, self.order_number, initial_metadata_name)
        cmd = f'cp {p_file} {initial_meatadata_file}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': 'initial metadata file 복사 완료',
                'false_meg': 'cp initial_metadata_file',
            }, p_stdout=False,
        )
        echo(initial_meatadata_file)

    def save_custom_sort_file(self, p_file):
        l_name = list()
        with open(p_file, 'r') as o_custom_name:
            for i in o_custom_name:
                if i.startswith('#'):
                    continue
                else:
                    l_name.append(i.strip().split('\t'))
        d_sample_list = dict()
        for tu in self.l_nt_sample_list:
            d_sample_list[tu.cur_name] = tu.new_name
        custom_sort_order_name = f'{self.order_number}_custom_sort.txt'
        custom_sort_order_file = os.path.join(self.analysis_base_path, self.order_number, custom_sort_order_name)
        with open(custom_sort_order_file, 'w') as o_custom_sort:
            for name in l_name:
                new_name = d_sample_list[name[0]]  # ex) SampleID New Group
                o_custom_sort.write(new_name)
                o_custom_sort.write('\n')
        secho('>>> Custom Sort 파일 생성', fg='cyan')
        echo(custom_sort_order_file)

    def check_custom_name(self):
        dk_custom_cur_name = self.__dic_custom_sample_name.keys()
        l_cur_name = [os.path.split(x)[1] for x in self.nt_rundata_path.sample_list]
        l_remained_name = l_cur_name.copy()  # 얕은 복사
        cur_max = max([len(x) for x in l_cur_name])
        new_max = max([len(x) for x in list(self.__dic_custom_sample_name.values())])
        secho('>>> 변경 파일 확인', fg='cyan')
        for name in list(dk_custom_cur_name):
            if name in l_cur_name:
                echo('{cur:>{di1}} --> {new:>{di2}} [OK]'.format(
                    cur=name,
                    di1=cur_max,
                    new=self.__dic_custom_sample_name.get(name),
                    di2=new_max)
                )
                l_remained_name.remove(name)
            else:
                secho('{cur:>{di1}} --> {new:>{di2}} [Not Found]'.format(
                        cur=name,
                        di1=cur_max,
                        new=self.__dic_custom_sample_name.get(name),
                        di2=new_max), fg='red', bold=True)
        if l_remained_name:
            secho('>>> 미변경 시료들', fg='bright_magenta')
            echo('{}'.format(', '.join(l_remained_name)))

    @property
    def order_number(self):
        return self.__order_number

    @property
    def rawdata_base_path(self):
        return self.__rawdata_base_path

    @property
    def analysis_base_path(self):
        return self.__analysis_base_path

    @property
    def target_dir_suffix(self):
        return self.__target_dir_suffix

    @property
    def r1_suffix(self):
        return self.__r1_suffix

    @property
    def r2_suffix(self):
        return self.__r2_suffix

    @property
    def index_kit(self):
        return self.__index_kit

    @property
    def copy_fastqc(self):
        return self.__copy_fastqc

    @property
    def sample_read(self):
        return self.__sample_read

    @property
    def q30(self):
        return self.__q30

    @property
    def n_base(self):
        return self.__n_base

    @property
    def n_read(self):
        return self.__n_read

    @property
    def sample_name_mode(self):
        return self.__sample_name_mode

    @property
    def copy_rawdata(self):
        return self.__copy_rawdata

    @property
    def my_story(self):
        return self.__my_story

    @property
    def microbe_and_me(self):
        return self.__microbe_and_me

    @property
    def dic_custom_sample_name(self):
        return self.__dic_custom_sample_name

    @property
    def dic_metadata_group(self):
        return self.__dic_metadata_group

    @property
    def adapter_trim_tool(self):
        return self.__adapter_trim_tool

    @property
    def trim_tail(self):
        return self.__trim_tail

    # --------------------------------------------

    @property
    def nt_rundata_path(self):
        return self.__nt_rundata_path

    @nt_rundata_path.setter
    def nt_rundata_path(self, p_arg):
        self.__nt_rundata_path = p_arg

    # @property
    # def dic_flowcell_id(self):
    #     return self.__dic_flowcell_id
    #
    # @dic_flowcell_id.setter
    # def dic_flowcell_id(self, p_arg):
    #     self.__dic_flowcell_id = p_arg

    @property
    def index_kit_info_from_rawdata(self):
        return self.__index_kit_info_form_rawdata

    @index_kit_info_from_rawdata.setter
    def index_kit_info_from_rawdata(self, p_arg):
        self.__index_kit_info_form_rawdata = p_arg

    @property
    def l_nt_sample_list(self):
        return self.__l_nt_sample_list

    @l_nt_sample_list.setter
    def l_nt_sample_list(self, p_arg):
        self.__l_nt_sample_list = p_arg

    @property
    def l_nt_sample_list_in_analysis(self):
        return self.__l_nt_sample_list_in_analysis

    @l_nt_sample_list_in_analysis.setter
    def l_nt_sample_list_in_analysis(self, p_arg):
        self.__l_nt_sample_list_in_analysis = p_arg

    @property
    def l_nt_sample_list_in_analysis_for_fastqc(self):
        return self.__l_nt_sample_list_in_analysis_for_fastqc

    @l_nt_sample_list_in_analysis_for_fastqc.setter
    def l_nt_sample_list_in_analysis_for_fastqc(self, p_arg):
        self.__l_nt_sample_list_in_analysis_for_fastqc = p_arg

    @property
    def l_summary_files(self):
        return self.__l_summary_files

    @l_summary_files.setter
    def l_summary_files(self, p_arg):
        self.__l_summary_files = p_arg

