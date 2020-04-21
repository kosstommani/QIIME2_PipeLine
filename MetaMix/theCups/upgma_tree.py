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

from SpoON.util import run_cmd, check_run_cmd
import os
from click import secho, echo


def get_upgma_tree_cmd(p_dm, p_out_file):
    cmd = \
        'upgma_cluster.py ' \
        '-i {dm} ' \
        '-o {out_tre}'.format(
            dm=p_dm,
            out_tre=p_out_file,
        )
    return cmd


def run_upgma_tree(p_dm, p_out_file):
    cmd = get_upgma_tree_cmd(p_dm, p_out_file)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'UPGMA Tree 완료',
            'false_meg': 'upgma_cluster',
        }, p_exit=False
    )
    if run.returncode == 0:
        return True
    else:
        return False


class TreeController(object):
    """
    FigTree_format_Controller_V2.1.py

    # 작성날짜 : 2016.11.01
    #------------------------------------------------------------------------
    # 수정내역
    # Ver. 1.1
    # 2017.01.06 : 6번(#FF5E00)과 12번(#ED4C00) 구분이 잘 안되어서 12번 제외.
    # -----------------------------------------------------------------------
    # Ver. 1.2
    # 수정일자 : 2017.06.12
    # - mappingfile.txt의 필드명으로 tre파일생성
    # - tre파일을 tre_ori으로 변경하지 않도록 수정
    # -----------------------------------------------------------------------
    # Ver. 2.0
    # 수정일자 : 2017.06.19
    # --grouping_field 옵션에서 그룹의 컬럼 인자를 여러개 입력받도록 변경
    # 그룹의 인자별로 .tre 파일 생성
    # upgma_cluster.html를 파싱하여 그룹별로 png를 출력하게끔 html 수정
    # upgma_cluster.html 출력시 그룹명 및 그룹 구성 이름 출력
    # -----------------------------------------------------------------------
    # Ver. 2.1
    # 수정일자 : 2017.06.29
    # mappingfile의 그룹명으로 데이터 정렬
    # HTML 페이지에서 그룹명 표시할시 그룹색으로 글자색 지정
    # -----------------------------------------------------------------------
    # 수정일자 : 2017.08.24
    # mappingfile 파싱 구분자 변경 : " " --> "\t"
    # -----------------------------------------------------------------------
    # 수정일자 : 2017.11.22
    # 색상 3개 추가
    # -----------------------------------------------------------------------
    # 수정일자 : 2018.06.15
    # format 설정 변경
    #  - Scale Bar 선택 해제
    #  - Scale Bar 선택 설정
    #  - Scale Bar : Show grid 선택 해제
    # -----------------------------------------------------------------------
    # 수정일자 : 2019.08.19
    # 기존 프로그램을 MetaMix 에 통합 및 변경

    """
    def __init__(self, p_metadata, p_tre, p_grouping_field):
        self.metadata = p_metadata
        self.tre = p_tre
        self.grouping_field = p_grouping_field
        self.group_data = None
        self.new_tre_file = None
        self.new_tre_name = None
        self.sample_group_dic = {}
        self.group_color_dic = {}
        self.color_list = ["#000000", "#FF0000", "#0100FF", "#1DDB16",
                           "#FF5E00", "#00D8FF", "#FF00DD", "#5F00FF", "#FFBB00",
                           "#8C8C8C",            "#F29661", "#997000", "#D5D5D5",
                           "#3DB7CC", "#FFB2F5", "#980000", "#86E57F", "#A3A2FF",
                           "#0000A5", "#DBC000", "#009300", "#FFA7A7", "#1993A8",
                           "#5E5E5E", "#FF007F", "#ABF200", "#A5587F", "#FFC6C6",
                           "#8A14B2", "#79B200", "#36A5B2", "#336633"]
        # yellow : #FFE400

    def read_metadata_file(self):
        echo(f'>>> Grouping Field : {self.grouping_field}')
        with open(self.metadata, 'r') as input_data:
            temp = []
            for i in input_data:
                if i.startswith('#'):
                    i_split = i.strip().split('\t')
                    field_index = i_split.index(self.grouping_field)
                else:
                    i_split = i.strip().split('\t')
                    temp_group = i_split[field_index]  # grouping field
                    self.sample_group_dic[i_split[0]] = temp_group
                    temp.append(temp_group)

        temp_set = set(temp)  # grouping field 의 중복 제거
        if len(temp_set) > len(self.color_list):
            secho(f'Error : 그룹의 개수가 {len(self.color_list)}개 이상입니다. 색상 추가가 필요합니다.')
            echo('프로그램 작성자에게 문의하십시오.')
            exit()
        self.group_data = list(temp_set)
        self.group_data.sort()
        for i in range(len(self.group_data)):
            self.group_color_dic[self.group_data[i]] = self.color_list[i]
        echo(f'>>> 샘플의 개수 : {len(self.sample_group_dic.keys())}')
        echo(f'>>> 그룹의 개수 : {len(temp_set)}')
        echo()

    def read_metadata_file_no_group(self):
        echo('>>> Grouping Field : No_Group')
        with open(self.metadata, 'r') as input_data:
            for i in input_data:
                if i.startswith('#'):
                    continue
                else:
                    i_split = i.strip().split('\t')
                    self.sample_group_dic[i_split[0]] = 'No_Group'

        self.group_data = ['No_Group']
        self.group_color_dic[self.group_data[0]] = self.color_list[0]

    def read_tree_file(self):
        with open(self.tre, 'r') as input_data:
            tre_input_data = input_data.readline()
            
        # tre파일명에 mappingfile.txt의 필드명 추가
        self.new_tre_file = self.tre.replace('.tre', f'.{self.grouping_field}.tre')
        self.new_tre_name = os.path.basename(self.new_tre_file)
        
        with open(self.new_tre_file, 'w') as output:
            output_text_1 = """\
#NEXUS
begin taxa;
    dimensions ntax={0};
    taxlabels
""".format(len(self.sample_group_dic.keys()))

            output_text_2 = """\
;
end;

begin trees;
    tree tree_1 = [&R] {0}
end;

begin figtree;
    set appearance.backgroundColorAttribute="Default";
    set appearance.backgroundColour=#ffffff;
    set appearance.branchColorAttribute="User selection";
    set appearance.branchColorGradient=false;
    set appearance.branchLineWidth=3.0;
    set appearance.branchMinLineWidth=0.0;
    set appearance.branchWidthAttribute="Fixed";
    set appearance.foregroundColour=#000000;
    set appearance.hilightingGradient=false;
    set appearance.selectionColour=#2d3680;
    set branchLabels.colorAttribute="User selection";
    set branchLabels.displayAttribute="Branch times";
    set branchLabels.fontName="Arial";
    set branchLabels.fontSize=8;
    set branchLabels.fontStyle=0;
    set branchLabels.isShown=false;
    set branchLabels.significantDigits=4;
    set layout.expansion=0;
    set layout.layoutType="RECTILINEAR";
    set layout.zoom=0;
    set legend.attribute=null;
    set legend.fontSize=10.0;
    set legend.isShown=false;
    set legend.significantDigits=4;
    set nodeBars.barWidth=4.0;
    set nodeBars.displayAttribute=null;
    set nodeBars.isShown=false;
    set nodeLabels.colorAttribute="User selection";
    set nodeLabels.displayAttribute="Node ages";
    set nodeLabels.fontName="Arial";
    set nodeLabels.fontSize=12;
    set nodeLabels.fontStyle=0;
    set nodeLabels.isShown=true;
    set nodeLabels.significantDigits=4;
    set nodeShape.colourAttribute="User selection";
    set nodeShape.isShown=false;
    set nodeShape.minSize=10.0;
    set nodeShape.scaleType=Width;
    set nodeShape.shapeType=Circle;
    set nodeShape.size=4.0;
    set nodeShape.sizeAttribute="Fixed";
    set polarLayout.alignTipLabels=false;
    set polarLayout.angularRange=0;
    set polarLayout.rootAngle=0;
    set polarLayout.rootLength=100;
    set polarLayout.showRoot=true;
    set radialLayout.spread=0.0;
    set rectilinearLayout.alignTipLabels=false;
    set rectilinearLayout.curvature=0;
    set rectilinearLayout.rootLength=100;
    set scale.offsetAge=0.0;
    set scale.rootAge=1.0;
    set scale.scaleFactor=1.0;
    set scale.scaleRoot=false;
    set scaleAxis.automaticScale=true;
    set scaleAxis.fontSize=8.0;
    set scaleAxis.isShown=true;
    set scaleAxis.lineWidth=1.0;
    set scaleAxis.majorTicks=1.0;
    set scaleAxis.origin=0.0;
    set scaleAxis.reverseAxis=false;
    set scaleAxis.showGrid=false;
    set scaleBar.automaticScale=true;
    set scaleBar.fontSize=10.0;
    set scaleBar.isShown=false;
    set scaleBar.lineWidth=1.0;
    set scaleBar.scaleRange=0.0;
    set tipLabels.colorAttribute="User selection";
    set tipLabels.displayAttribute="Names";
    set tipLabels.fontName="Arial";
    set tipLabels.fontSize=16;
    set tipLabels.fontStyle=1;
    set tipLabels.isShown=true;
    set tipLabels.significantDigits=4;
    set trees.order=false;
    set trees.orderType="increasing";
    set trees.rooting=false;
    set trees.rootingType="User Selection";
    set trees.transform=false;
    set trees.transformType="cladogram";
end;""".format(tre_input_data)

            output.write(output_text_1)
            temp_sample = list(self.sample_group_dic.keys())
            temp_sample.sort()
            for i in temp_sample:  # i <-- sample name
                sample_name = i
                color = self.group_color_dic[self.sample_group_dic[sample_name]]
                output.write(f'\t{sample_name}[&!color={color}]\n')
            output.write(output_text_2)
