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

from click import echo, secho, progressbar, style
from collections import defaultdict, namedtuple
from SpoON.util import read_yaml, make_dir_using_input_file, run_cmd, check_run_cmd, read_data
from time import sleep
import os


BiomTable_L2 = namedtuple('BiomTable_L2', 'L1 L2 ratio')
BiomTable_L5 = namedtuple('BiomTable_L4', 'L1 L2 L3 L4 L5 ratio')
BiomTable_L6 = namedtuple('BiomTable_L6', 'L1 L2 L3 L4 L5 L6 ratio')
BiomTable_L7 = namedtuple('BiomTable_L7', 'L1 L2 L3 L4 L5 L6 L7 ratio')


class Score:
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        self.d_all_bacteria = p_bacteria
        self.d_all_dist_db = p_dist_db
        self.d_genus_table = p_genus
        self.d_species_table = p_species
        self.main_item = None
        self.save_detail_file = False
        self.o_score_detail_file = None

    @property
    def d_bacteria(self):
        """
        항목별로 검출되어야할 미생물 목록
        :return:
        """
        return self.d_all_bacteria[self.main_item]

    @property
    def d_dist_db(self):
        """
        한국인 장내미생물 분포 데이터
        :return: 
        """
        return self.d_all_dist_db[self.main_item]

    def set_save_score_detail_file(self, p_file):
        self.save_detail_file = True
        self.o_score_detail_file = open(p_file, 'w')
        secho('>>> score_detail 파일 저장 설정 완료', fg='cyan')

    def make_effect_rank_dict_from_db(self, item: str) -> dict:
        """
        Bacteria DB의 index에 해당되는 effect 및 rank를 추출하고, table에서 해당 taxon의 정보를 추출한다.
        추출된 정보는 딕션너리로 반환한다.

        main items:
         - 장 건강         intestinal_health
         - 웰니스          wellness
         - 장 환경 더 보기  more_bowel_environments
         - 부록            supplement

        :param item: items_and_bacteria.yaml 파일의 main items 에 있는 sub items의 이름  
        :return: d_effect_rank
        """
        d_effect_rank = defaultdict(dict)
        for effect in self.d_bacteria[item].keys():
            for rank in self.d_bacteria[item][effect].keys():
                d_rank = defaultdict(dict)
                for taxon in self.d_bacteria[item][effect][rank]:
                    d_rank[taxon] = self.extract_taxon(rank, taxon)
                d_effect_rank[effect][rank] = d_rank
        return d_effect_rank

    def extract_taxon(self, rank: str, taxon: str) -> namedtuple or None:
        """
        Summarized_taxa 에 의해서 생성된 각 레벨별 abundace 결과에서 items_and_bacteria.yaml 파일에 기재된
        각 항목에 일치하는 미생물이 있는 경우 해당 결과값을 반환한다.
        일치하는 항목이 없는 경우 None을 반환한다.
        
        :param rank: Taxon Level
        :param taxon: 존재여부 확인 및 추출하고자하는 Taxon
        :return:
        """
        if rank == 'phylum':
            table = self.d_phylum_table
        elif rank == 'family':
            table = self.d_family_table
        elif rank == 'genus':
            table = self.d_genus_table
        elif rank == 'species':
            table = self.d_species_table
        else:
            raise ValueError(f'rank 값이 올바르지 않습니다. rank: {rank}')
        return table.get(taxon)

    @staticmethod
    def _transform(p_type: str, p_min: float, p_max: float, p_value: float) -> float:
        """
        최소값, 최대값을 이용하여 1~99 범위로 데이터의 값을 표준화한다.

        :param p_type: [positive | negative] 해당 항목의 영향
        :param p_min: 해당 항목의 한국인 장내미생물 분포에 대한 최소값
        :param p_max: 해당 항목의 한국인 장내미생물 분포에 대한 최대값
        :param p_value: 해당 항목의 시료의 Abundance 
        :return: 1-99 표준화 수치
        """
        score = (p_value - p_min) / (p_max - p_min) * 100
        if p_type == 'positive':
            if p_value <= p_min:
                return 1
            elif p_value >= p_max:
                return 99
            else:
                if score < 1:
                    return 1
                elif score > 99:
                    return 99
                else:
                    return score
        elif p_type == 'negative':
            if p_value <= p_min:
                return 99
            elif p_value >= p_max:
                return 1
            else:
                if score < 1:
                    return 99
                elif score > 99:
                    return 1
                else:
                    return 99 - score
        else:
            raise ValueError(f'p_type 옵션값이 올바르지 않습니다. p_type: {p_type}')

    def transform_data(self, p_item: str, p_d_effect_rank: dict) -> dict:
        """
        분포 레퍼런스의 최소값, 최대값을 이용하여 1~99 범위로 데이터의 값을 표준화한다.

        :param p_item: 2단계 분류 아이템 이름. ex) constipation, gas
        :param p_d_effect_rank: Effect[positive | negative] 항목을 Key로 가지는 딕션너리.
                                Value는 각 Effect에 해당되는 미생물[Phylum | Family | Gensus | Species]
        :return: 1-99 표준화 점수 데이터를 가지는 딕션너리
                 Key: Effect[positive | negative]
                 Value: 1-99 표준화 점수
        """
        d_sum = defaultdict(dict)
        for effect in p_d_effect_rank.keys():
            d_rank = defaultdict()
            for rank in p_d_effect_rank[effect].keys():
                rank_sum = sum(
                    [value.ratio for value in p_d_effect_rank[effect][rank].values() if value is not None])
                d_rank[rank] = rank_sum
            d_sum[effect] = d_rank

        d_transformed_data = defaultdict(float)
        for effect in d_sum.keys():
            sum_values = sum(d_sum[effect].values())
            db_min = self.d_dist_db[p_item][effect]['min']
            db_max = self.d_dist_db[p_item][effect]['max']
            d_transformed_data[effect] = self._transform(effect, db_min, db_max, sum_values)
            # 변환에 필요한 데이터 출력
            if self.save_detail_file is True:
                echo(file=self.o_score_detail_file)
                echo(f'item: {p_item}', file=self.o_score_detail_file)
                echo(f'effect: {effect}', file=self.o_score_detail_file)
                echo(f'db_min: {db_min}', file=self.o_score_detail_file)
                echo(f'db_max: {db_max}', file=self.o_score_detail_file)
                echo(f'sum_values: {sum_values}', file=self.o_score_detail_file)
                echo('='*20, file=self.o_score_detail_file)
                echo(f'trans_values: {self._transform(effect, db_min, db_max, sum_values)}', file=self.o_score_detail_file)
                echo(file=self.o_score_detail_file)
        return d_transformed_data

    def transform_data_one_effect(self, p_item: str, p_d_effect_rank: dict) -> float:
        d_transformed_data = self.transform_data(p_item, p_d_effect_rank)
        for index, effect in enumerate(p_d_effect_rank, 1):
            value = d_transformed_data[effect]
        if index > 1:
            raise RuntimeError(f'effect가 1개가 아닙니다.\n p_d_effect_rank: {p_d_effect_rank}')
        return value

    def transform_data_two_effect(self, p_item: str, p_d_effect_rank: dict) -> float:
        """
        장 건강, 웰니스, 부록(영양소 대사)
        positive 와 negative 가 모두 있는 경우에 대해서 점수 변환을 진행한다.
        
        :param p_item:
        :param p_d_effect_rank:
        :return: positive 와 negative
        """
        d_transformed_data = self.transform_data(p_item, p_d_effect_rank)
        return sum(d_transformed_data.values()) * 0.5


class IntestinalHealth(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        """
        장 건강 항목에 대한 점수를 계산한다.
        장 건강 항목에는 다음과 같은 세부 항목이 있다.
         - 변비
         - 방귀
         - 복부팽만감
         - 신경성 복부 불편감
         - 설사
         
        :param p_bacteria: 항목별로 검출되어야할 미생물 목록
        :param p_dist_db: 한국인 장내미생물 분포 데이터
        :param p_genus: 분석시료의 Genus Level의 Abundance 데이터
        :param p_species: 분석시료의 Species Level의 Abundance 데이터
        """
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'intestinal_health'

        # property
        self.__constipation = None
        self.__gas = None
        self.__abdominal_distension = None
        self.__neurotic_abdominal_discomfort = None
        self.__diarrhea = None

    @property
    def constipation(self):
        """
        한국인 분포값에 따라 변환된 변비 항목에 대한 점수
        :return: 
        """
        return self.__constipation

    @property
    def gas(self):
        """
        한국인 분포값에 따라 변환된 방귀 항목에 대한 점수
        :return: 
        """
        return self.__gas

    @property
    def abdominal_distension(self):
        """
        한국인 분포값에 따라 변환된 복부팽만감 항목에 대한 점수
        :return:
        """
        return self.__abdominal_distension

    @property
    def IBS(self):
        """
        한국인 분포값에 따라 변환된 신경성 복부 물편감 항목에 대한 점수
        :return:
        """
        return self.__neurotic_abdominal_discomfort

    @property
    def diarrhea(self):
        """
        한국인 분포값에 따라 변환된 설사 항목에 대한 점수
        :return:
        """
        return self.__diarrhea

    def compute_constipation(self):
        item_name = 'constipation'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__constipation = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_gas(self):
        """
        only negative
        :return:
        """
        item_name = 'gas'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__gas = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_abdominal_distension(self):
        item_name = 'abdominal_distension'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__abdominal_distension = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_ibs(self):
        item_name = 'neurotic_abdominal_discomfort'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__neurotic_abdominal_discomfort = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_diarrhea(self):
        item_name = 'diarrhea'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__diarrhea = self.transform_data_two_effect(item_name, d_effect_rank)


class Wellness(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict,
                 p_phylum: dict, p_family: dict,
                 p_genus: dict, p_species: dict):
        """
        웰니스 항목에 대한 점수를 계산한다.
        웰니스 항목에는 다음과 같은 세부 항목이 있다.
         - 행복 지수
         - 피로 지수
         - 면역 지수
         - 비만 지수
         - 수면 지수
         - 노화 지수
         
        :param p_bacteria: 항목별로 검출되어야할 미생물 목록
        :param p_dist_db: 한국인 장내미생물 분포 데이터
        :param p_phylum: 분석시료의 Phylum Level의 Abundance 데이터
        :param p_family: 분석시료의 Family Level의 Abundance 데이터
        :param p_genus: 분석시료의 Genus Level의 Abundance 데이터
        :param p_species: 분석시료의 Species Level의 Abundance 데이터
        """
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'wellness'
        self.d_phylum_table = p_phylum
        self.d_family_table = p_family

        # property
        self.__happiness_index = None
        self.__tiredness_index = None
        self.__immunity_index = None
        self.__obesity_index = None
        self.__sleep_index = None
        self.__aging_index = None

    @property
    def behavior(self):
        """
        한국인 분포값에 따라 변환된 행복 지수에 대한 점수
        :return:
        """
        return self.__happiness_index

    @property
    def fatigue(self):
        """
        한국인 분포값에 따라 변환된 피로 지수에 대한 점수
        :return:
        """
        return self.__tiredness_index

    @property
    def immunity(self):
        """
        한국인 분포값에 따라 변환된 면역 지수에 대한 점수
        :return: 
        """
        return self.__immunity_index

    @property
    def obesity(self):
        """
        한국인 분포값에 따라 변환된 비만 지수에 대한 점수
        :return:
        """
        return self.__obesity_index

    @property
    def sleep(self):
        """
        한국인 분포값에 따라 변환된 수면 지수에 대한 점수
        :return:
        """
        return self.__sleep_index

    @property
    def aging(self):
        """
        한국인 분포값에 따라 변환된 노화 지수에 대한 점수
        :return:
        """
        return self.__aging_index

    def compute_happiness_index(self):
        """
        only positive
        :return:
        """
        item_name = 'happiness_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__happiness_index = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_tiredness_index(self):
        item_name = 'tiredness_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__tiredness_index = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_immune_index(self):
        item_name = 'immunity_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__immunity_index = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_obesity_index(self):
        item_name = 'obesity_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__obesity_index = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_sleep_index(self):
        item_name = 'sleep_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__sleep_index = self.transform_data_two_effect(item_name, d_effect_rank)

    def compute_aging_index(self):
        item_name = 'aging_index'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__aging_index = self.transform_data_two_effect(item_name, d_effect_rank)


class MoreBowelEnvironments(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict,
                 p_score_path: str, p_sample_name: str):
        """
        장 환경 더 보기 항목에 대한 점수를 계산한다.
        장 환경 더 보기 항목에는 다음과 같은 세부 항목이 있다.
         - 유익균
         - 유해균
         - 장 유형
        장 유형의 경우, 한국인 장내미생물 분포 데이터를 이용하여 군집분석(K-means Clustering)을 통해 세 개 군집으로 분류하였다.
        분류된 군집 데이터는 머싱 러닝 기법인 Random Forest을 이용하여 학습하였다.
        학습된 정보를 바탕으로 분석대상 시료의 Enterotype을 분류한다.
        
        :param p_bacteria: 항목별로 검출되어야할 미생물 목록
        :param p_dist_db: 한국인 장내미생물 분포 데이터
        :param p_genus: 분석시료의 Genus Level의 Abundance 데이터
        :param p_species: 분석시료의 Species Level의 Abundance 데이터
        :param p_score_path: Score 디렉터리 경로
        :param p_sample_name: 시료명
        """
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'more_bowel_environments'
        self.score_path = p_score_path
        self.sample_name = p_sample_name
        self.__beneficial_bacteria = None
        self.__harmful_bacteria = None
        self.__enterotype = None
        self.__d_enterotype_value = None

    @property
    def benefit(self):
        """
        한국인 분포값에 따라 변환된 유익균에 대한 점수
        :return:
        """
        return self.__beneficial_bacteria

    @property
    def harmful(self):
        """
        한국인 분포값에 따라 변환된 유해균에 대한 점수
        :return:
        """
        return self.__harmful_bacteria

    @property
    def enterotype(self):
        """
        한국인 분포값에 따라 예측된 장환경(유형)
        :return:
        """
        return self.__enterotype

    @property
    def d_enterotype_value(self):
        """
        장환경(유형)예측에 필요한 3개 Genus의 비율
        :return:
        """
        return self.__d_enterotype_value

    def compute_beneficial_bacteria(self):
        item_name = 'beneficial_bacteria'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__beneficial_bacteria = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_harmful_bacteria(self):
        item_name = 'harmful_bacteria'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__harmful_bacteria = self.transform_data_one_effect(item_name, d_effect_rank)

    def predict_enterotype(self):
        item_name = 'enterotype'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        enterotype_name = 'enterotype_genus.txt'
        enterotype_score_file = os.path.join(self.score_path, enterotype_name)
        predicted_enterotype_file = os.path.join(self.score_path, "enterotype.predicted")
        self.__d_enterotype_value = dict()
        with open(enterotype_score_file, 'w') as o_enterotype:
            o_enterotype.write(f'Genus\t{self.sample_name}\n')
            for genus, nt_values in d_effect_rank['type']['genus'].items():
                if nt_values is None:
                    o_enterotype.write(f'{genus}\t0\n')
                    self.__d_enterotype_value[genus] = 0
                else:
                    o_enterotype.write(f'{genus}\t{nt_values.ratio}\n')
                    self.__d_enterotype_value[genus] = nt_values.ratio
        cmd = '/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/R_env_bin/Rscript ' \
              '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/MicrobeAndMe/enterotype.R ' \
              f'{enterotype_score_file} ' \
              f'{predicted_enterotype_file}'
        run = run_cmd(cmd)
        check_run_cmd(
            {
                'run': run,
                'true_meg': None,
                'false_meg': 'enterotype.R'
            }
        )
        l_predicted_enterotype = read_data(predicted_enterotype_file, ',', False)
        # R - 시료명이 숫자일 경우 컬럼명 설정시 X가 붙음.
        if l_predicted_enterotype[0][0].startswith('X'):
            predicted_sample_name = l_predicted_enterotype[0][0].lstrip('X')
            try:
                float(predicted_sample_name)
            except ValueError:
                predicted_sample_name = l_predicted_enterotype[0][0]
        else:
            predicted_sample_name = l_predicted_enterotype[0][0]
        if predicted_sample_name == self.sample_name:
            self.__enterotype = l_predicted_enterotype[0][1]
        else:
            secho('Error: enterotype.predicted 파일에 기재된 시료명이 분석대상 시료명과 다릅니다.',
                  fg='red', blink=True, err=True)
            echo(f'분석시료명: {self.sample_name}')
            echo(f'enterotype.predicted 파일: {l_predicted_enterotype[0][0]}')
            exit(1)


class Supplement(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'supplement'
        # property
        self.__dietary_fibre = None
        self.__lactose = None
        self.__resistant_starch = None
        self.__protein = None
        self.__essential_amino_acid = None
        self.__monoenoic_fatty_acid = None
        self.__b1 = None
        self.__b2 = None
        self.__b3 = None
        self.__b6 = None
        self.__b7 = None
        self.__b9 = None
        self.__b12 = None
        self.__k2 = None

    @property
    def fiber(self):
        """
        한국인 분포값에 따라 변환된 식이섬유에 대한 점수
        :return:
        """
        return self.__dietary_fibre

    @property
    def lactose(self):
        """
        한국인 분포값에 따라 변환된 유당에 대한 점수
        :return: 
        """
        return self.__lactose

    @property
    def starch(self):
        """
        한국인 분포값에 따라 변환된 저항성전분에 대한 점수
        :return:
        """
        return self.__resistant_starch

    @property
    def protein(self):
        return self.__protein

    @property
    def bcaa(self):
        """
        한국인 분포값에 따라 변환된 필수아미노산에 대한 점수
        :return:
        """
        return self.__essential_amino_acid

    @property
    def scfa(self):
        """
        한국인 분포값에 따라 변환된 단쇄지방산에 대한 점수
        :return: 
        """
        return self.__monoenoic_fatty_acid

    @property
    def vitaminb1(self):
        """
        한국인 분포값에 따라 변환된 비타민B1에 대한 점수
        :return:
        """
        return self.__b1

    @property
    def vitaminb2(self):
        """
         한국인 분포값에 따라 변환된 비타민B2에 대한 점수
        :return:
        """
        return self.__b2

    @property
    def vitaminb3(self):
        """
         한국인 분포값에 따라 변환된 비타민B3에 대한 점수
        :return:
        """
        return self.__b3

    @property
    def vitaminb6(self):
        """
         한국인 분포값에 따라 변환된 비타민B6에 대한 점수
        :return:
        """
        return self.__b6

    @property
    def vitaminb7(self):
        """
         한국인 분포값에 따라 변환된 비타민B7에 대한 점수
        :return:
        """
        return self.__b7

    @property
    def vitaminb9(self):
        """
         한국인 분포값에 따라 변환된 비타민B9에 대한 점수
        :return:
        """
        return self.__b9

    @property
    def vitaminb12(self):
        """
         한국인 분포값에 따라 변환된 비타민B12에 대한 점수
        :return:
        """
        return self.__b12

    @property
    def vitamink2(self):
        """
         한국인 분포값에 따라 변환된 비타민K2에 대한 점수
        :return:
        """
        return self.__k2

    def compute_dietary_fibre(self):
        item_name = 'dietary_fibre'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__dietary_fibre = self.transform_data_two_effect(item_name, d_effect_rank)
        pass

    def compute_lactose(self):
        """
        only positive
        :return:
        """
        item_name = 'lactose'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__lactose = self.transform_data_one_effect(item_name, d_effect_rank)
    
    def compute_resistant_starch(self):
        """
        only positive
        :return:
        """
        item_name = 'resistant_starch'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__resistant_starch = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_protein(self):
        """
        only positive
        :return:
        """
        item_name = 'protein'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__protein = self.transform_data_one_effect(item_name, d_effect_rank)
    
    def compute_essential_amino_acid(self):
        """
        only positive
        :return:
        """
        item_name = 'essential_amino_acid'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__essential_amino_acid = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_monoenoic_fatty_acid(self):
        """
        only positive
        :return:
        """
        item_name = 'monoenoic_fatty_acid'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__monoenoic_fatty_acid = self.transform_data_one_effect(item_name, d_effect_rank)
    
    def compute_b1(self):
        """
        only positive
        :return:
        """
        item_name = 'B1'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b1 = self.transform_data_one_effect(item_name, d_effect_rank)
    
    def compute_b2(self):
        """
        only positive
        :return:
        """
        item_name = 'B2'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b2 = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_b3(self):
        """
        only positive
        :return:
        """
        item_name = 'B3'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b3 = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_b6(self):
        """
        only positive
        :return:
        """
        item_name = 'B6'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b6 = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_b7(self):
        """
        only positive
        :return:
        """
        item_name = 'B7'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b7 = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_b9(self):
        """
        only positive
        :return:
        """
        item_name = 'B9'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b9 = self.transform_data_one_effect(item_name, d_effect_rank)
    
    def compute_b12(self):
        """
        only positive
        :return:
        """
        item_name = 'B12'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__b12 = self.transform_data_one_effect(item_name, d_effect_rank)

    def compute_k2(self):
        """
        only positive
        :return:
        """
        item_name = 'K2'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__k2 = self.transform_data_one_effect(item_name, d_effect_rank)


class SupplementProbiotics21(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_species: dict):
        """
        
        :param p_bacteria:
        :param p_dist_db:
        :param p_species: Probiotics19 DB를 사용한 Taxonomy Assignment 결과
        """
        super().__init__(p_bacteria, p_dist_db, None, p_species)
        self.main_item = 'supplement'
        # property
        self.__probiotics_21 = None

    @property
    def probiotics_21(self):
        return self.__probiotics_21

    def compute_probiotics_21(self):
        """
        only positive
        :return:
        """
        item_name = 'probiotics_21'
        d_effect_rank = self.make_effect_rank_dict_from_db(item_name)
        self.__probiotics_21 = d_effect_rank


class DiversityIndex(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_diversity_file: str):
        super().__init__(p_bacteria, p_dist_db, None, None)
        self.main_item = 'diversity_index'
        self.diversity_file = p_diversity_file
        self.__ori_shannon = None
        self.__ori_simpson = None
        self.__shannon = None
        self.__simpson = None

    @property
    def ori_shannon(self):
        return self.__ori_shannon

    @property
    def ori_simpson(self):
        return self.__ori_simpson

    @property
    def shannon(self):
        """
        표준화된 Shannon Index
        :return:
        """
        return self.__shannon

    @property
    def simpson(self):
        """
        계산되지 않음. 
        :return:
        """
        return self.__simpson

    def _read_alpha_diversity(self):
        """
        Alpha_Diversity.txt 파일을 읽고 데이터를 반환한다.

        :return: namedtuple - AlphaDiversity.shannon
                                            .simpson
        """
        data = None
        header = None
        with open(self.diversity_file, 'r') as o_text:
            for line, text in enumerate(o_text):
                if line == 0:
                    header = text.strip().split('\t')
                elif line == 1:
                    if text.strip():
                        data = text.strip().split()[1:]  # Sample ID 제거
                    else:
                        secho('Error: 데이터가 없습니다.', fg='red')
                        secho(f'---> {self.diversity_file}를 확인하세요.', fg='magenta')
                        exit(1)
                else:
                    if text.strip():
                        secho('Error: 데이터가 여러개 존재합니다. 시료 1개에 대한 데이터만 있어야 합니다.', fg='red')
                        secho(f'---> {self.diversity_file}를 확인하세요.(데이터 개수 확인, 개행 여부 확인)',
                              fg='magenta')
                        exit(1)
        shannon = float(data[header.index('shannon')])
        simpson = float(data[header.index('simpson')])
        self.__ori_shannon = shannon
        self.__ori_simpson = simpson
        # secho('>>> Alpha Diversity 파일 읽기 완료', fg='cyan')
        # echo(self.__diversity_file)

    def transform_data(self):
        """
        shannon index를 표준화시킨다.
        :return:
        """
        db_min = self.d_dist_db['shannon']['min']
        db_max = self.d_dist_db['shannon']['max']
        transformed_data = self._transform('positive', db_min, db_max, self.__ori_shannon)
        if self.save_detail_file is True:
            echo(file=self.o_score_detail_file)
            echo(f'item: {self.main_item}', file=self.o_score_detail_file)
            echo(f'effect: positive', file=self.o_score_detail_file)
            echo(f'db_min: {db_min}', file=self.o_score_detail_file)
            echo(f'db_max: {db_max}', file=self.o_score_detail_file)
            echo(f'shannon: {self.__ori_shannon}', file=self.o_score_detail_file)
            echo('='*20, file=self.o_score_detail_file)
            echo(f'trans_values: {transformed_data}', file=self.o_score_detail_file)
            echo(file=self.o_score_detail_file)
        return transformed_data

    def compute_shannon(self):
        self._read_alpha_diversity()
        self.__shannon = self.transform_data()
        

class AndMe:
    def __init__(self, p_kargs):
        self.bact_file = p_kargs['bact_file']
        self.dist_db_file = p_kargs['dist_db_file']
        self.tax_id_file = p_kargs['tax_id_file']
        self.phylum_file = p_kargs['phylum_file']
        self.family_file = p_kargs['family_file']
        self.genus_file = p_kargs['genus_file']
        self.species_file = p_kargs['species_file']
        self.probiotics19_species_file = p_kargs['probiotics19_species_file']
        self.alpha_diversity_file = p_kargs['alpha_diversity_file']
        self.sample_name = p_kargs['sample_name']
        self.out_path = p_kargs['out_path']
        self.__d_bacteria = None
        self.__d_dist_db = None
        self.__d_tax_id = None
        self.__d_phylum_table = None
        self.__d_family_table = None
        self.__d_genus_table = None
        self.__d_species_table = None
        self.__d_probiotics21_species_table = None
        self.__i_intestinal_health = None
        self.__i_wellness = None
        self.__i_more_bowel_environments = None
        self.__i_supplement = None
        self.__i_probiotics21 = None
        self.__i_diversity_index = None
        self.__score_path = None

    @property
    def d_bacteria(self):
        return self.__d_bacteria

    @property
    def d_dist_db(self):
        return self.__d_dist_db

    @property
    def d_tax_id(self):
        return self.__d_tax_id

    @property
    def report_version(self):
        item_report_version = self.d_bacteria['Items_and_Bacteria']['report_version']
        dist_ref_report_version = self.d_dist_db['MicrobeAndMe_Ref']['report_version']
        if item_report_version == dist_ref_report_version:
            return item_report_version
        else:
            secho('Error: items_and_bacteria.yaml 과 MicrobeAndMe_Ref_V2.1.yaml 에 기재된 report_version의 정보가'
                  '다릅니다.', fg='red', blink=True, err=True)
            echo(f'items_and_bacteria.yaml: {item_report_version}', err=True)
            echo(f'MicrobeAndMe_Ref_V2.1.yaml: {dist_ref_report_version}', err=True)
            exit()

    @property
    def d_phylum_table(self):
        return self.__d_phylum_table

    @property
    def d_family_table(self):
        return self.__d_family_table

    @property
    def d_genus_table(self):
        return self.__d_genus_table

    @property
    def d_species_table(self):
        return self.__d_species_table

    @property
    def d_probiotics19_species_table(self):
        return self.__d_probiotics21_species_table

    @property
    def i_intestinal_health(self):
        return self.__i_intestinal_health

    @property
    def i_wellness(self):
        return self.__i_wellness

    @property
    def i_more_bowel_environments(self):
        return self.__i_more_bowel_environments

    @property
    def i_supplement(self):
        return self.__i_supplement

    @property
    def i_probiotics21(self):
        return self.__i_probiotics21

    @property
    def i_diversity_index(self):
        return self.__i_diversity_index

    @property
    def score_path(self):
        return self.__score_path

    def make_score_dir(self):
        self.__score_path = make_dir_using_input_file(self.out_path, 'Score', 0, False)

    def _read_table(self, p_l: int):
        """
        otu_table_L?.txt 파일을 읽어 리스트에 저장한다.
        :param p_l: 읽어들일 otu table의 레벨.[2|4|6|7|19]
                    2: Phylum
                    4: Family
                    6: Genus
                    7: Species
                    19: Probiotics19 DB을 이용한 Species Level의  otu table.
        :return: data
        """
        d_data = defaultdict()
        sample_indicator = False
        if p_l == 2:
            table_file = self.phylum_file
        elif p_l == 5:
            table_file = self.family_file
        elif p_l == 6:
            table_file = self.genus_file
        elif p_l == 7:
            table_file = self.species_file
        elif p_l == 19:
            table_file = self.probiotics19_species_file
        else:
            raise ValueError(f'지원하지 않는 레벨입니다. p_l: {p_l}')
        with open(table_file, 'r') as o_table:
            for ele in o_table:
                if ele.startswith('#'):
                    if self.sample_name in ele:
                        sample_indicator = True
                    continue
                else:
                    l_ele = ele.strip().split('\t')
                    rank = [x.strip() for x in l_ele[0].split(';__')]
                    ratio = float(l_ele[1].strip())
                    if p_l == 2:
                        d_data[rank[1]] = BiomTable_L2(L1=rank[0], L2=rank[1], ratio=ratio)
                    elif p_l == 5:
                        d_data[rank[3]] = BiomTable_L5(L1=rank[0], L2=rank[1], L3=rank[3],
                                                       L4=rank[3], L5=rank[4], ratio=ratio)
                    elif p_l == 6:
                        d_data[rank[5]] = BiomTable_L6(L1=rank[0], L2=rank[1], L3=rank[2],
                                                       L4=rank[3], L5=rank[4], L6=rank[5],
                                                       ratio=ratio)
                    elif (p_l == 7) or (p_l == 19):
                        d_data[rank[6]] = BiomTable_L7(L1=rank[0], L2=rank[1], L3=rank[2],
                                                       L4=rank[3], L5=rank[4], L6=rank[5],
                                                       L7=rank[6], ratio=ratio)
                    else:
                        raise ValueError(f'지원하지 않는 레벨입니다. p_l: {p_l}')
        if sample_indicator is False:
            secho('Error: OTU Table 파일내에 존재하는 시료명과 결과 시료명이 일치하지 않습니다.', fg='red', err=True)
            secho('----> OTU Table 파일내에 기재된 시료명과 디렉터리명을 확인하세요.', fg='magenta', err=True)
            echo(self.__table_file, err=True)
            exit(1)
        # secho('>>> BIOM Table 파일 읽기 완료', fg='cyan')
        # echo(self.__table_file)
        return d_data

    def read_bacteria_file(self):
        self.__d_bacteria = read_yaml(self.bact_file)
        secho('>>> Items & Bacteria 읽기 완료', fg='cyan')
        echo(self.bact_file)

    def read_dist_db_file(self):
        self.__d_dist_db = read_yaml(self.dist_db_file)
        secho('>>> 분포 DB 읽기 완료', fg='cyan')
        echo(self.dist_db_file)

    def read_tax_id_file(self):
        self.__d_tax_id = read_yaml(self.tax_id_file)
        secho('>>> Tax ID 파일 읽기 완료', fg='cyan')
        echo(self.tax_id_file)

    def read_phylum_file(self):
        self.__d_phylum_table = self._read_table(2)

    def read_family_file(self):
        self.__d_family_table = self._read_table(5)

    def read_genus_file(self):
        self.__d_genus_table = self._read_table(6)

    def read_species_file(self):
        self.__d_species_table = self._read_table(7)

    def read_probiotics19_species_file(self):
        self.__d_probiotics21_species_table = self._read_table(19)

    def compute_intestinal_health(self):
        self.__i_intestinal_health = IntestinalHealth(self.d_bacteria, self.d_dist_db,
                                                      self.d_genus_table, self.d_species_table)
        score_detail_file = os.path.join(self.score_path, f'{self.__i_intestinal_health.main_item}.detail')
        self.__i_intestinal_health.set_save_score_detail_file(score_detail_file)
        l_function = [self.i_intestinal_health.compute_constipation, self.i_intestinal_health.compute_gas,
                      self.i_intestinal_health.compute_abdominal_distension, self.i_intestinal_health.compute_diarrhea,
                      self.i_intestinal_health.compute_ibs]
        with progressbar(l_function, label='장 건강'.center(12),
                         fill_char=style('#', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        self.__i_intestinal_health.o_score_detail_file.close()
        secho('>>> 장 건강: 점수 계산 완료', fg='cyan')

    def compute_wellness(self):
        self.__i_wellness = Wellness(self.d_bacteria, self.d_dist_db,
                                     self.d_phylum_table, self.d_family_table,
                                     self.d_genus_table, self.d_species_table)
        score_detail_file = os.path.join(self.score_path, f'{self.__i_wellness.main_item}.detail')
        self.__i_wellness.set_save_score_detail_file(score_detail_file)
        l_function = [self.i_wellness.compute_happiness_index, self.i_wellness.compute_tiredness_index,
                      self.i_wellness.compute_immune_index, self.i_wellness.compute_obesity_index,
                      self.i_wellness.compute_sleep_index, self.i_wellness.compute_aging_index]
        with progressbar(l_function, label='웰니스'.center(12), fill_char=style('#', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        self.__i_wellness.o_score_detail_file.close()
        secho('>>> 웰니스: 점수 계산 완료', fg='cyan')

    def compute_more_bowel_environments(self):
        self.__i_more_bowel_environments = MoreBowelEnvironments(self.d_bacteria, self.d_dist_db,
                                                                 self.d_genus_table, self.d_species_table,
                                                                 self.score_path, self.sample_name)
        score_detail_file = os.path.join(self.score_path, f'{self.__i_more_bowel_environments.main_item}.detail')
        self.__i_more_bowel_environments.set_save_score_detail_file(score_detail_file)
        l_function = [self.__i_more_bowel_environments.compute_beneficial_bacteria,
                      self.__i_more_bowel_environments.compute_harmful_bacteria,
                      self.__i_more_bowel_environments.predict_enterotype]
        with progressbar(l_function, label='장 환경 더 보기', fill_char=style('#', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        self.__i_more_bowel_environments.o_score_detail_file.close()
        secho('>>> 장 환경 더 보기: 점수 계산 완료', fg='cyan')

    def compute_supplement(self):
        self.__i_supplement = Supplement(self.d_bacteria, self.d_dist_db,
                                         self.d_genus_table, self.d_species_table)
        score_detail_file = os.path.join(self.score_path, f'{self.__i_supplement.main_item}.detail')
        self.__i_supplement.set_save_score_detail_file(score_detail_file)
        l_function = [self.i_supplement.compute_dietary_fibre, self.i_supplement.compute_lactose,
                      self.i_supplement.compute_resistant_starch, self.i_supplement.compute_protein,
                      self.i_supplement.compute_essential_amino_acid, self.i_supplement.compute_monoenoic_fatty_acid,
                      self.i_supplement.compute_b1, self.i_supplement.compute_b2,
                      self.i_supplement.compute_b3, self.i_supplement.compute_b6,
                      self.i_supplement.compute_b7, self.i_supplement.compute_b9,
                      self.i_supplement.compute_b12, self.i_supplement.compute_k2]
        with progressbar(l_function, label='부록'.center(13), fill_char=style('#', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        self.__i_supplement.o_score_detail_file.close()
        secho('>>> 부록: 점수 계산 완료', fg='cyan')

    def compute_supplement_probiotics21(self):
        self.__i_probiotics21 = SupplementProbiotics21(self.d_bacteria, self.d_dist_db,
                                                       self.d_probiotics19_species_table)
        self.__i_probiotics21.compute_probiotics_21()

    def compute_diversity_index(self):
        self.__i_diversity_index = DiversityIndex(self.d_bacteria, self.d_dist_db, self.alpha_diversity_file)
        score_detail_file = os.path.join(self.score_path, f'{self.__i_diversity_index.main_item}.detail')
        self.__i_diversity_index.set_save_score_detail_file(score_detail_file)
        self.__i_diversity_index.compute_shannon()
        self.__i_diversity_index.o_score_detail_file.close()

    def get_score(self):
        intestinal_health = ['constipation', 'gas', 'abdominal_distension', 'IBS', 'diarrhea']
        wellness = ['behavior', 'fatigue', 'immunity', 'obesity', 'sleep', 'aging']
        more_bowel_environments = ['benefit', 'harmful']
        supplement = ['fiber', 'lactose', 'starch', 'protein', 'bcaa', 'scfa', 'vitaminb1', 'vitaminb2', 'vitaminb3',
                      'vitaminb6', 'vitaminb7', 'vitaminb9', 'vitaminb12', 'vitamink2']
        diversity_index = ['shannon']
        main_items = [self.i_intestinal_health, self.i_wellness, self.i_more_bowel_environments,
                      self.i_supplement, self.i_diversity_index]
        sub_items = [intestinal_health, wellness, more_bowel_environments, supplement, diversity_index]
        for main_item, l_items in zip(main_items, sub_items):
            for item in l_items:
                yield item, getattr(main_item, item)

    def get_probiotics_score(self):
        for key, values in self.i_probiotics21.probiotics_21['positive']['species'].items():
            if values is None:
                yield key, values
            else:
                yield key, values.ratio

    def print_scores(self):
        for item, score in self.get_score():
            echo(f'{item}: {score}')

        echo('=' * 30)
        echo('=== Probiotics21 ===')
        for item, score in self.get_probiotics_score():
            echo(f'{item}: {score}')

    def save_scores(self):
        point_db_file = os.path.join(self.score_path, f'{self.sample_name}_points.csv')
        point_file = os.path.join(self.score_path, 'points.csv')
        probiotics_db_file = os.path.join(self.score_path, f'{self.sample_name}_probiotics.csv')
        probiotics_file = os.path.join(self.score_path, 'probiotics.csv')
        with open(point_db_file, 'w') as o_points_db:
            with open(point_file, 'w') as o_points:
                # 장 건강, 웰니스, 장 환경 더 보기, 부록, 다양성 지수 결과 저장
                for item, score in self.get_score():
                    o_points_db.write(f'0,{self.sample_name},{item},{score},{self.report_version}\n')
                    o_points.write(f'{item},{score}\n')
                # Enterotype 결과 저장
                if self.i_more_bowel_environments.enterotype == 'Bacteroides':
                    enterotype_code = 1
                elif self.i_more_bowel_environments.enterotype == 'Prevotella':
                    enterotype_code = 2
                elif self.i_more_bowel_environments.enterotype == 'Ruminococcus':
                    enterotype_code = 3
                o_points_db.write(f'0,{self.sample_name},enterotype,{enterotype_code},{self.report_version}\n')
                o_points.write(f'enterotype,{self.i_more_bowel_environments.enterotype}\n')
                # Enteortype Value 저장
                l_enterotype_value = [str(self.i_more_bowel_environments.d_enterotype_value['Bacteroides'] * 100),
                                      str(self.i_more_bowel_environments.d_enterotype_value['Prevotella'] * 100),
                                      str(self.i_more_bowel_environments.d_enterotype_value['Ruminococcus'] * 100)]
                o_points_db.write(f'0,{self.sample_name},enterotype_value,'
                                  f'{"|".join(l_enterotype_value)},{self.report_version}\n')
                for genus in ['Bacteroides', 'Prevotella', 'Ruminococcus']:
                    o_points.write(f'{genus},{self.i_more_bowel_environments.d_enterotype_value[genus]}\n')
        secho('>>> Points Score 저장', fg='cyan')
        echo(point_db_file)
        echo(point_file)
        with open(probiotics_db_file, 'w') as o_probiotics_db:
            with open(probiotics_file, 'w') as o_probiotics:
                for item, score in self.get_probiotics_score():
                    d_tax_id = self.d_tax_id.get(item)
                    if d_tax_id is None:
                        secho('Error: Probiotics에 해당 tax_id가 없습니다.', fg='red', blink=True, err=True)
                        echo(f'{item}: {d_tax_id}')
                        secho('\t--> tax_id_for_probiotics21.yaml 파일을 확인하세요.')
                        exit(1)
                    else:
                        o_probiotics.write(f'{item},{score}\n')
                        if score is None:
                            score = 0
                        o_probiotics_db.write(f'0,{self.sample_name},{d_tax_id["tax_id"]},'
                                              f'{score},probiotics,{self.report_version}\n')
        secho('>>> Probiotics Score 저장', fg='cyan')
        echo(probiotics_db_file)

