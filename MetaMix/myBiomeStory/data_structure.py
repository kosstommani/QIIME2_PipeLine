# ----------------------------------------------------------------------------------------------------------------------
#                        888b     d888          888             888b     d888 d8b
#                        8888b   d8888          888             8888b   d8888 Y8P
#                        88888b.d88888          888             88888b.d88888
# 88888b.d88b.  888  888 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b 888  888 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 Y88b 888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888  "Y88888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#                    888
#               Y8b d88P
#                "Y88P"
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '0.1.2'

import os
from collections import namedtuple, defaultdict
from click import echo, secho

BiomTable = namedtuple('BiomTable_L7', 'L1 L2 L3 L4 L5 L6 L7 ratio')
AlphaDiversity = namedtuple('AlphaDiversity', 'shannon simpson')
HealthIndexBoth = namedtuple('HealthIndexBoth', 'positive negative')
HealthIndexPosi = namedtuple('HealthIndexPositive', 'positive')
HealthIndexNega = namedtuple('HealthIndexNegative', 'negative')
SamplePoints = namedtuple('SamplePoints',
                          ['total_point', 'health_point', 'balance_beneficial_point',
                           'balance_harmful_point', 'enterotype_code', 'enterotype_div_bacteroides',
                           'enterotype_div_prevotella', 'enterotype_div_ruminococcus', 'constipation',
                           'diarrhea', 'bloating', 'obesity',
                           'blood_glucose', 'cholesterol', 'non_alcoholic',
                           'irritable_bowel', 'inflammatory_bowel', 'crohn',
                           'db_version'])
SampleMicrobiome = namedtuple('SampleMicrobiome', 'probiotics beneficial harmful db_version')


class StoryIndex:
    db_path = '/crystal/Tools/Amplicon_MetaGenome/QIIME2_PipeLine/MetaMix/myBiomeStory'
    tax_id = 'tax_id.yaml'
    tax_id_yaml = os.path.join(db_path, tax_id)
    d_tax_id = None
    d_db = None

    def __init__(self, p_sample_name):
        self.kitid = p_sample_name
        # Alpha Diversity
        self.alpha = None

        # main index
        self.balance_index1 = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.balance_index2 = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.gut_type_L6 = defaultdict(list)
        self.probiotics_L6 = defaultdict(list)
        self.probiotics_L7 = defaultdict(list)

        # health index
        self.constipation = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.diarrhea = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.bloating = HealthIndexNega(negative=defaultdict(dict))
        self.obesity = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.blood_glucose = HealthIndexPosi(positive=defaultdict(dict))
        self.cholesterol = HealthIndexPosi(positive=defaultdict(dict))
        self.non_alcoholic = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.irritable_bowel = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))
        self.inflammatory_bowel = HealthIndexNega(negative=defaultdict(dict))
        self.crohn = HealthIndexBoth(positive=defaultdict(dict), negative=defaultdict(dict))

        self.sample_points = None
        self.sample_biome = None

    @staticmethod
    def compute_trans_alpha(p_index, p_db_min, p_db_max):
        """
        Alpha Diversity에 해당되는 Index를 DB에 등록된 정보를 바탕으로 표준화하고 표준화된 점수를 반환한다.

        :type p_index: float
        :param p_index: s
        :param p_db_min:
        :param p_db_max:
        :rtype: float
        :return: 0 or float
        """
        if p_index > p_db_min:
            return round((p_index - p_db_min) / (p_db_max - p_db_min) * 100, 5)
        else:
            return 0

    def transform_alpha(self, p_shannon, p_simpson):
        """
        다양성 지수인 shannon과 simpson 값을 DB의 정보를 바탕으로 표준화한다.
        표준화한 점수를 바탕으로 health point를 산출한다.

        :type p_shannon: dict
        :param p_shannon: DB의 shannon 정보
        :type p_simpson: dict
        :param p_simpson: DB의 simpson 정보
        :rtype: float
        :return: health_point
        """
        trans_shannon = self.compute_trans_alpha(self.alpha.shannon, p_shannon['min'], p_shannon['max'])
        trans_simpson = self.compute_trans_alpha(self.alpha.simpson, p_simpson['min'], p_simpson['max'])
        health_point = round((trans_shannon + trans_simpson) / 2, 5)
        return health_point

    def transform_balance1(self):
        """
        balance_index1에 해당되는 값들을 변환하여 변환점수 및 변환비율을 반환한다.
        유익균 점수: (Bifidobacterium + Lactobacillus) / Total sum *100
        유해균 점수: 100 - 유익균 점수
        유익균 비율: 유익균 점수 * 0.1
        유해균 비율: 10 - 유익균 비율

        :return: beneficial_point, beneficial_relative_ratio, harmful_point, harmful_relative_ratio
        """
        positive_sum = 0
        for rank in self.balance_index1.positive.keys():
            for taxon in self.balance_index1.positive[rank].keys():
                positive_sum += sum(self.balance_index1.positive[rank][taxon])
        negative_sum = 0
        for rank in self.balance_index1.negative.keys():
            for taxon in self.balance_index1.negative[rank].keys():
                negative_sum += sum(self.balance_index1.negative[rank][taxon])
        beneficial_point = round(positive_sum / (positive_sum + negative_sum) * 100)
        # harmful_point = round(negative_sum / (positive_sum + negative_sum) * 100)
        harmful_point = 100 - beneficial_point
        beneficial_relative_ratio = round(beneficial_point * 0.1)
        # harmful_relative_ratio = round(harmful_point * 0.1)
        harmful_relative_ratio = 10 - beneficial_relative_ratio
        return beneficial_point, beneficial_relative_ratio, harmful_point, harmful_relative_ratio

    def max_gut_type(self):
        """
        장유형을 구분하는 지표 미생물의 분포 비율을 확인하여 해당되는 장유형(enterotype)을 코드로 반환한다.
        Bacteroides: 1
        Prevotella: 2
        Ruminococcus: 3

        :rtype: int or None
        :return: 장유형[1 | 2 | 3 | None]
        """
        max_name = None
        max_value = 0
        for taxon, ratio in self.gut_type_L6.items():
            if sum(ratio) > max_value:
                max_name = taxon
                max_value = sum(ratio)
        if max_name == 'Bacteroides':
            return 1
        elif max_name == 'Prevotella':
            return 2
        elif max_name == 'Ruminococcus':
            return 3
        else:
            secho(f'Warning :  max_gut_type : {max_name}', fg='yellow')
            echo(self.gut_type_L6.items())
            return None

    @staticmethod
    def compute_both_point(p_posi_sum, p_nega_sum, p_db_posi, p_db_nega):
        """
        Ref. DB의 항목(index)에서 positive effect 와 negative effect를 바탕으로 전달된 positive effect의 총합과
        negative effect의 총합을 표준화하고, 표준화된 점수를 반환한다.

        :type p_posi_sum: float
        :param p_posi_sum: 분석된 시료에서의 해당 항목의 positive effect 에 해당되는 미생물의 비율의 합
        :type p_nega_sum: float
        :param p_nega_sum: 분석된 시료에서의 해당 항목의 negative effect 에 해당되는 미생물의 비율의 합
        :type p_db_posi: dict
        :param p_db_posi: Ref. DB 에서 해당 항목의 positive effect 의 min, max
        :type p_db_nega: dict
        :param p_db_nega: Ref. DB 에서 해당 항목의 negative effect 의 min, max
        :return: point
        """
        if p_posi_sum < p_db_posi['min']:
            posi_sum = p_db_posi['min']
        elif p_posi_sum > p_db_posi['max']:
            posi_sum = p_db_posi['max']
        else:
            posi_sum = p_posi_sum
        if p_nega_sum < p_db_nega['min']:
            nega_sum = p_db_nega['min']
        elif p_nega_sum > p_db_nega['max']:
            nega_sum = p_db_nega['max']
        else:
            nega_sum = p_nega_sum
        positive_ratio = (posi_sum - p_db_posi['min']) / (p_db_posi['max'] - p_db_posi['min']) * 100
        negative_ratio = (nega_sum - p_db_nega['min']) / (p_db_nega['max'] - p_db_nega['min']) * 100
        point = round((positive_ratio + (100 - negative_ratio)) / 2, 5)
        return point

    @staticmethod
    def compute_one_point(p_sum, p_db):
        """
        Ref. DB의 항목(index)에서 positive effect 또는 negative effect를 바탕으로 전달된 positive effect의 총합 또는
        negative effect의 총합을 표준화하고, 표준화된 점수를 반환한다.

        :type p_sum: float
        :param p_sum: 분석된 시료에서의 해당 항목의 positive effect 또는 negative effect에 해당되는 미생물의 비율의 합
        :type p_db: dict
        :param p_db: Ref. DB의 해당 항목(index)의 positive 또는 negative 의 min, max
        :return:
        """
        if p_sum < p_db['min']:
            one_sum = p_db['min']
        elif p_sum > p_db['max']:
            one_sum = p_db['max']
        else:
            one_sum = p_sum
        point = round((one_sum - p_db['min']) / (p_db['max'] - p_db['min']) * 100, 5)
        return point

    def transform_both(self, p_nt_index, p_db_posi, p_db_nega):
        """
        입력값에서 positive effect 와 negative effect 에 해당 되는 미생물의 비율의 총합을 계산한다.
        Ref. DB의 positive effect 와 negative effect 의 점수를 이용하여 구해진 총합을 표준화한다.

        :type p_nt_index: namedtuple
        :param p_nt_index: 계산할 항목(index)의 attribute
        :type p_db_posi: dict
        :param p_db_posi: Ref. DB의 해당 항목(index)의 positive
        :type p_db_nega: dict
        :param p_db_nega: Ref. DB의 해당 항목(index)의 negative
        :return: point
        """
        positive_sum = 0
        for rank in p_nt_index.positive.keys():
            for taxon in p_nt_index.positive[rank].keys():
                positive_sum += sum(p_nt_index.positive[rank][taxon])
        negative_sum = 0
        for rank in p_nt_index.negative.keys():
            for taxon in p_nt_index.negative[rank].keys():
                negative_sum += sum(p_nt_index.negative[rank][taxon])
        point = self.compute_both_point(positive_sum, negative_sum, p_db_posi, p_db_nega)
        return point

    def transform_one(self, p_nt_index, p_db):
        """
        입력값에서 positive effect 또는 negative effect 에 해당 되는 미생물의 비율의 총합을 계산한다.
        Ref. DB의 positive effect 또는 negative effect 의 점수를 이용하여 구해진 총합을 표준화한다.

        :type p_nt_index: namedtuple
        :param p_nt_index: 계산할 항목(index)의 attribute
        :type p_db: dict
        :param p_db: Ref. DB의 해당 항목(index)의 positive 또는 negative
        :return: point
        """
        one_sum = 0
        for rank in p_nt_index.keys():
            for taxon in p_nt_index[rank].keys():
                one_sum += sum(p_nt_index[rank][taxon])
        point = self.compute_one_point(one_sum, p_db)
        return point

    def transform_constipation(self, p_positive, p_negative):
        return self.transform_both(self.constipation, p_positive, p_negative)

    def transform_diarrhea(self, p_positive, p_negative):
        return self.transform_both(self.diarrhea, p_positive, p_negative)

    def transform_bloating(self, p_negative):
        return self.transform_one(self.bloating.negative, p_negative)

    def transform_obesity(self, p_positive, p_negative):
        return self.transform_both(self.obesity, p_positive, p_negative)

    def transform_blood_glucose(self, p_positive):
        return self.transform_one(self.blood_glucose.positive, p_positive)

    def transform_cholesterol(self, p_positive):
        return self.transform_one(self.cholesterol.positive, p_positive)

    def transform_non_alcoholic(self, p_positive, p_negative):
        return self.transform_both(self.non_alcoholic, p_positive, p_negative)

    def transform_irritable_bowel(self, p_positive, p_negative):
        return self.transform_both(self.irritable_bowel, p_positive, p_negative)

    def transform_inflammatory_bowel(self, p_negative):
        return self.transform_one(self.inflammatory_bowel.negative, p_negative)

    def transform_crohn(self, p_positive, p_negative):
        return self.transform_both(self.crohn, p_positive, p_negative)

    @staticmethod
    def transform_total_point(p_health_point, p_beneficial_ratio):
        """
        beneficial relative ratio 의 균형 비율을 균형지수 환산 테이블을 이용하여 점수로 변환한다.
        health point 와 균형지수 환산점수를 이용하여 total point를 계산하고 반환한다.

        비율   균형지수 환산값
        (비율 - beneficial_relative_ratio : harmful_relative_ratio)
        7:3  - 100
        8:2  - 90
        9:1  - 80
        6:4  - 70
        5:5  - 60
        4:6  - 50
        3:7  - 40
        2:8  - 30
        1:9  - 20
        0:10 - 10
        10:0 - 60

        :type p_health_point: float
        :param p_health_point: health_point
        :type p_beneficial_ratio: int
        :param p_beneficial_ratio: balance_index1 에서 계산한 beneficial ratio
        :return: (p_health_point + d_transformed_table[p_beneficial_ratio]) / 2
        """
        d_transformed_table = {7: 100, 8: 90, 9: 80, 6: 70, 5: 60,
                               4: 50, 3: 40, 2: 30, 1: 20, 0: 10,
                               10: 60}
        return (p_health_point + d_transformed_table[p_beneficial_ratio]) / 2

    def sum_enterotype_genus(self, p_genus):
        """
        gut_type_L6 Attribute에서 매개변수로 전달된 genus에 해당되는 genus의 비율들을 모두 합하여 반환한다.

        :param p_genus: entero type 를 확인할 수 있는 속 이름.
        :rtype: float
        :return: 지표 속의 총 비율
        """
        ratio = self.gut_type_L6.get(p_genus)
        if ratio is None:
            return 0
        else:
            return round(sum(ratio), 5)

    def transform_index(self):
        """
        분석되어진 미생물의 분포 정보들을 Ref. DB를 이용하여 표준화 및 변환한다.
        변환된 점수들을 이용하여 MyBiomeStory DB에 입력할 Sample Points 를 생성한다.

        :return: None
        """
        health_point = self.transform_alpha(self.d_db['shannon'], self.d_db['simpson'])
        beneficial_point, beneficial_relative_ratio, harmful_point, harmful_relative_ratio = self.transform_balance1()
        balance_beneficial_point = f'{beneficial_point}|{beneficial_relative_ratio}'
        balance_harmful_point = f'{harmful_point}|{harmful_relative_ratio}'
        entero_type_code = self.max_gut_type()
        enterotype_div_bacteroides = self.sum_enterotype_genus('Bacteroides')
        enterotype_div_prevotella = self.sum_enterotype_genus('Prevotella')
        enterotype_div_ruminococcus = self.sum_enterotype_genus('Ruminococcus')
        constipation_point = self.transform_constipation(self.d_db['constipation']['positive'],
                                                         self.d_db['constipation']['negative'])
        diarrhea_point = self.transform_diarrhea(self.d_db['diarrhea']['positive'],
                                                 self.d_db['diarrhea']['negative'])
        bloating_point = self.transform_bloating(self.d_db['bloating']['negative'])
        obesity_point = self.transform_obesity(self.d_db['obesity']['positive'],
                                               self.d_db['obesity']['negative'])
        blood_glucose_point = self.transform_blood_glucose(self.d_db['blood_glucose']['positive'])
        cholesterol_point = self.transform_cholesterol(self.d_db['cholesterol']['positive'])
        non_alcoholic_point = self.transform_non_alcoholic(self.d_db['non_alcoholic']['positive'],
                                                           self.d_db['non_alcoholic']['negative'])
        irritable_bowel_point = self.transform_irritable_bowel(self.d_db['irritable_bowel']['positive'],
                                                               self.d_db['irritable_bowel']['positive'])
        inflammatory_bowel_point = self.transform_inflammatory_bowel(self.d_db['inflammatory_bowel']['negative'])
        crohn_point = self.transform_crohn(self.d_db['crohn']['positive'],
                                           self.d_db['crohn']['negative'])
        self.sample_points = SamplePoints(
            total_point=self.transform_total_point(health_point, beneficial_relative_ratio),
            health_point=health_point,
            balance_beneficial_point=balance_beneficial_point,
            balance_harmful_point=balance_harmful_point,
            enterotype_code=entero_type_code,
            enterotype_div_bacteroides=enterotype_div_bacteroides,
            enterotype_div_prevotella=enterotype_div_prevotella,
            enterotype_div_ruminococcus=enterotype_div_ruminococcus,
            constipation=constipation_point,
            diarrhea=diarrhea_point,
            bloating=bloating_point,
            obesity=obesity_point,
            blood_glucose=blood_glucose_point,
            cholesterol=cholesterol_point,
            non_alcoholic=non_alcoholic_point,
            irritable_bowel=irritable_bowel_point,
            inflammatory_bowel=inflammatory_bowel_point,
            crohn=crohn_point,
            db_version=self.d_db['myBiomeDB']['version'])

    @classmethod
    def read_tax_id(cls):
        """
        tax_id.yaml를 읽어 딕션너리로 저장한다.

        :return: None
        """
        data = MyBiomeStory.read_yaml(cls.tax_id_yaml)
        cls.d_tax_id = data

    def arrange_sample_microbiome(self, p_bact_list):
        """
        balance_index2 index 와 probiotics index 에 해당되는 값들을 계산 및 정리한다.
        SampleMicrobiome 정보를 생성하고 sample_biome Attribute에 저장한다.
        SampleMicrobiome.probiotics: dict
                        .beneficial: dict
                        .harmful: dict
                        .db_version: str

        :type p_bact_list: dict
        :param p_bact_list: bacteria list 정보
        :return: None
        """
        d_beneficial = dict()
        d_harmful = dict()
        d_probiotics = dict()
        for effect in p_bact_list['balance_index2'].keys():
            if effect == 'positive':
                d_save = d_beneficial
            elif effect == 'negative':
                d_save = d_harmful
            else:
                raise ValueError(f'effect: {effect}')
            for rank in p_bact_list['balance_index2'][effect].keys():
                for taxon in p_bact_list['balance_index2'][effect][rank]:
                    ratio = getattr(self.balance_index2, effect)[rank].get(taxon)
                    if ratio is None:
                        d_save.update({self.d_tax_id[rank][taxon]: 0})
                    else:
                        d_save.update({self.d_tax_id[rank][taxon]: sum(ratio)})

        for rank in p_bact_list['probiotics'].keys():
            for taxon in p_bact_list['probiotics'][rank]:
                ratio = getattr(self, f'probiotics_{rank}')[taxon]
                if ratio is None:
                    d_probiotics.update({self.d_tax_id[rank][taxon]: 0})
                else:
                    d_probiotics.update({self.d_tax_id[rank][taxon]: sum(ratio)})
        self.sample_biome = SampleMicrobiome(probiotics=d_probiotics,
                                             beneficial=d_beneficial,
                                             harmful=d_harmful,
                                             db_version=self.d_db['myBiomeDB']['version'])

    def save_points_to_csv(self, p_out_file):
        """
        Story Index에서 계산한 정보 중 sample points 에 해당되는 값들을 MyBiomeStroy DB에 삽입할 수 있도록
        CSV 형태로 저장한다.

        :param p_out_file: sample points 를 저장할 파일명
        :return: True or False
        """
        try:
            with open(p_out_file, 'w') as o_output:
                # o_output.write('idx,kitid,label_point,point,version\n')
                for point in dir(self.sample_points):
                    if point.startswith('_'):
                        continue
                    elif point in ['count', 'index', 'db_version']:
                        continue
                    else:
                        # FIELD: idx kitid label_point point version
                        ratio = getattr(self.sample_points, point)
                        o_output.write(f'0,{self.kitid},{point},{ratio},'
                                       f'{self.sample_points.db_version}\n')
        except Exception as err:
            secho('Error: Exception - 아래의 메시지를 확인하세요.', fg='red')
            echo(err)
            return False
        else:
            return True

    def save_biome_to_csv(self, p_out_file):
        """
        Story Index에서 계산한 정보 중 sample biome 에 해당되는 값들을 MyBiomeStroy DB에 삽입할 수 있도록
        CSV 형태로 저장한다.

        :param p_out_file: sample biome 을 저장할 파일명
        :return: True or False
        """
        try:
            with open(p_out_file, 'w') as o_output:
                # o_output.write('idx,kitid,taxid,organism_ratio,section,version\n')
                for biome_type in dir(self.sample_biome):
                    if biome_type.startswith('_'):
                        continue
                    elif biome_type in ['count', 'index', 'db_version']:
                        continue
                    else:
                        type_tax = getattr(self.sample_biome, biome_type)
                        for taxid in type_tax.keys():
                            # FIELD: idx kitid taxid organism_ratio section version
                            ratio = type_tax[taxid]
                            o_output.write(f'0,{self.kitid},{taxid},{ratio},{biome_type},'
                                           f'{self.sample_biome.db_version}\n')
        except Exception as err:
            secho('Error: Exception - 아래의 메시지를 확인하세요.', fg='red')
            echo(err)
            return False
        else:
            return True


class MyBiomeStory:
    db_path = '/garnet/Tools/Amplicon_MetaGenome/PipeLine_Dev/MetaMix/myBiomeStory'
    biom_db = 'myBiomDB.yaml'
    bacteria = 'bacteria_list.yaml'
    biom_db_yaml = os.path.join(db_path, biom_db)
    bacteria_yaml = os.path.join(db_path, bacteria)
    d_bact_list = None
    d_db = None

    def __init__(self, p_sample_name, p_table_file, p_diversity_file):
        self.__sample_name = p_sample_name
        self.__table_file = p_table_file
        self.__diversity_file = p_diversity_file
        self.__l_table = self._read_table()
        self.story_index = None

    @property
    def sample_name(self):
        return self.__sample_name

    @property
    def table_file(self):
        return self.__table_file

    @property
    def diversity_file(self):
        return self.__diversity_file

    @property
    def l_table(self):
        return self.__l_table

    def _read_table(self):
        """
        otu_table_L7.txt 파일을 읽어 리스트에 저장한다.

        :return: data
        """
        data = list()
        sample_indicator = False
        with open(self.__table_file, 'r') as o_table:
            for ele in o_table:
                if ele.startswith('#'):
                    if self.sample_name in ele:
                        sample_indicator = True
                    continue
                else:
                    l_ele = ele.strip().split('\t')
                    rank = [x.strip() for x in l_ele[0].split(';__')]
                    ratio = float(l_ele[1].strip()) * 100
                    data.append(BiomTable(L1=rank[0], L2=rank[1], L3=rank[2], L4=rank[3],
                                          L5=rank[4], L6=rank[5], L7=rank[6], ratio=ratio))
        if sample_indicator is False:
            secho('Error: OTU Table 파일내에 존재하는 시료명과 결과 시료명이 일치하지 않습니다.', fg='red', err=True)
            secho('----> OTU Table 파일내에 기재된 시료명과 디렉터리명을 확인하세요.', fg='magenta', err=True)
            echo(self.__table_file, err=True)
            exit(1)
        # secho('>>> BIOM Table 파일 읽기 완료', fg='cyan')
        # echo(self.__table_file)
        return data

    def _read_alpha_diversity(self):
        """
        Alpha_Diversity.txt 파일을 읽고 데이터를 반환한다.

        :return: namedtuple - AlphaDiversity.shannon
                                            .simpson
        """
        data = None
        header = None
        with open(self.__diversity_file, 'r') as o_text:
            for line, text in enumerate(o_text):
                if line == 0:
                    header = text.strip().split('\t')
                elif line == 1:
                    if text.strip():
                        data = text.strip().split()[1:]  # Sample ID 제거
                    else:
                        secho('Error: 데이터가 없습니다.', fg='red')
                        secho(f'---> {self.__diversity_file}를 확인하세요.', fg='magenta')
                        exit(1)
                else:
                    if text.strip():
                        secho('Error: 데이터가 여러개 존재합니다. 시료 1개에 대한 데이터만 있어야 합니다.', fg='red')
                        secho(f'---> {self.__diversity_file}를 확인하세요.(데이터 개수 확인, 개행 여부 확인)',
                              fg='magenta')
                        exit(1)
        shannon = float(data[header.index('shannon')])
        simpson = float(data[header.index('simpson')])
        alpha_diversity = AlphaDiversity(shannon=shannon, simpson=simpson)
        # secho('>>> Alpha Diversity 파일 읽기 완료', fg='cyan')
        # echo(self.__diversity_file)
        return alpha_diversity

    @staticmethod
    def read_yaml(p_file):
        """
        yaml 파일을 읽어 딕션너리 객체로 반환한다.

        :param p_file: yaml 파일
        :return: data
        """
        import yaml
        with open(p_file, 'r') as o_yaml:
            data = yaml.load(o_yaml)
        return data

    @classmethod
    def read_bacteria_list(cls):
        """
        bacteria_list.yaml 파일을 읽고, 해당 데이터를 딕션너리 객체로 반환환다.

        :return: None
        """
        cls.d_bact_list = cls.read_yaml(cls.bacteria_yaml)
        secho('>>> Bacteria List YAML 파일 일기 완료', fg='cyan')
        echo(cls.bacteria_yaml)

    @classmethod
    def read_mybiome_db(cls):
        """
        myBiomDB.yaml 파일을 읽고 데이터를 반환한다.

        :return: data
        """
        cls.d_db = cls.read_yaml(cls.biom_db_yaml)
        secho('>>> myBiomeStroy-DB YAML 파일 일기 완료', fg='cyan')
        echo(cls.biom_db_yaml)

    def _find_data(self, p_text, *args):
        """
        찾고자하는 미생물의 이름이 self.d_bact_list(bacteria_list.yaml)에 존재 여부를 확인한다.

        :type p_text: str
        :param p_text: 찾고자하는 미생물의 이름.
        :param args: 매개변수들.
                     2개일 경우 : index, rank
                     3개일 경우 : index, effect, rank
                     4개일 경우 : index1, index2, effect, rank

                     index: [gut_type | probiotics]
                     index1: [main_index | health_index]
                     index2:
                            main_index :
                                [balance_index1 | balance_index2]
                            health_index :
                                [constipation | diarrhea | bloating | obesity | blood_glucose | cholesterol |
                                 non_alcoholic | irritable_bowel | inflammatory_bowel | crohn]
                     effect: [positive | negative]
                     rank: [L1 | L2 | L3 | L4 | L5 | L6 | L7]
        :rtype: bool
        :return: True or False
        """
        if len(args) == 2:
            index, rank = args
            try:
                _ = self.d_bact_list
                _ = self.d_bact_list[index][rank].index(p_text)
            except ValueError:
                return False
            else:
                return True
        elif len(args) == 3:
            index, effect, rank = args
            try:
                _ = self.d_bact_list[index][effect][rank].index(p_text)
            except ValueError:
                return False
            else:
                return True
        elif len(args) == 4:
            index1, index2, effect, rank = args
            try:
                _ = self.d_bact_list[index1][index2][effect][rank].index(p_text)
            except ValueError:
                return False
            else:
                return True
        else:
            raise ValueError

    def _check_index(self, p_taxon, *args):
        """
        Bacteria DB의 Index에 전달된 미생물의 존재 여부를 확인한 후, StoryIndex 인스턴스의 Index Attribute 에 저장한다.

        :param p_taxon: 찾고자하는 미생물의 이름.
        :param args: 매개변수들.
                    2개일 경우: index_name, l_rank
                    3개일 경우: index_name1, index_name2, l_effect_rank
        :return:
        """
        if len(args) == 2:
            index_name, l_rank = args
            for rank in l_rank:
                taxon_rank = getattr(p_taxon, rank)
                if self._find_data(taxon_rank, index_name, rank):
                    getattr(self.story_index, f'{index_name}_{rank}')[taxon_rank].append(p_taxon.ratio)
        elif len(args) == 3:
            index_name1, index_name2, l_effect_rank = args
            for effect in l_effect_rank.keys():
                for rank in l_effect_rank[effect]:
                    taxon_rank = getattr(p_taxon, rank)
                    if index_name1 == 'main_index':
                        if self._find_data(taxon_rank, index_name2, effect, rank):
                            story_index = getattr(getattr(self.story_index, f'{index_name2}'), effect)
                            if story_index[rank].get(taxon_rank) is None:
                                story_index[rank][taxon_rank] = list()
                            story_index[rank][taxon_rank].append(p_taxon.ratio)
                    else:
                        if self._find_data(taxon_rank, index_name1, index_name2, effect, rank):
                            story_index = getattr(getattr(self.story_index, f'{index_name2}'), effect)
                            if story_index[rank].get(taxon_rank) is None:
                                story_index[rank][taxon_rank] = list()
                            story_index[rank][taxon_rank].append(p_taxon.ratio)
        else:
            raise ValueError

    def _make_effect_rank_dict_from_db(self, *args):
        """
        Bacteria DB의 index에 해당되는 effect 및 rank를 추출하여 딕션너리로 반환한다.

        :param args: 매개변수들.
                    1개일 경우: index - [balance_index1 | balance_index2]
                    2개일 경우: index1 - [health_index]
                               index2 - [constipation | diarrhea | bloating | obesity | blood_glucose | cholesterol |
                                         non_alcoholic | irritable_bowel | inflammatory_bowel | crohn]

        :return: d_effect_rank
        """
        d_effect_rank = defaultdict(list)
        if len(args) == 1:
            index, = args
            for effect in self.d_bact_list[index].keys():
                for rank in self.d_bact_list[index][effect].keys():
                    d_effect_rank[effect].append(rank)
        elif len(args) == 2:
            index1, index2 = args
            for effect in self.d_bact_list[index1][index2].keys():
                for rank in self.d_bact_list[index1][index2][effect].keys():
                    d_effect_rank[effect].append(rank)
        else:
            raise ValueError
        return d_effect_rank

    def _make_rank_list_from_db(self, p_index):
        """
        Bacteria DB의 index에 해당되는 rank를 추출하여 리스트로 반환한다.

        :param p_index: [gut_type | probiotics]
        :return: list
        """
        return list(self.d_bact_list[p_index].keys())

    # 균형지수1
    def _check_balance_index1(self, p_taxon):
        index_name = 'balance_index1'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name)
        return self._check_index(p_taxon, 'main_index', index_name, d_effect_rank)

    # 균형지수2
    def _check_balance_index2(self, p_taxon):
        index_name = 'balance_index2'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name)
        return self._check_index(p_taxon, 'main_index', index_name, d_effect_rank)

    # 장유형
    def _check_gut_type(self, p_taxon):
        index_name = 'gut_type'
        l_rank = self._make_rank_list_from_db(index_name)
        return self._check_index(p_taxon, index_name, l_rank)

    # 유산균
    def _check_probiotics(self, p_taxon):
        index_name = 'probiotics'
        l_rank = self._make_rank_list_from_db(index_name)
        return self._check_index(p_taxon, index_name, l_rank)

    # 변비
    def _check_constipation(self, p_taxon):
        index_name1 = 'health_index'
        index_naem2 = 'constipation'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_naem2)
        return self._check_index(p_taxon, index_name1, index_naem2, d_effect_rank)

    # 설사
    def _check_diarrhea(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'diarrhea'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 복부팽만감
    def _check_bloating(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'bloating'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 비만
    def _check_obesity(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'obesity'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 혈당조절
    def _check_blood_glucose(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'blood_glucose'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 콜레스테롤
    def _check_cholesterol(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'cholesterol'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 비알콜성 지방간
    def _check_non_alcoholic(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'non_alcoholic'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 과민성대장증후군
    def _check_irritable_bowel(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'irritable_bowel'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 염증성 장질환
    def _check_inflammatory_bowel(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'inflammatory_bowel'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    # 크론병
    def _check_crohn(self, p_taxon):
        index_name1 = 'health_index'
        index_name2 = 'crohn'
        d_effect_rank = self._make_effect_rank_dict_from_db(index_name1, index_name2)
        return self._check_index(p_taxon, index_name1, index_name2, d_effect_rank)

    def build_data_structure(self):
        """
        StoryIndex 객체를 이용하여 시료에 대한 인스턴스을 생성한다.
        StoryIndex 인스턴스에 myBiomDB.yaml 파일을 읽은 DB 정보를 할당한다.

        :return: None
        """
        self.story_index = StoryIndex(self.sample_name)
        self.story_index.d_db = self.d_db

    def compute_data(self):
        """
        Alpha_Diversity.txt 파일을 읽고 StoryIndex 인스턴스에 저장한다.
        otu_table_L7.txt 파일의 데이터로부터 각 항목에 해당하는 미생물의 존재여부를 확인한다.

        :return: None
        """
        self.story_index.alpha = self._read_alpha_diversity()
        for taxon in self.l_table:
            self._check_balance_index1(taxon)
            self._check_balance_index2(taxon)
            self._check_gut_type(taxon)
            self._check_probiotics(taxon)
            self._check_constipation(taxon)
            self._check_diarrhea(taxon)
            self._check_bloating(taxon)
            self._check_obesity(taxon)
            self._check_blood_glucose(taxon)
            self._check_cholesterol(taxon)
            self._check_non_alcoholic(taxon)
            self._check_irritable_bowel(taxon)
            self._check_inflammatory_bowel(taxon)
            self._check_crohn(taxon)
