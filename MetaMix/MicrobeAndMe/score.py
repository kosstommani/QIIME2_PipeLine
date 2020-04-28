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
from SpoON.util import read_yaml
from time import sleep

BiomTable_L6 = namedtuple('BiomTable_L6', 'L1 L2 L3 L4 L5 L6 ratio')
BiomTable_L7 = namedtuple('BiomTable_L7', 'L1 L2 L3 L4 L5 L6 L7 ratio')


class DiversityIndex:
    def __init__(self, p_file):
        self.diversity_file = p_file
        self._shannon = None
        self._simpson = None

    @property
    def shannon(self):
        return self._shannon

    @property
    def simpson(self):
        return self._simpson

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
        self._shannon = shannon
        self._simpson = simpson
        # secho('>>> Alpha Diversity 파일 읽기 완료', fg='cyan')
        # echo(self.__diversity_file)


class Score:
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        self.d_all_bacteria = p_bacteria
        self.d_all_dist_db = p_dist_db
        self.d_genus_table = p_genus
        self.d_species_table = p_species
        self.main_item = None

    @property
    def d_bacteria(self):
        return self.d_all_bacteria[self.main_item]

    @property
    def d_dist_db(self):
        return self.d_all_dist_db[self.main_item]

    def make_effect_rank_dict_from_db1(self, item: str) -> dict:
        """
        Bacteria DB의 index에 해당되는 effect 및 rank를 추출하고, table에서 해당 taxon의 정보를 추출한다.
        딕션너리로 반환한다.

        :param item: 1개일 경우: index - [balance_index1 | balance_index2]
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

    # def make_effect_rank_dict_from_db2(self, index1: str, index2: str) -> dict:
    #     """
    #     Bacteria DB의 index1, index2에 해당되는 effect 및 rank를 추출하여 딕션너리로 반환한다.
    #
    #     :param index1: index1 - [health_index]
    #     :param index2: index2 - [constipation | diarrhea | bloating | obesity | blood_glucose | cholesterol |
    #                    non_alcoholic | irritable_bowel | inflammatory_bowel | crohn]
    #
    #     :return: d_effect_rank
    #     """
    #     d_effect_rank = defaultdict(list)
    #     for effect in self.d_bacteria[index1][index2].keys():
    #         for rank in self.d_bacteria[index1][index2][effect].keys():
    #             d_effect_rank[effect].append(rank)
    #     return d_effect_rank

    def extract_taxon(self, rank: str, taxon: str) -> namedtuple:
        if rank == 'phylum':
            table = None
        elif rank == 'family':
            table = None
        elif rank == 'genus':
            table = self.d_genus_table
        elif rank == 'species':
            table = self.d_species_table
        else:
            raise ValueError(f'rank 값이 올바르지 않습니다. rank: {rank}')
        return table.get(taxon)

    def _transform(self, p_type: str, p_min: float, p_max: float, p_value: float) -> float:
        if p_type == 'positive':
            if p_value <= p_min:
                return 0
            elif p_value >= p_max:
                return 100
            else:
                return (p_value - p_min) / (p_max - p_min) * 100
        elif p_type == 'negative':
            if p_value <= p_min:
                return 100
            elif p_value >= p_max:
                return 0
            else:
                return 100 - ((p_value - p_min) / (p_max - p_min) * 100)
        else:
            raise ValueError(f'p_type 옵션값이 올바르지 않습니다. p_type: {p_type}')

    def transform_data(self, p_item: str, p_d_effect_rank: dict) -> float:
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
        return sum(d_transformed_data.values()) * 0.5


class IntestinalHealth(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'intestinal_health'

        # property
        self.__constipation = None
        self.__gas = None
        self.__bloating = None
        self.__neurotic_abdominal_discomfort = None
        self.__diarrhea = None

    @property
    def constipation(self):
        return self.__constipation

    @property
    def gas(self):
        return self.__gas

    @property
    def bloating(self):
        return self.__bloating

    @property
    def neurotic_abdominal_discomfort(self):
        return self.__neurotic_abdominal_discomfort

    @property
    def diarrhea(self):
        return self.__diarrhea

    def compute_constipation(self):
        item_name = 'constipation'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__constipation = self.transform_data(item_name, d_effect_rank)

    def compute_gas(self):
        item_name = 'gas'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__gas = self.transform_data(item_name, d_effect_rank)

    def compute_bloating(self):
        item_name = 'bloating'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__bloating = self.transform_data(item_name, d_effect_rank)

    def compute_neurotic_abdominal_discomfort(self):
        item_name = 'neurotic_abdominal_discomfort'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__neurotic_abdominal_discomfort = self.transform_data(item_name, d_effect_rank)

    def compute_diarrhea(self):
        item_name = 'diarrhea'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__diarrhea = self.transform_data(item_name, d_effect_rank)


class Wellness(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict,
                 p_phylum: dict, p_family: dict,
                 p_genus: dict, p_species: dict):
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.d_phylum_table = p_phylum
        self.d_family_table = p_family
        self.main_item = 'wellness'
        # property
        self.__happiness_index = None
        self.__tiredness_index = None
        self.__immune_index = None
        self.__obesity_index = None
        self.__sleep_index = None
        self.__aging_index = None

    @property
    def happiness_index(self):
        return self.__happiness_index

    @property
    def tiredness_index(self):
        return self.__tiredness_index

    @property
    def immune_index(self):
        return self.__immune_index

    @property
    def obesity_index(self):
        return self.__obesity_index

    @property
    def sleep_index(self):
        return self.__sleep_index

    @property
    def aging_index(self):
        return self.__aging_index

    def compute_happiness_index(self):
        item_name = 'happiness_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__happiness_index = self.transform_data(item_name, d_effect_rank)

    def compute_tiredness_index(self):
        item_name = 'tiredness_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__tiredness_index = self.transform_data(item_name, d_effect_rank)

    def compute_immune_index(self):
        item_name = 'immune_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__immune_index = self.transform_data(item_name, d_effect_rank)

    def compute_obesity_index(self):
        item_name = 'obesity_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__obesity_index = self.transform_data(item_name, d_effect_rank)

    def compute_sleep_index(self):
        item_name = 'sleep_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__sleep_index = self.transform_data(item_name, d_effect_rank)

    def compute_aging_index(self):
        item_name = 'aging_index'
        d_effect_rank = self.make_effect_rank_dict_from_db1(item_name)
        self.__aging_index = self.transform_data(item_name, d_effect_rank)


class MoreBowelEnvironments(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        super().__init__(p_bacteria, p_dist_db, p_genus, p_species)
        self.main_item = 'MoreBowelEnvironments'
        self.__beneficial_bacteria = None
        self.__harmful_bacteria = None
        self.__enterotype = None

    @property
    def beneficial_bacteria(self):
        return self.__beneficial_bacteria

    @property
    def harmful_bacteria(self):
        return self.__harmful_bacteria

    @property
    def enterotype(self):
        return self.__enterotype

    def compute_beneficial_bacteria(self):
        pass

    def compute_harmful_bacteria(self):
        pass

    def predict_enterotype(self):
        pass


class Supplement(Score):
    def __init__(self, p_bacteria: dict, p_dist_db: dict, p_genus: dict, p_species: dict):
        super().__init__()
        # property
        self.__dietary_fibre = None
        self.__lactose = None
        self.__resistant_starch = None
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
        self.__probiotics_19 = None

    @property
    def dietary_fibre(self):
        return self.__dietary_fibre

    @property
    def lactose(self):
        return self.__lactose

    @property
    def resistant_starch(self):
        return self.__resistant_starch

    @property
    def protein(self):
        return self.__protein

    @property
    def essential_amino_acid(self):
        return self.__essential_amino_acid

    @property
    def monoenoic_fatty_acid(self):
        return self.__monoenoic_fatty_acid

    @property
    def b1(self):
        return self.__b1

    @property
    def b2(self):
        return self.__b2

    @property
    def b3(self):
        return self.__b3

    @property
    def b6(self):
        return self.__b6

    @property
    def b7(self):
        return self.__b7

    @property
    def b9(self):
        return self.__b9

    @property
    def b12(self):
        return self.__b12

    @property
    def k2(self):
        return self.__k2

    @property
    def probiotics_19(self):
        return self.__probiotics_19

    def compute_dietary_fibre(self):
        pass

    def compute_lactose(self):
        pass
    
    def compute_resistant_starch(self):
        pass

    def compute_protein(self):
        pass
    
    def compute_essential_amino_acid(self):
        pass
    
    def compute_monoenoic_fatty_acid(self):
        pass
    
    def compute_b1(self):
        pass
    
    def compute_b2(self):
        pass
    
    def compute_b3(self):
        pass

    def compute_b6(self):
        pass

    def compute_b7(self):
        pass
    
    def compute_b9(self):
        pass
    
    def compute_b12(self):
        pass
    
    def compute_k2(self):
        pass
    
    def compute_probiotics_19(self):
        pass


class AndMe:
    # def __init__(self, p_bact_file, p_distribution_db_file, p_genus_file, p_species_file, p_sample_name):
    def __init__(self, p_kargs):
        self.bact_file = p_kargs['bact_file']
        self.dist_db_file = p_kargs['dist_db_file']
        self.phylum_file = p_kargs['phylum_file']
        self.family_file = p_kargs['family_file']
        self.genus_file = p_kargs['genus_file']
        self.species_file = p_kargs['species_file']
        self.sample_name = p_kargs['sample_name']
        self.__d_bacteria = None
        self.__d_dist_db = None
        self.__d_genus_table = None
        self.__d_species_table = None

    @property
    def d_bacteria(self):
        return self.__d_bacteria

    @property
    def d_dist_db(self):
        return self.__d_dist_db

    @property
    def d_genus_table(self):
        return self.__d_genus_table

    @property
    def d_species_table(self):
        return self.__d_species_table

    def _read_table(self, p_l: int):
        """
        otu_table_L?.txt 파일을 읽어 리스트에 저장한다.

        :return: data
        """
        d_data = defaultdict()
        sample_indicator = False
        if p_l == 6:
            table_file = self.genus_file
        elif p_l == 7:
            table_file = self.species_file
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
                    ratio = float(l_ele[1].strip()) * 100
                    if p_l == 6:
                        d_data[rank[5]] = BiomTable_L6(L1=rank[0], L2=rank[1], L3=rank[2],
                                                       L4=rank[3], L5=rank[4], L6=rank[5],
                                                       ratio=ratio)
                    elif p_l == 7:
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

    def read_genus_file(self):
        self.__d_genus_table = self._read_table(6)

    def read_species_file(self):
        self.__d_species_table = self._read_table(7)

    def compute_intestinal_health(self):
        i_intestinal_health = IntestinalHealth(self.d_bacteria, self.d_dist_db,
                                               self.d_genus_table, self.d_species_table)
        l_function = [i_intestinal_health.compute_constipation, i_intestinal_health.compute_gas,
                      i_intestinal_health.compute_bloating, i_intestinal_health.compute_diarrhea,
                      i_intestinal_health.compute_neurotic_abdominal_discomfort]
        with progressbar(l_function, label='장 건강',
                         fill_char=style('#', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        secho('>>> 장 건강: 점수 계산 완료', fg='cyan')

    def compute_wellness(self):
        i_wellness = Wellness(self.d_bacteria, self.d_dist_db,
                              self.d_genus_table, self.d_species_table)
        l_function = [i_wellness.compute_happiness_index, i_wellness.compute_tiredness_index,
                      i_wellness.compute_immune_index, i_wellness.compute_obesity_index,
                      i_wellness.compute_sleep_index, i_wellness.compute_aging_index]
        with progressbar(l_function, label='웰니스',
                         fill_char=style('$', fg='green')) as bar_l_function:
            sleep(1)
            for function in bar_l_function:
                function()
        secho('>>> 웰니스: 점수 계산 완료', fg='cyan')

    def compute_more_bowel_environments(self):
        i_more = MoreBowelEnvironments()

    def compute_supplement(self):
        i_supplement = Supplement()

    def compute_diversity_index(self):
        i_diversity = DiversityIndex()
