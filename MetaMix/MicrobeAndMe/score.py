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


from click import echo, secho
from collections import defaultdict, namedtuple
from SpoON.util import read_yaml


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
    def __init__(self, p_bacteria: dict):
        self.d_all_bacteria = p_bacteria
        self.items_name = None

    @property
    def d_bacteria(self):
        return self.d_all_bacteria[self.items_name]

    def make_effect_rank_dict_from_db1(self, index: str) -> dict:
        """
        Bacteria DB의 index에 해당되는 effect 및 rank를 추출하여 딕션너리로 반환한다.

        :param index: 매개변수들.
                    1개일 경우: index - [balance_index1 | balance_index2]
        :return: d_effect_rank
        """
        d_effect_rank = defaultdict(dict)
        for effect in self.d_bacteria[index].keys():
            for rank in self.d_bacteria[index][effect].keys():
                # d_effect_rank[effect].append(rank)
                d_effect_rank[effect] = {rank: self.d_bacteria[index][effect][rank]}
        return d_effect_rank

    def make_effect_rank_dict_from_db2(self, index1: str, index2: str) -> dict:
        """
        Bacteria DB의 index1, index2에 해당되는 effect 및 rank를 추출하여 딕션너리로 반환한다.

        :param index1: index1 - [health_index]
        :param index2: index2 - [constipation | diarrhea | bloating | obesity | blood_glucose | cholesterol |
                       non_alcoholic | irritable_bowel | inflammatory_bowel | crohn]

        :return: d_effect_rank
        """
        d_effect_rank = defaultdict(list)
        for effect in self.d_bacteria[index1][index2].keys():
            for rank in self.d_bacteria[index1][index2][effect].keys():
                d_effect_rank[effect].append(rank)
        return d_effect_rank


class IntestinalHealth(Score):
    def __init__(self, p_bacteria: dict, p_genus: dict, p_species: dict):
        super().__init__(p_bacteria)
        self.items_name = 'intestinal_health'
        self.d_genus_table = None
        self.d_species_table = None

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
        index_name = 'constipation'
        d_effect_rank = self.make_effect_rank_dict_from_db1(index_name)
        


class Wellness(Score):
    def __init__(self):
        # property
        self.__happiness_index = None
        self.__tireness_index = None
        self.__immune_index = None
        self.__obesity_index = None
        self.__sleep_index = None
        self.__aging_index = None

    @property
    def happiness_index(self):
        return self.__happiness_index

    @property
    def tireness_index(self):
        return self.__tireness_index

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


class MoreBowelEnvironments(Score):
    def __init__(self):
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


class Supplement(Score):
    def __init__(self):
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


class AndMe:
    def __init__(self, p_bact_file, p_genus_file, p_species_file, p_sample_name):
        self.bact_file = p_bact_file
        self.genus_file = p_genus_file
        self.species_file = p_species_file
        self.sample_name = p_sample_name
        self._d_bacteria = None
        self.__d_genus_table = None
        self.__d_species_table = None

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
        self._d_bacteria = read_yaml(self.bact_file)
        secho('>>> Items & Bacteria 읽기 완료', fg='cyan')
        echo(self.bact_file)

    def read_genus_file(self):
        self.__d_genus_table = self._read_table(6)

    def read_species_file(self):
        self.__d_species_table = self._read_table(7)

    def compute_intestinal_health(self):
        i_intestinal_health = IntestinalHealth(self._d_bacteria, self.d_genus_table, self.d_species_table)

    def compute_wellness(self):
        i_wellness = Wellness()

    def compute_more_bowel_environments(self):
        i_more = MoreBowelEnvironments()

    def compute_supplement(self):
        i_supplement = Supplement()

    def compute_diversity_index(self):
        i_diversity = DiversityIndex()


if __name__ == '__main__':
    and_me = AndMe('./MicrobeAndMe/items_and_bacteria.yaml',
                   '/garnet/Analysis/BI/AmpliconMetaGenome_MAM/HN00122513/Analysis_1/35.8/Summarized_Taxa/ASVs.NCBI_16S_L6.txt',
                   '/garnet/Analysis/BI/AmpliconMetaGenome_MAM/HN00122513/Analysis_1/35.8/Summarized_Taxa/ASVs.NCBI_16S_L7.txt',
                   '35.8')
    and_me.read_bacteria_file()
    and_me.read_genus_file()
    and_me.read_species_file()
    
    # test = DiversityIndex('/garnet/Analysis/BI/AmpliconMetaGenome_MAM/HN00122513/Analysis_1/35.8/Alpha_Diversity/adiv.txt')
    