import sys
sys.path.append('/crystal/Tools/Amplicon_MetaGenome/QIIME2_PipeLine/MetaMix/myBiomeStory')
from data_struct import MyBiomeStory
import os
from pprint import pprint
import pdb

# path = '/crystal/Analysis/BI/myBiomeStory/HN00103524/Analysis_1/018375'
path = '/crystal/Analysis/BI/myBiomeStory/HN00103524/Analysis_1'
table = 'Taxonomy_Assignment/otu_table_L7.txt'
alpha = 'Alpha_Diversity.txt'

path2 = '/crystal/Tools/Amplicon_MetaGenome/QIIME2_PipeLine/MetaMix/myBiomeStory'
biom_db = 'myBiomDB.yaml'
bacteria = 'bacteria_list.yaml'
tax_id = 'tax_id.yaml'

for sample_name in ['018375', '040689', '083660', '099176', '099665', '130766',
                    '140834', '206621', '219880', '239260', '260446', '261284',
                    '275912', '278832', '316679', '324565', '331920', '332666',
                    '355349', '356727', '399840', '400348', '400800', '409490',
                    '412649', '432756', '435491', '445577', '455181', '462165',
                    '468234', '477203', '514896', '558934', '559320', '559321',
                    '559322', '559331', '559354', '559950', '559988', '559993',
                    '560168', '560176', '560178', '560203', '560394']:
    print(sample_name)
    table_file = os.path.join(path, sample_name, table)
    alpha_file = os.path.join(path, sample_name, alpha)
    biom_db_file = os.path.join(path2, biom_db)
    bacteria_file = os.path.join(path2, bacteria)
    tax_id_file = os.path.join(path2, tax_id)
    mybiom = MyBiomeStory(table_file, alpha_file, bacteria_file, biom_db_file)
    mybiom.build_data_structure()
    mybiom.compute_data()
    sample_points = mybiom.story_index.transform_index(mybiom.d_db)
    print("- SamplePoint -")
    for i in dir(sample_points):
        if i.endswith('point'):
            pprint(f'{i}: {getattr(sample_points, i)}')
        elif i.startswith('_'):
            pass
        else:
            pprint(f'{i}: {getattr(sample_points, i)}')
    sample_biome = mybiom.story_index.arrange_sample_microbiome(tax_id_file, mybiom.d_bact_list)
    print()
    print('- Sample Microbiome - ')
    for i in ['beneficial', 'harmful', 'probiotics']:
            print(f'{i}:')
            pprint(getattr(sample_biome, i))
    # pdb.set_trace()
    print()

import pymysql
conn = pymysql.connect(host='211.192.85.228',
                       user='mybiome',
                       password='Qnstjr3qn!',
                       db='my-biomestory')