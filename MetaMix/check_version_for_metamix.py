#!/garnet/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
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
# |   __|___ ___|     |   | |
# |__   | . | . |  |  | | | |
# |_____|  _|___|_____|_|___|
#       |_|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------
__author__ = 'JungWon Park(KOSST)'
__version__ = '1.1'

import os


def print_version(p_modules):
    print('-'*20)
    for module_text in p_modules:
        try:
            exec(f'import {module_text}')
            module = eval(f'{module_text}')
            print('{:<34} : {}'.format(module.__name__, module.__version__))
        except ModuleNotFoundError as err:
            print(f'{module_text:<34} : {err}')
    print('-'*20)
    print()

print('-'*20)
os.system('./PreMA.py --version')
os.system('./CoFI.py --version')
os.system('./MetaMix.py --version')
os.system('./theCups.py --version')
print('-'*20)


print('PreMA')
l_prema_modules = ['PreMA.adapter_trim', 'PreMA.core', ' PreMA.data_structure', 'PreMA.rawdata_handler']
print_version(l_prema_modules)


print('CoFI')
l_cofi_modules = ['CoFI.alignment', 'CoFI.ASVs', 'CoFI.core', 'CoFI.data_handler', 'CoFI.data_structure',
                  'CoFI.diversity', 'CoFI.OTU', 'CoFI.phylogeny', 'CoFI.read_assembly', 'CoFI.taxonomy']
print_version(l_cofi_modules)

print('SpoON')
l_spoon_modules = ['SpoON.fastq_handler', 'SpoON.util',  'SpoON.run_time']
print_version(l_spoon_modules)

# print('DBcontrol')
# l_dbcontrol_modules = [DBcontrol.dbBase, DBcontrol.dbMyBiomeStory]

# print('myBiomeStory')
# l_mybiomestory = [myBiomeStory.data_struct,  myBiomeStory.my_util, myBiomeStory.qiime1_cmd]

print('theCups')
l_thecups_modules = ['theCups.beta_diversity', 'theCups.clustering_and_otu_picking',
                     'theCups.core', 'theCups.data_structure',
                     'theCups.diversity_index', 'theCups.pcoa',
                     'theCups.probiotics', 'theCups.rarefaction',
                     'theCups.read_assembly', 'theCups.report',
                     'theCups.report_miseq_v1', 'theCups.report_miseq_v2',
                     'theCups.taxonomy', 'theCups.upgma_tree',
                     ]
print_version(l_thecups_modules)


