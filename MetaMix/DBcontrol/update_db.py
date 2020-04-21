#!/crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
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
__version__ = '0.1'

import sys
from dbMyBiomeStory import DBmyBiomeStory

if len(sys.argv) != 2:
    print('update_db.py [data]')
    print('ex) data')
    print('\tbalance_harmful_point	26|3	MBS190425000112')
    print('\tbalance_beneficial_point	74|26	MBS190425000112')
    exit()
with open(sys.argv[1], 'r') as o_txt:
    data = list()
    for text in o_txt:
        temp = text.strip().split('\t')
        data.append((temp[1], temp[2], temp[0]))

# print(data)

db_story = DBmyBiomeStory({'ip': '172.19.87.50',
                           'user': 'mybiome',
                           'password': 'Qnstjr3qn!',
                           'db': 'my-biomestory'
                           })

db_story.start_transaction()

sql = 'UPDATE sample_point SET point = %s WHERE kitid = %s AND label_point = %s'
try:
    for ele in data:
        db_story.execute_sql(sql, ele)
except:
    db_story.rollback()
else:
    db_story.commit()
    db_story.close()
