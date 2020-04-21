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
__version__ = '0.1.3'

from DBcontrol.dbBase import DBbase
from click import secho, echo
import pymysql


class DBmyBiomeStory(DBbase):
    def __init__(self, kargs):
        """

        :type kargs: dict
        :param kargs: db 연결에 필요한 정보를 담은 딕션너리
                [ MySQL ]
                    host - MySQL 서버의 주소.
                    user - 계정.
                    password - 비번.
                    db - DB 이름.
        """
        try:
            self.conn = pymysql.connect(host=kargs['ip'],
                                        user=kargs['user'],
                                        password=kargs['password'],
                                        db=kargs['db'],
                                        charset='utf8')

        except pymysql.err.OperationalError as err:
            secho('Error: MyBiomeStroy DB에 접근할 수 없습니다.', fg='red', blink=True)
            echo(err)
            secho('----> CSV 파일을 이용하여 DB에 입력하세요.', fg='megenta')
        else:
            secho('>>> MyBiomeStory DB 연결', fg='cyan')
        try:
            self.cursor = self.conn.cursor()
        except Exception as err:
            secho('Error: (MyBiomeStroy DB) Cursor 생성에 문제가 있습니다.', fg='red', blink=True, err=True)
            echo(err, err=True)
            secho('----> 개발자에게 문의하세요.', fg='megenta', err=True)

    def close(self):
        self.conn.close()
        secho('>>> MyBiomStroy DB 연결 해제', fg='cyan')

    def read_csv_excutemany(self, p_sql, p_csv_file):
        return super()._read_csv_excutemany_for_mysql(p_sql, p_csv_file)

    def execute_sql(self, p_sql, p_data):
        return super()._execute_sql_for_mysql(p_sql, p_data)

    @property
    def sql_point(self):
        sql = 'INSERT INTO sample_point(idx, kitid, label_point, point, version) ' \
              'VALUES(%s, %s, %s, %s, %s)'
        return sql

    @property
    def sql_biome(self):
        sql = 'INSERT INTO sample_microbiome(idx, kitid, taxid, organism_ratio, section, version) ' \
              'VALUES(%s, %s, %s, %s, %s, %s)'
        return sql

    @property
    def sql_sample(self):
        # Table Filed
        # idx	client_id_number	client_name	client_birth	client_sex	client_note	client_race
        sql = '''\
INSERT INTO sample(idx, client_id_number, client_name, client_birth, client_sex, client_note, client_race)
VALUES(%s, %s, %s, %s, %s, %s, %s)'''
        return sql

    @property
    def sql_sample_analysis(self):
        # Table Filed
        # idx	client_id_number	order_number	KitId	code_reseller	sample_source	sample_type	sample_kit
        # sample_application	date_collected	date_received	sample_note	del
        sql = '''\
INSERT INTO sample_analysis(idx, client_id_number, order_number, KitId, code_reseller, sample_source, sample_type, 
                            sample_kit, sample_application, date_collected, date_received, sample_note, del)
VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''
        return sql

    @property
    def sql_analysis_version(self):
        # Table Filed
        # idx	KitId	version	date_analysis	date_reported	del
        sql = 'INSERT INTO analysis_version(idx, KitId,	version, date_analysis, del) ' \
              'VALUES(%s, %s, %s, %s, %s)'
        return sql

    def insert_analysis_version(self, p_data):
        self.execute_sql(self.sql_analysis_version, p_data)

    def insert_point(self, p_csv_file):
        self.read_csv_excutemany(self.sql_point, p_csv_file)
        # secho('>>> Sample Points 입력 완료', fg='cyan')

    def insert_biome(self, p_csv_file):
        self.read_csv_excutemany(self.sql_biome, p_csv_file)
        # secho('>>> Sample Microbiome 입력 완료', fg='cyan')

    def insert_sample_info(self, p_info, p_order_number, p_time, p_race='Mongolid'):
        if p_info['client_id_number'] is None:
            client_id_number = p_time
        elif p_info['client_id_number'] == ' ':
            client_id_number = p_time
        elif p_info['client_id_number'] == '':
            client_id_number = p_time
        elif p_info['client_id_number']:
            client_id_number = p_info['client_id_number']
        else:
            raise ValueError
        sample = (0, client_id_number, p_info['client_name'], p_info['Birth'], p_info['Sex'], 'NULL', p_race)
        sample_analysis = (0, client_id_number, p_order_number, p_info['KitId'], p_info['code_reseller'],
                           p_info['sample_source'], p_info['sample_type'], p_info['sample_kit'], 'AmpliconMetagenome',
                           p_info['date_collected'], p_info['date_received'], p_info['sample_note'], 0)
        self.execute_sql(self.sql_sample, sample)
        self.execute_sql(self.sql_sample_analysis, sample_analysis)
