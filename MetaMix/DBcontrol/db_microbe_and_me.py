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


from DBcontrol.dbBase import DBbase
from click import secho, echo
import pymysql


class DBMicrobeAndMe(DBbase):
    def __init__(self, kargs):
        """

        :type kargs: dict
        :param kargs: db 연결에 필요한 정보를 담은 딕션너리
                [ mariaDB ]
                    host - maria 서버의 주소.
                    user - 계정.
                    password - 비번.
                    db - DB 이름.
        """
        try:
            self.conn = pymysql.connect(host=kargs['ip'],
                                        user=kargs['user'],
                                        password=kargs['password'],
                                        db=kargs['db'],
                                        port=kargs['port'],
                                        charset='utf8')

        except pymysql.err.OperationalError as err:
            secho('Error: Microbe&Me DB에 접근할 수 없습니다.', fg='red', blink=True)
            echo(err)
        else:
            secho('>>> Microbe&Me DB 연결', fg='cyan')
        try:
            self.cursor = self.conn.cursor()
        except Exception as err:
            secho('Error: (Microbe&Me DB) Cursor 생성에 문제가 있습니다.', fg='red', blink=True, err=True)
            echo(err, err=True)
            secho('----> 개발자에게 문의하세요.', fg='magenta', err=True)

    def close(self):
        self.conn.close()
        secho('>>> Microbe&Me DB 연결 해제', fg='cyan')

    def read_csv_excutemany(self, p_sql, p_csv_file):
        return self._read_csv_excutemany_for_mysql(p_sql, p_csv_file)

    def execute_sql(self, p_sql, p_data):
        return self._execute_sql_for_mysql(p_sql, p_data)

    def execute_sql_for_select(self, p_sql, p_data):
        if 'insert' in p_sql.lower():
            raise ValueError('SQL에 INSERT문이 있을 경우 메서드를 실행할 수 없습니다.')
        elif 'delete' in p_sql.lower():
            raise ValueError('SQL에 DELETE문이 있을 경우 메서드를 실행할 수 없습니다.')
        self.cursor.execute(p_sql, p_data)
        db_data = self.cursor.fetchall()
        return db_data

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

    @property
    def sql_del_sample(self):
        sql = 'DELETE FROM sample WHERE client_name = %s AND client_birth = %s AND client_sex = %s'
        return sql

    @property
    def sql_del_sample_analysis(self):
        sql = 'DELETE FROM sample_analysis WHERE kitid = %s'
        return sql

    @property
    def sql_del_sample_point(self):
        sql = 'DELETE FROM sample_point WHERE kitid = %s AND version = %s'
        return sql

    @property
    def sql_del_sample_microbiome(self):
        sql = 'DELETE FROM sample_microbiome WHERE kitid = %s AND version = %s'
        return sql
    
    @property
    def sql_del_analysis_version(self):
        sql = 'DELETE FROM analysis_version WHERE kitid = %s AND version = %s'
        return sql

    def del_sample(self, p_name, p_birth, p_gender):
        self.execute_sql(self.sql_del_sample, (p_name, p_birth, p_gender))

    def del_sample_analysis(self, p_kitid):
        self.execute_sql(self.sql_del_sample_analysis, p_kitid)

    def del_sample_point(self, p_kitid, p_version):
        self.execute_sql(self.sql_del_sample_point, (p_kitid, p_version))

    def del_sample_microbiome(self, p_kitid, p_version):
        self.execute_sql(self.sql_del_sample_microbiome, (p_kitid, p_version))

    def del_analysis_version(self, p_kitid, p_version):
        self.execute_sql(self.sql_del_analysis_version, (p_kitid, p_version))

    @property
    def sql_select_sample(self):
        sql = 'SELECT * FROM sample WHERE client_name = %s AND client_birth = %s AND client_sex = %s'
        return sql

    @property
    def sql_select_sample_analysis(self):
        sql = 'SELECT * FROM sample_analysis WHERE kitid = %s'
        return sql

    @property
    def sql_select_sample_point(self):
        sql = 'SELECT * FROM sample_point WHERE kitid = %s AND version = %s'
        return sql

    @property
    def sql_select_sample_microbiome(self):
        sql = 'SELECT * FROM sample_microbiome WHERE kitid = %s AND version = %s'
        return sql

    @property
    def sql_select_analysis_version(self):
        sql = 'SELECT * FROM analysis_version WHERE kitid = %s AND version = %s'
        return sql

    def select_sample(self, p_name, p_birth, p_gender):
        """
        
        :param p_name: client_name
        :param p_birth: Birth
        :param p_gender: Sex
        :return:
        """
        data = self.execute_sql_for_select(self.sql_select_sample, (p_name, p_birth, p_gender))
        return data

    def select_sample_analysis(self, p_kitid):
        """
        
        :param p_kitid:
        :param p_version:
        :return:
        """
        data = self.execute_sql_for_select(self.sql_select_sample_analysis, p_kitid)
        return data

    def select_sample_point(self, p_kitid, p_version):
        """
        
        :param p_kitid:
        :param p_version:
        :return:
        """
        data = self.execute_sql_for_select(self.sql_select_sample_point, (p_kitid, p_version))
        return data

    def select_sample_microbiome(self, p_kitid, p_version):
        """
        
        :param p_kitid:
        :param p_version:
        :return:
        """
        data = self.execute_sql_for_select(self.sql_select_sample_microbiome, (p_kitid, p_version))
        return data

    def select_analysis_version(self, p_kitid, p_version):
        """
        
        :param p_kitid:
        :param p_version:
        :return:
        """
        data = self.execute_sql_for_select(self.sql_select_analysis_version, (p_kitid, p_version))
        return data


class DBDownloadSite(DBbase):
    def __init__(self, kargs):
        try:
            self.conn = pymysql.connect(
                host=kargs['ip'],
                user=kargs['user'],
                password=kargs['password'],
                db=kargs['db'],
                port=kargs['port'],
                charset='utf8'
            )
        except pymysql.err.OperationalError as err:
            secho('Error: Download Site DB에 접근할 수 없습니다.', fg='red', blink=True)
            echo(err)
        else:
            secho('>>> Download Site DB 연결', fg='cyan')
        try:
            self.cursor = self.conn.cursor()
        except Exception as err:
            secho('Error: (Download Site DB) Cursor 생성에 문제가 있습니다.', fg='red', blink=True, err=True)
            echo(err, err=True)
            secho('----> 개발자에게 문의하세요.', fg='magenta', err=True)

    def close(self):
        self.conn.close()
        secho('>>> Microbe&Me DB 연결 해제', fg='cyan')

    def read_csv_excutemany(self, p_sql, p_csv_file):
        return self._read_csv_excutemany_for_mysql(p_sql, p_csv_file)

    def execute_sql(self, p_sql, p_data):
        return self._execute_sql_for_mysql(p_sql, p_data)