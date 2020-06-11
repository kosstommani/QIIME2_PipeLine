# ----------------------------------------------------------------------------------------------------------------------
# 888b     d888          888             888b     d888 d8b
# 8888b   d8888          888             8888b   d8888 Y8P
# 88888b.d88888          888             88888b.d88888
# 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#  ____  _____             _           _
# |    \| __  |___ ___ ___| |_ ___ ___| |
# |  |  | __ -|  _| . |   |  _|  _| . | |
# |____/|_____|___|___|_|_|_| |_| |___|_|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '0.1'

from click import secho, echo, style
import csv
from pprint import pprint


class DBbase:
    def rollback(self):
        secho('>>> RollBack 진행', fg='red', blink=True, bold=True)
        self.cursor.close()
        echo('Cursor Close [ {} ]'.format(style('OK', fg='yellow')))
        self.conn.rollback()
        echo('RollBack [ {} ]'.format(style('OK', fg='yellow')))
        self.conn.close()
        echo('DB Connection Close [ {} ]'.format(style('OK', fg='yellow')))
        exit(1)

    def commit(self):
        self.conn.commit()
        secho('>>> DB COMMIT 완료', fg='cyan')

    def start_transaction(self):
        self.cursor.execute('START TRANSACTION')

    def _read_csv_excutemany_for_mysql(self, p_sql, p_csv_file):
        """
        CSV 파일을 읽은 후, 데이터를 DB에 삽입한다.
        트랜잭션을 실행하고, 삽입 중 에러가 발생하면 롤백한다.
        
        :param p_sql:
        :param p_csv_file:
        :return:
        """
        import pymysql
        with open(p_csv_file, 'r') as o_csv:
            o_csv_data = csv.reader(o_csv)
            l_csv_data = [tuple(x) for x in o_csv_data]
            try:
                self.cursor.executemany(p_sql, l_csv_data)
            except pymysql.err.IntegrityError as err:
                secho('Error: SQL IntegrityError 발생', fg='red')
                echo(err)
                echo('=' * 20)
                echo(p_sql)
                pprint(l_csv_data)
                echo('=' * 20)
                self.rollback()
            except pymysql.err.ProgrammingError as err:
                secho('Error: ProgrammingError 발생', fg='red')
                echo(err)
                echo('=' * 20)
                echo(p_sql)
                pprint(l_csv_data)
                echo('=' * 20)
                self.rollback(self.cursor)
            except Exception as err:
                secho('Error: Exception 발생', fg='red')
                echo(err)
                echo('=' * 20)
                echo(p_sql)
                pprint(l_csv_data)
                echo('=' * 20)
                self.rollback()

    def _execute_sql_for_mysql(self, p_sql, p_data):
        """

        :param p_sql:
        :param p_data:
        :return:
        """
        import pymysql
        try:
            self.cursor.execute(p_sql, p_data)
        except pymysql.err.IntegrityError as err:
            secho('Error: SQL IntegrityError 발생', fg='red')
            echo(err)
            echo('=' * 20)
            echo(p_sql)
            pprint(p_data)
            echo('=' * 20)
            self.rollback()
        except pymysql.err.ProgrammingError as err:
            secho('Error: ProgrammingError 발생', fg='red')
            echo(err)
            echo('=' * 20)
            echo(p_sql)
            pprint(p_data)
            echo('=' * 20)
            self.rollback()
        except Exception as err:
            secho('Error: Exception 발생', fg='red')
            echo(err)
            echo('=' * 20)
            echo(p_sql)
            pprint(p_data)
            echo('=' * 20)
            self.rollback()
