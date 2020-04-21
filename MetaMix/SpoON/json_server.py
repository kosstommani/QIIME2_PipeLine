#!/crystal/Tools/Amplicon_MetaGenome/QIIME2/Miniconda3/bin/python3
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
__version__ = '1.0.0'

import socket
import sys
from click import echo
from urllib import request
import ssl

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server_address = ('172.19.85.10', 8877)  # denovo07
echo('staring up on %s port %s' % server_address)
server.bind(server_address)
server.listen(10)

while True:
    echo('waiting for a connection', file=sys.stderr)
    connection, client_address = server.accept()
    try:
        echo('connection from {}'.format(client_address))

        while True:
            try:
                lims_json = connection.recv(1024)
            except ConnectionResetError:
                echo('ConnectionResetError 발생')
                continue
            except Exception:
                echo('Exception Error 발생')
                continue
            else:
                if lims_json:
                    echo('received {}'.format(lims_json))
                    response = request.urlopen(lims_json.decode('utf-8'), context=ssl._create_unverified_context())
                    json_text = response.read()

                    if json_text:
                        echo('sending data back to the client')
                        connection.sendall(json_text, socket.MSG_EOR)
                        break
                else:
                    echo('no data from {}'.format(client_address))
                    break
    finally:
        connection.close()
