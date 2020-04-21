# ----------------------------------------------------------------------------------------------------------------------
# 888b     d888          888             888b     d888 d8b
# 8888b   d8888          888             8888b   d8888 Y8P
# 88888b.d88888          888             88888b.d88888
# 888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#  _   _       _____
# | |_| |_ ___|     |_ _ ___ ___
# |  _|   | -_|   --| | | . |_ -|
# |_| |_|_|___|_____|___|  _|___|
#                       |_|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------

__author__ = 'JungWon Park(KOSST)'
__version__ = '1.0.0'

from SpoON.util import run_cmd, check_run_cmd


def get_make_2d_pcoa_cmd(p_pc, p_metadata, p_out_path):
    cmd = \
        'make_2d_plots.py ' \
        '-i {pc} ' \
        '-m {metadata} ' \
        '-o {out_path}'.format(
            pc=p_pc,
            metadata=p_metadata,
            out_path=p_out_path,
        )
    return cmd


def run_make_2d_pcoa(p_pc, p_metadata, p_out_path):
    cmd = get_make_2d_pcoa_cmd(p_pc, p_metadata, p_out_path)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': '2D PCoA 완료',
            'false_meg': 'make_2d_plots',
        }, p_exit=False, p_stdout=False,
    )
    if run.returncode == 0:
        return True
    else:
        return False
