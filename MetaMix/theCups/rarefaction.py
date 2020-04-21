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


def get_alpha_rarefaction_cmd(p_biom, p_metadata, p_out_dir, p_tre, p_max_depth=None):
    """

    :param p_biom:
    :param p_metadata:
    :param p_out_dir:
    :param p_tre:
    :param p_max_depth:
    :return:
    """
    if p_max_depth is None:
        max_rare_depth = ''
    else:
        max_rare_depth = f'-e {p_max_depth}'
    cmd = \
        'alpha_rarefaction.py ' \
        '-i {biom} ' \
        '-m {metadata} ' \
        '-o {out_dir} ' \
        '-t {tre} ' \
        '-f ' \
        '{max_rare_depth}'.format(
            biom=p_biom,
            metadata=p_metadata,
            out_dir=p_out_dir,
            tre=p_tre,
            max_rare_depth=max_rare_depth,
        )
    return cmd


def run_alpha_rarefaction(p_kargs):
    """

    :param p_kargs:
                biom:
                out_dir:
                metadata:
                tre:
                max_rare_depth:
    :return:
    """
    cmd = get_alpha_rarefaction_cmd(p_kargs['biom'], p_kargs['metadata'],
                                    p_kargs['out_dir'], p_kargs['tre'],
                                    p_kargs['max_rare_depth'])
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'Alpha Rarefaction 완료',
            'false_meg': 'alpha_rarefaction',
        }, p_exit=False, p_stdout=False)
    if run.returncode == 0:
        return True
    else:
        return False
