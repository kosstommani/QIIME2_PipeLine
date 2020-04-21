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


def get_beta_diversity_through_plots_cmd(p_biom, p_metadata, p_out_dir, p_tre, p_no_3d):
    if p_no_3d is True:
        plot_3d = '--suppress_emperor_plots'
    else:
        plot_3d = ''
    cmd = \
        'beta_diversity_through_plots.py ' \
        '-i {biom} ' \
        '-m {metadata} ' \
        '-o {out_dir} ' \
        '-t {tre} ' \
        '{plot_3d}'.format(
            biom=p_biom,
            metadata=p_metadata,
            out_dir=p_out_dir,
            tre=p_tre,
            plot_3d=plot_3d,
        )
    return cmd


def run_beta_diversity_through_plots(p_biom, p_metadata, p_out_dir, p_tre, p_no_3d):
    cmd = get_beta_diversity_through_plots_cmd(p_biom, p_metadata, p_out_dir, p_tre, p_no_3d)
    run = run_cmd(cmd)
    check_run_cmd(
        {
            'run': run,
            'true_meg': 'Beta Diversity & 3D PCoA 완료',
            'false_meg': 'beta_diversity_through_plots'
        }, p_exit=False, p_stdout=False,
    )
    if run.returncode == 0:
        return True
    else:
        return False
