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


def get_closed_reference_otu_cmd(p_fasta, p_out_dir):
    parameter = '/crystal/Tools/Amplicon_MetaGenome/QIIME2_PipeLine/MetaMix/CoFI/closed_OTU_parameter.qiime_1.9.txt'
    cmd = \
        '/crystal/FLX/Tools/Python-2.7.10/bin/pick_closed_reference_otus.py ' \
        '-i {fasta} ' \
        '-o {out_dir} ' \
        '--parallel ' \
        '--jobs_to_start 10 ' \
        '--assign_taxonomy ' \
        '--parameter_fp {para}'.format(
            fasta=p_fasta,
            out_dir=p_out_dir,
            para=parameter)
    return cmd


def get_cmd_alpha_diversity_qiime1(p_biom, p_output):
    cmd = \
        '/crystal/FLX/Tools/Python-2.7.10/bin/alpha_diversity.py ' \
        '-i {biom} ' \
        '-m shannon,simpson ' \
        '-o {output}'.format(
            biom=p_biom,
            output=p_output)
    return cmd


def get_cmd_summarize_taxa_qiime1(p_biom, p_output):
    cmd = \
        '/crystal/FLX/Tools/Python-2.7.10/bin/summarize_taxa.py ' \
        '-i {biom} ' \
        '-o {output} ' \
        '-L 7'.format(
            biom=p_biom,
            output=p_output)
    return cmd
