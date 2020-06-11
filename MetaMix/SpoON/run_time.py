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
__version__ = '1.0.1'

from click import echo
import time
import string
import random

DATE_FMT = '%Y-%m-%d %a %H:%M:%S'
STORY = random.choice(['cafe', 'spiderman'])


def compute_run_time(p_start, p_end):
    run_time = p_end - p_start
    hour = int(run_time // 60 // 60)
    minute = int(run_time // 60 % 60)
    second = int(run_time % 60)
    return '{0}h {1}m {2}s'.format(hour, minute, second)


def print_prema_run_time(p_start, p_main, p_integ, p_multiqc, p_end):
    global DATE_FMT, STORY
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    main_copy_time = p_main['copy_run_time']
    main_trim_time = p_main['trim_run_time']
    main_fastqc_time = p_main['fastqc_run_time']
    main_excel_time = p_main['excel_run_time']
    total_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    if p_integ != 'None':
        integ_copy_time = list()
        integ_trim_time = list()
        integ_fastqc_time = list()
        integ_excel_time = list()
        for ele in p_integ:
            integ_copy_time.append(ele['copy_run_time'])
            integ_trim_time.append(ele['trim_run_time'])
            integ_fastqc_time.append(ele['fastqc_run_time'])
            integ_excel_time.append(ele['excel_run_time'])
    else:
        pass

    if p_integ != 'None':
        integ_meg = """- other order -
 copy   : {copy} 
 fastqc : {fastqc}
 trim   : {trim}
 excel  : {excel}""".format(copy=','.join(integ_copy_time),
                            fastqc=','.join(integ_fastqc_time),
                            trim=','.join(integ_trim_time),
                            excel=','.join(integ_excel_time))
    else:
        integ_meg = '- other order -\n None'

    if STORY == 'cafe':
        msg = """
             -/        
            y/         주문시간 : {start}
           -h`
           o/
       ```.h.          - MetaMix 제조 시간 -
    ``` ``++````       
  ``  ```-h-.....       PreMA - main order
  .` `...os-...`..      copy    : {copy} 
 .--:./oshs+/::---.     fastqc  : {fastqc}
 `:/::::////:////-`     trim    : {trim}   
  -hhdddhhdhhhddy       excel   : {excel}   
   ydddmddmdddmm/       MultiQC : {multiqc} 
   odmmmmmmmmmmm.       integ. : 별도(other order 참조)
   :mmmmmmNNmmmd        other  : 생략
   `mmmNNNNNNNNs        ==================   
    ymNNNNNNNNN/        total  : {total}
    +mNNNNNNNNN.
    `dNNNNNNNNs
     `/yhddh+-         완성시간 : {end}

{integ}     
""".format(start=start_time,
           copy=main_copy_time,
           fastqc=main_fastqc_time,
           trim=main_trim_time,
           excel=main_excel_time,
           multiqc=p_multiqc,
           total=total_time,
           end=end_time,
           integ=integ_meg)
        echo(msg)
    elif STORY == 'spiderman':
        spiderman = '''
    __                                                      
     /  l                                                     
   .'   :               __.....__..._  ____                   
  /  /   \\          _.-" $$SSSSSS$$SSSSSSSSSp.                
 (`-: .qqp:    .--.'  .p.S$$$$SSSSS$$$$$$$$SSSSp.             
  """yS$SSSb,.'.g._\\.SSSSS^^""       `S""^^$$$SSSb.
    :SS\\$S\\$\\$\\$\\$SSSSS^"""-. _.        `.   "^$$$SSb._.._     
    SSS$$S$$SSP^/       `.               \\     "^$SSS$$SSb.   
    :SSSS$SP^" :          \\  `-.          `-      "^TSS$$SSb  
     $$$$S'    ;          db               ."        TSSSSSS$,
     :$$P      ;$b        $ ;    (        /   .-"   .S$$$$$$$;
       ;-"     :$ ^s.    d' $           .g. .'    .SSdSSSS$P" 
      /     .-  T.  `b. 't._$ .- / __.-j$'.'   .sSSSdP^^^'    
     /  /      `,T._.dP   "":'  /-"   .'       TSSP'          
    :  :         ,""       ; .'    .'      .-""              
   _J  ;         ; `.      /.'    _/    \\.-"                  
  /  "-:        /"--.b-..-'     .'       ;                    
 /     /  ""-..'            .--'.-'/  ,  :                    
:S.   :     dS$ bug         `-i" ,',_:  _ \\                   
:S$b  '._  dS$;             .'.-"; ; ; j `.l                  
 TS$b          "-._         `"  :_/ :_/                       
  `T$b             "-._                   스파이더맨 등장 : {start}                    
    :S$p._             "-.                                    
     `TSSS$ "-.     )     `.              - 스파이더맨 변신 시간 -                      
        ""^--""^-. :        \\                                 
                  ";         \\             PreMA                    
                  :           `._           copy    : {copy}                 
                  ; /    \\ `._   ""---.     fastqc  : {fastqc}                     
                 / /   _      `.--.__.'     trim    : {trim}                  
                : :   / ;  :".  \\           excel   : {excel}                     
                ; ;  :  :  ;  `. `.         MultiQC : {multiqc}                  
               /  ;  :   ; :    `. `.       integ. : 별도(other order 참조)                  
              /  /:  ;   :  ;     "-'       other  : 생략                  
             :_.' ;  ;    ; :               ==================                     
                 /  /     :_l               total  : {total}                  
                 `-'                     
                                          스파이더맨 퇴장 : {end}
{integ}
'''.format(start=start_time,
           copy=main_copy_time,
           fastqc=main_fastqc_time,
           trim=main_trim_time,
           excel=main_excel_time,
           multiqc=p_multiqc,
           total=total_time,
           end=end_time,
           integ=integ_meg)
        echo(spiderman)


def print_legnth_trim_time(p_start, p_end):
    global DATE_FMT
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    trim_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    msg = """
            .                        .
            Wm      - Bug Time -    ~#`
            )#`                     Hd
             Q\\                    `#~
             I#`                   o#`
             -g*   -E_       hv   `Bv
 uv`         `I#_  `gw      -@!   H#`          =y_
 'Tgy'        :@k   i#`     y0   :@o         )0P:
    ~$0v`      M#~  i@L    -@D` `R#.      !GQ]`
      :y$Ir`   ~@#' Z@*.mI^ #@. m@X    "lDM*`
         )D#z-  <Q#y^B##@@@B@Tx$#V  `vQ#y`
           -c@#M]^rkB@@@@@@@@#Zx;vUQ@d~
             `!]UM68#@@@##@@@@BDb3Vr.
              `*TIZg@@@@@@@@@@BO3Vv,
           _TB@Qmr_~M@#@@@@@#Q)_!c$@#G~
        `T8#w=. `X#6y0@@@@@@#KI#E!  !vg#K:
      :3$y*    :B@]i@@@@@@@@@@D=Q@]    _uZ0)`
    !gg)      _QD!.@@@@@@@@@@@@i`y@L      "WBY
  *EU.       "B0` _@@@@@@@@@@@@]  T@]``     `vgu`
 P3_        -dV`   ]@@@@@@@@@@0`   *Q^`       `Y9:
           -05      =O@@@@@@BY`     v@*
           d0'        '*RMy_         T#-
          V$                          YQ'
         x#.                           UE
         ,`                             :
         나타난 시간 : {start}
         
         - 움직인 시간
           Length Trim : {trim}
         
         사라진 시간 : {end}
""".format(start=start_time, trim=trim_time, end=end_time)
    echo(msg)


def print_cofi_flash_run_time(p_start, p_assembly, p_end):
    global DATE_FMT, STORY
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    flash_time = p_assembly['flash_run_time']
    filtering_time = p_assembly['filter_run_time']
    total_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    if STORY == 'cafe':
        msg = """
             -/        
            y/         주문시간 : {start}
           -h`
           o/
       ```.h.          
    ``` ``++````       - MetaMix 제조 시간 -
  ``  ```-h-.....       
  .` `...os-...`..      CoFI - FLASH
 .--:./oshs+/::---.     
 `:/::::////:////-`     FLASH     : {flash} 
  -hhdddhhdhhhddy       Length   
   ydddmddmdddmm/       Filtering : {filter}   
   odmmmmmmmmmmm.       other     : 생략
   :mmmmmmNNmmmd        =====================   
   `mmmNNNNNNNNs        total     : {total}
    ymNNNNNNNNN/        
    +mNNNNNNNNN.
    `dNNNNNNNNs
     `/yhddh+-         완성시간 : {end}
""".format(start=start_time,
           flash=flash_time,
           filter=filtering_time,
           total=total_time,
           end=end_time)
        echo(msg)
    elif STORY == 'spiderman':
        spiderman = """
                   ,,,,
             ,;) .';;;;',
 ;;,,_,-.-.,;;'_,|I\\;;;/),,_
  `';;/:|:);[ ;;;|| \\;/ /;;;\\__                  스파이더맨 등장 : {start} 
      L;/-';/ \\;;\\',/;\\/;;;.') \\
      .:`''` - \\;;'.__/;;;/  . _'-._             - 스파이더맨 변신 시간 - 
    .'/   \\     \\;;;;;;/.'_7:.  '). \\_        
  .''/     | '._ );][;//.'    '-:  '.,L            CoFI - FLASH
.'. /       \\  ( |;;;/_/         \\._./;\\   _,   
 . /        |\\ ( /;;/_/             ';;;\\,;;_,     FLASH     : {flash}
. /         )__(/;;/_/                (;;'''''     Length
 /        _;:':;;;;:';-._             );           Filtering : {filter}
/        /   \\  `'`   --.'-._         \\/\\          other     : 생략
       .'     '.  ,'         '-,                   =====================
      /    /   r--,..__       '.\\                  total     : {total}
    .'    '  .'        '--._     ]
    (     :.(;>        _ .' '- ;/                스파이더맨 퇴장 : {end}
    |      /:;(    ,_.';(   __.'
     '- -'"|;:/    (;;;;-'--'
           |;/      ;;(
KOSST      ''      /;;|          Amplicon Meta Team
                   \\;;|          윤선미, 박정원, 임세은, 김종범, 김상겸
                    \\/
""".format(start=start_time,
           flash=flash_time,
           filter=filtering_time,
           total=total_time,
           end=end_time)
        echo(spiderman)


def print_cofi_legnth_filter_run_time(p_start, p_run_time, p_end):
    global DATE_FMT
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    flash_excel_time = p_run_time['flash_excel_run_time']
    filtering_time = p_run_time['filter_run_time']
    total_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    msg = """
             -/        
            y/         주문시간 : {start}
           -h`
           o/
       ```.h.          
    ``` ``++````       - MetaMix 제조 시간 -
  ``  ```-h-.....       
  .` `...os-...`..      CoFI - FLASH
 .--:./oshs+/::---.     
 `:/::::////:////-`     FLASH Excel : {flash} 
  -hhdddhhdhhhddy       Length   
   ydddmddmdddmm/       Filtering   : {filter}   
   odmmmmmmmmmmm.       
   :mmmmmmNNmmmd        =====================   
   `mmmNNNNNNNNs        total     : {total}
    ymNNNNNNNNN/        
    +mNNNNNNNNN.
    `dNNNNNNNNs
     `/yhddh+-         완성시간 : {end}
""".format(start=start_time,
           flash=flash_excel_time,
           filter=filtering_time,
           total=total_time,
           end=end_time)
    echo(msg)


def print_cofi_cd_hit_run_time(p_start, p_cd_hit, p_end):
    global DATE_FMT, STORY
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    flash_time = p_cd_hit['cd-hit_run_time']
    filtering_time = p_cd_hit['otus_file_run_time']
    total_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    if STORY == 'cafe':
        msg = """
             -/        
            y/         주문시간 : {start}
           -h`
           o/
       ```.h.          
    ``` ``++````       - MetaMix 제조 시간 -
  ``  ```-h-.....       
  .` `...os-...`..      CoFI - CD-HIT-OTU
 .--:./oshs+/::---.     
 `:/::::////:////-`     CD-HIT-OTU     : {cd_hit} 
  -hhdddhhdhhhddy       pick_otus.txt &    
   ydddmddmdddmm/       otus_rep.fasta : {otus}   
   odmmmmmmmmmmm.       other     : 생략
   :mmmmmmNNmmmd        ======================   
   `mmmNNNNNNNNs        total     : {total}
    ymNNNNNNNNN/        
    +mNNNNNNNNN.
    `dNNNNNNNNs
     `/yhddh+-         완성시간 : {end}
""".format(start=start_time,
           cd_hit=flash_time,
           otus=filtering_time,
           total=total_time,
           end=end_time)
        echo(msg)
    elif STORY == 'spiderman':
        spiderman = '''
                                   z11
                                   z@d
                                    0@@.                z@.
                z1jd01j".          j11j
                                  110jd@z@000@jjjjjz@1jd.
                                zj@j @jd    "zzzzzj@0dz0
jdd@j1""                    "zj@d@.   00z       .jjd1 @@
  "z11j@00jjjzz   .zzzz"zjj00@1z.   .111@01jz".zj10.  110
         """j1j1jzzz10j11z""      .zj@0@0@j1j@1@1@j   .@1@
             .j@@""j1000@11111111@0@11@1j@00@@@01@0     1d0"
              "@d.     jjd11@0@@dd@1@dd1jd0@1d0@zj10z    z@@0z
               z00      00"   @1d110@d0@0@1@@j@j0@1@d@11111@1@@111@@@
                1@j     j0.   "@d   "j0@01j0jjjd1@@dddjzzzzz0@dj"
                j@@     j0.   11d"j10@@@1d@1dz10?zz00.     z@0j
                j@d"   z?0. .j@100@j   j11 jj0?@1"@d"     .@@"
                z01   z1@jj@@d@1.     "@0z  "z1@0d1j      011
               1z0. "@0d@@@110z"      d@d    .@j@d0j     j1@
              0@0 .j@1jjj"  "zj@11"  .d?zzj10@@@0@jd     0@z
            "??d1j@@0@1.       .z@@010@00@1".    @@@0.   @d
           1@000@11j100dd0@11z.    z1@d0111@01111j0@0d   @d
        .j@0@?0@zzz.      "zz@d1"   1ddd0@""""""zzj@00d. d?
     zj0d011111jzj0000@11z.   z0d0" 0z10.           .100"0dj
 "z10@@j.             .zz10@1   "0d?zj1               j@dd@@
@@1jz                     z10dj   00@1    ..zj"1@zzzj".j@000z
                            zdd@  11@  "z001jjjjjjzzj1100111d
CD-HIT-OTU  : {cd_hit}        .1??j0@1z0d@.               zzzj@
                               j?d1010@.                   .1@@"
                                jz00j                        jd?"
                                jj@                           .d01
                               .11j                             j@0"
                               j@0                               zj1
                               j0j

    거미줄 설치 시작: {start}
     - 스파이더맨 작업 시간 -
       CoFI - CD-HIT-OTU
       CD-HIT-OTU     : {cd_hit}
       pick_otus.txt &
       otus_rep.fasta : {otus}
       other          : 생략
       ======================
       total      : {total}
    거미줄 완성시간 : {end}
'''.format(start=start_time,
           cd_hit=flash_time,
           otus=filtering_time,
           total=total_time,
           end=end_time)
        echo(spiderman)


def print_flash_hit_run_time(p_start, p_flash, p_cd_hit, p_closed):
    global DATE_FMT, STORY
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    total_run_time = compute_run_time(p_start, p_closed)
    end_time = time.strftime(DATE_FMT, time.localtime(p_closed))
    if STORY == 'cafe':
        msg = """
             -/        
            y/         주문시간 : {start}
           -h`
           o/
       ```.h.          - MetaMix 제조 시간 -
    ``` ``++````       
  ``  ```-h-.....        CoFI - FLASH_HIT
  .` `...os-...`..       FLASH          : {flash}
 .--:./oshs+/::---.      Length 
 `:/::::////:////-`      Filtering      : {filter}
  -hhdddhhdhhhddy        CD-HIT-OTU     : {cd_hit}
   ydddmddmdddmm/        pick_otus.txt &
   odmmmmmmmmmmm.        otus_rep.fasta : {otus}
   :mmmmmmNNmmmd         other          : 생략
   `mmmNNNNNNNNs         =========================
    ymNNNNNNNNN/         total          : {total}
    +mNNNNNNNNN.         
    `dNNNNNNNNs          
     `/yhddh+-         완성시간 : {end}
""".format(start=start_time,
           flash=p_flash['flash_run_time'],
           filter=p_flash['filter_run_time'],
           cd_hit=p_cd_hit['cd-hit_run_time'],
           otus=p_cd_hit['otus_file_run_time'],
           total=total_run_time,
           end=end_time)
        echo(msg)

    elif STORY == 'spiderman':
        spiderman = '''
                             .-"""-.    __
                            /       \\.q$$$$p.
                         __:  db     $$$$$$$$b.
                  _._.-""  :  $"b.   :$$$S$$$$$b
                .'   "-.  "   T. `b d$$$S$S$$$$P^.              .-,
    .qp.       :        `.     TsP' TSSS$P'S$P'   `.            `dP
 ,q$$$$$b      ;b         \\  '.     /"T$P  :P       `.      __ dP_,
 :$$$$SS$b_    $$b.  __.   ;_  `-._/   Y    \\         `.   ( dP".';
 $$$$$$$S$$b.  :S$$$p._    $$$$p./      ;    "-._       `--dP .'  ;
:$$$$P^^TS$$$b  SS$$$$$$   'T$$$$b.     ;        ""--.   dP  /   /
$$$$P    :$$$$bd$SSS$$$;\\  . "^$$$$b___/            __\\dP .-"_.-"
$$$P     $$$$$$b`T$SSS$  "-.J   "^T$/  ""-._       ( dP\\ /   /
:$$      ; T$$$$b.`T$$;     d$+.     ""-.   ""--.._dP-, `._."
 T;     :   T$$$$$b.`^'   _d$P .$p.___   "         \\`-'
  `.    ;    T$$$$$$b._.dS$$$ :$$$b"--..__..---g,   \\
    `. :      $$$$$$$$$$S$$$P\\ TP^"\\       ,-dP ;    ;
      \\;   .-'$$$$$$$SSSP^^"  \\     `._,-.-dP-' |    ;
       :     :"^^"" """        `.   `._:'.`.\\   :    ;\\
        \\  , :              bug  "-. (,j\\ ` /   ;\\(// \\\\
         `:   \\                     "dP__.-"    '-\\\\   \\;
           \\   :                .--dP,             \\;
            `--'                `dP`-\'
                              .-j
                              `-:_
                                 \\)
                                  `--\'        
    스파이더맨 등장 : {start}

    - 스파이더맨 변신 시간 -
      CoFI - FLASH_HIT
      FLASH          : {flash}
      Length 
      Filtering      : {filter}
      CD-HIT-OTU     : {cd_hit}
      pick_otus.txt &
      otus_rep.fasta : {otus}
      other          : 생략
    =========================
    total          : {total}
    
    스파이더맨 퇴장 : {end}                      
'''.format(start=start_time,
           flash=p_flash['flash_run_time'],
           filter=p_flash['filter_run_time'],
           cd_hit=p_cd_hit['cd-hit_run_time'],
           otus=p_cd_hit['otus_file_run_time'],
           total=total_run_time,
           end=end_time)
        echo(spiderman)


def print_cofi_closed_otu_run_time(p_start, p_closed, p_end):
    global DATE_FMT
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    total_run_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    msg = """
                 -/        
                y/         주문시간 : {start}
               -h`
               o/
           ```.h.
        ``` ``++````       - MetaMix 제조 시간 -
      ``  ```-h-.....      
      .` `...os-...`..       CoFI - CLOSED
     .--:./oshs+/::---.      
     `:/::::////:////-`      
      -hhdddhhdhhhddy        Closed-Ref. OTU Picking    
       ydddmddmdddmm/        : {closed}
       odmmmmmmmmmmm.            
       :mmmmmmNNmmmd         Total
       `mmmNNNNNNNNs         : {total}
        ymNNNNNNNNN/
        +mNNNNNNNNN.
        `dNNNNNNNNs
         `/yhddh+-         완성시간 : {end}
    """.format(start=start_time,
               closed=p_closed,
               total=total_run_time,
               end=end_time)
    echo(msg)


def print_metamix_run_time(p_start, p_prema_end, p_pipeline, p_end):
    """

    :param p_start:
    :param p_prema_end:
    :param p_pipeline: kargs
    :param p_end:
    :return:
    """
    global DATE_FMT
    msg = string.Template('''
               ,        ,          주문시간: $start 
              /(        )`              
              \\ \\___   / |         - MetaMix 제조 시간 -
              /- _  `-/  '                        
             (/\\/ \\ \\   /\\           PreMA 
             / /   | `    \\            : $prema
             O O   ) /    |                     
             `-^--'`<     '          CoFI          
            (_.)  _  )   /             FLASH_HIT: $flash_hit
             `.___/`    /              ALIGNMENT: $alignment
               `-----' /               PHYLOGENY: $phylogeny   
    <----.     __ / __   \\              TAXONOMY: $taxonomy         
    <----|====O)))==) \\) /====              BIOM: $biom
    <----'    `--' `.__,' \\              theCups: $thecups
               |        |             Total     
                \\       /              : $total             
           ______( (_  / \\______        
         ,'  ,-----'   |        \\   완성시간: $end              
         `--{__________)        \\/    
      악마의 유혹 - MetaMix 카페
    ''')
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    prema_run_time = compute_run_time(p_start, p_prema_end)
    if isinstance(p_pipeline, dict):
        flash_hit_time = compute_run_time(p_pipeline['flash_hit_start_time'], p_pipeline['flash_hit_end_time'])
        alignment_time = compute_run_time(p_pipeline['flash_hit_end_time'], p_pipeline['alignment_end_time'])
        phylogeny_time = compute_run_time(p_pipeline['alignment_end_time'], p_pipeline['phylogeny_end_time'])
        taxonomy_time = compute_run_time(p_pipeline['phylogeny_end_time'], p_pipeline['taxonomy_end_time'])
        biom_time = compute_run_time(p_pipeline['taxonomy_end_time'], p_pipeline['biom_end_time'])
        thecups_time = compute_run_time(p_pipeline['biom_end_time'], p_pipeline['theCups_end_time'])
    elif isinstance(p_pipeline, list):
        flash_hit_time = list()
        alignment_time = list()
        phylogeny_time = list()
        taxonomy_time = list()
        biom_time = list()
        thecups_time = list()
        for ele in p_pipeline:
            flash_hit_time.append(compute_run_time(ele['flash_hit_start_time'], ele['flash_hit_end_time']))
            alignment_time.append(compute_run_time(ele['flash_hit_end_time'], ele['alignment_end_time']))
            phylogeny_time.append(compute_run_time(ele['alignment_end_time'], ele['phylogeny_end_time']))
            taxonomy_time.append(compute_run_time(ele['phylogeny_end_time'], ele['taxonomy_end_time']))
            biom_time.append(compute_run_time(ele['taxonomy_end_time'], ele['biom_end_time']))
            thecups_time.append(compute_run_time(ele['biom_end_time'], ele['theCups_end_time']))
    total_run_time = compute_run_time(p_start, p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))

    if isinstance(p_pipeline, dict):
        echo(msg.substitute({
            'start': start_time,
            'prema': prema_run_time,
            'flash_hit': flash_hit_time,
            'alignment': alignment_time,
            'phylogeny': phylogeny_time,
            'taxonomy': taxonomy_time,
            'biom': biom_time,
            'thecups': thecups_time,
            'total': total_run_time,
            'end': end_time
        }))
    elif isinstance(p_pipeline, list):
        echo(msg.substitute({
            'start': start_time,
            'prema': prema_run_time,
            'flash_hit': ', '.join(flash_hit_time),
            'alignment': ', '.join(alignment_time),
            'phylogeny': ', '.join(phylogeny_time),
            'taxonomy': ', '.join(taxonomy_time),
            'biom': ', '.join(biom_time),
            'thecups': ', '.join(thecups_time),
            'total': total_run_time,
            'end': end_time
        }))


def print_mam_metamix_run_time(p_start, p_pipeline, p_end):
    global DATE_FMT
    start_time = time.strftime(DATE_FMT, time.localtime(p_start))
    dada2_time = compute_run_time(p_pipeline['make_sample'], p_pipeline['dada2'])
    assign_taxonomy_time = compute_run_time(p_pipeline['dada2'], p_pipeline['assign_taxonomy'])
    biom_time = compute_run_time(p_pipeline['assign_taxonomy'], p_pipeline['biom'])
    alpha_diversity_time = compute_run_time(p_pipeline['biom'], p_pipeline['summarize_taxa'])
    score_time = compute_run_time(p_pipeline['summarize_taxa'], p_pipeline['score'])
    info_time = compute_run_time(p_pipeline['score'], p_pipeline['info'])
    insert_score_time = compute_run_time(p_pipeline['info'], p_pipeline['insert_score'])
    find_score_time = compute_run_time(p_pipeline['insert_score'], p_end)
    end_time = time.strftime(DATE_FMT, time.localtime(p_end))
    total_time = compute_run_time(p_start, p_end)
    K = '\033[1;32;40mK\033[m'
    O = '\033[1;32;40mO\033[m'
    S = '\033[1;32;40mS\033[m'
    T = '\033[1;32;40mT\033[m'
    mgs = f'''
,---,---,---,---,---,---,---,---,---,---,---,---,---,-------,
| ~ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0 | [ | ] | <-    |
|---'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-----|
| ->| | " | , | . | P | Y | F | G | C | R | L | / | = |  \\  |
|-----',--',--',--',--',--',--',--',--',--',--',--',--'-----|
| Caps | A | {O} | E | U | I | D | H | {T} | N | {S} | - |  Enter |
|------'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'-,-'--------|
|        | ; | Q | J | {K} | X | B | M | W | V | Z |          |
|------,-',--'--,'---'---'---'---'---'---'-,-'---',--,------|
| ctrl |  | alt |           \033[1;32;40mKOSST\033[m          | alt  |  | ctrl |
'------'  '-----'--------------------------'------'  '------'
키보드 입력 시작 시간: {start_time}
    DADA2          : {dada2_time}
    Assign Taxonomy: {assign_taxonomy_time}
    BIOM           : {biom_time}
    Alpha Diversity: {alpha_diversity_time}
    Score          : {score_time}
    Sample Info    : {info_time}
    Insert Score   : {insert_score_time}
    Find SCore     : {find_score_time}
    ===============================
    Total          : {total_time} 
입력 완료 시간: {end_time}
'''
    echo(mgs)

