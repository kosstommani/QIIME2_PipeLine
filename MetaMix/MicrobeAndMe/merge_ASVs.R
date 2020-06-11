# ----------------------------------------------------------------------------------------------------------------------
#                                      888b     d888          888             888b     d888 d8b
#                                      8888b   d8888          888             8888b   d8888 Y8P
#                                      88888b.d88888          888             88888b.d88888
# 88888b.d88b.   8888b.  88888b.d88b.  888Y88888P888  .d88b.  888888  8888b.  888Y88888P888 888 888  888
# 888 "888 "88b     "88b 888 "888 "88b 888 Y888P 888 d8P  Y8b 888        "88b 888 Y888P 888 888 `Y8bd8P'
# 888  888  888 .d888888 888  888  888 888  Y8P  888 88888888 888    .d888888 888  Y8P  888 888   X88K
# 888  888  888 888  888 888  888  888 888   "   888 Y8b.     Y88b.  888  888 888   "   888 888 .d8""8b.
# 888  888  888 "Y888888 888  888  888 888       888  "Y8888   "Y888 "Y888888 888       888 888 888  888
#  ____  _____ ____  _____ ___
# |    \|  _  |    \|  _  |_  |
# |  |  |     |  |  |     |  _|
# |____/|__|__|____/|__|__|___|
#
# 작성자 : 박정원(JungWon Park, KOSST)
# ----------------------------------------------------------------------------------------------------------------------
# Title     : Merge ASVs
# Objective : DADA2 의 결과(ASVs)들을 하나의 테이블로 통합
# Created by: 박정원
# Created on: 2020-03-31


library(argparse)
library(crayon)


parser <- ArgumentParser(description='Microbe & Me - Merge DADA2 ASVs')
parser$add_argument('rds_dir', action='store', help='RDS파일들이 있는 디렉터리')
parser$add_argument('--rds_name',
                    default='seqtab_noChimera.rds',
                    choices=c('seqtab_noChimera.rds', 'seqtab.rds'))
args <- parser$parse_args(commandArgs(TRUE))


### 함수 ###
errQuit <- function(mesg, status=1) { message('Error: ', mesg); q(status=status) }
errMesg <- function(mesg) message(red('Error: ', mesg))
successText <- function(mesg) cyan('>>> ', mesg)
infoText <- function(mesg) green('Info: ', mesg)
getN <- function(x) sum(getUniques(x))
fileCat <- function(..., file='', sep='', fill=FALSE, labels=NULL, append=FALSE) {
  cat(..., file=file, sep=sep, fill=fill, labels=labels, append=append)
  cat(..., sep=sep, fill=fill, labels=labels, append=append)
}

out_path <- dirname(args$rds_dir)
log.name <- 'merge_ASVs.log'
log.file <- file.path(out_path, log.name)

# RDS파일 목록
rds.files <- sort(list.files(args$rds_dir, args$rds_name, full.names=TRUE))
samples.name <- sapply(strsplit(basename(rds.files), paste0('_', args$rds_name)), `[`, 1)
names(rds.files) <- samples.name
fileCat(infoText('파일 개수\n'), file=log.file, append=TRUE)
fileCat(length(rds.files), file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat(infoText('파일명 개수\n'), file=log.file, append=TRUE)
fileCat(length(samples.name), file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat(successText('RDS 파일 목록화 완료\n'), file=log.file, append=TRUE)

# RDS 파일 읽기
rds <- vector('list', length(samples.name))
names(rds) <- samples.name
for (sample in samples.name) {
  seqtab.noChimera <- readRDS(rds.files[[sample]])
  rownames(seqtab.noChimera) <- sample
  rds[[sample]] <- seqtab.noChimera
}
fileCat(successText('RDS 파일 읽기 완료\n'), file=log.file, append=TRUE)

# DADA2 라이브러리
library(dada2)
fileCat('--- 라이브러리 버전 ---', '\n', file=log.file, append=TRUE)
fileCat('DADA2       : ', as.character(packageVersion('dada2')), '\n', file=log.file, append=TRUE)
fileCat('Rcpp        : ', as.character(packageVersion('Rcpp')), '\n', file=log.file, append=TRUE)
fileCat('RcppParallel: ', as.character(packageVersion('RcppParallel')), '\n', file=log.file, append=TRUE)
fileCat('-----------------------', '\n', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)

# DADA2 디렉터리 생성
dada2.dir <- file.path(out_path, 'R_DADA2')
dir.create(dada2.dir, mode='0770')

# ASVs 통합
seqtab_all <- mergeSequenceTables(tables=rds)
fileCat(successText('mergeSequenceTables 완료\n'), file=log.file, append=TRUE)
seqtab_all.rds_name <- paste0('all', '_', args$rds_name)
seqtab_all.rds_file <- file.path(dada2.dir, seqtab_all.rds_name)
saveRDS(seqtab_all, seqtab_all.rds_file)
fileCat(successText(paste(seqtab_all.rds_name, '저장 완료\n')))
trans_st.all <- t(seqtab_all)

# FASTA 파일 생성
dada2_output.fasta <- file.path(dada2.dir, 'all_ASVs.fasta')
asv_id <- paste0('ASV', seq(1, length(colnames(seqtab_all))))
uniquesToFasta(seqtab_all, dada2_output.fasta, id=asv_id)
fileCat(successText(paste(dada2_output.fasta, '저장 완료\n')), file=log.file, append=TRUE)

# Tables 파일 생성
col.names <- samples.name
col.names[[1]] <- paste0('#ASVs_ID\t', col.names[[1]])
row.names <-  paste0('ASV', seq(1, length(rownames(trans_st.all))))
dada2_output.tsv <- file.path(dada2.dir, 'all_ASVs.tsv')
write.table(trans_st.all, dada2_output.tsv,
            sep='\t', row.names=row.names, col.names=col.names, quote=FALSE)
fileCat(successText(paste(dada2_output.tsv, '저장 완료\n')), file=log.file, append=TRUE)
fileCat(successText('merge_ASVs.R 완료\n'), file=log.file, append=TRUE)
