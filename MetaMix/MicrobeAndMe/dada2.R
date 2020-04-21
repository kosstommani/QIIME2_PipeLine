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
# Version : 1.0
# ----------------------------------------------------------------------------------------------------------------------
# Title     : DADA2 (Microbe & Me)
# Objective : DADA2 실행
# Created by: 박정원
# Created on: 2020-03-31


library(argparse)
library(crayon)

parser <- ArgumentParser(description='DTC - Microbe&Me.')
parser$add_argument('Order_Number', action='store', help='수주번호.')
parser$add_argument('Analysis_Number', action='store', help='분석 작업 디렉터러리.')
parser$add_argument('Run_Dir', action='store', help='시료별 FASTQ 파일이 있는 Run Dir. ex) RawData/200120_CPWV8_2')
parser$add_argument('--analysis_base_path',
                    default='/garnet/Analysis/BI/AmpliconMetaGenome_MAM',
                    help='분석 작업 경로. 기본값: /garnet/Analysis/BI/AmpliconMetaGenome_MAM')
                    # default='/garnet/Tools/Amplicon_MetaGenome/MicrobeAndMe/',
                    # help='분석 작업 경로. 기본값: /garnet/Tools/Amplicon_MetaGenome/MicrobeAndMe')
parser$add_argument('--sample_name',
                    help='시료명. 미입력시 RawData/Run 디렉터리에 있는 시료명 디렉터리를 검색하여 검색되는 모든 시료들을 분석함.')
parser$add_argument('--out_path',
                    help='기본값: 자동생성(HN000XXX/Analysis_X). 결과 생성 경로.')
parser$add_argument('--multithread',
                    default=TRUE,
                    help=paste('기본값: TRUE.',
                               '병렬처리 여부를 선탬함.',
                               'TRUE: 사용가능한 threads 자동 선택.',
                               'FLASE: 미적용.',
                               '숫자: 해당 숫자만큼 thread 사용.', sep=' '))
parser$add_argument('--truncLenF', type='integer',
                    default=250,
                    help='기본값: 250.')
parser$add_argument('--truncLenR', type='integer',
                    default=200,
                    help='기본값: 200.')
parser$add_argument('--trimLeftF', type='integer',
                    default=0,
                    help='기본값: 0.')
parser$add_argument('--trimLeftR', type='integer',
                    default=0,
                    help='기본값: 0.')
parser$add_argument('--maxN', type='integer',
                    default=0,
                    help='기본값: 0.')
parser$add_argument('--maxEEF', type='integer',
                    default=5,
                    help='기본값: 5.')
parser$add_argument('--maxEER', type='integer',
                    default=5,
                    help='기본값: 5.')
parser$add_argument('--truncQ', type='integer',
                    default=2,
                    help='기본값: 2.')
parser$add_argument('--minOverlap', type='integer',
                    default=15,
                    help='기본값: 15.')
parser$add_argument('--chimera_method',
                    default='consensus', choices=c('consensus', 'pooled', 'per-sample'),
                    help='기본값: consensus.')
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

### Set Path ###
order_number.dir <- file.path(args$analysis_base_path, args$Order_Number)
rawdata_run.dir <- file.path(order_number.dir, args$Run_Dir)
analysis.dir <- file.path(order_number.dir, args$Analysis_Number)

### 분석 수주디렉터리: 수주번호 --> RawData --> run_dir 확인 ###
if (!dir.exists(rawdata_run.dir)) {
    errMesg('디렉터리가 존재하지 않습니다.')
    message(rawdata_run.dir)
    quit('no', 1)
}

### Analysis -> Sample Dir 생성
analysis_sample.dir <- file.path(analysis.dir, args$sample_name)
dir.create(analysis_sample.dir, mode='0770')
log.name <- 'run_dada2.log'
log.file <- file.path(analysis_sample.dir, log.name)
samples.name <- args$sample_name
run_sample.dir <- file.path(rawdata_run.dir, args$sample_name)

if (!dir.exists(run_sample.dir)) {
    errMesg('시료명 디렉터리가 존재하지 않습니다.')
    message(run_sample.dir)
    quit('no', 1)
}

### 분석 수주디렉터리: 수주번호 --> RawData --> run_dir --> 시료명 디렉터리 확인 ###
if (length(run_sample.dir) != 0) {
    read1.files <- list.files(run_sample.dir, '_1.fastp.fastq.gz$', full.names=TRUE)
    read2.files <- list.files(run_sample.dir, '_2.fastp.fastq.gz$', full.names=TRUE)
} else {
    errMesg('시료 디렉터리가 존재하지 않습니다.')
    message(run_sample.dir)
    quit('no', 1)
}

### 파일 확인 ###
if (length(read1.files) == 0) {
    errMesg('Read1 파일이 없습니다.')
    quit('no', 1)
}
if (length(read2.files) == 0) {
    errMesg('Read2 파일이 없습니다.')
    quit('no', 1)
}
if (length(read1.files) != length(read2.files)) {
    errMesg('Read1 과 Read2 의 파일 개수가 다릅니다.')
    message('Read1: ', length(read1.files))
    message('Read2: ', length(read2.files))
    quit('no', 1)
}


### Analysis_? 디렉터리 확인 ###
if (!dir.exists(analysis.dir)) {
    errMesg('Analysis_? 디렉터리가 없습니다.')
    message(analysis.dir)
    quit('no', 1)
}


#--------------------#
### DADA2 라이브러리 ###
#--------------------#
suppressWarnings(library(dada2))
#library(RcppParallel)
#setThreadOptions(numThreads=16)

### Set option ###
if (is.logical(args$multithread) == TRUE) {
    multithread <- args$multithread
} else {
    multithread <- as.numeric(args$multithread)
}
verbose <- TRUE
fileCat('--- 라이브러리 버전 ---', '\n', file=log.file, append=TRUE)
fileCat('DADA2       : ', as.character(packageVersion('dada2')), '\n', file=log.file, append=TRUE)
fileCat('Rcpp        : ', as.character(packageVersion('Rcpp')), '\n', file=log.file, append=TRUE)
fileCat('RcppParallel: ', as.character(packageVersion('RcppParallel')), '\n', file=log.file, append=TRUE)
fileCat('-----------------------', '\n', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
### DADA2 Option 출력 ###
fileCat(R.version$version.string, '\n', file=log.file)
fileCat(infoText('DADA2 Options'),'\n', file=log.file, append=TRUE)
fileCat('multithread: ', multithread, '\n', sep='', file=log.file, append=TRUE)
fileCat('verbose: ', verbose, '\n', sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat('--- filterAnd Trim Options ---\n', file=log.file, append=TRUE)
fileCat('truncLenF: ', args$truncLenF, '\n', sep='', file=log.file, append=TRUE)
fileCat('truncLenR: ', args$truncLenR, '\n', sep='', file=log.file, append=TRUE)
fileCat('trimLeftF: ', args$trimLeftF, '\n', sep='', file=log.file, append=TRUE)
fileCat('trimLeftR: ', args$trimLeftR, '\n', sep='', file=log.file, append=TRUE)
fileCat('maxN: ', args$maxN, '\n', sep='', file=log.file, append=TRUE)
fileCat('maxEEF: ', args$maxEEF, '\n', sep='', file=log.file, append=TRUE)
fileCat('maxEER: ', args$maxEER, '\n', sep='', file=log.file, append=TRUE)
fileCat('truncQ: ', args$truncQ, '\n', sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat('--- mergePairs ---\n', file=log.file, append=TRUE)
fileCat('minOverlap: ', args$minOverlap, '\n', sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat('--- Length Filter ---\n', file=log.file, append=TRUE)
fileCat('no_lenFilter: ', args$no_lenFilter, '\n', sep='', file=log.file, append=TRUE)
fileCat('lenFilterF: ', args$lenFilterF, '\n', sep='', file=log.file, append=TRUE)
fileCat('lenFilterR: ', args$lenFilterR, '\n', sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat('--- removeBimeraDenovo Options ---\n', file=log.file, append=TRUE)
fileCat('chimera_method: ', args$chimera_method, '\n', sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)


### DADA2 디렉터리 생성 ###
dada2.dir <- file.path(analysis_sample.dir, 'R_DADA2')
dir.create(dada2.dir, mode='0770')
dada2.filter_dir <- file.path(dada2.dir, 'Filter')
dir.create(dada2.filter_dir, mode='0770')


#----------------------------#
### Plot - Quality Profile ###
#----------------------------#
qPlot.name <- 'Quality_Profile_'
qPlot_R1.file <- file.path(dada2.dir, paste0(qPlot.name, 'R1.pdf'))
qPlot_R2.file <- file.path(dada2.dir, paste0(qPlot.name, 'R2.pdf'))
# Read1
pdf(qPlot_R1.file)
plotQualityProfile(read1.files)
dev.off()
fileCat(successText('Quality Profile Plot - Read1 완료\n'), file=log.file, append=TRUE)
# Read2
pdf(qPlot_R2.file)
plotQualityProfile(read2.files)
dev.off()
fileCat(successText('Quality Profile Plot - Read2 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)


#---------------------#
### Trim And Filter ###
#---------------------#
fileCat('1) Filtering\n', file=log.file, append=TRUE)
filtered_read1.files <- file.path(dada2.filter_dir, basename(read1.files))
filtered_read2.files <- file.path(dada2.filter_dir, basename(read2.files))
names(filtered_read1.files) <- samples.name
names(filtered_read2.files) <- samples.name
filterAndTrim.out <- filterAndTrim(read1.files, filtered_read1.files,
                                   read2.files, filtered_read2.files,
                                   maxN=args$maxN,
                                   maxEE=c(args$maxEEF, args$maxEER),
                                   truncLen=c(args$truncLenF, args$truncLenR),
                                   truncQ=args$truncQ,
                                   rm.phix=FALSE,
                                   compress=TRUE,
                                   multithread=multithread,
                                   verbose=verbose)
# Filter & Trim 출력 & Log
fileCat(infoText('Results of filtering and trimming\n'), file=log.file, append=TRUE)
print(filterAndTrim.out)
write.table(filterAndTrim.out, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# Filter & Trim : Table 저장
filterAndTrim.out_file <- file.path(dada2.dir, 'dada2.filterAndTrim')
write.table(filterAndTrim.out, filterAndTrim.out_file,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
fileCat(successText('FilterAndTrim.txt 생성 완료\n'), file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
# Filtered Files 존재 확인
fileCat('Filtered Read1 파일 확인 : ', file=log.file, append=TRUE)
fileCat(ifelse(file.exists(filtered_read1.files), '.', 'x'), sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat('Filtered Read2 파일 확인 : ', file=log.file, append=TRUE)
fileCat(ifelse(file.exists(filtered_read2.files), '.', 'x'), sep='', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
# 파일명 검출
fileCat('\n', file=log.file, append=TRUE)
filtered_read1.files <- list.files(dada2.filter_dir, pattern='_1.fastp.fastq.gz$', full.names=TRUE)
filtered_read2.files <- list.files(dada2.filter_dir, pattern='_2.fastp.fastq.gz$', full.names=TRUE)
names(filtered_read1.files) <- samples.name
names(filtered_read2.files) <- samples.name
# 파일명 출력
fileCat('------------------\n', file=log.file, append=TRUE)
fileCat(names(filtered_read1.files), sep=', ', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat(filtered_read1.files, sep='\n', file=log.file, append=TRUE)
fileCat('------------------\n', file=log.file, append=TRUE)
fileCat(names(filtered_read2.files), sep=', ', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
fileCat(filtered_read2.files, sep='\n', file=log.file, append=TRUE)
fileCat('------------------\n', file=log.file, append=TRUE)

if (length(filtered_read1.files) == 0) {
    errMesg('Filter And Trim 단계에서 모든 read들이 제거되었습니다.')
    message(magenta('\t--> truncLenF/R 이 read length 보다 더 긴지 확인하세요.'))
    quit('no', 2)
}
fileCat(successText('Filtering 완료\n\n\n'), file=log.file, append=TRUE)


#-----------------------#
### Learn Error Rates ###
#-----------------------#
fileCat('2) Learning Error Rates\n', file=log.file, append=TRUE)
fileCat(length(filtered_read1.files), file=log.file, append=TRUE)
fileCat(filtered_read1.files, sep='\n', file=log.file, append=TRUE)
err_read1 <- learnErrors(filtered_read1.files, multithread=multithread, verbose=verbose)
err_read2 <- learnErrors(filtered_read2.files, multithread=multithread, verbose=verbose)
fileCat(successText('Learning Error Rates 완료\n\n\n'), file=log.file, append=TRUE)


#----------------------#
### Error Rates Plot ###
#----------------------#
errPlot.name = 'Estimated_Error_Rates_'
errPlot_R1.file = file.path(dada2.dir, paste0(errPlot.name, 'R1.pdf'))
errPlot_R2.file = file.path(dada2.dir, paste0(errPlot.name, 'R2.pdf'))
pdf(errPlot_R1.file)
plotErrors(err_read1, nominalQ=TRUE, err_in=TRUE, err_out=TRUE)
dev.off()
fileCat(successText('Estimated Error Rates Plot - Read1 완료\n'), file=log.file, append=TRUE)
pdf(errPlot_R2.file)
plotErrors(err_read2, nominalQ=TRUE, err_in=TRUE, err_out=TRUE)
dev.off()
fileCat(successText('Estimated Error Rates Plot - Read2 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)

#-------------------#
### Dereplication ###  Ver. 1.10 에서는 별도로 실행해야 함.
#-------------------#
fileCat('3) Dereplicate FASTQ\n', file=log.file, append=TRUE)
derep_R1 <- derepFastq(filtered_read1.files, verbose=verbose)
derep_R2 <- derepFastq(filtered_read2.files, verbose=verbose)
fileCat(successText('Dereplicate FASTQ 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)


#----------------------#
### Sample Inference ###
#----------------------#
fileCat('4) Denoise -  Sample Inference\n', file=log.file, append=TRUE)
dada_R1 <- dada(derep_R1, err=err_read1, multithread=multithread, verbose=verbose)
dada_R2 <- dada(derep_R2, err=err_read2, multithread=multithread, verbose=verbose)
# DADA2 1.41.1
#dada_R1 <- dada(filtered_read1.files, err=err_read1, multithread=multithread, verbose=verbose)
#dada_R2 <- dada(filtered_read2.files, err=err_read2, multithread=multithread, verbose=verbose)
if (length(samples.name) == 1) {
    denoise_track <- cbind(getN(dada_R1), getN(dada_R2))
} else {
    denoise_track <- cbind(sapply(dada_R1, getN), sapply(dada_R2, getN))
}
colnames(denoise_track) <- c('Denoised_R1', 'Denoised_R2')
rownames(denoise_track) <- samples.name
# Denoise Track 출력 & Log
fileCat(infoText('Results of denoising\n'), file=log.file, append=TRUE)
print(denoise_track)
write.table(denoise_track, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# Denoise Track 저장
denoise.out_file <- file.path(dada2.dir, 'dada2.denoise')
write.table(denoise_track, denoise.out_file,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
fileCat(successText('dada2.denoise 생성 완료\n'), file=log.file, append=TRUE)
fileCat(successText('Denoise - Sample Inference 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)

#-------------------------#
#### Merge Paired Reads ###
#-------------------------#
fileCat('5) Merge Paired Reads\n', file=log.file, append=TRUE)
mergers <- mergePairs(dada_R1, derep_R1,
                      dada_R2, derep_R2,
                      verbose=verbose, minOverlap=args$minOverlap)
#mergers <- mergePairs(dada_R1, filtered_read1.files, dada_R2, filtered_read2.files, verbose=verbose, minOverlap=15)
if (length(samples.name) == 1) {
    mergers_track <- cbind(getN(mergers))
} else {
    mergers_track <- cbind(sapply(mergers, getN))
}
## Track 데이터 생성
colnames(mergers_track) <- 'Merged'
rownames(mergers_track) <- samples.name
# Mergers Track 출력 & Log
fileCat(infoText('Results of merging\n'), file=log.file, append=TRUE)
print(mergers_track)
write.table(mergers_track, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# Mergers Track 저장
merge.out_file <- file.path(dada2.dir, 'dada2.merge')
write.table(mergers_track, merge.out_file,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
fileCat(successText('dada2.merge 생성 완료\n'), file=log.file, append=TRUE)
fileCat(successText('Merge Paired Reads 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)


#------------------------------#
### Construct Sequence Table ###
#------------------------------#
fileCat('6) Construct Sequence Table\n', file=log.file, append=TRUE)
seqtab <- makeSequenceTable(mergers)
# ASVs 개수 출력
fileCat(infoText('시료개수\tASVs개수\n'), file=log.file, append=TRUE)
fileCat(dim(seqtab), sep='\t', file=log.file, append=TRUE)
fileCat('\n', file=log.file, append=TRUE)
# 길이별 분포 출력
fileCat(infoText('길이별 분포\n'), file=log.file, append=TRUE)
dist_length <- table(nchar(getSequences(seqtab)))
print(dist_length)
write.table(dist_length, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# Seqtab.Length Track 저장
seqtab_length.out_file <- file.path(dada2.dir, 'dada2.seqtab.length')
table.names <- c('Length', 'Freq')
write.table(dist_length, seqtab_length.out_file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=table.names, quote=FALSE)
fileCat('\n\n', file=log.file, append=TRUE)
# RDS 파일 저장
seqtab_output.rds <- file.path(dada2.dir, 'seqtab.rds')
saveRDS(seqtab, seqtab_output.rds)
fileCat(successText('seqtab 저장 완료\n'))
fileCat(successText('Construct Sequence Table 완료\n'), file=log.file, append=TRUE)
fileCat('\n\n', file=log.file, append=TRUE)


#--------------------#
### Remove Chimera ###
#--------------------#
fileCat('7) Remove Chimeras (method: ', args$chimera_method, ')\n', sep='', file=log.file, append=TRUE)
seqtab.noChimera <- removeBimeraDenovo(seqtab, method=args$chimera_method,
                                       multithread=multithread, verbose=verbose)
# Chimera 제거 후 ASVs개수 출력
fileCat(infoText('시료개수\tASVs개수(noChimera)\n'), file=log.file, append=TRUE)
if (is.null(args$sample_name)) {
    fileCat(dim(seqtab.noChimera), sep='\t', file=log.file, append=TRUE)
} else {
    fileCat('1', length(colnames(seqtab.noChimera)), sep='\t', file=log.file, append=TRUE)
}
fileCat('\n', file=log.file, append=TRUE)
# Chimera 제거 비율 출력  & Track 데이터 생성
fileCat('Chimera 제거 후 남은 비율 : ', file=log.file, append=TRUE)
if (length(samples.name) == 1) {
    fileCat(sum(seqtab.noChimera)/sum(seqtab), sep='\n', file=log.file, append=TRUE)
    noChimera_track <- cbind(sum(seqtab.noChimera))
} else {
    fileCat(rowSums(seqtab.noChimera) / rowSums(seqtab))
    noChimera_track <- cbind(rowSums(seqtab.noChimera))
}
# Chimera Track 출력
colnames(noChimera_track) <- 'noChimera'
rownames(noChimera_track) <- samples.name
# 출력 & Log
fileCat(infoText('Results of removing chimera\n'), file=log.file, append=TRUE)
print(noChimera_track)
write.table(noChimera_track, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# noChimera Track 파일 저장
dada2_noChimera.out_file <- file.path(dada2.dir, 'dada2.noChimera')
write.table(noChimera_track, dada2_noChimera.out_file,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
fileCat(successText('dada2.noChimera 생성 완료\n'), file=log.file, append=TRUE)
# RDS 저장
seqtab_noChimera_output.rds <- file.path(dada2.dir, 'seqtab_noChimera.rds')
saveRDS(seqtab.noChimera, seqtab_noChimera_output.rds)
fileCat(successText('seqtab_noChimera.rds 저장 완료\n'))
fileCat(successText('Remove Chimera 완료\n\n\n'), file=log.file, append=TRUE)

#--------------------------------------#
### Track Reads Through The Pipeline ###
#--------------------------------------#
if (length(samples.name) == 1) {
    track <- cbind(filterAndTrim.out,
                   getN(dada_R1),
                   getN(dada_R2),
                   getN(mergers),
                   sum(seqtab.noChimera))
} else {
    track <- cbind(filterAndTrim.out,
                   sapply(dada_R1, getN),
                   sapply(dada_R2, getN),
                   sapply(mergers, getN),
                   rowSums(seqtab.noChimera))
}
colnames(track) <- c('Input', 'Filtered', 'Denoised_R1', 'Denoised_R2', 'Merged', 'noChimera')
rownames(track) <- samples.name
# 출력 & Log
fileCat(infoText('DADA2 Summary\n'), file=log.file, append=TRUE)
print(track)
write.table(track, log.file, append=TRUE,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# 파일 저장
dada2_summary.out_file <- file.path(dada2.dir, 'dada2.summary')
write.table(track, dada2_summary.out_file,
            sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
fileCat(successText('dada2.summary 생성 완료\n\n\n'), file=log.file, append=TRUE)


#-----------------------------#
### Save file(fasta, table) ###
#-----------------------------#
fileCat('8) Write Output\n', file=log.file, append=TRUE)
# Table 저장 - ASVs.tsv
trans_seqtab.noChimera <- t(seqtab.noChimera)  # QIIME has OTUs as rows
if (length(args$sample_name) == 1) {
    col.names <- samples.name
    col.names[[1]] <- paste0('#ASVs_ID\t', col.names[[1]])
} else {
    # TODO: 시료가 여러개일 경우
}
row.names <-  paste0('ASV', seq(1, length(trans_seqtab.noChimera)))
dada2_output.tsv <- file.path(dada2.dir, 'ASVs.tsv')
write.table(trans_seqtab.noChimera, dada2_output.tsv,
            sep='\t', row.names=row.names, col.names=col.names, quote=FALSE)
fileCat(successText('ASVs.tsv 생성 완료\n'), file=log.file, append=TRUE)
# Table 저장 - seqASvs.tsv
# dada2_output.seq_tsv <- file.path(dada2.dir, 'seqASVs.tsv')
# write.table(trans_seqtab.noChimera, dada)

# FASTA 저장 - ASVs.fasta
dada2_output.fasta <- file.path(dada2.dir, 'ASVs.fasta')
# asv_id <- paste0('ASV', seq(1, length(seqtab.noChimera)), ';size=', unname(seqtab.noChimera),';')
asv_id <- paste0('ASV', seq(1, length(seqtab.noChimera)))
uniquesToFasta(seqtab.noChimera, dada2_output.fasta, id=asv_id)
fileCat(successText('ASVs.fasta 생성 완료\n'), file=log.file, append=TRUE)
quit(status=0)
