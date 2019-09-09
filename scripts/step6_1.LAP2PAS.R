#01/27/2018
#Ruijia Wang
rm(list = ls())
suppressMessages(library(GenomicRanges))
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(GenomicAlignments))
args = commandArgs(trailingOnly=TRUE)
#set parameter#
flank=24

### setting ###
rootdir=args[1]
project=args[2]
refdir=args[3]
geno=args[4]

LAPdir=paste(rootdir,project,"rawsam/LAP_raw/",sep="/")
setwd(paste(rootdir,project,"tbl/",sep="/"))


if(geno=="hg19")
{
  PASfile=paste(refdir,'human.cental.decompact.txt',sep='/')
}

if(geno=="mm9")
{
  PASfile=paste(refdir,'mouse.cental.decompact.txt',sep='/')
}

filelist = list.files(path = LAPdir, full.names = TRUE)
filelist = grep(".LAP.txt",filelist,value = TRUE)
namelist = basename(filelist)
namelist = gsub(".LAP.txt", "", namelist)

for (k in 1:length(filelist)){
if(k==1)
{
dftmp=read.table(filelist[k], header = TRUE, quote = "", sep = "\t", stringsAsFactors=FALSE)
colnames(dftmp)=c('PoII_ID','Chr','Pos','Strand',paste0('Reads_count_',namelist[k]))
}
if(k>1)
{
data=read.table(filelist[k], header = TRUE, quote = "", sep = "\t", stringsAsFactors=FALSE)
colnames(data)=c('PoII_ID','Chr','Pos','Strand',paste0('Reads_count_',namelist[k]))
dftmp=merge(dftmp,data,by=c('PoII_ID','Chr','Pos','Strand'),all=TRUE)
}
}

dftmp[is.na(dftmp)] <- 0
readcols=grep("Reads_count_", colnames(dftmp),value = TRUE)
dftmp$total_reads=rowSums(dftmp[,readcols], na.rm = TRUE)

dftmp=dftmp[!duplicated(dftmp), ]
dftmp$start=dftmp$Pos
dftmp$end=dftmp$Pos

ENDdb<-makeGRangesFromDataFrame(dftmp, keep.extra.columns = TRUE)
ENDdb$Gene_name='Not'

##build pAs database##
PAS = read.table(PASfile, header = TRUE, quote = "", sep = "\t", stringsAsFactors=FALSE)
PAS=PAS[,c("PASid","Chrom","Strand","Pos")]
PAS=PAS[!duplicated(PAS$PASid), ]

PAS$start=PAS$Pos
PAS$end=PAS$Pos
po2db<-makeGRangesFromDataFrame(PAS, keep.extra.columns = TRUE)
po2db_flank <- resize(po2db, width = flank*2 + width(po2db), fix = "center")
names(po2db_flank)=po2db_flank$PASid

####annotate LAP####
olp = findOverlaps(ENDdb, po2db_flank)
ENDdb[queryHits(olp), ]$Gene_name = po2db_flank[subjectHits(olp), ]$PASid
dfENDdb=as.data.frame(ENDdb)
dfENDdb=dfENDdb[!grepl("_random", dfENDdb$seqnames) & !grepl("BK000964", dfENDdb$seqnames)& !grepl("chrM", dfENDdb$seqnames),]
dfENDdb=dfENDdb[,-c(2:4)]
colnames(dfENDdb)[1:3]=c('Chrom','Strand','LAPID')
dfENDdb$LAPID=paste0(dfENDdb$Chrom, ':', dfENDdb$Pos, ':', dfENDdb$Strand)
colnames(dfENDdb)[length(colnames(dfENDdb))]='Match'
write.table(dfENDdb,"All_LAP2PAS.tbl", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



