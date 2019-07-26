rm(list = ls())
suppressMessages(library("DESeq"))
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("Rsamtools"))
suppressMessages(library("GenomicFeatures"))
args = commandArgs(trailingOnly=TRUE)

### setting ###
rootdir=args[1]
project=args[2]
refdir=args[3]
geno=args[4]


####seting####
setwd(paste(rootdir,project,"tbl/",sep="/"))
outputfile=paste0(geno,".pA2gene_usage.DRPM.fix.tbl")
inputfile='All_LAP2PAS.tbl'

if(geno=="hg19")
{
  PASfile=paste(refdir,'human.cental.decompact.txt',sep='/')
}

if(geno=="mm9")
{
  PASfile=paste(refdir,'mouse.cental.decompact.txt',sep='/')
}

dfinput = read.table(inputfile, , header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
dfdb = read.table(PASfile, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

dfinput=dfinput[which(dfinput$Match!='Not'),]
readscols=grep("Reads_count_", colnames(dfinput), value = TRUE)

dfinput=dfinput[,c('Chrom','Pos','Strand','Match',readscols)]
colnames(dfinput)=c('Chrom','Pos','Strand','pAid',readscols)

dfin=dfinput

dfquant=aggregate(dfin[,readscols], by=list(PASid=dfin$pAid), FUN=sum)

dfdb=dfdb[,c('PASid','Chrom','Pos','Strand','LOCATION','region','GENEID','gene_symbol','gene_name','Detailed_gene_type','ext','pA_type')]
dfquant=merge(dfquant,dfdb,by='PASid')
dfout=dfquant[,c('PASid','Chrom','Pos','Strand',readscols,'LOCATION','region','GENEID','gene_symbol','gene_name','Detailed_gene_type','ext','pA_type')]
names(dfout)<-c('pAid','chromosome','pA_pos','strand',readscols,'LOCATION','region','GENEID','gene_symbol','gene_desc','gene_Biotype','ext','pAtype_1')
names(dfout)<-gsub("Reads_count_","num_",names(dfout))
dfout=dfout[which(dfout$pAtype_1!='intergenic'),]

allFCcol=grep("num_", colnames(dfout), value = TRUE)
dfquery=dfout[which(dfout$pAtype_1!='UA' & dfout$pAtype_1!=""),]
dfquery=dfquery[which(dfquery$chromosome!='chrM'),]
dfquery=dfquery[!duplicated(dfquery$pAid), ]

dfquery=dfquery[,c('pAid', allFCcol)]
rownames(dfquery)=dfquery$pAid
readscount=dfquery[,c(-1)]

### build your condition
condition = factor(names(readscount)) ## without replicates
### get count table
dds = newCountDataSet(readscount, condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

DRPMcols=gsub('num_','DRPM_',allFCcol)
dfout[,DRPMcols]=t( t(dfout[,allFCcol])/sizeFactors(dds))

###save###
write.table(dfout, outputfile, sep="\t", row.names=FALSE,quote=FALSE)





