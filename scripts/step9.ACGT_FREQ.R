#09/20/2019
#Ruijia Wang

rm(list = ls())
suppressMessages(require(GenomicRanges))
suppressMessages(require(Biostrings))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(require(BSgenome.Mmusculus.UCSC.mm9))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg19))
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)


### setting ###
setwd(args[1])
geno=args[2]
mypath = paste0(getwd(),"/3UTR/")

##find all files in mypath folder, using 3'reads pipeline 3most results
filelist = list.files(path = mypath, full.names = TRUE)

##if all the file names start with "3most"
namelist = sub(".*\\/(3mostAPA.+)\\.DRPM.fix.tbl", "\\1", filelist)
namelist = gsub("3mostAPA.", "", namelist)

##define 3 simple functions##
absmax <- function(x) { x[which.max( abs(x) )]}
cald<-function(kmer1, kmer2){
	d=0
	k1=strsplit(kmer1, "")[[1]]
	k2=strsplit(kmer2, "")[[1]]
	for (p in 1:length(k1)){
	  if (k1[p] != k2[p]){
	  d=d+1
	  }
	}
	return(d)
}

calFREQ_prx <- function(datain,geno,i) {
    gr_region_prx = GRanges(seqnames = datain$chromosome, 
    						ranges = IRanges(start = datain[,paste(prx_region[i],"from",sep = "_")], 
    						end = datain[,paste(prx_region[i],"to",sep = "_")]),
    						strand = datain$strand)
	if(geno=="mm9")
	{
	gr_prx = getSeq(Mmusculus, gr_region_prx)
	}
	if(geno=="hg19")
	{
	gr_prx = getSeq(Hsapiens, gr_region_prx)
	}
	counts_ACGT = t(as.data.frame(consensusMatrix(gr_prx)))[,c('A','C','G','T')]
	dfsingle=colSums(counts_ACGT)
	allbp=sum(dfsingle)
	dfsingle_FREQ=as.data.frame(dfsingle/allbp*100)
	dfsingle_FREQ$bp=rownames(dfsingle_FREQ)
	colnames(dfsingle_FREQ)[1]='Freq'
	dfsingle_FREQ=dfsingle_FREQ[,c('bp','Freq')]
	dfsingle_FREQ$bp[4]='U'
	return(dfsingle_FREQ)
}

calFREQ_dis <- function(datain,geno,i) {
    gr_region_prx = GRanges(seqnames = datain$chromosome, 
    						ranges = IRanges(start = datain[,paste(dis_region[i],"from",sep = "_")], 
    						end = datain[,paste(dis_region[i],"to",sep = "_")]),
    						strand = datain$strand)
	if(geno=="mm9")
	{
	gr_prx = getSeq(Mmusculus, gr_region_prx)
	}
	if(geno=="hg19")
	{
	gr_prx = getSeq(Hsapiens, gr_region_prx)
	}
	counts_ACGT = t(as.data.frame(consensusMatrix(gr_prx)))[,c('A','C','G','T')]
	dfsingle=colSums(counts_ACGT)
	allbp=sum(dfsingle)
	dfsingle_FREQ=as.data.frame(dfsingle/allbp*100)
	dfsingle_FREQ$bp=rownames(dfsingle_FREQ)
	colnames(dfsingle_FREQ)[1]='Freq'
	dfsingle_FREQ=dfsingle_FREQ[,c('bp','Freq')]
	dfsingle_FREQ$bp[4]='U'	
	return(dfsingle_FREQ)
}	
	
##major loop start##
wb <- createWorkbook()
for (k in 1:length(filelist)){

  data = read.table(filelist[k], header = TRUE, quote = "", sep = "\t", stringsAsFactors=FALSE)
  type = grep("^pA1.pAutype_", colnames(data), value = TRUE)
  df_NC=data[data[,type]=='NC',]
  rownames(df_NC)=1:nrow(df_NC)
  
  df_diff=data[data[,type]!='NC',]
  countdiff=nrow(df_diff)
  
  df_NCsub=df_NC[sample(nrow(df_NC), countdiff), ]
  
  data=rbind(df_diff,df_NCsub)
  rownames(data)=1:nrow(data)
  data=data[,c('gene_symbol','pA_pos_pA1','pA_pos_pA2','chromosome','strand',type)]
  
  prx_region = c("prx_region_1", "prx_region_2", "prx_region_3", "prx_region_4")
  dis_region = c("dis_region_1", "dis_region_2", "dis_region_3", "dis_region_4")  
 
  ####build the region around PAS including -100 to -41 nt,-40 to -1 nt and +1 to +100 nt, region is based on Wencheng's Plos Genetics paper
  data$prx_region_1_from <- NA
  data$prx_region_1_to <- NA
  data$prx_region_2_from <- NA
  data$prx_region_2_to <- NA
  data$prx_region_3_from <- NA
  data$prx_region_3_to <- NA
  data$prx_region_4_from <- NA
  data$prx_region_4_to <- NA  
  
  data$dis_region_1_from <- NA
  data$dis_region_1_to <- NA
  data$dis_region_2_from <- NA
  data$dis_region_2_to <- NA
  data$dis_region_3_from <- NA
  data$dis_region_3_to <- NA
  data$dis_region_4_from <- NA
  data$dis_region_4_to <- NA 
  
  for (i in 1:nrow(data)){
  	if (data$strand[i] == '+'){
  		data[,"prx_region_1_from"][i] = data$pA_pos_pA1[i] - 150
  		data[,"prx_region_1_to"][i] = data$pA_pos_pA1[i] - 51
  		data[,"prx_region_2_from"][i] = data$pA_pos_pA1[i] - 50
  		data[,"prx_region_2_to"][i] = data$pA_pos_pA1[i] - 1
  		data[,"prx_region_3_from"][i] = data$pA_pos_pA1[i] + 1
  		data[,"prx_region_3_to"][i] = data$pA_pos_pA1[i] + 50
  		data[,"prx_region_4_from"][i] = data$pA_pos_pA1[i] + 51
  		data[,"prx_region_4_to"][i] = data$pA_pos_pA1[i] + 150		
		
  		data[,"dis_region_1_from"][i] = data$pA_pos_pA2[i] - 150
  		data[,"dis_region_1_to"][i] = data$pA_pos_pA2[i] - 51
  		data[,"dis_region_2_from"][i] = data$pA_pos_pA2[i] - 50
  		data[,"dis_region_2_to"][i] = data$pA_pos_pA2[i] - 1
  		data[,"dis_region_3_from"][i] = data$pA_pos_pA2[i] + 1
  		data[,"dis_region_3_to"][i] = data$pA_pos_pA2[i] + 50
  		data[,"dis_region_4_from"][i] = data$pA_pos_pA2[i] + 51
  		data[,"dis_region_4_to"][i] = data$pA_pos_pA2[i] + 150		
  	}
  		if (data$strand[i] == '-'){
  			data[,"prx_region_1_from"][i] = data$pA_pos_pA1[i] + 51
  			data[,"prx_region_1_to"][i] = data$pA_pos_pA1[i] + 150
  			data[,"prx_region_2_from"][i] = data$pA_pos_pA1[i] + 1
  			data[,"prx_region_2_to"][i] = data$pA_pos_pA1[i] + 50
  			data[,"prx_region_3_from"][i] = data$pA_pos_pA1[i] - 50
  			data[,"prx_region_3_to"][i] = data$pA_pos_pA1[i] - 1
  			data[,"prx_region_4_from"][i] = data$pA_pos_pA1[i] - 150
  			data[,"prx_region_4_to"][i] = data$pA_pos_pA1[i] - 51			
			
			
  			data[,"dis_region_1_from"][i] = data$pA_pos_pA2[i] + 51
  			data[,"dis_region_1_to"][i] = data$pA_pos_pA2[i] + 150
  			data[,"dis_region_2_from"][i] = data$pA_pos_pA2[i] + 1
  			data[,"dis_region_2_to"][i] = data$pA_pos_pA2[i] + 50
  			data[,"dis_region_3_from"][i] = data$pA_pos_pA2[i] - 50
  			data[,"dis_region_3_to"][i] = data$pA_pos_pA2[i] - 1
  			data[,"dis_region_4_from"][i] = data$pA_pos_pA2[i] - 150
  			data[,"dis_region_4_to"][i] = data$pA_pos_pA2[i] - 51			
			
  		}
  }

  for (i in 1:4){
########proximal#######
	df_diff=data[data[,type]!='NC',]
	df_NC=data[data[,type]=='NC',]

	dfATGT_PRX_DIFF=calFREQ_prx(df_diff,geno,i)
	dfATGT_PRX_NC=calFREQ_prx(df_NC,geno,i)
	colnames(dfATGT_PRX_DIFF)[2]=paste(prx_region[i],'diff_PAS',sep='_')	
	colnames(dfATGT_PRX_NC)[2]=paste(prx_region[i],'NC_PAS',sep='_')
	
	
	
	dfATGT_DIS_DIFF=calFREQ_dis(df_diff,geno,i)
	dfATGT_DIS_NC=calFREQ_dis(df_NC,geno,i)	
	colnames(dfATGT_DIS_DIFF)[2]=paste(dis_region[i],'diff_PAS',sep='_')	
	colnames(dfATGT_DIS_NC)[2]=paste(dis_region[i],'NC_PAS',sep='_')

	dfFREQ_ALL=merge(dfATGT_PRX_DIFF,dfATGT_PRX_NC,by='bp')
	dfFREQ_ALL=merge(dfFREQ_ALL,dfATGT_DIS_DIFF,by='bp')
	dfFREQ_ALL=merge(dfFREQ_ALL,dfATGT_DIS_NC,by='bp')
	dfFREQ_ALL=dfFREQ_ALL[!duplicated(dfFREQ_ALL), ]
	
	if(i==1)
	{
	dfout = dfFREQ_ALL
	} else {
	dfout=merge(dfout,dfFREQ_ALL,by='bp')
	}
}	

print(namelist[k])
addWorksheet(wb, namelist[k])
writeData(wb, sheet = namelist[k], x = dfout)	
}
saveWorkbook(wb, "PAS.nuc_freq.diff_vs_NC.xlsx", overwrite = TRUE)  