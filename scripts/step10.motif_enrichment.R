#01/27/2017
#Ruijia Wang

rm(list = ls())
suppressMessages(require(GenomicRanges))
suppressMessages(require(Biostrings))
# suppressMessages(library(xlsx))
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
# project=args[2]
# refdir=args[3]
# 

# setwd("/home/soft/QSr/APAlyzer_qrev/tbl")
# mypath = paste0(getwd(),"/3UTR/test")
# geno='mm9'

# setwd("D:/R/kmer")
# mypath = paste0(getwd(),"/3UTR/test")
# geno='mm9'

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

##major loop start##

for (k in 1:length(filelist)){
wb <- createWorkbook()
  # jgc()
  ##test for file 1 in the folder
  data = read.table(filelist[k], header = TRUE, quote = "", sep = "\t", stringsAsFactors=FALSE)
  type = grep("^pA1.pAutype_", colnames(data), value = TRUE)
  df_NC=data[data[,type]=='NC',]
  rownames(df_NC)=1:nrow(df_NC)
  
  df_diff=data[data[,type]!='NC',]
  countdiff=nrow(df_diff)
  
  df_NCsub=df_NC[sample(nrow(df_NC), countdiff), ]
  
  data=rbind(df_diff,df_NCsub)
  # data = data[which(abs(data$pA_pos_pA2-data$pA_pos_pA1)<20000),]
  # data = data[which(data$type!='NC'),]
  rownames(data)=1:nrow(data)
  data=data[,c('gene_symbol','pA_pos_pA1','pA_pos_pA2','chromosome','strand',type)]
  # print(paste(nrow(df_NC),nrow(df_diff),nrow(df_NCsub),nrow(data),sep="::"))
  # type = grep("^type", colnames(data), value = TRUE)
  
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

  
  ###mm9
  # require(BSgenome.Mmusculus.UCSC.mm9)
  
  ##hg19
  # require(BSgenome.Hsapiens.UCSC.hg19)
  for (j in c(4,6)){
  for (i in 1:4){
########proximal#######
    gr_region_prx = GRanges(seqnames = data$chromosome, 
    						ranges = IRanges(start = data[,paste(prx_region[i],"from",sep = "_")], 
    						end = data[,paste(prx_region[i],"to",sep = "_")]),
    						strand = data$strand)
 
	if(geno=="mm9")
	{
	gr_prx = getSeq(Mmusculus, gr_region_prx)
	}
	if(geno=="hg19")
	{
	gr_prx = getSeq(Hsapiens, gr_region_prx)
	}	
	
	counts_prx = oligonucleotideFrequency(gr_prx, width = j, step = 1)
	lower_index = as.numeric(rownames(data[data[,type]=='NC',]))
	upper_index = as.numeric(rownames(data[data[,type] != 'NC',]))
	###calculate each k-mer frequency
	counts_prx = t(rbind(colSums(counts_prx[upper_index,]), colSums(counts_prx[lower_index,])))
	colnames(counts_prx) = c("proximal_UP", "proximal_NOTUP")
	###generate fisher exact test table 
	sum_up = rep(sum(counts_prx[,"proximal_UP"]),nrow(counts_prx))
	sum_noup = rep(sum(counts_prx[,"proximal_NOTUP"]),nrow(counts_prx))
	df_tmp <- data.frame(counts_prx, sum_up, sum_noup)
	df_tmp$other_up = df_tmp[,"sum_up"] - df_tmp[,"proximal_UP"]
	df_tmp$other_noup = df_tmp[,"sum_noup"] - df_tmp[,"proximal_NOTUP"]
	df_prx = df_tmp[,-c(3,4)]
	###fisher exact test p-value and odds_ratio
	df_prx$p.value = apply(df_prx, 1, function(x) fisher.test(matrix(x,nr=2))$p.value)
	df_prx$odds_ratio = apply(df_prx[,-5], 1, function(x) fisher.test(matrix(x,nr=2))$estimate)
	###calculate S for S = 1 -- up-regulated, S = -1 -- down-regulated
	df_prx$S <- NA
	for (h in 1:nrow(df_prx)){
		if (df_prx$odds_ratio[h] >= 1){
			df_prx$S[h] = 1
		}
			if (df_prx$odds_ratio[h] < 1){
				df_prx$S[h] = -1
			}
	}


	###calculate SS score
	df_prx$SS = (-log10(df_prx$p.value))*df_prx$S
	df_prx = df_prx[order(df_prx$p.value, decreasing = FALSE),]

	df_prx$kmer=rownames(df_prx)
	df_prx=df_prx[,-c(6:7)]
	df_prx<-df_prx[c(7,6,5,1:4)]

	df_prx_export=df_prx[,c(1,3)]
	df_prx_export=df_prx
	df_prx_export[,c(1:2)]=df_prx_export[,c(1:2)] %>% mutate_if(is.numeric, funs(round(., 2)))
	df_prx_export$kmer <- gsub('T', 'U', df_prx_export$kmer)
	df_prx_Sh=df_prx_export[which(df_prx_export$SS>=0),]
	df_prx_Ln=df_prx_export[which(df_prx_export$SS<0),]

	addWorksheet(wb, paste0('enrich_',prx_region[i],"_",j,"-mer"))
	addWorksheet(wb, paste0('depleted_',prx_region[i],"_",j,"-mer"))
	writeData(wb, sheet = paste0('enrich_',prx_region[i],"_",j,"-mer"), x = df_prx_Sh[,c("kmer","SS","p.value")])	
	writeData(wb, sheet = paste0('depleted_',prx_region[i],"_",j,"-mer"), x = df_prx_Ln[,c("kmer","SS","p.value")])
	
########Distal#######
    gr_region_dis = GRanges(seqnames = data$chromosome, 
    						ranges = IRanges(start = data[,paste(dis_region[i],"from",sep = "_")], 
    						end = data[,paste(dis_region[i],"to",sep = "_")]),
    						strand = data$strand)
	if(geno=="mm9")
	{
	gr_dis = getSeq(Mmusculus, gr_region_dis)
	}
	if(geno=="hg19")
	{
	gr_dis = getSeq(Hsapiens, gr_region_dis)
	}	
	
	counts_dis = oligonucleotideFrequency(gr_dis, width = j, step = 1)
	lower_index = as.numeric(rownames(data[data[,type]=='NC',]))
	upper_index = as.numeric(rownames(data[data[,type] != 'NC',]))	
	
	###calculate each k-mer frequency
	counts_dis = t(rbind(colSums(counts_dis[upper_index,]), colSums(counts_dis[lower_index,])))
	colnames(counts_dis) = c("distal_UP", "distal_NOTUP")
	###generate fisher exact test table 
	sum_up = rep(sum(counts_dis[,"distal_UP"]),nrow(counts_dis))
	sum_noup = rep(sum(counts_dis[,"distal_NOTUP"]),nrow(counts_dis))
	df_tmp <- data.frame(counts_dis, sum_up, sum_noup)
	df_tmp$other_up = df_tmp[,"sum_up"] - df_tmp[,"distal_UP"]
	df_tmp$other_noup = df_tmp[,"sum_noup"] - df_tmp[,"distal_NOTUP"]
	df_dis = df_tmp[,-c(3,4)]
	###fisher exact test p-value and odds_ratio
	df_dis$p.value = apply(df_dis, 1, function(x) fisher.test(matrix(x,nr=2))$p.value)
	df_dis$odds_ratio = apply(df_dis[,-5], 1, function(x) fisher.test(matrix(x,nr=2))$estimate)
	###calculate S for S = 1 -- up-regulated, S = -1 -- down-regulated
	df_dis$S <- NA
	for (h in 1:nrow(df_dis)){
		if (df_dis$odds_ratio[h] >= 1){
			df_dis$S[h] = 1
		}
			if (df_dis$odds_ratio[h] < 1){
				df_dis$S[h] = -1
			}
	}


	###calculate SS score
	df_dis$SS = (-log10(df_dis$p.value))*df_dis$S
	df_dis = df_dis[order(df_dis$p.value, decreasing = FALSE),]

	df_dis$kmer=rownames(df_dis)
	df_dis=df_dis[,-c(6:7)]
	df_dis<-df_dis[c(7,6,5,1:4)]

	df_dis_export=df_dis[,c(1,3)]
	df_dis_export=df_dis
	df_dis_export[,c(1:2)]=df_dis_export[,c(1:2)] %>% mutate_if(is.numeric, funs(round(., 2)))
	df_dis_export$kmer <- gsub('T', 'U', df_dis_export$kmer)
	df_dis_Sh=df_dis_export[which(df_dis_export$SS>=0),]
	df_dis_Ln=df_dis_export[which(df_dis_export$SS<0),]

	addWorksheet(wb, paste0('enrich_',dis_region[i],"_",j,"-mer"))
	addWorksheet(wb, paste0('depleted_',dis_region[i],"_",j,"-mer"))		
	writeData(wb, sheet = paste0('enrich_',dis_region[i],"_",j,"-mer"), x = df_dis_Sh[,c("kmer","SS","p.value")])	
	writeData(wb, sheet = paste0('depleted_',dis_region[i],"_",j,"-mer"), x = df_dis_Ln[,c("kmer","SS","p.value")])		
	

	}
}
saveWorkbook(wb, paste0(namelist[k],".motif_enrichment.diff_vs_NC.xlsx"), overwrite = TRUE)  
}
 