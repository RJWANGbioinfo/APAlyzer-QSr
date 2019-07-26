rm(list = ls())
suppressMessages(library("DEXSeq"))
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("Rsamtools"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
args = commandArgs(trailingOnly=TRUE)

### setting ###
rootdir=args[1]
project=args[2]
samplefile=args[3]
geno=args[4]
sample_t=unlist(strsplit(as.character(args[5])," "))
sample_c=unlist(strsplit(as.character(args[6])," "))

setwd(rootdir)
PASfile=paste0(rootdir,'/',project,"/tbl/",geno,".pA2gene_usage.DRPM.fix.tbl")

cutoff=2
dfsample = read.table(samplefile, header = TRUE, sep = "\t", quote = "", stringsAsFactor = FALSE) #args[1]
dataraw = read.table(PASfile, header = TRUE, sep = "\t", quote = "", stringsAsFactor = FALSE)

addUPDN2<-function(data_dex_in, control_samples, treat_samples, cutoff, dfrpm){
  colnames(dfrpm)[1]='groupID'
  dfPRO=data_dex_in[which(data_dex_in$type=='UPS'),]
  dfDIS=data_dex_in[which(data_dex_in$type=='3UTR'),]

  dfALL=merge(dfPRO, dfDIS, by = "groupID", sort = FALSE)  
  names(dfALL) <- sub(".x$", ".pA1", sub(".y$", ".pA2", names(dfALL)))
  dfALL=merge(dfALL, dfrpm[dfrpm$type=='UPS',], by ='groupID', sort = FALSE)
  dfALL=merge(dfALL, dfrpm[dfrpm$type=='3UTR',], by ='groupID', sort = FALSE)
  names(dfALL) <- sub(".x$", ".pA1", sub(".y$", ".pA2", names(dfALL)))  
  
  ###calculate DeltaRA###
  control_col1=c(paste0("countData.num_",control_samples,".pA1"))
  treat_col1=c(paste0("countData.num_",treat_samples,".pA1"))
  
  if (length(control_samples)!=1){
    dfALL$avgControl_reads_pA1=rowMeans(dfALL[,control_col1], na.rm = TRUE)
  }
  if (length(control_samples)==1){
    dfALL$avgControl_reads_pA1=dfALL[,control_col1]
  }
  
  if (length(treat_samples)!=1){
    dfALL$avgTreat_reads_pA1=rowMeans(dfALL[,treat_col1], na.rm = TRUE)
  }
  
  
  if (length(treat_samples)==1){
    dfALL$avgTreat_reads_pA1=dfALL[,treat_col1]
  }  
  
  control_col2=c(paste0("countData.num_",control_samples,".pA2"))
  treat_col2=c(paste0("countData.num_",treat_samples,".pA2")) 
  
  if (length(control_samples)!=1){
    dfALL$avgControl_reads_pA2=rowMeans(dfALL[,control_col2], na.rm = TRUE)
  }
  if (length(control_samples)==1){
    dfALL$avgControl_reads_pA2=dfALL[,control_col2]
  }
  
  if (length(treat_samples)!=1){
    dfALL$avgTreat_reads_pA2=rowMeans(dfALL[,treat_col2], na.rm = TRUE)
  }
  if (length(treat_samples)==1){
    dfALL$avgTreat_reads_pA2=dfALL[,treat_col2]
  }
  
  dfALL$Delta_RA=(dfALL$avgTreat_reads_pA1/(dfALL$avgTreat_reads_pA1+dfALL$avgTreat_reads_pA2))-(dfALL$avgControl_reads_pA1/(dfALL$avgControl_reads_pA1+dfALL$avgControl_reads_pA2))
  
  ##Add 5 reads cutoff###
  index5Reads=which(dfALL$avgControl_reads_pA1>=cutoff & dfALL$avgTreat_reads_pA1>=cutoff & dfALL$avgControl_reads_pA2>=cutoff & dfALL$avgTreat_reads_pA2>=cutoff)  
  dfALL=dfALL[index5Reads,]
  
  ###Add UP, DN, NC Based on pvalue and fold###
  dfALL$type='NC'
  indexUP=which(dfALL$Delta_RA>0.05 & dfALL$pvalue.pA1<0.05)
  indexDN=which(dfALL$Delta_RA<(-0.05) & dfALL$pvalue.pA1<0.05)
  if(length(indexUP>0)){
    dfALL[indexUP,]$type='UP'
  }
  if(length(indexDN>0)){
    dfALL[indexDN,]$type='DN'
  }

  dfALL=dfALL[!duplicated(dfALL$groupID), ]
  return(dfALL)
}


for(i in 1:length(sample_c))
{
	#file setting
	filename=paste(sample_t[i], sample_c[i], sep='_')
	outputfile1=paste0(rootdir,'/',project,"/tbl/UPS/UPS_APA.",filename,'.cut2.tbl')

	###grab samples for each group###
	treat_samples=dfsample[which(dfsample$Group==sample_t[i]),]$Sample
	control_samples=dfsample[which(dfsample$Group==sample_c[i]),]$Sample

	###make a trim, remove pAs with total reads <=5###  
	numbercol=c(paste0("num_",treat_samples),paste0("num_",control_samples))  
	RPMcol= c(paste0("DRPM_",treat_samples),paste0("DRPM_",control_samples)) 
	dataraw$total_reads=rowSums(dataraw[,numbercol], na.rm = TRUE)
	datall=dataraw[which(dataraw$total_reads>5),]
	treatRPMcols=paste0("DRPM_",treat_samples)
	controlRPMcols=paste0("DRPM_",control_samples)
	treatNUMcols=paste0("num_",treat_samples)
	controlNUMcols=paste0("num_",control_samples)

	###calculate RPM###	
	datall$totalRPM=rowSums(datall[,c(treatRPMcols,controlRPMcols)], na.rm=TRUE)
	datall=datall[with(datall, order(gene_symbol, -totalRPM)),]

	###3'UTR###
	index3UTR=which(datall$pAtype_1=="F" | datall$pAtype_1=="L" | datall$pAtype_1=="M" | datall$pAtype_1=="S")
	df3UTR=datall[index3UTR,]
	df3UTR=df3UTR[,c('gene_symbol',numbercol,RPMcol)]
	df3UTR_2=df3UTR %>% group_by(gene_symbol) %>% summarise_all(list(sum))
	df3UTR_2=data.frame(df3UTR_2)
	df3UTR_2$type='3UTR'
	df3UTR_2$col=2
	
	###upstream###
	index_UPS=which(datall$pAtype_1=="intron" | datall$pAtype_1=="CDS")
	dfUPS=datall[index_UPS,]
	dfUPS=dfUPS[,c('gene_symbol',numbercol,RPMcol)]
	dfUPS_2=dfUPS %>% group_by(gene_symbol) %>% summarise_all(list(sum))
	dfUPS_2=data.frame(dfUPS_2)
	dfUPS_2$type='UPS'
	dfUPS_2$col=1
	
	dfallnew=rbind(dfUPS_2,df3UTR_2)	

	###prepare input for DEXseq###
	data_rename = dfallnew %>%
	  group_by(gene_symbol) %>%
	  mutate(count = n()) %>%
	  filter(count > 1) ##remove gene without APA
	data_rename_2 = data.frame(data_rename);dim(data_rename_2)

	data_rename_2$gene = paste(data_rename_2$gene_symbol, data_rename_2$col, sep = ":")
	data_rename_2$genomicData = data_rename_2$type
	data_rename_3 = data_rename_2[,c("gene","type")]
	data_subset = data_rename_2[,c("gene","gene_symbol","col",grep("^num_", colnames(dfallnew),value = TRUE))]
	data_subset=data_subset[with(data_subset, order(gene_symbol, col)),]
	
	if (length(control_samples)!=1){
		data_rename_2$controlRPM=rowMeans(data_rename_2[,controlRPMcols])
	}
	if (length(control_samples)==1){
		data_rename_2$controlRPM=data_rename_2[,controlRPMcols]
	}

	if (length(treat_samples)!=1){
		data_rename_2$treatRPM=rowMeans(data_rename_2[,treatRPMcols])
	}

	if (length(treat_samples)==1){
		data_rename_2$treatRPM=data_rename_2[,treatRPMcols]
	}
		
	
	dfrpm=data_rename_2
	dfrpm$Log2Ratio=log2(dfrpm$treatRPM/dfrpm$controlRPM)
	dfrpm=dfrpm[,c('gene_symbol', 'Log2Ratio','type', RPMcol)]	

	###Run DEXseq###
	data_in = data_subset[,c("gene","gene_symbol","col",treatNUMcols,controlNUMcols)]
	rownames(data_in) = data_in[,1]
	data_in_2 = data_in[,-c(1,2,3)]
	featureID <- as.factor(data_in$col)
	groupID <- as.factor(data_in$gene_symbol)
	design1 <- formula( ~ sample + exon + condition:exon )
	sampleTable = data.frame(row.names = c(controlNUMcols,
										   treatNUMcols),
							 condition = c(rep("Control",length(control_samples)),
										   rep("Treatment",length(treat_samples))),
							 replicate = c(1:length(control_samples), 1:length(treat_samples)),
							 type = c(rep("single",(length(control_samples)+length(treat_samples)))))			

	data_run = DEXSeqDataSet(data_in_2, sampleTable, design1, featureID, groupID)
	data_run = estimateSizeFactors(data_run)
	data_run = estimateDispersions(data_run)
	data_run = testForDEU(data_run)
	data_run = estimateExonFoldChanges(data_run, fitExpToVar = "condition")
	data_run_result = DEXSeqResults(data_run)

	data_run_result_2 = data.frame(data_run_result)
	data_run_result_2$gene = rownames(data_run_result_2)
	data_run_final = merge(data_run_result_2, data_rename_3, by = "gene", sort = FALSE)
	data_run_f=data_run_final[ , -which(names(data_run_final) %in% c("genomicData"))]
	data_run_f=data_run_f[,-c(3:6)]
	data_run_part_out = addUPDN2(data_run_f, control_samples, treat_samples, cutoff, dfrpm)
	
	names(data_run_part_out) = gsub('countData.num_','num_',names(data_run_part_out))
	numbercol_PD=grep("^num_", colnames(data_run_part_out),value = TRUE)
	RPMcol_PD= grep("^DRPM_", colnames(data_run_part_out),value = TRUE)
	dftblout=data_run_part_out[,c("groupID",
	"pvalue.pA1",
	numbercol_PD, RPMcol_PD, "Log2Ratio.pA1", "Log2Ratio.pA2",
	"Delta_RA", "type")]	
	dftblout$RED=dftblout$Log2Ratio.pA1-dftblout$Log2Ratio.pA2
	names(dftblout)[1:2]=c("gene_symbol","pvalue")
	names(dftblout)=gsub('type','status',names(dftblout))
	names(dftblout)=gsub('.pA1','_UPS',names(dftblout))
	names(dftblout)=gsub('.pA2','_exon3',names(dftblout))
	write.table(dftblout, outputfile1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


