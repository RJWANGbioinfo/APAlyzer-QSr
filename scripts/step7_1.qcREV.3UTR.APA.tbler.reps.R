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

cutoff=5
dfsample = read.table(samplefile, header = TRUE, sep = "\t", quote = "", stringsAsFactor = FALSE) #args[1]
dataraw = read.table(PASfile, header = TRUE, sep = "\t", quote = "", stringsAsFactor = FALSE)

addUPDN<-function(data_dex_in, control_samples, treat_samples, cutoff, dfrpm){
  ###Separate dataset###
  data_dex_in$sepcol=data_dex_in$pAid
  data_dex_in=separate(data = data_dex_in, col = sepcol, into = c("Chrom", "Pos","Strand"), sep = ":")
  data_run_sub=data_dex_in[,c('groupID', 'Strand')]
  data_run_sub=data_run_sub[!duplicated(data_run_sub$groupID), ]
  dfMinMAX=data_dex_in  %>%  group_by(groupID)  %>%  summarise(minPos = min(Pos), maxPos = max(Pos))
  dfMinMAX=merge(dfMinMAX, data_run_sub, by = "groupID", sort = FALSE)
  indexN=which(dfMinMAX$Strand=='-')
  indexP=which(dfMinMAX$Strand=='+')
  dfMinMAX$promixal=''
  dfMinMAX$distal=''
  
  dfMinMAX[indexP,]$promixal=dfMinMAX[indexP,]$minPos
  dfMinMAX[indexP,]$distal=dfMinMAX[indexP,]$maxPos
  dfMinMAX[indexN,]$promixal=dfMinMAX[indexN,]$maxPos
  dfMinMAX[indexN,]$distal=dfMinMAX[indexN,]$minPos
  dfPRO=dfMinMAX[,c('groupID', 'promixal')]
  colnames(dfPRO)[2]='Pos'
  dfDIS=dfMinMAX[,c('groupID', 'distal')]
  colnames(dfDIS)[2]='Pos'
  data_slim=data_dex_in

  dfPRO=merge(data_slim, dfPRO, by = c("groupID", "Pos"), sort = FALSE)
  dfDIS=merge(data_slim, dfDIS, by = c("groupID", "Pos"), sort = FALSE)
  dfALL=merge(dfPRO, dfDIS, by = "groupID", sort = FALSE)
  
  names(dfALL) <- sub(".x$", ".pA1", sub(".y$", ".pA2", names(dfALL)))
  dfALL=merge(dfALL, dfrpm, by.x='pAid.pA1',by.y='pAid', sort = FALSE)
  dfALL=merge(dfALL, dfrpm, by.x='pAid.pA2',by.y='pAid', sort = FALSE)
  names(dfALL) <- sub(".x$", ".pA1", sub(".y$", ".pA2", names(dfALL)))
  
  subsamples=unique(c(treat_samples,control_samples))
  for(subsample in subsamples)
  {

	subsamplePA1= c(paste0("DRPM_",subsample,".pA1"))
	subsamplePA2= c(paste0("DRPM_",subsample,".pA2"))
    subsampleRE= c(paste0("RE_",subsample))
	dfALL[,subsampleRE]=log2(dfALL[,subsamplePA2]/dfALL[,subsamplePA1])
  }
  
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
  
  ##reads cutoff###
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
	outputfile1=paste0(rootdir,'/',project,"/tbl/3UTR/3mostAPA.",filename,'.DRPM.fix.tbl')

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
	dfrpm=datall[,c('pAid',RPMcol)]

	if (length(control_samples)!=1){
		dfrpm$controlRPM=rowMeans(dfrpm[,controlRPMcols])
	}
	if (length(control_samples)==1){
		dfrpm$controlRPM=dfrpm[,controlRPMcols]
	}

	if (length(treat_samples)!=1){
		dfrpm$treatRPM=rowMeans(dfrpm[,treatRPMcols])
	}

	if (length(treat_samples)==1){
		dfrpm$treatRPM=dfrpm[,treatRPMcols]
	}
		
	dfrpm$Log2Ratio=log2(dfrpm$treatRPM/dfrpm$controlRPM)
	dfrpm=dfrpm[,c('pAid', 'Log2Ratio', RPMcol)]

	###3'UTR###
	index3UTR=which(datall$pAtype_1=="F" | datall$pAtype_1=="L" | datall$pAtype_1=="M" | datall$pAtype_1=="S")
	datall=datall[index3UTR,]
	datall=datall[with(datall, order(gene_symbol, -totalRPM)),]

	###prepare input for DEXseq###
	data_rename = datall %>%
	  group_by(gene_symbol) %>%
	  mutate(col = 1:n()) %>% ##create featureID
	  mutate(count = n()) %>%
	  filter(count > 1) ##remove gene without APA
	data_rename_2 = data.frame(data_rename);dim(data_rename_2)

	###Pick up Top2 pAs###
	top2index = which(data_rename_2$col==1 | data_rename_2$col==2)
	data_rename_2 = data_rename_2[top2index,]
	data_rename_2$gene = paste(data_rename_2$gene_symbol, data_rename_2$col, sep = ":")
	data_rename_2$genomicData = data_rename_2$pAid
	data_rename_3 = data_rename_2[,c("gene","pAid")]
	data_subset = data_rename_2[,c("gene","gene_symbol","col",grep("^num_", colnames(datall),value = TRUE))]

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
	data_run_part_out = addUPDN(data_run_f, control_samples, treat_samples, cutoff, dfrpm)	
	names(data_run_part_out) = gsub('countData.num_','num_',names(data_run_part_out))
	numbercol_PD=grep("^num_", colnames(data_run_part_out),value = TRUE)
	RPMcol_PD= grep("^DRPM_", colnames(data_run_part_out),value = TRUE)
	REcol_PD= grep("^RE_", colnames(data_run_part_out),value = TRUE)
	dftblout=data_run_part_out[,c("groupID", "Pos.pA1", "Pos.pA2",
	"pvalue.pA1",
	numbercol_PD, RPMcol_PD,REcol_PD, "Log2Ratio.pA1", "Log2Ratio.pA2",
	"Delta_RA", "type")]
	dftblout$RED=dftblout$Log2Ratio.pA2-dftblout$Log2Ratio.pA1
	names(dftblout)[1]="gene_symbol"
	names(dftblout)[4]="pvalue"
	names(dftblout)=gsub('type',paste0('pA1.pAutype_',sample_t[i],'_',sample_c[i]),names(dftblout))
	names(dftblout)=gsub('.pA1','_pA1',names(dftblout))
	names(dftblout)=gsub('.pA2','_pA2',names(dftblout))
	write.table(dftblout, outputfile1, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

