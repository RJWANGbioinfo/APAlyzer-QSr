# APAlyzer-QSr: APA analysis toolkit for Quant-Seq REV

- [Overview](#Overview)
  * [Authors](#Authors)
  * [License](#License)
  * [Languages](#Languages)  
- [Requirements](#Requirements)
- [Quick Start](#Quick-Start)
- [Path setting example](#Path-setting-example)
- [Output](#Output)
- [Run the Toolkit Step by Step](#Run-the-Toolkit-Step-by-Step)



## Overview
APAlyzer-QSr is a automatic toolkit performing APA analysis using Quant-Seq REV data. In general it uses PolyA DB3 to clean the PASs indetified from the sequencing data.

### Authors
* **Ruijia Wang**
* **Bin Tian**

### License
If you plan to use the APAlyzer-QSr for-profit, you will need to purchase a license. Please contact rjwang.bioinfo@gmail.com for more information.

### Languages
Python 3 and R

### Supported Genome
Genomes that covered by PolyA DB3 including mm9, hg19 and rn5.


## Requirements
#### Python Modules
numpy
pandas
pylab
scipy
HTSeq
matplotlib

All these can be install through
```
pip install packagename
```

#### R library
GenomicRanges
Biostrings
dplyr
GenomicAlignments

All these can be install through
```
install.packages(packagename)
```
or
```
BiocManager::install(packagename)
```


### Optional Requirements
Although the PASs indetification and clean are using bam files as input, this toolkit also provide scripts to handel the raw fastq files from Quan-Seq REV.
The following packages are used for QC, trimming and mapping:

#### FastQC for quality control
https://www.bioinformatics.babraham.ac.uk/projects/fastqc

#### BBmap for trimming
https://sourceforge.net/projects/bbmap/

#### STAR as RNAseq mapper
https://github.com/alexdobin/STAR

#### Sambamba for bam converting
https://lomereiter.github.io/sambamba


### Reference file
All the reference files used for fastq trimming and PAS clean are stored in REF/ folder

## Quick-Start
The toolkit contains a shell scripts can simply run all steps in one shot aftering setting the path. You can first put all the fastq files under rootdir/project/rawfastq/, then run:
```
./raw.pip.example.sh
```

## Path-setting-example
And a expample of path setting in the shell:

##### 1)Define the number of threads will be used
```
THREADS=24
```

##### 2)Define the path of the 'scripts' folder in APAlyzer_qrev
```
scrdir=/xxx/APAlyzer_qrev/scripts/
```	

##### 3)Define project name
```				
project=project1	
```

##### 4)Define path of the rootdir, usually your project-dir will be rootdir/project/
``` 				
rootdir=/xxx/my_rootdir/
```	

##### 5)Define the genmoe version; mm9, hg19, or rn5	
```			
geno=mm9		
```

##### 6)Define path of the gemome folder used for STAR mapping 
```			
genodir=/xxx/genome/mm9_star/
```

##### 7)Define path of the REF folder containing reference files 
```					
refdir=/xxx/REF/		
```			

##### 8)Define the analysis design; 'YES' or 'NO'; setting whether the analysis design is single sample or multiple replicates
```
Reps='YES' 											
samplefile=/xxx/Samples/sample_list.txt	  # Only need when Reps='YES'	
```

##### 9)Define the treat groups
```			
treats="COM1_NT COM2_NT"
```

##### 10)Define the control groups
```
controls="COM1_TRT COM2_TRT"
```
Using the setting above, we can analysis a mouse mm9 dataset, and compare the APA in "COM1_TRT" vs "COM1_NT", and "COM2_TRT" vs "COM2_NT".


### Sample file
Sample file are only need when Reps='YES'. This file is a two-column table contain the information of the sample and group. For instance, in the above case, the sample file should looks like:

| Sample | Group |
| --- | --- |
| NT1_rep1 | COM1_NT |
| NT1_rep2 | COM1_NT |
| TRT1_rep1 | COM1_TRT |
| TRT1_rep2 | COM1_TRT |
| NT2_rep1 | COM2_NT |
| NT2_rep2 | COM2_NT |
| TRT2_rep1 | COM2_TRT |
| TRT2_rep2 | COM2_TRT |


## Output
#### 1. The output of the toolkits convering different files in different folders:

| File | Folder | Note |
| --- | --- | --- |
| QC file | rootdir/project/qccheck/ | fastq QC results |
| *.clipped.fastq | rootdir/project/fastq/ | trimmed fq file |
| *.sorted.bam | rootdir/project/rawsam/ | mapping bam file |
| mapping.summary.txt | rootdir/project/tbl/ | mapping summary |
| *.pA2gene_usage.DRPM.fix.tbl | rootdir/project/tbl/ | PAS expression profile |
| 3mostAPA.*.DRPM.fix.tbl | rootdir/project/tbl/3UTR/ | 3'UTR APA results |
| UPS.*.cut2.tbl | rootdir/project/tbl/UPS/ | Upstream APA results |
| *.png | rootdir/project/plot/ | Scatter plots of 3'UTR and UPS APA |
| PAS.nuc_freq.diff_vs_NC.xlsx | rootdir/project/tbl/ | nucleotide frequency (+/- 150 nt from PAS) of the sites that are changing vs same number of sites that do not change |
| *.motif_enrichment.diff_vs_NC.xlsx | rootdir/project/plot/ | motif enrichment (4-mers and 6-mers) within +/- 150 nts of the PAS of the sites that are changing vs sites that do not change |


#### 2. Columns of *.pA2gene_usage.DRPM.fix.tbl:
| Column | Description |
| --- | --- |
| pAid | ID of each PAS, shown as chromosome:Position:Strand |
| chromosome | chromosome ID of PAS |
| pA_pos | genomic position of PAS |
| strand | strand information of PAS |
| GENEID | Entrez Gene ID |
| gene_symbol | gene symbol |
| gene_desc | gene name description |
| gene_Biotype | type of the gene, protein-coding or various types of ncRNAs |
| LOCATION | Specific PAS annotation for ncRNAs, including 5'exon, 3'exon, singel(S) exon, other exon and intron  |
| region | PAS location in annotated genes, including 5'UTR, CDS, intron, 3'UTR and various types of ncRNAs. |
| pAtype_1 | PAS location in annotated genes, including 5'UTR, CDS, intron, 3'UTR and various types of ncRNAs. For PASs in 3'UTRs, they are further divided into First (F), Middle (M), and Last (L). If there is only one PAS in 3'UTR, it is called S. |
| ext | Whether PAS is located on an extended 3' end region beyond RefSeq/Ensembl annotations, the extrension region is annotated by polyA_DB3, Yes/No |
| num_* | reads count columns of each sample |
| DRPM_* | normalized expression columns count columns of each sample |


#### 3. Columns of 3mostAPA.*.DRPM.fix.tbl:
| Column | Description |
| --- | --- |
| gene_symbol | gene symbol |
| chromosome | chromosome ID of PAS |
| strand | strand information of PAS |
| pA_pos_pA1 | genomic position of proximal PAS |
| pA_pos_pA2 | genomic position of distal PAS |
| pvalue | p-value of alternative polyadenylation |
| num_* | reads count columns of each sample |
| DRPM_* | normalized expression columns count columns of each sample |
| RE_* | Relative expression of each gene in each sample, RE=distal/proximal |
| Log2Ratio_pA1 | log2FC of proximal PAS between two groups |
| Log2Ratio_pA2 | log2FC of distal PAS between two groups |
| Delta_RA | Delta Relative abundance |
| RED | Delta Relative expression, RED=Log2Ratio_pA2 - Log2Ratio_pA1  |
| pA1.pAutype | alternative polyadenylation pattern, 'UP' for shortening, 'DN' for lengthening, 'NC' for no change |


#### 4. PAS Region definitions in 'PAS.nuc_freq.diff_vs_NC.xlsx' and '*.motif_enrichment.diff_vs_NC.xlsx':
| Region | Range | PAS |
| --- | --- | --- |
| prx_region_1 | -150 ~ -51 | proximal PAS |
| prx_region_2 | -50 ~ -1 | proximal PAS |
| prx_region_3 | +1 ~ +50 | proximal PAS |
| prx_region_4 | +51 ~ +150 | proximal PAS |
| dis_region_1 | -150 ~ -51 | distal PAS |
| dis_region_2 | -50 ~ -1 | distal PAS |
| dis_region_3 | +1 ~ +50 | distal PAS |
| dis_region_4 | +51 ~ +150 | distal PAS |


## Run-the-Toolkit-Step-by-Step
#### 1. QC check, trimming and mapping Quan-Seq REV
```
python scripts/step1_2_3.QC_and_mapping.qcREV.py \
				--rootdir ROOTPATH  \
				--project PRJNAME  \
				--genodir GENOMEPATH  \
				--refdir REFPATH  \
				--threads NUM
```

#### 2. Summary mapping results
```
python scripts/step4.STAR_log_summarizer.qcREV.py --rootdir ROOTPATH --project PRJNAME
```

#### 3. Identify the 3' End Alignment Position
```
python scripts/step5.LAP.hunter.qcREV.py --rootdir ROOTPATH --project PRJNAME
```

#### 4. Clean PASs using PolyA DB3
```
Rscript scripts/step6_1.LAP2PAS.R $rootdir $project $refdir $geno
```

#### 5. Annotation of cleaned PAS
```
Rscript scripts/step6_2.Quant_R.pas2gene.builder.R $rootdir $project $refdir $geno
```

#### 6. Analysis of 3'UTR APA and upstream APA (no replicate design)
```
python scripts/step7_1.qcREV.3UTR.APA.tbler.norep.py --project $project --rootdir $rootdir --genome $geno --control $controls --treatment $treats
python scripts/step7_2.qcREV.UR.APA.tbler.norep.py --project $project --rootdir $rootdir --genome $geno --control $controls --treatment $treats

```


#### 6. Analysis of 3'UTR APA and upstream APA (replicate design)
```
Rscript scripts/step7_1.qcREV.3UTR.APA.tbler.reps.R $rootdir $project $samplefile $geno "$treats" "$controls"
Rscript scripts/step7_2.qcREV.UR.APA.tbler.reps.R $rootdir $project $samplefile $geno "$treats" "$controls"

```

#### 7. Plot 3'UTR APA and UPS APA pattern using scatter plots
```
python scripts/step8_1.qcREV.3UTR.scatter.plotter.py --project $project --rootdir $rootdir --control $controls --treatment $treats
python scripts/step8_2.qcREV.UPS.scatter.plotter.py --project $project --rootdir $rootdir --control $controls --treatment $treats
```

#### 8. calculate nucleotide frequency of PAS regions
```
Rscript $scrdir/step9.ACGT_FREQ.R $rootdir/$project/tbl/ $geno
```

#### 9. motif enrichment of PAS regions
```
Rscript $scrdir/step10.motif_enrichment.R $rootdir/$project/tbl/ $geno
```
