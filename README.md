# APAlyzer-QSr: APA analysis toolkit for Quant-Seq REV

- [Overview](#Overview)
  * [Authors](#Authors)
  * [License](#License)
  * [Languages](#Languages)  
- [Requirements](#Requirements)
- [Quick Start](#Quick-Start)
- [Run-the-pipline](#Run-the-pipline)
- [Output](#Output) 
- [Quick-start-example](#Quick-start-example) 


## Overview
APAlyzer-QSr is a automatic toolkit performing APA analysis using Quant-Seq REV data. In general it uses PolyA DB3 to clean the PASs indetified from the sequencing data.

### Authors
* **Ruijia Wang**
* **Bin Tian**

### License
This project is licensed under the LGPL-3 License.

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
And a expample of path setting in the shell:

##### 1)Define the number of threads will be used
```
THREADS=24
```

##### 2)Define the path of the APAlyzer_qrev
```
scrdir=/xxx/APAlyzer_qrev/
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
The output of the toolkits convering different files in different folders:

| File | Folder | Note |
| --- | --- | --- |
| QC file | rootdir/project/qccheck/ | fastq QC results |
| trimmed fq file | rootdir/project/fastq/ | mapping file |
| bam file | rootdir/project/rawsam/ | mapping file |
| mapping summary | rootdir/project/tbl/ | mapping summary |
| 3'UTR APA results | rootdir/project/tbl/3UTR | tables of 3'UTR APA |
| Upstream APA results | rootdir/project/tbl/UPS | tables of UPS APA |
| Scatter plots of 3'UTR and UPS APA | rootdir/project/plot/ | scatter plots (.png) |




A file named as "sample_list.txt", which is a tab separated file containing 4 colums: 'Run', 'LibraryLayout', 'samplename', 'condition'.
'Run' is SRA ID for each sample usually start with "SRR". 
'LibraryLayout' define the sequencing library type: 'SINGLE' is single-end; PAIRED' is pair-end. 
'samplename' define the really sample name of this file, e.g., 'KD_GENE_REP1' is the samplename of SRR12345.
'condition' define group name of the sample, for example, 'KD' or 'NT'. The information of each sample can be easily obtained from EBI (https://www.ebi.ac.uk/ena/browse).

## Run-the-pipline
To run the pipline a) put .py and .R file into sample folder, and b) put all the reference file in the same folder (e.g., the REF folder)
```
python all_steps.RNAseq_APA.pip.py --rootdir ROOTPATH  \
				--project PRJNAME  \
				--sradir SRAPATH  \
				--genodir GENOMEPATH  \
				--refdir REFPATH  \
				--genome GENO  \
				--KEEPRAW NO  \
				--threads NUM
```

About the options:
```
--project(required); define the project name, e.g., GEO123456
--rootdir(required); define the rootdir, e.g., /scratch/user/  
--sradir(required); define the sra dir (~/sra/sra/), e.g., /scratch/user/soft/sra/srt
--genodir(required); define the gemome dir for mapping, e.g., /scratch/xxxx/genomes/mm9_star
--refdir(required); define the ref dir; this dir should contain all the reference files, correspoding files of mm9 and hg19 have been store in REF/                
--threads(required); define # of threads, e.g., 24
--genome(Optional); define genome version, e.g., hg19, default is mm9
--KEEPRAW(Optional); whether keep the downloaded fastq and bam files; default is NO, if set to YES, then the fastq and bam files will be store at rootdir/project/KEEPRAW/
```

Note: 
"sample_list.txt" MUST be put into rootdir/project/ 

## Output
In general, rootdir/project/ is the overall output folder.
The combine results will be stored in "report" folder, including a expression file ".CDS.allsample.txt", a 3'UTR APA file ".UTR.allsample.txt", and a IPA file ".IPA_LE.allsample.txt"
The correspoding files for each sample are stored in "rawout" folder.

## Quick-start-example
Let's try a really mouse RNA-seq case in GEO(GSE112698). First creat the project folder and put the sample_list.txt in to the project folder:
```
cd /scratch/your_account/
mkdir GSE112698
```

Then start to run the pipline use Rutgers-Perceval with 24 threads:
```
python all_steps.RNAseq_APA.pip.py --rootdir /scratch/your_account/  \
				--project GSE112698  \
				--sradir /scratch/your_account/sra/sra/  \
				--genodir /scratch/xxxx/genomes/mm9_star/  \
				--refdir /scratch/your_account/REF/  \
				--genome mm9  \
				--KEEPRAW NO  \
				--threads 24
```

The combine results can be found in "/scratch/your_account/GSE112698/report" folder, including:

"GSE112698.CDS.allsample.txt", "GSE112698.UTR.allsample.txt", "GSE112698.IPA_LE.allsample.txt".

The correspoding files for each sample are stored in "/scratch/your_account/GSE112698/rawout" folder.
