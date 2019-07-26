# GEO RNA-seq APA profiler: Automatic pipline profiling APA,IPA and gene expression using RNAseq data from GEO

- [Overview](#Overview)
  * [Feature](#Version-feature)
  * [Authors](#Authors)
  * [License](#License)
  * [Languages](#Languages)  
- [Requirements](#Requirements)
- [Run-the-pipline](#Run-the-pipline)
- [Output](#Output) 
- [Quick-start-example](#Quick-start-example) 


## Overview
GEO RNA-seq APA profiler is a automatic pipline performing APA,IPA and gene expression profiling using RNAseq data from GEO. In general it uses SRA toolkit, STAR, and APAanalizer to finished the entire work.

### Version feature
1) RNA-seq sample is processed one by one to address the issues of sequencing library diversity. E.g, within a GEO project, some samples are pair-end, some are single-end, some are forward sequencing, some are reverse sequencing. 
2) Added a extra check section to address the downloading failed issues.
3) User now can select whether to keep the raw fastq and bam files.

### Authors
* **Ruijia Wang**
* **Bin Tian**

### License
This project is licensed under the LGPL-3 License.

### Languages
Python 2.7 and R 3.5.0

## Requirements

### Pythod_Modules as 'import section' 
The required pythod modules can be found in in all_steps.RNAseq_APA.pip.py

### APAlyzer in bioconductor 
```
R CMD INSTALL APAlyzer.tar.gz
```

### STAR as RNAseq mapper
https://github.com/alexdobin/STAR

### SRA toolkit
https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/


Add it to the .bashrc
```
# added path of sra toolkits
export PATH="/scratch/PATH_TO_THE/sratoolkit.2.8.2-1-centos_linux64/bin:$PATH"
```
Config the tools:
```
./vdb-config -i
```



### RSeQ
http://rseqc.sourceforge.net

### Sample file
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
