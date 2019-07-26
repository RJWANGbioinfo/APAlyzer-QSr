#!/bin/bash
###user setting section###
THREADS=24
scrdir=/xxx/APAlyzer_qrev/
project=project1
rootdir=/xxx/my_rootdir/
genodir=/xxx/genome/mm9_star/
refdir=/xxx/REF/
geno=mm9
samplefile=/xxx/Samples/sample_list.txt
Reps='YES'  # 'YES' or 'NO'
treats="COM1_NT COM2_NT"
controls="COM1_TRT COM2_TRT"

############################### Main steps ##################################################

###Step1, 2, and 3####
echo "--------------Step1, 2, and 3 ----------------"
python $scrdir/step1_2_3.QC_and_mapping.qcREV.py \
		--project $project \
		--rootdir $rootdir \
		--refdir $refdir \
		--genodir $genodir \
		--threads $THREADS		

###Step4 Summary mapping stat####
echo "--------------Step4: Summary mapping stat----------------"	
python $scrdir/step4.STAR_log_summarizer.qcREV.py --project $project --rootdir $rootdir

###step5 capture LAP####
echo "--------------step5 capture LAP----------------"	
python $scrdir/step5.LAP.hunter.qcREV.py --project $project --rootdir $rootdir

###step6 clean and annotation using ploA_DB3####
echo "--------------step6 clean and annotation using ploA_DB3----------------"	
Rscript $scrdir/step6_1.LAP2PAS.R $rootdir $project $refdir $geno
Rscript $scrdir/step6_2.Quant_R.pas2gene.builder.R $rootdir $project $refdir $geno

###Step7 call 3'UTR and UR APA ####
echo "--------------step7 call 3'UTR and UR APA----------------"
if [ $Reps == 'NO' ]
then
	python $scrdir/step7_1.qcREV.3UTR.APA.tbler.norep.py --project $project --rootdir $rootdir --genome $geno --control $controls --treatment $treats
	python $scrdir/step7_2.qcREV.UR.APA.tbler.norep.py --project $project --rootdir $rootdir --genome $geno --control $controls --treatment $treats
else
	Rscript $scrdir/step7_1.qcREV.3UTR.APA.tbler.reps.R $rootdir $project $samplefile $geno "$treats" "$controls"
	Rscript $scrdir/step7_2.qcREV.UR.APA.tbler.reps.R $rootdir $project $samplefile $geno "$treats" "$controls"
fi
	
###Step8 Plot call 3'UTR and UR APA ####
echo "--------------step8 Plot call 3'UTR and UR APA----------------"
python $scrdir/step8_1.qcREV.3UTR.scatter.plotter.py --project $project --rootdir $rootdir --control $controls --treatment $treats
python $scrdir/step8_2.qcREV.UPS.scatter.plotter.py --project $project --rootdir $rootdir --control $controls --treatment $treats


