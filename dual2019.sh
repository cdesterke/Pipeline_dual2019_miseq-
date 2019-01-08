#!/bin/bash
##miseq pipeline with GATK 3.8 and dual calling HC and Vardict
##Christophe Desterke_january 2019
## sh dual2019.sh R1.fastq R2.fastq 

variable=${1}
if [ -z "${variable}" ]
then 
	echo "TEXTTAB file not passed as parameter"
	exit 1
fi

variable=${2}
if [ -z "${variable}" ]
then 
	echo "TEXTTAB file not passed as parameter"
	exit 1
fi

nom_fichier=$(echo $1 | sed -re 's/(.*).fastq/\1/')
echo "--------------------"


log_file="log_file.log" 

./fastqc $1
./fastqc $2

echo "analyse du ficher = $nom_fichier" >> $log_file

date >> $log_file

echo "step1 : bowtie2 alignement"
UPTIME1=`uptime`
## aligned fastq with bowti2 and gap parameters
bowtie2 -x genome -1 $1 -2 $2 -S aligned.sam --rdg 5,1 --rfg 5,1 --dpad 100


## bam compression and sorted with samtools
samtools view -bSu aligned.sam > aligned.bam 
samtools sort aligned.bam -o aligned.sorted.bam

echo "--------------------"
echo "step2 : bam compression and sorted samtools"
UPTIME2=`uptime`

##picardtools remove duplicates
PicardCommandLine MarkDuplicates INPUT=aligned.sorted.bam  OUTPUT=marked.bam METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
echo "--------------------"
echo "step3 : picard remove duplicates"
UPTIME3=`uptime`


## picard ajouter un ReadGroup pour rendre compatible GATK
PicardCommandLine AddOrReplaceReadGroups I=marked.bam O=markedRG.bam LB=whatever PL=illumina PU=whatever SM=whatever
echo "--------------------"
echo "step4 : picard add read group"
UPTIME4=`uptime`

## restriction of the bam to the myeloid panel
bedtools intersect -abam markedRG.bam -b sophia.bed > myeloid.bam
PicardCommandLine ReorderSam I=markedRG.bam O=markedRGorder.bam REFERENCE=genome.fa
echo "--------------------"
echo "step5 : bedtools myeloid restriction"
UPTIME5=`uptime`

## samtools filtration Q30 quality
samtools view -h -b -q 30 myeloid.bam -o myeloidQ30.bam 
echo "--------------------"
echo "step6 : quality Q30 filtration"
UPTIME6=`uptime`

## samtools index of the bam file
samtools index myeloidQ30.bam
echo "--------------------"
echo "step7 : index of the bam"
UPTIME7=`uptime`


rm aligned.sam
rm aligned.bam
rm aligned.sorted.bam
rm markedRG.bam
rm marked.bai
rm marked.bam
rm myeloid.bam


## GATK "local realignment around indels"

# creation d'une table des possibles indels
java -Xmx4g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -known Mills_and_1000G_gold_standard.indels.hg19.vcf -R genomeHG19.fa -o marked.bam.list -I myeloidQ30.bam --filter_reads_with_N_cigar

# realigne les reads autour de ces possibles indels
java -Xmx4g -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar -I myeloidQ30.bam -R genomeHG19.fa -T IndelRealigner -known Mills_and_1000G_gold_standard.indels.hg19.vcf -targetIntervals marked.bam.list -o marked.realigned.bam --filter_reads_with_N_cigar

echo "--------------------"
echo "step8 : GATK realignement around indels"
UPTIME8=`uptime`

## picard fixmate (fixation de l'information mate comme l'alignement peut etre modifié pendant le processus de realignement)
PicardCommandLine FixMateInformation I=marked.realigned.bam O=fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

echo "--------------------"
echo "step9 : Picard commandline fixmate"
UPTIME9=`uptime`

## base recalibrator
java -Xmx4G -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I fixed.bam -R genomeHG19.fa -knownSites dbsnp.vcf -o recaldata.table

java -Xmx4G -jar GenomeAnalysisTK.jar -T PrintReads -R genomeHG19.fa -I fixed.bam -BQSR recaldata.table -o recal.bam

echo "--------------------"
echo "step10 : GATK base recalibration"
UPTIME10=`uptime`

## samtools 
samtools sort recal.bam -o recal.srt.bam

samtools index recal.srt.bam

echo "--------------------"
echo "step11 : Samtools sort and index"
UPTIME11=`uptime`




#java -jar GenomeAnalysisTK.jar -T MuTect2 -R genomeHG19.fa -I:tumor recal.srt.bam -o mutect2.vcf
#python Platypus.py callVariants --bamFiles=recal.srt.bam --refFile=genomeHG19.fa --regions=myeloid.txt --output=platypus.vcf
#./lofreq call -f genomeHG19.fa -o lofreqsnp.vcf recal.srt.bam -l sophia.bed
#./lofreq call -f genomeHG19.fa -o lofreqindels.vcf recal.srt.bam --call-indels -l sophia.bed
#samtools mpileup -uf genomeHG19.fa recal.961809040944_S8_L001_R1_001.fastq.srt.bam | bcftools call -mv -Oz > call.vcf
#samtools mpileup -B -f genomeHG19.fa recal.srt.bam | java -jar VarScan.v2.3.9.jar pileup2snp --output-vcf 1 > snp.vcf
#samtools mpileup -B -f genomeHG19.fa recal.srt.bam | java -jar VarScan.v2.3.9.jar pileup2indel --output-vcf 1 > indel.vcf
#samtools mpileup -B -f genomeHG19.fa recal.srt.bam | java -jar VarScan.v2.3.9.jar pileup2cnv --output-vcf 1 > cnv.vcf


java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R genomeHG19.fa -I recal.srt.bam -o hc.vcf 

echo "--------------------"
echo "step12 : Haplotype caller variant calling"
UPTIME12=`uptime`



AF_THR="0.01" # minimum allele frequency

perl vardict.pl -G genomeHG19.fa -f $AF_THR -N patient -b recal.srt.bam -c 1 -S 2 -E 3 -g 4 sophia.bed | perl var2vcf_valid.pl -N patient -E -f $AF_THR > vardict.vcf
echo "--------------------"
echo "step13 : Vardict variant calling"
UPTIME13=`uptime`




rm marked.realigned.bai
rm myeloidQ30.bam
rm myeloidQ30.bam.bai
rm marked.bam.list
rm marked.realigned.bam
rm fixed.bam
rm fixed.bai
rm recaldata.table
rm recal.bam
rm recal.bai



##conversion du VCF hc en AVinput ANNOVAR
perl convert2annovar.pl -format vcf4  hc.vcf -outfile outav

UPTIME14=`uptime`

perl table_annovar.pl outav humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,dbnsfp30a,exac03,kaviar_20150923,esp6500siv2_all,avsnp147,clinvar_20160302,cosmic70,icgc21 -operation g,r,f,f,f,f,f,f,f,f -nastring . -csvout -polish -xreffile example/gene_fullxref.txt

rm outav

UPTIME15=`uptime`



##conversion du VCF vardict en AVinput ANNOVAR
perl convert2annovar.pl -format vcf4  vardict.vcf -outfile outav2

UPTIME16=`uptime`

perl table_annovar.pl outav2 humandb/ -buildver hg19 -out myanno2 -remove -protocol refGene,cytoBand,dbnsfp30a,exac03,kaviar_20150923,esp6500siv2_all,avsnp147,clinvar_20160302,cosmic70,icgc21 -operation g,r,f,f,f,f,f,f,f,f -nastring . -csvout -polish -xreffile example/gene_fullxref.txt

rm outav2

UPTIME17=`uptime`






##log

echo "---------------------------------">> $log_file
echo "log of the pipeline preGATK">> $log_file
echo "analyse du ficher = $nom_fichier">> $log_file
echo "--------------------">> $log_file
cal
echo "--------------------"
echo "step1 bowtie2 alignement">> $log_file
echo "Uptime = $UPTIME1">> $log_file
echo "-------------"
echo "step2 bam compression and sorted samtools">> $log_file
echo "Uptime = $UPTIME2">> $log_file
echo "------------"
echo "step3 picard remove duplicates">> $log_file
echo "Uptime = $UPTIME3">> $log_file
echo "------------"
echo "step4 picard add read group">> $log_file
echo "Uptime = $UPTIME4">> $log_file
echo "------------"
echo "step5 bedtools myeloid restriction">> $log_file
echo "Uptime = $UPTIME5">> $log_file
echo "------------"
echo "step6 quality Q30 filtration">> $log_file
echo "Uptime = $UPTIME6">> $log_file
echo "------------"
echo "step7 index of the bam">> $log_file
echo "Uptime = $UPTIME7">> $log_file
echo "------------"
echo "step8 GATK 3.8 realignement around indels">> $log_file
echo "Uptime = $UPTIME8">> $log_file
echo "------------"
echo "step9 Picard commandline fixmate">> $log_file
echo "Uptime = $UPTIME9">> $log_file
echo "------------"
echo "step10 GATK 3.8 base recalibration">> $log_file
echo "Uptime = $UPTIME10">> $log_file
echo "------------"
echo "step11 Samtools sort and index">> $log_file
echo "Uptime = $UPTIME11">> $log_file
echo "------------"
echo "step12 Haplotype caller variant calling">> $log_file
echo "Uptime = $UPTIME12">> $log_file
echo "------------"
echo "step13 VarDict variant calling">> $log_file
echo "Uptime = $UPTIME13">> $log_file
echo "------------"

echo "---------------------------------">> $log_file
echo "log of the pipeline in ANNOVAR">> $log_file
echo "--------------------">> $log_file
echo "step14 VCF convertion ok for haplotype caller">> $log_file
echo "Uptime = $UPTIME14">> $log_file
echo "------------"
echo "step15 ANNOVAR annotation ok for haplotype calller">> $log_file
echo "Uptime = $UPTIME15">> $log_file
echo "------------">> $log_file
echo "step16 VCF convertion ok for VarDict">> $log_file
echo "Uptime = $UPTIME16">> $log_file
echo "------------"
echo "step17 ANNOVAR annotation ok for VarDict">> $log_file
echo "Uptime = $UPTIME17">> $log_file
echo "------------">> $log_file



cat myanno.hg19_multianno.csv myanno2.hg19_multianno.csv > mytotal.csv



mkdir RESULTS
mv log_file.log RESULTS
mv recal.srt.bam.bai RESULTS
mv recal.srt.bam RESULTS
mv myanno.hg19_multianno.csv RESULTS
mv hc.vcf RESULTS
mv hc.vcf.idx RESULTS

mv myanno2.hg19_multianno.csv RESULTS
mv vardict.vcf RESULTS

mv mytotal.csv RESULTS

cd RESULTS


mv log_file.log $(echo log_file.log | sed "s/\./".$nom_fichier"\./")
mv recal.srt.bam $(echo recal.srt.bam | sed "s/\./".$nom_fichier"\./")
mv recal.srt.bam.bai $(echo recal.srt.bam.bai | sed "s/\./".$nom_fichier"\./")
mv myanno.hg19_multianno.csv $(echo myanno.hg19_multianno.csv | sed "s/\./".$nom_fichier"\./")
mv myanno2.hg19_multianno.csv $(echo myanno2.hg19_multianno.csv | sed "s/\./".$nom_fichier"\./")


mv mytotal.csv $(echo mytotal.csv | sed "s/\./".$nom_fichier"\./")

mv hc.vcf $(echo hc.vcf | sed "s/\./".$nom_fichier"\./")
mv hc.vcf.idx $(echo hc.vcf.idx | sed "s/\./".$nom_fichier"\./")
mv vardict.vcf $(echo vardict.vcf | sed "s/\./".$nom_fichier"\./")


cd ..

mv RESULTS RESULTS_$nom_fichier

exit 0

