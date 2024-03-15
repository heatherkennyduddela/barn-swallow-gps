##########################################################################################
#																						 #
#																						 #	
#																						 #
#			  	Workflow for analysis of fertilization networks for Colorado			 #
#		 		 	  and Ohio barn swallows using low coverage whole					 #
#							  	  genome sequencing data							 	 #
#																						 #
#																						 #
#																						 #
##########################################################################################

Created 11.28.2022 by Drew Schield

This README contains details on processing and analysis steps taken in our study of barn
swallow fertilization networks for populations from Colorado and Ohio. The workflow will
focus on low coverage (i.e., ~2x) whole-genome resequencing data from these populations.

General overview of sections in this README:
i. WGS data processing, mapping, and variant calling
1. Parentage analysis using AlphaAssign
	


##########################################################################################
#																						 #
#																						 #
#		  		   i. WGS data processing, mapping, and variant calling					 #
#																						 #
#																						 #
##########################################################################################

This section details the steps taken to process and filter raw fastq data, map filtered
data to the reference genome, and call variant sites among samples. These steps will include
previously generated data (i.e., used in Schield et al. 2021 Mol Ecol) with our newly-
generated WGS data.

General specifics:
	i.i Processing raw data with Trimmomatic
	i.ii Mapping filtered data to reference genome with BWA
	i.iii Variant calling and filtering using GATK and bcftools


##########################################################################################
#																						 #
#		  				 i.i Processing raw data with Trimmomatic						 #
#																						 #
##########################################################################################

We'll use Trimmomatic to filter and quality trim our raw read data.

We will impose these filters to trim reads:
	* Remove 5' end bases if quality is below 20
	* Remove 3' end bases if quality is below 20
	* Minimum read length = 32
	* Remove reads if average quality is < 30

Raw input data are in `/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/01.RawData/`.
Individual read files are in sample-specific subdirectories.

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd /data2/
$mkdir hirundo_fertilization_network
$cd hirundo_fertilization_network
$mkdir fastq_filtered
$mkdir log

------------------------------------------------------------------------------------------
Generated sample list:
------------------------------------------------------------------------------------------

$ls -lh /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/01.RawData/ | tail -n+2 | cut -d" " -f9 > sample.all.list

------------------------------------------------------------------------------------------
Ran time test on first sample:
------------------------------------------------------------------------------------------

time trimmomatic PE -phred33 -threads 4 /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/01.RawData/CO_00632/CO_00632_ZKDN220001168-1A_H233CDSX5_L4_1.fq.gz /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/01.RawData/CO_00632/CO_00632_ZKDN220001168-1A_H233CDSX5_L4_2.fq.gz ./fastq_filtered/test_1_P.trim.fq.gz ./fastq_filtered/test_1_U.trim.fq.gz ./fastq_filtered/test_2_P.trim.fq.gz ./fastq_filtered/test_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30

Time output:
real	6m20.788s
user	18m24.133s
sys	0m28.491s

This would take several days to run on a single sample list. We'll parallelize across a
number of separate lists.

How many parallel runs to have going? We Could do 12 x 4 threads = 48, leaving some
headroom for other processes. Doing so should complete the job in 7-8 hours:

380 seconds * 878 samples / 60 seconds per minute = 5560 minutes / 60 minutes per hour = 92 hours / 12 processes = 7.7 hours.

------------------------------------------------------------------------------------------
Formatted 12 sample subset lists:
------------------------------------------------------------------------------------------

sample.*.list

* = 1-12

------------------------------------------------------------------------------------------
Wrote script to run Trimmomatic on each sample per list:
------------------------------------------------------------------------------------------

As noted above, some individuals have four input files due to being sequenced on two
lanes. We'll run if statements to deal with cases where there are two versus four input
files.

runTrimmomatic.sh <sample-list>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
samples=$1
rawpath="/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/01.RawData"
for indv in `cat $samples`; do 
	inputs=`ls -lh $rawpath/$indv/*.fq.gz | wc -l`
	if [ "$inputs" -eq 2 ]; then
		trimmomatic PE -phred33 -threads 4 $rawpath/$indv/${indv}_*_1.fq.gz $rawpath/$indv/${indv}_*_2.fq.gz ./fastq_filtered/${indv}_1_P.trim.fq.gz ./fastq_filtered/${indv}_1_U.trim.fq.gz ./fastq_filtered/${indv}_2_P.trim.fq.gz ./fastq_filtered/${indv}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
    fi
    # Write temporary concatenated input if individual was sequenced on two lanes, then remove
	if [ "$inputs" -eq 4 ]; then
		cat $rawpath/$indv/${indv}_*_1.fq.gz > ./tmp.${indv}_1.fq.gz
		cat $rawpath/$indv/${indv}_*_2.fq.gz > ./tmp.${indv}_2.fq.gz
		trimmomatic PE -phred33 -threads 4 ./tmp.${indv}_1.fq.gz ./tmp.${indv}_2.fq.gz ./fastq_filtered/${indv}_1_P.trim.fq.gz ./fastq_filtered/${indv}_1_U.trim.fq.gz ./fastq_filtered/${indv}_2_P.trim.fq.gz ./fastq_filtered/${indv}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
		rm ./tmp.${indv}_*.fq.gz
	fi
done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

------------------------------------------------------------------------------------------
Ran script on sample subset lists:
------------------------------------------------------------------------------------------

$nohup sh runTrimmomatic.sh sample.1.list > ./log/runTrimmomatic.1.log &
$nohup sh runTrimmomatic.sh sample.2.list > ./log/runTrimmomatic.2.log &
$nohup sh runTrimmomatic.sh sample.3.list > ./log/runTrimmomatic.3.log &
$nohup sh runTrimmomatic.sh sample.4.list > ./log/runTrimmomatic.4.log &
$nohup sh runTrimmomatic.sh sample.5.list > ./log/runTrimmomatic.5.log &
$nohup sh runTrimmomatic.sh sample.6.list > ./log/runTrimmomatic.6.log &
$nohup sh runTrimmomatic.sh sample.7.list > ./log/runTrimmomatic.7.log &
$nohup sh runTrimmomatic.sh sample.8.list > ./log/runTrimmomatic.8.log &
$nohup sh runTrimmomatic.sh sample.9.list > ./log/runTrimmomatic.9.log &
$nohup sh runTrimmomatic.sh sample.10.list > ./log/runTrimmomatic.10.log &
$nohup sh runTrimmomatic.sh sample.11.list > ./log/runTrimmomatic.11.log &
$nohup sh runTrimmomatic.sh sample.12.list > ./log/runTrimmomatic.12.log &

------------------------------------------------------------------------------------------
Remove unpaired reads:
------------------------------------------------------------------------------------------

rm ./fastq_filtered/*U.trim.fq.gz


##########################################################################################
#																						 #
#		  		  i.ii Mapping filtered read data to reference using BWA				 #
#																						 #
##########################################################################################

We'll use bwa mem to map the filtered read data to the barn swallow reference genome.

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd /data2/hirundo_fertilization_network
$mkdir bam

------------------------------------------------------------------------------------------
Retrieve reference genome:
------------------------------------------------------------------------------------------

$cp /media/drewschield/VernalBucket/hirundo/Hirundo_rustica_bHirRus1.final.* .

------------------------------------------------------------------------------------------
Wrote script to run bwa mem on samples in lists:
------------------------------------------------------------------------------------------

To save disk space on `/data2/`, the script will move the filtered fastq data for each
sample to `/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/fastq_filtered/` after
completing the analysis.

We'll run the mapping script on the sample subset lists above, with the number of analysis
threads set to 8, which means we can only run 6 parallel analyses on Terminator at a time.
A pilot analysis with 8 threads took 16 minutes and 27 seconds, meaning that the whole
analysis will take ~2 days to complete (in two batches, lists 1-6 then 7-12).

runBwaMem.sh <sample-list>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
filpath="/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/fastq_filtered/"
samples=$1
for indv in `cat $samples`; do
	name=$indv
	echo "Mapping filtered $name data to reference"
	bwa mem -t 8 -R "@RG\tID:$name\tLB:Hirundo\tPL:illumina\tPU:NovaSeq6000\tSM:$name" Hirundo_rustica_bHirRus1.final.fasta ./fastq_filtered/${name}_1_P.trim.fq.gz ./fastq_filtered/${name}_2_P.trim.fq.gz | samtools sort -@ 8 -O bam -T $name.temp -o ./bam/$name.bam -
	echo "indexing..."
	samtools index -@ 8 ./bam/$name.bam
	mv ./fastq_filtered/${name}_*_P.trim.fq.gz $filpath
done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

------------------------------------------------------------------------------------------
Ran script on sample subset lists:
------------------------------------------------------------------------------------------

$nohup sh runBwaMem.sh sample.1.list > ./log/runBwaMem.1.log &
$nohup sh runBwaMem.sh sample.2.list > ./log/runBwaMem.2.log &
$nohup sh runBwaMem.sh sample.3.list > ./log/runBwaMem.3.log &
$nohup sh runBwaMem.sh sample.4.list > ./log/runBwaMem.4.log &
$nohup sh runBwaMem.sh sample.5.list > ./log/runBwaMem.5.log &
$nohup sh runBwaMem.sh sample.6.list > ./log/runBwaMem.6.log &
$nohup sh runBwaMem.sh sample.7.list > ./log/runBwaMem.7.log &
$nohup sh runBwaMem.sh sample.8.list > ./log/runBwaMem.8.log &
$nohup sh runBwaMem.sh sample.9.list > ./log/runBwaMem.9.log &
$nohup sh runBwaMem.sh sample.10.list > ./log/runBwaMem.10.log &
$nohup sh runBwaMem.sh sample.11.list > ./log/runBwaMem.11.log &
$nohup sh runBwaMem.sh sample.12.list > ./log/runBwaMem.12.log &

Reminder: the filtered fastq data have been moved to
/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/fastq_filtered/


##########################################################################################
#																						 #
#		  		  		i.iii Individual variant calling using GATK						 #
#																						 #
##########################################################################################

We'll use GATK HaplotypeCaller to call individual variant sites across the genome.

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd /data2/hirundo_fertilization_network
$mkdir gvcf

Retrieved GATK distributions for analysis:

$cp -r /data3/hirundo/gatk-* .

------------------------------------------------------------------------------------------
Ran pilot analysis to see how many threads GATK is going to use per instance:
------------------------------------------------------------------------------------------

./gatk-4.0.8.1/gatk HaplotypeCaller -R Hirundo_rustica_bHirRus1.final.fasta --ERC GVCF -I ./bam/CO_00632.bam -O ./gvcf/CO_00632.raw.snps.indels.g.vcf

------------------------------------------------------------------------------------------
Wrote script to run GATK HaplotypeCaller on sample subset lists:
------------------------------------------------------------------------------------------

To save disk space on `/data2/`, the script will move the bam file for each sample to
`/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/bam/` after completing the analysis.

runGatkHaplotypeCaller.sh <sample-list>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
filpath="/media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/bam/"
samples=$1
for indv in `cat $samples`; do
	echo "Calling variants for $indv using GATK HaplotypeCaller"
	./gatk-4.0.8.1/gatk HaplotypeCaller -R Hirundo_rustica_bHirRus1.final.fasta --ERC GVCF -I ./bam/$indv.bam -O ./gvcf/$indv.raw.snps.indels.g.vcf
	bgzip ./gvcf/$indv.raw.snps.indels.g.vcf
	mv ./bam/$indv.bam $filpath 
done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

------------------------------------------------------------------------------------------
Ran script on sample subset lists:
------------------------------------------------------------------------------------------

$nohup sh runGatkHaplotypeCaller.sh sample.1.list > ./log/runGatkHaplotypeCaller.1.log &
$nohup sh runGatkHaplotypeCaller.sh sample.2.list > ./log/runGatkHaplotypeCaller.2.log &
$nohup sh runGatkHaplotypeCaller.sh sample.3.list > ./log/runGatkHaplotypeCaller.3.log &
$nohup sh runGatkHaplotypeCaller.sh sample.4.list > ./log/runGatkHaplotypeCaller.4.log &
$nohup sh runGatkHaplotypeCaller.sh sample.5.list > ./log/runGatkHaplotypeCaller.5.log &
$nohup sh runGatkHaplotypeCaller.sh sample.6.list > ./log/runGatkHaplotypeCaller.6.log &
$nohup sh runGatkHaplotypeCaller.sh sample.7.list > ./log/runGatkHaplotypeCaller.7.log &
$nohup sh runGatkHaplotypeCaller.sh sample.8.list > ./log/runGatkHaplotypeCaller.8.log &
$nohup sh runGatkHaplotypeCaller.sh sample.9.list > ./log/runGatkHaplotypeCaller.9.log &
$nohup sh runGatkHaplotypeCaller.sh sample.10.list > ./log/runGatkHaplotypeCaller.10.log &
$nohup sh runGatkHaplotypeCaller.sh sample.11.list > ./log/runGatkHaplotypeCaller.11.log &
$nohup sh runGatkHaplotypeCaller.sh sample.12.list > ./log/runGatkHaplotypeCaller.12.log &

[Update]: These processes were taking too long on Terminator alone, so I ran another 12
parallel processes on Nostromo, taking the last 37 samples from each list into a new list
and moving the bam files for these on Terminator to a temporary directory to avoid 
duplicating effort.

On Nostromo:

$nohup sh runGatkHaplotypeCaller.sh sample.1-tmp.list > ./log/runGatkHaplotypeCaller.1-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.2-tmp.list > ./log/runGatkHaplotypeCaller.2-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.3-tmp.list > ./log/runGatkHaplotypeCaller.3-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.4-tmp.list > ./log/runGatkHaplotypeCaller.4-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.5-tmp.list > ./log/runGatkHaplotypeCaller.5-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.6-tmp.list > ./log/runGatkHaplotypeCaller.6-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.7-tmp.list > ./log/runGatkHaplotypeCaller.7-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.8-tmp.list > ./log/runGatkHaplotypeCaller.8-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.9-tmp.list > ./log/runGatkHaplotypeCaller.9-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.10-tmp.list > ./log/runGatkHaplotypeCaller.10-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.11-tmp.list > ./log/runGatkHaplotypeCaller.11-tmp.log &
$nohup sh runGatkHaplotypeCaller.sh sample.12-tmp.list > ./log/runGatkHaplotypeCaller.12-tmp.log &

------------------------------------------------------------------------------------------
Rejoice:
------------------------------------------------------------------------------------------

It is a joyous day, for the individual variant calls are now complete and ready for cohort
variant calling among individuals.

To speed up this step, I'll use scaffold lists in `.intervals` files to run a set of parallel
analyses.


##########################################################################################
#																						 #
#		  		  		  i.iv Cohort variant calling using GATK						 #
#																						 #
##########################################################################################

We'll use GATK GenotypeGVCFs to call variant sites among individuals across the genome.

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

cd /data2/hirundo_fertilization_network/
mkdir vcf

------------------------------------------------------------------------------------------
Wrote sample list with path to gVCF files for GATK:
------------------------------------------------------------------------------------------

`sample.all.gvcf.list`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
./gvcf/CO_00632.raw.snps.indels.g.vcf.gz
./gvcf/CO_06912.raw.snps.indels.g.vcf.gz
./gvcf/CO_21085.raw.snps.indels.g.vcf.gz
./gvcf/CO_33601.raw.snps.indels.g.vcf.gz
./gvcf/CO_33602.raw.snps.indels.g.vcf.gz
./gvcf/CO_33603.raw.snps.indels.g.vcf.gz
./gvcf/CO_33604.raw.snps.indels.g.vcf.gz
./gvcf/CO_33605.raw.snps.indels.g.vcf.gz
./gvcf/CO_33606.raw.snps.indels.g.vcf.gz
./gvcf/CO_33607.raw.snps.indels.g.vcf.gz
...
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

------------------------------------------------------------------------------------------
Indexed input GVCFs with tabix:
------------------------------------------------------------------------------------------

runTabix.sh <list>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
list=$1
for gvcf in `cat $list`; do
	echo indexing $gvcf
	tabix -p vcf $gvcf
done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

$nohup sh runTabix.sh sample.all.gvcf.list > ./log/runTabix.log &

------------------------------------------------------------------------------------------
Wrote .intervals files with lists of scaffolds to run analyses on:
------------------------------------------------------------------------------------------

These are split out based on the largest scaffold, which represents nearly 14% of the
genome sequence. I've arranged sets of scaffolds to not exceed this proportion.

Running everything in sequence would take several weeks, as estimated by GATK, so this
should allow for faster processing.

`scaffold.list1.intervals`
`scaffold.list2.intervals`
`scaffold.list3.intervals`
`scaffold.list4.intervals`
`scaffold.list5.intervals`
`scaffold.list6.intervals`
`scaffold.list7.intervals`
`scaffold.list8.intervals`
`scaffold.list9.intervals`
`scaffold.list10.intervals`
`scaffold.list11.intervals`
`scaffold.list12.intervals`
`scaffold.list13.intervals`

Scaffold interval lists are in `./log/`.

------------------------------------------------------------------------------------------
Ran GenotypeGVCFs to call 'raw' variants across scaffold intervals:
------------------------------------------------------------------------------------------

Here, we're specifying `-stand_call_conf 15` to permit lower phred quality scores.

$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list1.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list1.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list1.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list2.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list2.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list2.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list3.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list3.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list3.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list4.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list4.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list4.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list5.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list5.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list5.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list6.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list6.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list6.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list7.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list7.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list7.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list8.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list8.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list8.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list9.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list9.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list9.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list10.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list10.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list10.vcf.recalc.log &
$nohup java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ./Hirundo_rustica_bHirRus1.final.fasta -stand_call_conf 15 -L ./log/scaffold.list11.intervals -V ./sample.all.gvcf.list -o ./vcf/hirundo_rustica.raw.scaffold.list11.vcf.gz > ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.list11.vcf.recalc.log &

We're not running intervals in lists 12 and 13 - these are the sex chromosomes.

To check in on running progress:

$for i in ./log/GenotypeGVCFs.hirundo_rustica.raw.scaffold.*.recalc.log; do echo $i; grep 'INFO' $i | tail -n 1; done

------------------------------------------------------------------------------------------
Wrote list with paths to input VCF files:
------------------------------------------------------------------------------------------

`vcf.recalc.list`

This omits the list 12 and 13 scaffolds (sex linked).

------------------------------------------------------------------------------------------
Ran MergeVCFs to output combined VCF file with all scaffolds:
------------------------------------------------------------------------------------------

$nohup java -jar picard.jar MergeVcfs -I vcf.recalc.list -O ./vcf/hirundo_rustica.recalc.vcf.gz > ./log/mergeVcfs.recalc.picard.log &

------------------------------------------------------------------------------------------
Removed 'raw' scaffold list VCFs taking up disk space:
------------------------------------------------------------------------------------------

$rm ./vcf/hirundo_rustica.raw.scaffold.*.vcf.*


##########################################################################################
#																						 #
#																						 #	
#																						 #
#			  1. Workflow for running lcMLkin, an algorithm for parentage				 #
#		 		 		  assignment using genome-wide SNP data							 #
#																						 #
#																						 #
#																						 #
##########################################################################################

This sections contains details on the installation of lcMLkin, testing the software
on test data, and preparing inputs and running the program on our data based on the best-
practices recommendations for the method.

This section also includes critical details on recalculation of raw genotype likelihoods
for sets of target SNPs, based on the unfiltered VCF produced in the previous section.


##########################################################################################
#				  		  	 1.1 Installation of lcMLkin								 #
##########################################################################################

Arielle Fogel suggested using this new method, lcMLkin, to Heather, which is available on
GitHub:

https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation

And best practices are described here:

https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation/wiki/Best-Practices

Note, the clone from the main repository did not work on Ubuntu. Another user fixed this
issue in a forked repo:

https://github.com/didillysquat/maximum-likelihood-relatedness-estimation

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd ~/tmp/

------------------------------------------------------------------------------------------
Install dependency clang on the system:
------------------------------------------------------------------------------------------

$sudo apt install clang

------------------------------------------------------------------------------------------
Install lcMLkin:
------------------------------------------------------------------------------------------

$git clone https://github.com/didillysquat/maximum-likelihood-relatedness-estimation.git
$cd maximum-likelihood-relatedness-estimation
$make

This compiles the executable `lcmlkin`.

Copied to analysis directory:

$cp lcmlkin /data2/hirundo_fertilization_network/

------------------------------------------------------------------------------------------
Install raw genotype likelihood variant caller for lcMLkin:
------------------------------------------------------------------------------------------

The authors of lcMLkin have written a companion python script `SNPbam2vcf.py` to calculate
genotype likelihoods for SNP sites given an input list of BAM files and SNP positions,
then outputting a VCF file.

The script requires pysam and numpy, both of which are already in our default Conda env.

Retrieved script from https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation/blob/master/src_python/SNPbam2vcf/SNPbam2vcf.py

`SNPbam2vcf.py`


##########################################################################################
#																						 #
#			 1.2 Produce target SNP sets for lcMLkin variant calling,					 #
#			  perform raw GL recalculation, and run lcMLkin analysis					 #
#																						 #
##########################################################################################

Pilot analyses (see Appendix) indicate two interesting features of analysis of the 'raw'
GLs from the in-house variant caller for lcMLkin:
	1. Relatedness estimates for known relatives increase
	2. Relatedness estimates for known non-relatives decrease
	
This is good, and argues for using the most inclusive set of 'target' SNPs in our dataset.

To do this, we've relaxed the phred quality filter in the above GATK genotypeGVCFs command
and will only filter for SNPs that are not in repeats, followed by appropriate MAF and
LD filters (i.e., we won't impose any further depth/quality filters).

------------------------------------------------------------------------------------------
Set up environment:
------------------------------------------------------------------------------------------

$cd /data2/hirundo_fertilization_network/

------------------------------------------------------------------------------------------
Removed repeats, indels, and non-biallelic SNPs using bcftools:
------------------------------------------------------------------------------------------

$bcftools view --threads 24 -T ^/data3/hirundo/genome_annotation/repeatmasker/Hirundo_rustica_bHirRus1.final.fasta.repeat.bed -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.mask.vcf.gz ./vcf/hirundo_rustica.recalc.vcf.gz

------------------------------------------------------------------------------------------
Formatted sample lists:
------------------------------------------------------------------------------------------

`list.CO-2021.indv`
`list.CO-2022.indv`
`list.OH-2022.indv`

------------------------------------------------------------------------------------------
Parsed samples by population/year, MAF >= 0.05:
------------------------------------------------------------------------------------------

$bcftools view --threads 16 -S list.CO-2021.indv -i 'MAF >= 0.05' -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.vcf.gz ./vcf/hirundo_rustica.recalc.mask.vcf.gz
$bcftools view --threads 16 -S list.CO-2022.indv -i 'MAF >= 0.05' -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.vcf.gz ./vcf/hirundo_rustica.recalc.mask.vcf.gz
$bcftools view --threads 16 -S list.OH-2022.indv -i 'MAF >= 0.05' -m2 -M2 -U -v snps -O z -o ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.vcf.gz ./vcf/hirundo_rustica.recalc.mask.vcf.gz

-----------------------------------------------------------------------------------------
Ran vcftools to thin SNPs to prune for LD:
------------------------------------------------------------------------------------------

$vcftools --gzvcf ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.vcf.gz --thin 10000 --recode --stdout | bgzip -c > ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.thin10k.vcf.gz
$vcftools --gzvcf ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.vcf.gz --thin 10000 --recode --stdout | bgzip -c > ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.thin10k.vcf.gz
$vcftools --gzvcf ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.vcf.gz --thin 10000 --recode --stdout | bgzip -c > ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.thin10k.vcf.gz

------------------------------------------------------------------------------------------
Formatted bam lists:
------------------------------------------------------------------------------------------

$touch list.CO-2021.indv-bam; for i in `cat list.CO-2021.indv`; do echo /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/bam/${i}.bam >> list.CO-2021.indv-bam; done
$touch list.CO-2022.indv-bam; for i in `cat list.CO-2022.indv`; do echo /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/bam/${i}.bam >> list.CO-2022.indv-bam; done
$touch list.OH-2022.indv-bam; for i in `cat list.OH-2022.indv`; do echo /media/drewschield/VernalBucket/hirundo_lcwgs_CO-OH/bam/${i}.bam >> list.OH-2022.indv-bam; done

`list.CO-2021.indv-bam`
`list.CO-2022.indv-bam`
`list.OH-2022.indv-bam`

------------------------------------------------------------------------------------------
Formatted SNP positions file:
------------------------------------------------------------------------------------------

$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.thin10k.vcf.gz > snp.recalc.CO-2021.maf05.thin10kb.list
$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.thin10k.vcf.gz > snp.recalc.CO-2022.maf05.thin10kb.list
$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.thin10k.vcf.gz > snp.recalc.OH-2022.maf05.thin10kb.list

------------------------------------------------------------------------------------------
Run raw GL variant caller script on input data:
------------------------------------------------------------------------------------------

Usage: ./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file>

$./SNPbam2vcf.py list.CO-2021.indv-bam ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.thin10k.rawGL.vcf snp.recalc.CO-2021.maf05.thin10kb.list
$./SNPbam2vcf.py list.CO-2022.indv-bam ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.thin10k.rawGL.vcf snp.recalc.CO-2022.maf05.thin10kb.list
$./SNPbam2vcf.py list.OH-2022.indv-bam ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.thin10k.rawGL.vcf snp.recalc.OH-2022.maf05.thin10kb.list

------------------------------------------------------------------------------------------
Run lcMLkin on recalculated GL data:
------------------------------------------------------------------------------------------

$./lcmlkin -i ./vcf/hirundo_rustica.recalc.mask.CO-2021.maf05.thin10k.rawGL.vcf -o ./lcmlkin_analysis/CO-2021.thin10k.GL-recalc.relate -g all -t 12
$./lcmlkin -i ./vcf/hirundo_rustica.recalc.mask.CO-2022.maf05.thin10k.rawGL.vcf -o ./lcmlkin_analysis/CO-2022.thin10k.GL-recalc.relate -g all -t 12
$./lcmlkin -i ./vcf/hirundo_rustica.recalc.mask.OH-2022.maf05.thin10k.rawGL.vcf -o ./lcmlkin_analysis/OH-2022.thin10k.GL-recalc.relate -g all -t 12
