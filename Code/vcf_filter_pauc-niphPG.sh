#!/bin/bash

#PBS -l ncpus=1
#PBS -l mem=20GB
#PBS -l jobfs=6GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=03:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

## Run SNP filtering for pauc/niph samples that have phenotype data
## March 28 2024

## run it like:
# qsub -v keeplist=lowmiss_geo_allspp -v CHROM=Chr01 Code/vcf_filter2.sh

module use /g/data/if89/apps/modulefiles
module load vcftools

echo ${CHROM}
echo ${keeplist}.keeplist

newfile=${keeplist}_${CHROM}_biasnps_DP5-25_mac4_mis70 # stub of the outputs

cd /scratch/xe2/jb5097/tmp/

vcftools --bcf ${CHROM}.bcf \
	--keep ${keeplist}.keeplist \
	--remove-indels \
	--min-alleles 2 --max-alleles 2 \
	--minDP 5 \
	--maxDP 25 \
	--mac 4 \
	--max-missing 0.7 \
	--recode --stdout | gzip -c > ${newfile}.vcf.gz

## Now run FST in windows for a quick look at divergence:

i=niph33
j=pauc233
vcftools --gzvcf ${newfile}.vcf.gz \
	--weir-fst-pop $i.popfile \
	--weir-fst-pop $j.popfile \
	--fst-window-size 5000 \
	--fst-window-step 2500 \
	--out ${newfile}_fst_${i}_${j}_5kb


vcftools --gzvcf ${newfile}.vcf.gz \
        --weir-fst-pop $i.popfile \
        --weir-fst-pop $j.popfile \
	--out ${newfile}_fst_${i}_${j}


