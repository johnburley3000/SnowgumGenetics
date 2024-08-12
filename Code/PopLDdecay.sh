#!/bin/bash
#PBS -l ncpus=8
#PBS -l mem=64GB
#PBS -l jobfs=16GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=6:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

# Caculate LD decay for meliodora and sideroxylon using high qual phased SNPs (short read seq data) from Scott's GigaScience paper
# This script uses the program called PopLDDecay https://academic.oup.com/bioinformatics/article/35/10/1786/5132693

export PATH=/g/data/xe2/John/Software/bcftools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/PopLDdecay/bin:$PATH

cd /scratch/xe2/jb5097/tmp/Scotts_snps

# convert BCF to VCF

stub=E_meliodora

# convert the bcf to vcf bcf file must be indexed

#bcftools index $stub.bcf
#bcftools view -r Chr01 -O z -o ${stub}_Chr01.vcf.gz $stub.bcf

PopLDdecay -InVCF ${stub}_Chr01.vcf.gz -OutStat LDdecay_${stub}_Chr01

Plot_OnePop.pl -inFile LDdecay_${stub}_Chr01.stat.gz -output ${stub}_PopLDdecay
