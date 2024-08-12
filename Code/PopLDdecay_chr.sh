#!/bin/bash
#PBS -l ncpus=8
#PBS -l mem=7GB
#PBS -l jobfs=16GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

# If re-doing, 8GB mem is plenty. Last go timed out at 48 hrs, only used 6GB MEM. Got to Chr 9 (incomplete)

# Caculate LD decay for meliodora and sideroxylon using high qual phased SNPs (short read seq data) from Scott's GigaScience paper
# This script uses the program called PopLDDecay https://academic.oup.com/bioinformatics/article/35/10/1786/5132693

export PATH=/g/data/xe2/John/Software/bcftools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/PopLDdecay/bin:$PATH

cd /scratch/xe2/jb5097/tmp/Scotts_snps

stub=E_meliodora # expects $stub.bcf

# convert the bcf to vcf bcf file must be indexed
bcftools index $stub.bcf

for i in 09 10 11
do
	echo 'Chr'$i
	bcftools view -r Chr${i} -O z -o ${stub}_Chr${i}.vcf.gz $stub.bcf
	PopLDdecay -InVCF ${stub}_Chr${i}.vcf.gz -OutStat ${stub}_Chr${i}_LD
done

ls *_LD.stat.gz > LD_chroms.list

Plot_OnePop.pl -inList LD_chroms.list \
        -bin1 10 \
        -bin2 100 \
        -maxX 500 \
        -output ${stub}_PopLDdecay_allchroms

# This generates an R script and a dataset for plotting.
# I have modified the R script for log x-axis plot and to show the dist at which LD decays to half the max.

# do -inList if doing multiple chrs (list the .stat.gz files

