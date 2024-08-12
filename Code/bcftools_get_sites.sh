#!/bin/bash
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l jobfs=64GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

# Purpose:
# We have to re-run variant calling using some samples that were accidentally not included. 
# To make this process run faster, get the coordinates of SNPs in the large dataset, to restrict the next round of variant calling. 

## run it like: 
# qsub -v CHROM=01 Code/bcftools_get_sites.sh

# Genetics virtual environment
module load python3/3.12.1
source /g/data/xe2/John/Nalu/bin/activate
# New installs of samtools,bcftools,htslib (https://www.htslib.org/download/)
export PATH=/g/data/xe2/John/Software/samtools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/bcftools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/htslib-1.20/bin:$PATH

cd /scratch/xe2/jb5097/tmp/

in=/g/data/xe2/projects/paneuc2/final-variants/2024-02-29/mpileup~bwa~Emelliodora_sf2~melsider~filtered-default.vcf.gz

# set the output file stub:
file=adnat_20240229
# set the chromosome number (as string e.g. 01)
i=$CHROM
#i=01

keeplist=adnataria_in_vcf_all_no_gorillas.csv
keeplist_desc=nogorilla


echo "starting Chr$i"

bcftools view -r Chr$i -S $keeplist -m2 -M2 -v snps $in --threads 48 | \
bcftools filter -e 'MAF < 0.01 && QUAL < 30' -Ou | \
bcftools query -f '%CHROM\t%POS\n' > ${file}_S${keeplist_desc}_Chr${i}_MAF01_Q30.bed

echo "finished Chr${i} (made bed file).."


