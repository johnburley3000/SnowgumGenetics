#!/bin/bash
#PBS -l ncpus=48
#PBS -l mem=164GB
#PBS -l jobfs=12GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

## run it like:
# qsub Code/vcf_filter_20240229.sh

module use /g/data/if89/apps/modulefiles
module load vcftools
module load htslib/1.16
module load bcftools

cd /scratch/xe2/jb5097/tmp/

in=/g/data/xe2/projects/paneuc2/final-variants/2024-02-29/mpileup~bwa~Emelliodora_sf2~melsider~filtered-default.vcf.gz
keeplist=mel_sider_list_filt3b.txt

# set the output file stub:
file=20240229_Chr10
out=${file}_filt3b

bcftools view -r Chr10 $in --samples-file $keeplist -r Chr10:1000000-5000000 --threads 48 | \
bcftools view -m2 -M2 -v snps --threads 48 | \
bcftools filter -e 'QUAL < 20 || MAF < 0.01' --threads 24 -Oz -o ${out}.vcf.gz

vcftools --gzvcf ${out}.vcf.gz --depth --out $out
vcftools --gzvcf ${out}.vcf.gz --site-mean-depth --out $out
vcftools --gzvcf ${out}.vcf.gz --het --out $out
vcftools --gzvcf ${out}.vcf.gz --relatedness2 --out $out
vcftools --gzvcf ${out}.vcf.gz --missing-indv --out $out

echo "finished $out"
