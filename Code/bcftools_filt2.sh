#!/bin/bash
#PBS -l ncpus=24
#PBS -l mem=164GB
#PBS -l jobfs=24GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=14:00:00
#PBS -l storage=scratch/xe2+gdata/xe2+gdata/if89
#PBS -l wd

## run it like:
# qsub Code/vcf_filter_20240229.sh

module use /g/data/if89/apps/modulefiles
module load vcftools
module load bcftools htslib/1.16

cd /scratch/xe2/jb5097/tmp/

in=/g/data/xe2/projects/paneuc2/final-variants/2024-02-29/mpileup~bwa~Emelliodora_sf2~melsider~filtered-default.vcf.gz
keeplist=Adnataria_after_lowqual_individs_dropped_miss_gt_20pct.csv

# set the output file stub:
file=20240229_Chr10
out=${file}_filt2b

bcftools view -r Chr10 $in --samples-file $keeplist --threads 24 | \
bcftools view -m2 -M2 -v snps --threads 24 | \
bcftools filter -e 'F_MISSING > 0.2 || MAF<0.01' --threads 24 | \
bcftools +prune --window 200bp --nsites-per-win 1 --nsites-per-win-mode "rand" -Oz -o ${out}.vcf.gz

vcftools --gzvcf ${out}.vcf.gz --depth --out $out
vcftools --gzvcf ${out}.vcf.gz --site-mean-depth --out $out
vcftools --gzvcf ${out}.vcf.gz --het --out $out
vcftools --gzvcf ${out}.vcf.gz --relatedness2 --out $out
vcftools --gzvcf ${out}.vcf.gz --missing-indv --out $out

echo "finished $out"
