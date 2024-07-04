#!/bin/bash
#PBS -l ncpus=48
#PBS -l mem=164GB
#PBS -l jobfs=12GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

# Activate python virtual env for genetics stuff:
module load python3/3.12.1
source /g/data/xe2/John/Nalu/bin/activate

export PATH=/g/data/xe2/John/Software/samtools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/bcftools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/htslib-1.20/bin:$PATH

cd /scratch/xe2/jb5097/tmp/

## Set the input VCF and sample list
vcf_in=/g/data/xe2/projects/paneuc2/final-variants/2024-02-29/mpileup~bwa~Emelliodora_sf2~melsider~filtered-default.vcf.gz

keeplist=adnataria_in_vcf_all_no_gorillas.csv
keeplist_desc=nogorilla

filtercmds=(
    "bcftools +fill-tags - -Ou -- -d  -t all,F_MISSING |"
    "bcftools view -Ob1 -'"
)

blsl vcfparallel \
    -o adna20240229_${keeplist_desc}_snpsbi_maf01.vcf.gz \
    -O z \
    -c "${filtercmds[*]}" \
    -f "-S $keeplist -m2 -M2 -v snps -i 'QUAL>30 && INFO/MAF >=0.01" \
    -t 48 \
    $vcf_in

