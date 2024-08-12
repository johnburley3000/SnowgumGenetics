#!/bin/bash
#PBS -l ncpus=8
#PBS -l mem=64GB
#PBS -l jobfs=4GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=6:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

# Caculate LD decay for meliodora and sideroxylon using high qual phased SNPs (short read seq data) from Scott's GigaScience paper
# note, because its doing pairwise SNP stuff, the memory required can be huge if not careful...
# It seems you only need small %  of SNPs (randomly sampled) to get the LD decay curve right.. 

cd /scratch/xe2/jb5097/tmp/Scotts_snps

stub=E_meliodora

# Activate python virtual env for genetics stuff:
module load python3/3.12.1
source /g/data/xe2/John/Nalu/bin/activate

module load plink2/v2.00a3.7

export PATH=/g/data/xe2/John/Software/samtools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/bcftools-1.20/bin:$PATH
export PATH=/g/data/xe2/John/Software/htslib-1.20/bin:$PATH

plink=/g/data/xe2/John/Software/plink

plotLD=/home/106/jb5097/Projects/SnowgumGenetics/Code/plotLD.py

$plink --bcf $stub.bcf \
	--double-id \
	--allow-extra-chr \
	--set-missing-var-ids @:# \
	--chr 1 \
	--thin 0.1 \
	--r2 gz yes-really \
	--ld-window 9999999 \
	--ld-window-kb 300 \
	--ld-window-r2 0 \
	--make-bed \
	--out ${stub}

# key flag decisions:
# --ld-window <ct+1> 99999 ensures that nearly all variant pairs are considered, assuming other window restrictions are not more limiting​
# --ld-window-kb <x> max dist between SNP pairs for LD calculation
# --ld-window-r2 <x> min reported R2

python $plotLD --input ${stub}.ld.gz --bin_size_bp 10 --output_prefix E_mel_Chr01_thin10_win300kb_bin10

