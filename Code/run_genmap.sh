#!/bin/bash
#PBS -l ncpus=1
#PBS -l mem=12GB
#PBS -l jobfs=6GB
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -l wd

cd /scratch/xe2/jb5097/tmp
source /g/data/xe2/John/Nalu/bin/activate
export PATH=/g/data/xe2/John/Software/genmap-build/bin/:$PATH

ref=/g/data/xe2/datasets/NCBI_genomes/GCA_004368105.4_ASM436810v4_genomic.fna
tmp_dir=genmap_index
out_dir=/scratch/xe2/jb5097/tmp #<where to save the index and outputs>

# Set params
K=30
E=2

# print settings:
# Print settings
echo "Reference: $ref"
echo "Temporary Directory: $tmp_dir"
echo "Output Directory: $out_dir"
echo "K: $K"
echo "E: $E"

# Create index -- Only need to run this this once
if [ ! -d "$tmp_dir" ]; then
  echo "Index does not exist. Creating index..."
  genmap index -F $ref -I $tmp_dir
else
  echo "Index already exists."
fi

# Compute mappability
genmap map -K $K -E $E -I $tmp_dir -O $out_dir -t -w -bg


