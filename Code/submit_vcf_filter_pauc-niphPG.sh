for i in 01 02 03 04 05 06 07 08 09 10 11
do
echo Chr${i}
#qsub -v keeplist=pauc_niph_PG,CHROM=Chr${i} Code/vcf_filter_pauc-niphPG.sh
qsub -v keeplist=niph33_pauc233,CHROM=Chr${i} Code/vcf_filter_pauc-niphPG.sh
done

# in this implementation, have to change the popfile manually in the .sh script. Should update this... 
