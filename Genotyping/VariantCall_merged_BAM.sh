#!/usr/bin/env bash

#SBATCH -c 4                  # number of CPUs (here 4)
#SBATCH --mem=32000          # memory pool for all cores (here 32GB)
#SBATCH --partition=medium      # Partition to submit to
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH -J "Variant"             # Job name


#Script for variant calling
path=$(pwd)


#mkdir $path/Variant_Call

#cd $path/Variant_Call

for l in $path/*.dedup.bam; do
/opt/tools.old/bamutil/1.0.15/bin/bam trimBam $l ${l/.bam/}.trim2.bam 2

done

#cp $path/*trim2.bam $path/Variant_Call 




ls $path/*trim2.bam | xargs -n 1 basename |cut -d '_' -f1 > $path/list_Name.txt
BAM=$(ls $path/*trim2.bam)

#Genotypes
samtools mpileup -q30 -Q30 -B -f /opt/reference/human/hg19/hg19_complete.fasta -l /opt/reference/human/hg19/SNPCapBEDs/1240K.pos.list_hg19.0based.bed $BAM | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chr//' | /opt/tools/sequencetools/1.5.2/bin/pileupCaller --randomHaploid --sampleNameFile list_Name.txt  --samplePopName AP  -f /opt/genotypes/eigensoft/reich_lab/v44.3_1240K_public.snp -e Calls.Comp

rm slurm*
