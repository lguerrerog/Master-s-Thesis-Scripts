#!/usr/bin/env bash

#SBATCH -c 6                  # number of CPUs (here 4)
#SBATCH --mem=32000          # memory pool for all cores (here 32GB)
#SBATCH --partition=medium      # Partition to submit to
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH -J "Schmutzi"             # Job name

path=$(pwd) #This should be /nfs/geo/projects/south_america/Banks_Island
#path2=/nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/MT_DNA_sg/config

#mkdir -p $path/schmutzi_Lucia
#cd $path2

# Iterates into each folder in config and copies the bam file in the new schmutzi folder
#for i in $path2/*; do
#pp=$(pwd)
#cd $i/5-DeDup
#cd 5-DeDup
#cp *_rmdup.sorted.bam $path/schmutzi_Lucia/Shotgun/ #Change this for shotgun and capture
#cd $pp
#done

cd $path/schmutzi_Lucia/Capture2

#Indexes each bam file
for i in *bam; do
samtools index $i
done

mkdir -p $path/schmutzi_Lucia/Capture2/contdeam

#Takes the first part of the name of the bam file until "_" (UMIN018.A.2.1.1.TFA.1)
for i in *bam; do
name=$(echo "${i}" | cut -d "_" -f 1)

#Extracts mitochondrial reads from .bam and calculates MD (map of differences) and NM (number of missmatches) as .bam and stores it in contdeam folder
samtools view -b $i chrMT |samtools calmd -b -  /nfs/geo/projects/south_america/Banks_Island/references/hg19_MT.fasta> $path/schmutzi_Lucia/Capture2/contdeam/$name.MD.bam
done

cd $path/schmutzi_Lucia/Capture2/contdeam

#Indexes again
for i in *MD.bam; do
samtools index $i 
done

# estimate DNA damage patterns and prepare files for contamination estimation
#library double: Specifies the library type as double-stranded DNA
for i in *MD.bam; do
eager1 /projects1/tools/schmutzi/1.0/contDeam.pl --library double --out $path/schmutzi_Lucia/Capture2/contdeam/$i $i 
done

#(add —-length 5 paramether if nonUDG data)

#Estimate mitochondrial contamination and reconstruct a likely endogenous (authentic) mitochondrial genome
#notusepredC = not to use predicted contamination from external predictors
#uselength:Tells Schmutzi to incorporate fragment length distribution in contamination estimation
#197/freqs: Reference allele frequencies (population panel) for contamination estimation. 197: worldwide individuals dataset
for i in *MD.bam; do
eager1 /projects1/tools/schmutzi/1.0/schmutzi.pl --notusepredC -t 15 --uselength --ref /nfs/geo/projects/south_america/Banks_Island/references/hg19_MT.fasta $path/schmutzi_Lucia/Capture2/contdeam/$i /projects1/tools/schmutzi/1.0/alleleFreqMT/197/freqs $i
done



#samtools view -b -s 0.2 input.bam -o output.bam


# Fasta q filtering
mkdir fasta_Schmutzi

for i in *final_endo.log; do
eager1 /projects1/tools/schmutzi/1.0/log2fasta -q 0 $i > $path/schmutzi_Lucia/Capture2/contdeam/fasta_Schmutzi/${i}_q0.fasta
eager1 /projects1/tools/schmutzi/1.0/log2fasta -q 10 $i > $path/schmutzi_Lucia/Capture2/contdeam/fasta_Schmutzi/${i}_q10.fasta
eager1 /projects1/tools/schmutzi/1.0/log2fasta -q 20 $i > $path/schmutzi_Lucia/Capture2/contdeam/fasta_Schmutzi/${i}_q20.fasta
eager1 /projects1/tools/schmutzi/1.0/log2fasta -q 30 $i > $path/schmutzi_Lucia/Capture2/contdeam/fasta_Schmutzi/${i}_q30.fasta

done

cd fasta_Schmutzi

for i in *fasta; do
path2=$(pwd)
name=$(echo "${i}" | cut -d "M" -f 1)
find * | grep "final_endo.log_q0.fasta"$ | awk '{split($1,nm,"/");print "sed -e \"s/>MT/>"nm[1]"/\" "$1;}' | bash > $path2/consensus_mtDNAq0.fasta
find * | grep "final_endo.log_q10.fasta"$ | awk '{split($1,nm,"/");print "sed -e \"s/>MT/>"nm[1]"/\" "$1;}' | bash > $path2/consensus_mtDNAq10.fasta
find * | grep "final_endo.log_q20.fasta"$ | awk '{split($1,nm,"/");print "sed -e \"s/>MT/>"nm[1]"/\" "$1;}' | bash > $path2/consensus_mtDNAq20.fasta
find * | grep "final_endo.log_q30.fasta"$ | awk '{split($1,nm,"/");print "sed -e \"s/>MT/>"nm[1]"/\" "$1;}' | bash > $path2/consensus_mtDNAq30.fasta

done

# compile all contamination results
cd ..
find * | grep "final.cont.est"$ | awk '{print "paste <(echo \""$1"\") "$1;}' | bash >$path2/mtDNA_contEst.txt

#clean up
