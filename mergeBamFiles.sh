#!/bin/bash

##Merge BAM files
# For non-bleached
samtools merge -@ 16 -o /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/UMIN018.A.1.1_merged.bam /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/UMIN018.A.1.1.TFA.1__XXX001.A2262.SG1.1/5-DeDup/UMIN018.A.1.1.TFA.1__XXX001.A2262.SG1.1_S0_L004_R1_001.fastq.truncated.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned_rmdup.sorted.bam /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/UMIN018.A.1.1.TFA.2__XXX001.A2262.SG1.2/5-DeDup/UMIN018.A.1.1.TFA.2__XXX001.A2262.SG1.2_S0_L007_R1_001.fastq.truncated.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned_rmdup.sorted.bam -f

#For bleached
samtools merge -@ 16 -o /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/configUMIN018.A.2.1_merged.bam /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config/UMIN018.A.2.1.1.TFA.1/5-DeDup/UMIN018.A.2.1.1.TFA.1_L001_R1.fastq.truncated.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned_rmdup.sorted.bam /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config/UMIN018.A.2.1.TFA.2/5-DeDup/UMIN018.A.2.1.TFA.2_L002_R1.fastq.truncated.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned_rmdup.sorted.bam /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config/UMIN018.A.2.1.TFA.3/5-DeDup/UMIN018.A.2.1.TFA.3_L002_R1.fastq.truncated.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned_rmdup.sorted.bam  -f

#Sort
for i in *_merged.bam; do samtools sort -o "${i/.bam/.sorted.bam}" "$i"; done

#De-Duplicate
for i in *_merged.sorted.bam; do eager1 dedup -m -i $i -o .; done

#Sort
for i in *merged.sorted_rmdup.bam; do samtools sort -o ${i/merged.sorted_rmdup.bam/merged_rmdup}.sorted.bam $i; done

#Sanity check
#Count reads before de-duplication 
samtools view -c *.sorted.bam

#Count reads after de-duplication
samtools view -c *.dedup.bam

#In one line:
echo "Before:" && samtools view -c *.sorted.bam
echo "After:" && samtools view -c *.dedup.bam