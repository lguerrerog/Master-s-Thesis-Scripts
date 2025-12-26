#!/usr/bin/env bash

#SBATCH -c 4                  # number of CPUs (here 4)
#SBATCH --mem=32000          # memory pool for all cores (here 32GB)
#SBATCH --partition=medium      # Partition to submit to
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH -J "ChrMT"             # Job name

#Sex determination and x chromosome contamination 
path=$(pwd) #/nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Capture or /nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Shotgun

mkdir -p $path/SexDet
mkdir -p $path/XContamination
echo "SexDet and Xcontamination folders done"

#Capture
#SexDet
#cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/merged_bam/UMIN018.A.1.1_merged_rmdup.sorted.trim2.bam /nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Capture/SexDet/
#cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/merged_bam/UMIN018.A.2.1_merged_rmdup.sorted.trim2.bam /nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Capture/SexDet/
#X contamination
#cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/merged_bam/UMIN018.A.1.1_merged_rmdup.sorted.trim2.bam /nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Capture/XContamination/
#cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Capture/config_2/merged_bam/UMIN018.A.2.1_merged_rmdup.sorted.trim2.bam /nfs/geo/projects/south_america/Banks_Island/SexDet_XContamination_Lucia/Capture/XContamination/

#Shotgun
#SexDet
cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Shotgun_high/config/merged_bam/*.trim2.bam $path/SexDet/
#X contamination
cp /nfs/geo/projects/south_america/Banks_Island/EAGER_Lucia/Shotgun_high/config/merged_bam/*.trim2.bam $path/XContamination/

#copy BAM files
cd $path/SexDet
cp $path2/*trim2.bam .
cd ..
cd $path/XContamination
cp $path2/*trim2.bam $path/XContamination

#Indexing
for Index in $path/XContamination/*.trim2.bam; do
	samtools index $Index
done
echo "Indexing done"

#Sex Determination
cd $path/SexDet
for i in *trim2.bam; do
    BaseName=$(basename "$i" .trim2.bam)
    SampleID=$(echo "$BaseName" | cut -d '_' -f1)

    CMD="samtools depth -q30 -Q37 -aa -b /opt/reference/human/hg19/SNPCapBEDs/1240K.pos.list_hg19.0based.bed $i | \
awk -f /nfs/geo/projects/south_america/Banks_Island/scripts/sexDetermination.hg19.awk > $SampleID.sexDetermination.txt"

    sbatch --mem=32000 -c 2 -p short -o "$i.sexDetermination.log" --wrap="$CMD"
done
echo "SexDet done"

#X contamination
cd $path/XContamination
for Bam2 in  $(ls *trim2.bam); do
	BaseName=$(basename "$Bam2" .trim2.bam)
  SampleID=$(echo "$BaseName" | cut -d "_" -f1)
	sbatch --mem 32000 -c 6 -p short -o "${SampleID}.angsdCounts.log" --wrap="angsd -i $Bam2 -r chrX:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${SampleID}.angsdCounts"
done
echo "X contamination done"

#Final count for X contamination
mkdir -p FinalCont
for X in $(ls *.trim2.bam); do
	BaseName=$(basename "$X" .trim2.bam)
  SampleID=$(echo "$BaseName" | cut -d "_" -f1)
	XCont=$(/opt/sources/angsd-0.940-67b6b3b/misc/contamination contamination  -a $SampleID.angsdCounts.icnts.gz -h /opt/sources/angsd-0.940-67b6b3b/RES/HapMapChrX.gz 2> $path/XContamination/FinalCont/$SampleID.xcontamination.out)

	sbatch --mem=32000 -o $path/XContamination/FinalCont/$XCont.xCont.log --wrap="$XCont"
done
echo "Final count done"
