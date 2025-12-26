#!/bin/bash

path=$(pwd)

#mkdir Admixture

#cp $path/Variant_Call/Calls.Comp.geno $path/Admixture
#cp $path/Variant_Call/Calls.Comp.snp $path/Admixture
#cp $path/Variant_Call/Calls.Comp.ind $path/Admixture

#cd $path/Admixture

#subset Eigestrat to desired sampleset
sbatch --partition=long --mem=64000  -o $path/references/Admixture/subset.log --wrap="convertf -p /nfs/geo/projects/south_america/Banks_Island/references/Admixture/admixture_param.txt"

#convert eigenstrat to plink
sbatch --mem=64000 -o $path/references/Admixture/plink.log --wrap="convertf -p /nfs/geo/projects/south_america/Banks_Island/scripts/Admixture/EigenstratToPED.txt"

#create map
#cd $path/references/Admixture

#create needed input for Admixture
plink --file subset_ADMIXTURE --make-bed --out subset_ADMIXTURE
plink --bfile subset_ADMIXTURE --make-bed --out subset_ADMIXTURE --indep-pairwise 200 25 0.4 --allow-no-sex
plink --bfile subset_ADMIXTURE --extract subset_ADMIXTURE.prune.in --make-bed --out subset_ADMIXTURE.pruned --allow-no-sex


BED_FILE=$path/pop_gen_Lucia/Admixture/subset_ADMIXTURE.pruned.bed

#Run admixture with K from 3 to 15 and cross validation of 5
for K in {3..15}; do
  sbatch -c 8 --mem=64G -p long \
    --wrap="admixture --cv=5 -j8 -s time \
      $BED_FILE $K | tee $path/pop_gen_Lucia/Admixture/admixture.K${K}.log"
done


#Supervised Admixture set up

#for K in {3..3}; do
 #   CMD="cd $OUTDIR2; admixture --cv=5 -j8 --supervised -s time $BED_FILE2 $K" #Normally you should give --cv as first option to admixture
   # echo $CMD
  #   sbatch -c 8 --mem 64000 -p long --wrap="$CMD"
#done

#plot admixture; be aware that the script requires a different version than the newest installed one. You will have to change the input file manually as well as the output and  K-parameter to properly plot it. It will output a pdf file you can download. You will also have to change the ind file to tab separated before running the script or it will throw a fit
r-40-default Rscript Plot_admixturer
