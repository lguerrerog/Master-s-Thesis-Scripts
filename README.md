# Master-s-Thesis-Scripts
The scripts used for my Master's thesis are grouped by topic or test.

- *Population_Selection* Contains the script for defining the populations to include in the analyses (_population_selection.ipynb_).
- *Admixture* Contains the script to convert genotypes to Plink format (_EigenstratToPED.txt_), the parameters to run the program (_admixture_param.txt_), the script to run ADMIXTURE (_Admixture.sh_), and one to visualise cross-validation and admixture results (_CVerror_plot.R_ and _Admixture_plot.R_ respectively).
- *F_tests* Contains scripts to run F-statistics with qp3Pop (_f3.sh_) and qpDstat(_f4.sh_) and the visualisation of F4 results (_f4test_plots.ipynb_).
- *Genotyping* Contains the script to run variant calling (_VariantCall_merged_BAM.sh_) and to count the SNPs obtained (_CountSNPsGenotyping.sh_).
- *PCA* Contains the scripts to run worldwide PCA with smartpca (_Smart_PCA_full.txt_) and for plotting PCA results (_Worldwide_PCA_visualisation.R_).
- *Merge* Contains the scripts to merge Umingmak results from the same treatment (_mergeBamFiles.sh_) and to merge the Umingmak genotypes with the 1240K dataset (_parameters_merge1240K.txt_).
- *Schmutzi* Contains the script to run Schmutzi on the mitochondrial DNA (_Schmutzi_Pipeline.sh_).
- *SexDet_XCont* Contains the script to run sex determination and X contamination (_SexDet_XCont_Lucia.sh_) and the code to calculate the chromosome's rate (_sexDetermination.hg19.awk_).

