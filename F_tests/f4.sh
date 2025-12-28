#!/usr/bin/env bash

## F4 test
#modify parameter here
sbatch --mem 64000 --partition=long -o /nfs/geo/projects/south_america/Banks_Island/pop_gen_Lucia/f4/f4_contamination_test.log --wrap="qpDstat -p /nfs/geo/projects/south_america/Banks_Island/pop_gen_Lucia/f4/parameters/f4_contamination_test.txt"

#Get results table 
echo -e "Pop1\tPop2\tPop3\tPop4\tD-stat\tZ\tBABA\tABBA\tSNPs" > f4_contamination_test.out
grep 'result:' f4_contamination_test.log | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10}' >> f4_contamination_test.out
