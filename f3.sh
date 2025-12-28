#!/usr/bin/env bash

## F3 test
#modify parameter here
sbatch --mem=32000 -o /nfs/geo/projects/south_america/Banks_Island/pop_gen_Lucia/f3/qp3Pop_site_cap_b.log --wrap="qp3Pop -p /nfs/geo/projects/south_america/Banks_Island/pop_gen_Lucia/f3/parameters/f3_param_cap_bleached.txt"

#Get results as table
echo -e "Source1\tSource2\tTarget\tf3\tstderr\tZ\tSNPs" > f3_cap_b_results_table.out
grep 'result:' qp3Pop_site_cap_b.log | awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8}' >> f3_cap_b_results_table.out
