less consurf_scores_flat.txt.gz | awk -v OFS="\t" '{print $1,$2-1,$2,$3}' > consurf_scores_flat.bed
