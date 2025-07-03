cat genespace.orthogroups.txt | awk -F'\t' '{split($11, a, "|"); print $1"\t"a[3]}'  > groupid2transid.MHK.txt
