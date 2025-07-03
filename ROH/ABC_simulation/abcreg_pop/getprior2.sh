cat ../abc.reps.txt | tail -n+2 | grep -v "NA" > priors.txt
#/data/software/ABCreg/src/reg -p priors.txt -d data.txt -P 3 -S 3  -b reg -t 0.001 -T #get posterior by regression
/public/software/ABCreg/src/reg -p priors.txt -d data_resq_inds.txt -P 3 -S 3  -b reg.reseq -t 0.001 -T
