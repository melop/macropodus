MSMC2=/data/software/msmc2-2.1.3/build/release/msmc2

arrPops=( MOPld-5 MOPdxs-1 GL-4 MOPhn-1 YN-mop1 HK-1 MHKmls-1 MHKmls-6 JX-4.bam )
sIn=formsmc2_in
sOut=msmc2ret
recOverMu=1

mkdir -p $sOut


for sPop in "${arrPops[@]}"; do
	sFiles="${sIn}/${sPop}*/*/formsmc2.multihetsep.txt"
	sCmd="$MSMC2 -o $sOut/$sPop -r $recOverMu -i 50 -t 24 $sFiles > /dev/null "
	eval $sCmd > /dev/null 2>&1 &
done

wait

