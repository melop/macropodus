clade=MOP
#clade=Bifurcation
#clade=Macropodus
#clade=MHK
#clade=BTP
#clade=Callopanchax
#clade=Aphyosemion
#clade=Scriptaphyosemion
#clade=Archiaphyosemion
#clade=Epiplatys
#clade=Fundulopanchax
#clade=PronothoNothos
#clade=AKN
suffix=
clade=${clade}${suffix}

cat Relax_${clade}/ret*.txt | grep -v "RELAX error" | sort -t $'\t' -k 5n,5n  > sum_${clade}.txt
cat Relax_${clade}/ret*.txt | grep "RELAX error" | sort -t $'\t' -k 1,1  > sum_${clade}_err.txt

echo $clade `wc -l sum_${clade}.txt`
echo $clade error `cat Relax_${clade}/ret*.txt | grep "RELAX error" | wc -l`
