sPop=MOPdsy
cat *.het.bed | sort -k1,1 -k2,2n > temp.bed
bedtools merge -i temp.bed > $sPop.het.union.bed
rm temp.bed
