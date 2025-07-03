bedtools intersect -a DSY-1.het.bed -b DSY-2.het.bed > temp1.bed
bedtools intersect -a temp1.bed -b DSY-3.het.bed > temp2.bed
bedtools intersect -a temp2.bed -b DSY-4.het.bed > temp3.bed
bedtools intersect -a temp3.bed -b DSY-5.het.bed > temp4.bed
bedtools intersect -a temp4.bed -b DSY-6.het.bed > temp5.bed
bedtools intersect -a temp5.bed -b DSY-7.het.bed > MOPdsy.het.intersected.bed
rm temp*
