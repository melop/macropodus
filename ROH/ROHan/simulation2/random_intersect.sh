sPath1="/data2/projects/zwang/m.hk/ROH/ROHan/all_het"
#sPath2="/data/projects/zwang/m.hk/GATK/final/merge_samples/ROHan/simulation2.2/"
# SZ-1 SZ-2 SZ-3 SZ-4 JX-1 JX-2 JX-3 JX-4 FJ-1 PT-1 PT-2
for s1 in JX-1 JX-2 JX-3 JX-4; do
	for s2 in SZ-1 SZ-2 SZ-3 SZ-4 MHKqns-1 MHKqns-2 MHKqns-3; do
		for s3 in FJ-1 PT-1 PT-2; do
#			echo $s1 $s2 $s3
#-----------------------------------------------------------------------------------------------------------------------------------合并连续杂和片段
{(			mkdir ${s1}_${s2}_${s3} && cd ${s1}_${s2}_${s3}; \
			sPath2="/data2/projects/zwang/m.hk/ROH/simulation2/${s1}_${s2}_${s3}"
			for i in $s1 $s2 $s3; do
				less $sPath1/$i/*.hEst.gz | awk -v OFS="\t" '{if(substr($1,1,6)=="mhkscf"){split($1,arr,"_");if (arr[2]<24 && ($3 % 50000 == 0) && $4>=25000){print $1,$2,$3,$5}}}' >> $i.scf.het.txt
				less $i.scf.het.txt | awk -v OFS="\t" '{if($4>0.00082){print $0}}' >> $i.scf.Het.bed
				bedtools merge -i $i.scf.Het.bed > tmp.$i.scf.union.Het.bed
				less tmp.$i.scf.union.Het.bed | awk -v OFS="\t" 'BEGIN {h=0.01};{print $0, h}' >> $i.scf.union.Het.bed
				less $i.scf.het.txt | awk -v OFS="\t" '{if($4<=0.00082){print $0}}' >> $i.scf.Hom.bed
				cat $i.scf.union.Het.bed $i.scf.Hom.bed > $i.scf.unionHet.het.txt
			done
#-------------------------------------------------------------------------------------------------------------------------------------统计个体间实际杂和片段交集长度

			Rfiles=($(ls $sPath2/*.scf.union.Het.bed))
			Routput="$sPath2/${s1}_${s2}_${s3}.real.intersect.bed"
			cp ${Rfiles[0]} $Routput
			for Rfile in "${Rfiles[@]:1}"; do
				bedtools intersect -a ${s1}_${s2}_${s3}.real.intersect.bed -b $Rfile > temp.real.intersect.bed
				mv temp.real.intersect.bed ${s1}_${s2}_${s3}.real.intersect.bed
			done
			less ${s1}_${s2}_${s3}.real.intersect.bed | awk 'BEGIN {len=0};{len+=($3-$2)}; END {print len}' >> ${s1}_${s2}_${s3}.real.intersect_len.txt
			rm *scf.Het.bed *.scf.union.Het.bed *.scf.Hom.bed 

#--------------------------------------------------------------------------------------------------------------------------------------将杂和片段随机放回染色体中，循环10000次
			for n in $(seq 1 10000); do
				for i in $sPath2/*.scf.unionHet.het.txt; do
					sF=$(basename $i)
					sStem=${sF/.scf.unionHet.het.txt/}
					echo $sStem
					# 获取唯一的染色体名列表
					chromosomes=$(awk '{print $1}' $i | sort -u)
	
					# 对每个染色体进行处理
					for chrom in $chromosomes; do
				# 只选择当前染色体的行并随机排序，然后更新位置
						awk -v OFS="\t" -v chrom=$chrom '$1 == chrom' $i | shuf | 
						awk -v OFS="\t" 'BEGIN {start=0} 
							{len = $3 - $2; end = start + len; print $1, start, end, $4; start = end;}' >> $sStem.$n.tem.het.txt											
					done
				done
#---------------------------------------------------------------------------------------------------------------------------------------统计重新排序后的个体间杂合片段交集长度
				files=($(ls $sPath2/*.$n.tem.het.txt))
				output="$sPath2/$n.intersect.bed"
				cp ${files[0]} $output
				less $output | awk -v OFS="\t" '{if($4>0.00082){print $0}}' >> final.$n.intersect.bed
				for file in "${files[@]:1}"; do
					sF=$(basename $file)
					sStem=${sF/.tem.het.txt/}
					less $file | awk -v OFS="\t" '{if($4>0.00082){print $0}}' >> $sStem.het.bed
					bedtools intersect -a final.$n.intersect.bed -b $sStem.het.bed > temp.intersect.$n.bed
					mv temp.intersect.$n.bed final.$n.intersect.bed
				done
				less final.$n.intersect.bed | awk 'BEGIN {len=0};{len+=($3-$2)}; END {print len}' >> ${s1}_${s2}_${s3}.intersect_len.txt
				rm *.tem.het.txt $n.intersect.bed final.$n.intersect.bed *.het.bed 
			done
			cd ..
		)} &	
		done
	done
done
wait
