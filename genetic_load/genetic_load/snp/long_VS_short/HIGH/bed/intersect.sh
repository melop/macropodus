sPath=/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/snp/ROH
GL=./allpops.refMHK.NONSYN.polarized.out.filtered.withAD.dedup.withConsurf.dedup.bed
for i in $sPath/*.ROH.bed; do
    sF=$(basename $i)
    sStem=${sF%%.*}
    # 合并文件
    bedtools merge -i $i > $sStem.ROH.merged.bed
    short=$sStem.shortROH.bed
    long=$sStem.longROH.bed
    # 根据长度分割文件
    awk -v OFS="\t" -v short="$short" -v long="$long" '{if($3-$2<1000000){print $0 > short}else{print $0 > long}}' $sStem.ROH.merged.bed

    # 判断 short 文件是否存在
    if [ -f "$short" ]; then
        bedtools intersect -a $GL -b $short > allpops.refMHK.NONSYN.filtered.$sStem.shortROH.bed
    else
        echo "$sStem 不存在 short 文件"
    fi

    # 判断 long 文件是否存在
    if [ -f "$long" ]; then
        bedtools intersect -a $GL -b $long > allpops.refMHK.NONSYN.filtered.$sStem.longROH.bed
    else
        echo "$sStem 不存在 long 文件"
    fi
done

