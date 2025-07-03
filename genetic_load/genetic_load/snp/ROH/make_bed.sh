sPath=/data2/projects/zwang/m.hk/ROH/ROHan/all_het

for i in $(ls -d $sPath/* | grep -v "WL\|YN") $sPath/YN-mop*; do
    sH=$i/*.hEst.gz
    sStem=$(basename $i)

    # 输出文件名
    het_file="$sStem.HET.bed"
    roh_file="$sStem.ROH.bed"

    # 清空输出文件
    > $het_file
    > $roh_file

    # 处理文件
    {    less $sH | awk -v OFS="\t" -v het_file="$het_file" -v roh_file="$roh_file" '
        {
            if (substr($1,1,6) == "mhkscf") {
                split($1, arr, "_");
                if (arr[2] < 24 && $4 >= 25000) {
                    if ($5 > 0.00082) {
                        print $1, $2, $3 > het_file
                    } else {
                        print $1, $2, $3 > roh_file
                    }
                }
            }
        }'
} &
done
wait
