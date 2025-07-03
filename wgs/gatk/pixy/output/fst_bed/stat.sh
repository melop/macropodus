#caculate ave dxy of whole genome
#cat pixy_output_dxy.txt | cut -f1,2 |sort|uniq > pops_dxy.list
ln -s ../pops_dxy.list ./
ln -s ../pixy_output_fst.txt ./
cat pops_dxy.list | while read line
do
    a=$(echo $line | awk '{print $1}');
    b=$(echo $line | awk '{print $2}');
    touch $a.$b.fst.bed
    awk -v OFS="\t" '$1 == "'$a'" && $2 == "'$b'" {print $3,$4-1,$5,$6}' pixy_output_fst.txt >> $a.$b.fst.bed
done

# 使用awk直接处理pops_dxy.list和pixy_output_dxy.txt文件，避免在循环中多次使用cat命令
#awk '
#    BEGIN {
        # 设置输出字段分隔符为制表符
#        OFS="\t";
#    }
#    NR==FNR {
#        # 处理pops_dxy.list文件，将每行的前两列存储为数组
#        pop1[$1] = $2;
#        next;
#    }
#    {
#        # 处理pixy_output_dxy.txt文件
#        a = $1;
#        b = $2;
#        if (a in pop1 && b == pop1[a]) {
#            # 根据pop1和pop2的值构建文件名
#            file_name = a "." b ".dxy.bed";
#            print $3, $4 - 1, $5, $6 >> file_name;
#        }
#    }
#' pops_dxy.list pixy_output_dxy.txt
