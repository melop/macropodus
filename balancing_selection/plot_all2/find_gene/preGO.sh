# 输出文件名
output_file="balanced.gene.txt"

# 首先清空或创建文件
echo -e "#Population\tGeneList" > $output_file

# 遍历每个文件并提取基因列表
for i in *.intersected.HET.balanced.gff3; do
    pop=${i/.intersected.HET.balanced.gff3/}
    
    # 初始化变量用于存储所有基因名
    gene_list=""
    
    # 使用 awk 提取 gene 名并拼接到 gene_list
    gene_list=$(awk -v OFS="\t" '{split($9, a, ";")  
        for (i in a) {
            if (a[i] ~ /^gene=/) {
                gene_name = substr(a[i], index(a[i], "=") + 1)
                # 检查 gene_name 是否为 "unknown"
                if (gene_name != "unknown") {
                    if (gene_list == "") {
                        gene_list = gene_name
                    } else {
                        gene_list = gene_list","gene_name
                    }
                }
            }
        }
    } END {print gene_list}' $i)
    
    # 输出种群名和基因列表到文件
    echo -e "$pop\t$gene_list" >> $output_file
done

