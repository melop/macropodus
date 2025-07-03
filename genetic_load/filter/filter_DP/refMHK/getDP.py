import gzip
import re

# 定义文件路径
input_file = 'allpops.g.vcf.gz'
output_file = 'allpops.gvcf.depth.bed'

# 打开输出文件
with open(output_file, 'w') as out_f:
    # 打开压缩的 VCF 文件
    with gzip.open(input_file, 'rt') as vcf_f:
        for line in vcf_f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            
            # 切分字段
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1]) - 1
            info = fields[7]
            
            # 查找 DP 信息
            match = re.search(r'\bDP=(\d+)', info)
            if match:
                dp_value = match.group(1)
            else:
                dp_value = '0'
            
            # 写入结果
            out_f.write(f"{chrom}\t{pos}\t{pos+1}\t{dp_value}\n")

print(f"Finished processing {input_file}, results saved in {output_file}.")

