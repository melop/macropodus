def parse_ad(ad):
    """解析AD信息，计算实际的AD值。"""
    if ',' in ad:
        # 格式为 "*/", "*/*" 或 "*/." 或 "./*"
        parts = ad.split(',')
        # 将 "." 转换为 0
        ad1 = int(parts[0]) if parts[0] != '.' else 0
        ad2 = int(parts[1]) if parts[1] != '.' else 0
        return ad1 + ad2
    elif ad == '.':
        # 处理单独的 "."
        return 0
    else:
        # 单独的数字
        return int(ad)

def process_bed_file(input_file, output_file):
    """处理BED文件，保留AD值大于等于2的位点。"""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # 去除行末的换行符，并用制表符分隔
            row = line.strip().split('\t')
            
            chrom = row[0]
            start = row[1]
            end = row[2]
            ad_values = row[3:]
            
            # 计算所有个体的AD值，并检查是否有一个大于等于2
            ad_numeric_values = [parse_ad(ad) for ad in ad_values]
            if any(ad >= 2 for ad in ad_numeric_values):
                # 如果有任意AD值 >= 2，则保留该位点
                outfile.write('\t'.join([chrom, start, end] + ad_values) + '\n')

# 输入文件名和输出文件名
input_file = 'allpops.refMOP.cds.AD.bed'
output_file = 'allpops.refMOP.cds.AD.filtered.bed'

# 处理文件
process_bed_file(input_file, output_file)

