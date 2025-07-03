import re

# 读取 GFF3 文件到内存中
with open('/data/projects/rcui/mhk/annotations/funannotate/train/mhk.longest_isoform.genesymbol.gff3', 'r') as gff_file:
    gff_lines = gff_file.readlines()

# 打开输入文件和输出文件
with open('/data2/projects/zwang/test/DeepVatiant/formed_SingleCopyOrthologs.MHK_MOP.txt', 'r') as orthologs_file, \
     open('refMHK.cds.bed', 'w') as output_file:
    
    # 遍历每个基因名称
    for line in orthologs_file:
        # 提取 mhk 名称
        mhk = line.split()[1]
        
        # 计算索引名称中 "mhk" 的数量
        mhk_count_mhk = mhk.count('mhk')

        # 遍历 GFF3 文件的每一行
        for gff_line in gff_lines:
            # 检查 gff3 的第3列是否为"CDS"
            fields = gff_line.split('\t')
            if len(fields) > 8 and fields[2] == "CDS":
                attributes = fields[8]
                
                # 提取 ID 字段并计算 ID 字段中包含的 "mhk" 数量
                id_match = re.search(r'ID=([^;]+)', attributes)
                if id_match:
                    id_value = id_match.group(1)
                    mhk_count_id = id_value.count('mhk')
                    
                    # 检查 mhk 是否在 ID 字段中，并且索引名称中 "mop" 的数量与 ID 中 "mop" 的数量一致
                    if mhk in id_value and mhk_count_mhk == mhk_count_id:
                        # 输出符合条件的 CDS 区间到输出文件
                        output_file.write(f"{fields[0]}\t{int(fields[3])-1}\t{fields[4]}\t{mhk}\n")

