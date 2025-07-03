import gzip

def process_vcf_to_bed(file_path, output_file):
    with gzip.open(file_path, 'rt') as f, open(output_file, 'w') as fout:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            columns = line.strip().split('\t')
            if len(columns) < 10:
                continue  # Ensure there are enough columns to process
            chrom = columns[0]  # 第一列：染色体信息
            start = int(columns[1])-1    # 第二列：位置信息
            end = columns[1]
            format_col = columns[8].split(':')
            ad_index = format_col.index('AD')  # Find index of AD in format column
            sample_data = columns[9:]  # Get sample data columns
            ad_values = [sample.split(':')[ad_index] for sample in sample_data]
            bed_format = f"{chrom}\t{start}\t{end}\t" + '\t'.join(ad_values)
            fout.write(bed_format + '\n')

# 示例用法
vcf_file_path = '/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/allpops.merged.snpeff.vcf.gz'
output_file = 'allpops.merged.snpeff.vcf.AD.bed'
process_vcf_to_bed(vcf_file_path, output_file) 
