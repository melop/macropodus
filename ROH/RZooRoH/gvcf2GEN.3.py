import gzip


def process_site(site_info, sample_info):
    """
    处理单个位点，直接使用 PL 字段
    :param site_info: 位点信息，包含 REF、ALT 等
    :param sample_info: 样本信息，包含 PL 和 RGQ 字段
    :return: 处理后的位点信息，包含 PL 字段值
    """
    chrom = site_info[0]
    pos = int(site_info[1])
    ref = site_info[3]
    alt = site_info[4]

    result = [chrom, "-", str(pos), ref, alt]

    format_fields = site_info[8].split(':')
    for sample in sample_info:
        sample_fields = dict(zip(format_fields, sample.split(':')))
        if 'PL' in sample_fields and sample_fields['PL'] != '.':
            pl_values = sample_fields['PL'].split(',')
        elif 'RGQ' in sample_fields:
            rgq = sample_fields['RGQ']
            if rgq == '.':
                pl_values = ['999', '999', '999']
            else:
                rgq = int(rgq)
                if rgq == 99:
                    pl_values = ['0', '99', '999']
                else:
                    pl_values = ['999', '999', '999']
        else:
            pl_values = ['999', '999', '999']

        result.extend(pl_values)

    return result


def extract_sites(gvcf_file, output_file):
    """
    提取符合条件的位点并保存到输出文件
    :param gvcf_file: 输入的 GVCF 文件路径
    :param output_file: 输出文件路径
    """
    if gvcf_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    # 以写入模式打开 gz 文件
    with open_func(gvcf_file, mode) as file, gzip.open(output_file, 'wt') as out:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # 判断是否为二等位基因的 SNP 位点
            is_snp = len(ref) == 1 and len(alt) == 1 and alt != '.'

            if is_snp:
                processed_site = process_site(fields, fields[9:])
                # 使用空格分隔列
                out.write(' '.join(processed_site) + '\n')


if __name__ == "__main__":
    gvcf_file = "MHKmlh_qns_MOPdsy.gvcf.gz"
    output_file = "MHKmlh_qns_MOPdsy.gvcf.PL_SNP.gen.gz"
    extract_sites(gvcf_file, output_file)

