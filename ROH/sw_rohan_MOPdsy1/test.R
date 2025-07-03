# 定义染色体名称
arrChr <- paste0("mhkscf_", 1:23)
# 读取数据
datHet <- read.table(gzfile('DSY-1.50000.100.het.tsv.gz', 'rt'), header = F)
# 设置纯合判定阈值
nHomCutoff <- 10^(-3.9)
# 定义输出文件名
sOut <- "test.1e-3.9.tsv"

# 写入表头
write.table(data.frame(chr = 'chr', ambstart = 'ambstart', homstart = 'homstart', homend = 'homend', ambend = 'ambend', stringsAsFactors = F), file = sOut, row.names = F, col.names = F, sep = "\t", quote = F)

# 遍历每条染色体
for (sChr in arrChr) {
  # 提取当前染色体的数据
  datHetChr <- datHet[datHet$V1 == sChr, ]
  # 初始化当前纯合片段的起始和结束位置
  hom_start <- -1
  hom_end <- -1
  
  # 遍历当前染色体的每个窗口
  for (i in 1:nrow(datHetChr)) {
    # 获取窗口的中间位置
    nPosMid <- datHetChr[i, 4]
    # 获取窗口的杂合度
    nTheta <- datHetChr[i, 9]
    # 判断窗口状态，默认杂合
    nState <- ifelse(nTheta <= nHomCutoff, 0, 1)
    
    if (nState == 0) {
      if (hom_start == -1) {
        # 如果纯合片段还未开始，记录起始位置
        hom_start <- nPosMid
      }
      # 更新纯合片段的结束位置
      hom_end <- nPosMid
    } else {
      if (hom_start != -1) {
        # 如果当前纯合片段结束，写入文件
        write.table(data.frame(chr = sChr, ambstart = hom_start, homstart = hom_start, homend = hom_end, ambend = hom_end, stringsAsFactors = F), file = sOut, append = T, row.names = F, col.names = F, sep = "\t", quote = F)
        # 重置纯合片段的起始和结束位置
        hom_start <- -1
        hom_end <- -1
      }
    }
  }
  
  # 处理染色体末尾的纯合片段
  if (hom_start != -1) {
    write.table(data.frame(chr = sChr, ambstart = hom_start, homstart = hom_start, homend = hom_end, ambend = hom_end, stringsAsFactors = F), file = sOut, append = T, row.names = F, col.names = F, sep = "\t", quote = F)
  }
}

