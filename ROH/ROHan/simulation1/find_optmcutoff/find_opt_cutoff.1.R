# 设置工作目录，根据实际情况调整路径
setwd("/data2/projects/zwang/m.hk/ROH/simulation1.2/find_optmcutoff")

# 定义不同测序深度文件的文件名前缀列表，根据实际情况补充完整
depth_prefixes <- c("X10", "X12", "X14", "X16", "X18", "X2", "X20", "X22", "X24", "X4", "X6", "X8") 

# 用于存储每个文件对应的结果（最大y值和对应的x值）
all_results <- data.frame(Depth = character(), Max_y = numeric(), Cutoff_at_max_y = numeric(), stringsAsFactors = FALSE)

# 循环处理每个测序深度的文件
for(depth_prefix in depth_prefixes) {
  # 构建文件名
  het_file_name <- paste0(depth_prefix, ".rohan.het.bed")
  win_file_name <- paste0(depth_prefix, ".roh.wins.bed")
  
  # 读取杂合位点文件
  datROHAN <- read.table(het_file_name, sep = "\t", header = F, stringsAsFactors = F)
  datROHAN <- datROHAN[complete.cases(datROHAN), ]
  datROHAN <- datROHAN[datROHAN$V1 %in% paste0("mhkscf_", 1:23), ]
  
  colnames(datROHAN) <- c('chr', 'winstart', 'winend', 'esthet')
  datROHAN$coord <- paste(datROHAN$chr, datROHAN$winstart, datROHAN$winend, sep = ":")
  
  # 读取纯合区域交集文件
  datROHWin <- read.table(win_file_name, sep = "\t", header = F)
  datROHWin$coord <- paste(datROHWin$V1, datROHWin$V2, datROHWin$V3, sep = ":")
  
  datROHWin$ROH <- TRUE
  datROHAN$ROH <- FALSE
  datROHAN$ROH[datROHAN$coord %in% datROHWin$coord ] <- TRUE
  
  nFROH <- sum(datROHAN$ROH) / nrow(datROHAN)
  
  # 定义Cutoff值的范围
  arrCutoffs <- seq(1e-5, 0.1, 1e-5)
  datSimRets <- NULL
  
  for(nCutoff in arrCutoffs) {
    nROH_on_ROH <- sum((datROHAN$esthet < nCutoff) & (datROHAN$ROH))
    nROH_on_HET <- sum((datROHAN$esthet < nCutoff) & (datROHAN$ROH == F))
    nHET_on_ROH <- sum((datROHAN$esthet >= nCutoff) & (datROHAN$ROH))
    nHET_on_HET <- sum((datROHAN$esthet >= nCutoff) & (datROHAN$ROH == F))
    nEstFROH <- (nROH_on_ROH + nROH_on_HET) / (nROH_on_ROH + nROH_on_HET + nHET_on_ROH + nHET_on_HET)
    datSimRets <- rbind(datSimRets, data.frame(Cutoff = nCutoff, nROH_on_ROH = nROH_on_ROH, nROH_on_HET = nROH_on_HET, nHET_on_ROH = nHET_on_ROH, nHET_on_HET = nHET_on_HET, estroh = nEstFROH ))
  }
  
  datSimRets$correct <- (datSimRets$nROH_on_ROH + datSimRets$nHET_on_HET) / rowSums(datSimRets[, 2:5])
  
  # 找到最大y值对应的索引
  max_correct_idx <- which.max(datSimRets$correct)
  # 获取对应的x值（Cutoff）
  x_value_at_max_y <- datSimRets$Cutoff[max_correct_idx]
  max_y_value <- datSimRets$correct[max_correct_idx]
  
  # 将当前文件的结果添加到总结果数据框中
  all_results <- rbind(all_results, data.frame(Depth = depth_prefix, Max_y = max_y_value, Cutoff_at_max_y = x_value_at_max_y))
}

print(all_results)

write.table(all_results, "DP_map_Cutoff.txt", sep = "\t", row.names = FALSE)


#-------------------------------------------------------------------------------------
all_results <- read.table("DP_map_Cutoff.txt", sep = "\t", header = T)
all_results$Depth <- as.numeric(gsub("X", "", all_results$Depth))
all_results <- all_results[order(all_results$Depth), ]
# 获取Cutoff_at_max_y列的范围（最小值和最大值）
cutoff_range <- range(all_results$Cutoff_at_max_y)
# 获取Max_y列的范围（最小值和最大值）
# 直接指定Max_y的范围
max_y_range <- c(0.8521, 0.8584)

# 指定保存的PDF文件名，可根据实际需求修改
pdf_file_name <- "Cutoff_map_accuracy_differentDP.pdf"

# 开启PDF图形设备，设置文件名、宽度和高度等参数（单位为英寸，可根据需求调整尺寸）
pdf(file = pdf_file_name, width = 10, height = 8)


# 通过par函数设置字体大小相关参数
par(
  cex.main = 2,  # 设置主标题字体大小，这里设置为1.5倍默认大小，可根据需求调整
  cex.lab = 2,  # 设置坐标轴标签字体大小，为1.2倍默认大小
  cex.axis = 2,  # 设置坐标轴刻度标签字体大小，为1.2倍默认大小
  mar = c(8, 8, 8, 8)  # 设置页边距，依次为下、左、上、右的空白行数，可根据需求调整
)

# 绘制基础图形，先绘制Cutoff_at_max_y对应的蓝色曲线（使用左侧y轴）
plot(
  all_results$Depth,
  all_results$Cutoff_at_max_y,
  type = "l",
  col = "blue",
  xlab = "Depth",
  ylab = "Cutoff_at_max_y",
  main = "Cutoff_at_max_y and Max_y vs Depth",
  ylim = cutoff_range  # 设置左侧y轴范围
)

# 使用par函数添加新的坐标轴（右侧y轴，对应Max_y）
par(new = TRUE)
# 绘制Max_y对应的红色曲线（使用右侧y轴）
plot(
  all_results$Depth,
  all_results$Max_y,
  type = "l",
  col = "red",
  xlab = "",  # 清空x轴标签，避免重复显示
  ylab = "",  # 清空y轴标签，避免重复显示
  axes = FALSE,  # 不绘制默认坐标轴，后续手动添加合适的右侧坐标轴
  ylim = max_y_range  # 设置右侧y轴范围
)

# 添加右侧y轴（对应Max_y）
axis(
  side = 4,  # 右侧坐标轴
  at = seq(from = max_y_range[1], to = max_y_range[2], length.out = 8),  # 设置合适的刻度
  labels = seq(from = max_y_range[1], to = max_y_range[2], length.out = 8),  # 设置刻度对应的标签
  col.axis = "red",  # 设置坐标轴颜色与对应曲线颜色一致（红色），方便区分
  las = 1,  # 让刻度标签平行于坐标轴
  cex.axis = 2  # 设置右侧坐标轴刻度标签字体大小，保持与之前设置一致
)

# 添加图例
legend(
  "topright",  # 图例位置在右上角
  legend = c("Cutoff_at_max_y", "Max_y"),  # 图例内容对应两条曲线
  col = c("blue", "red"),  # 对应曲线的颜色
  lty = 1,  # 线条类型，这里都是实线
  cex = 1.5  # 图例字体大小，可根据需求调整
)

# 关闭PDF图形设备，完成图形保存
dev.off()
