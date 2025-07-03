setwd("/data2/projects/zwang/m.hk/ROH/ROHan/runs_of_hom")

# 加载必要的包
library(ggplot2)
library(ggridges)

# 获取当前工作目录
current_dir <- getwd()

# 获取当前工作目录下的所有目录名（即种群名）
population_dirs <- list.dirs(current_dir, full.names = TRUE, recursive = FALSE)

# 初始化一个空数据框，用于存储所有种群的区间长度信息
all_data <- data.frame()

# 遍历每个种群目录
for (population in population_dirs) {
  
  # 构建文件路径，文件名为种群名.ROH.bed
  file_path <- file.path(population, paste0(basename(population), ".ROH.bed"))
  
  # 读取文件，假设是以空格或制表符分隔的
  if (file.exists(file_path)) {
    bed_data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    
    # 计算区间长度 (终止位置 - 起始位置)
    bed_data$length <- bed_data$V3 - bed_data$V2
    
    # 添加种群列
    bed_data$population <- basename(population)  # 只获取目录名
    
    # 合并到总数据中
    all_data <- rbind(all_data, bed_data[, c("population", "length")])
  } else {
    message(paste("文件未找到:", file_path))
  }
}

# 绘制频率分布的山脊图
p <- ggplot(all_data, aes(x = length, y = population, fill = population)) +
  geom_density_ridges(aes(y = population, height = ..count..), 
                      scale = 1, 
                      alpha = 0.7, 
                      color = "white", 
                      stat = "binline", 
                      bins = 30) +  # 指定bin的数量
  theme_ridges() + 
  theme(legend.position = "none") + 
  labs(x = "区间长度", y = "种群", title = "区间长度的频率分布山脊图")

# 保存到 PDF 文件中
pdf("ridgeline_plot.pdf", width = 10, height = 8)
print(p)
dev.off()

