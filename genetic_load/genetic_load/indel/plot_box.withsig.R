# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel")

# 加载所需的包
library(ggplot2)
library(tidyr)
library(scales)
library(ggsignif) # 用于添加显著性标记

dat <- read.table("genetic_load_of_individual.allpops.refMHK.txt", header = TRUE) 
datSamples <- read.table("header_samples.csv", header=F)
# 添加种群信息
dat$population <- rep(c("MOPdsy", "MOPfq_pt", "MOPgl", "MHKhk", "MHKhn", "MHKjx", "MHKmls_pop1", "MHKmls_pop2", "MHKmls_pop1", "MHKmls_pop2", "MHKmls_pop1", "MHKmls_pop2", "MHKmlh_qns", "MOPdxs_sg", "MOPhn", "MOPld", "MOPdxs_sg", "MOPfq_pt", "MHKmlh_qns", "MOPyn"), 
                      times = c(7, 1, 4, 3, 8, 4, 4, 4, 3, 1, 2, 8, 3, 2, 7, 5, 8, 2, 4, 5))

# 将种群列转换为有序因子
dat$population <- factor(dat$population, 
                         levels = c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKhk", 
                                    "MHKmls_pop1", "MHKmls_pop2", "MOPyn", "MOPhn", 
                                    "MOPgl", "MOPdsy", "MOPld", 
                                    "MOPdxs_sg", "MOPfq_pt"))

# 计算原始 y 轴范围
y_min <- min(dat$standardized_genetic_load)
y_max <- 2*max(dat$standardized_genetic_load)

# 绘制箱线图
p <- ggplot(dat , aes(x = population, y = standardized_genetic_load, fill = population)) +
  geom_boxplot() +
  labs(title = "Standardized_Genetic_Load by Population",
       x = "Population",
       y = "Standardized_Genetic_Load") +
  scale_y_continuous(labels = scientific, expand = c(0, 0)) +  # 纵轴使用科学计数法
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 24),  # x轴文本字体大小
        axis.text.y = element_text(size = 24),  # y轴文本字体大小
        axis.title.x = element_text(size = 24),  # x轴标题字体大小
        axis.title.y = element_text(size = 25),  # y轴标题字体大小
        plot.title = element_text(size = 24),  # 标题字体大小和加粗
        legend.position = "none")  # 隐藏图例

# 获取所有种群
all_pops <- levels(dat$population)

# 生成所有种群的两两组合
pop_combinations <- combn(all_pops, 2, simplify = FALSE)

# 存储显著性标记信息
signif_labels <- list()
# 用于记录 MHK 和 MOP 各自的计数
MHK_count <- 0
MOP_count <- 0


# 进行 Wilcoxon 秩和检验并添加显著性标记
for (comb in pop_combinations) {
  i <- comb[1]
  j <- comb[2]
  
  # 仅对 MHK 内部或 MOP 内部的种群组合进行检验
  if ((startsWith(i, "MHK") && startsWith(j, "MHK")) || (startsWith(i, "MOP") && startsWith(j, "MOP"))) {
    data_i <- dat$standardized_genetic_load[dat$population == i]
    data_j <- dat$standardized_genetic_load[dat$population == j]
    test_result <- t.test(data_i, data_j)
    p_value <- test_result$p.value
    
    # 根据 p 值添加显著性标记
    if (p_value < 0.01) {
      sig_mark <- "**"
    } else if (p_value < 0.05) {
      sig_mark <- "*"
    } else {
      sig_mark <- ""
    }
    if (startsWith(i, "MHK") && startsWith(j, "MHK")) {
      MHK_count <- MHK_count + 1
      y_pos <- max(dat$standardized_genetic_load) + 0.1 * diff(range(dat$standardized_genetic_load)) * MHK_count
    } else {
      MOP_count <- MOP_count + 1
      # 调整 MOP 种群的 y 位置，这里乘以 2 让它们的位置更高
      y_pos <- max(dat$standardized_genetic_load) + 0.1 * diff(range(dat$standardized_genetic_load)) * MOP_count
    }
    
    # 存储标记信息
    signif_labels[[length(signif_labels) + 1]] <- list(
      annotations = sig_mark,
      y_position = y_pos,
      xmin = which(all_pops == i),
      xmax = which(all_pops == j)
    )
  }
}

# 将显著性标记信息转换为数据框
signif_df <- do.call(rbind, lapply(signif_labels, as.data.frame))

# 过滤掉无显著性标记的数据
signif_df <- signif_df[signif_df$annotations != "", ]

# 添加显著性标记到图中
p <- p + geom_signif(data = signif_df, 
                     aes(annotations = annotations, xmin = xmin, xmax = xmax, y_position = y_position),
                     manual = TRUE, 
                     textsize = 5, 
                     inherit.aes = FALSE)


# 显示图形
print(p)
# 保存图形为 PDF
ggsave("genetic_load_boxplot.indel.with_signinf.pdf", plot = p, width = 12, height = 8, units = "in")

# 提取 MOPgl 和 MOPfq_pt 的数据
pop1 <- "MOPgl"
pop2 <- "MOPld"
data_pop1 <- dat$standardized_genetic_load[dat$population == pop1]
data_pop2 <- dat$standardized_genetic_load[dat$population == pop2]
# 进行 t 检验
test_result <- t.test(data_pop1, data_pop2)
p_value <- test_result$p.value

