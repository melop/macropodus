# 设置工作目录
setwd("/data2/projects/zwang/m.hk/Snpeff/all_MHK_MOP/rho_V2.0/filtered_SOL_3.5/Individual/indel")

# 加载所需的包
library(ggplot2)
library(tidyr)
library(scales)

dat <- read.table("genetic_load_of_individual.allpops.refMHK.txt", header=TRUE) 
# 添加种群信息
#dat$population <- rep(c("MHKhn", "MHKmlh_qns", "MHKmls_pop1", "MHKmls_pop2", "MOPdsy", "MOPdxs_sg", "MOPhn"), times = c(8, 7, 9, 13, 7, 10, 7)) 上次错排的顺序
dat$population <- rep(c("MOPdsy", "MOPfq_pt", "MOPgl", "MHKhk", "MHKhn", "MHKjx", "MHKmls_pop1", "MHKmls_pop2", "MHKmls_pop1", "MHKmls_pop2", "MHKmls_pop1", "MHKmls_pop2", "MHKmlh_qns", "MOPdxs_sg", "MOPhn", "MOPld", "MOPdxs_sg", "MOPfq_pt", "MHKmlh_qns", "MOPyn"), 
                      times = c(7, 1, 4, 3, 8, 4, 4, 4, 3, 1, 2, 8, 3, 2, 7, 5, 8, 2, 4, 5))

# 将种群列转换为有序因子
dat$population <- factor(dat$population, 
                         levels = c("MHKhn", "MHKjx", "MHKmlh_qns", "MHKhk", 
                                    "MHKmls_pop1", "MHKmls_pop2", "MOPyn", "MOPhn", 
                                    "MOPgl", "MOPdsy", "MOPld", 
                                    "MOPdxs_sg", "MOPfq_pt"))

# 绘制箱线图
p <- ggplot(dat , aes(x = population, y = standardized_genetic_load, fill = population)) +
  geom_boxplot() +
  labs(title = "Standardized_Genetic_Load by Population",
       x = "Population",
       y = "Standardized_Genetic_Load") +
  scale_y_continuous(labels = scientific) +  # 纵轴使用科学计数法
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # x轴文本字体大小
        axis.text.y = element_text(size = 16),  # y轴文本字体大小
        axis.title.x = element_text(size = 20),  # x轴标题字体大小
        axis.title.y = element_text(size = 20),  # y轴标题字体大小
        plot.title = element_text(size = 24, face = "bold"),  # 标题字体大小和加粗
        legend.position = "none")  # 隐藏图例

# 显示图形
print(p)

# 保存图形为PDF
pdf_file <- "allpops.refMHK.SV.Individuals_Genetic_Load_Boxplot.pdf"
ggsave(pdf_file, plot = p, width = 10, height = 6)
