library(RZooRoH)
file <- "./MHKmlh_qns_MOPdsy.gvcf.PL_SNP.gen.gz"
data <- zoodata(genofile = file, zformat = "gl", min_maf = 0.05)

mix10r <- zoomodel() #2 4 8 16 32 64 128 256 512 最后一列为非HBD

# 运行分析
results <- zoorun(mix10r, data,localhbd = TRUE, nT = 12, ids = c(1))

# 打开PDF文件，指定文件名
pdf("MHKmlh_qns_MOPdsy.DSY-1.onlySNP.MAF0.05.inbreeding_coeff_hists.pdf", height=8, width=10)

# 计算近交系数（Inbreeding coefficient）
# 这里计算包含所有HBD类的近交系数
inbreeding_coeff <- 1 - results@realized[, 10]
print(inbreeding_coeff)
hist(inbreeding_coeff,nc=20,main="all HBD",xlab="Inbreeding coefficient",xlim=c(0,1),col='tomato')

rounded_table <- round(results@realized[1:14],3)
write.table(rounded_table, file = "realized.txt", sep = "\t", na = "nan")


local_hbdp <- t(results@hbdp[[1]])
write.table(local_hbdp, file = "local_hbdp.txt", sep = "\t", na = "nan")

hbdseg <- results@hbdseg
write.table(hbdseg, file = "hbdseg.txt", sep = "\t", na = "nan")

#T1 <- t(apply(results@realized[,1:7],1,cumsum)) #T<128
#print(paste("T=128", T1))
#hist(T1,nc=60,main="T=128",xlab="Inbreeding coefficient",xlim=c(0,1),col='tomato')

#T2 <- t(apply(results@realized[,1:6],1,cumsum)) #T<64
#print(paste("T=64", T2))
#hist(T2,nc=60,main="T=64",xlab="Inbreeding coefficient",xlim=c(0,1),col='tomato')

#T3 <- t(apply(results@realized[,1:5],1,cumsum)) #T<32
#print(paste("T=32", T3))
#hist(T3,nc=60,main="T=32",xlab="Inbreeding coefficient",xlim=c(0,1),col='tomato')

#T4 <- t(apply(results@realized[,1:4],1,cumsum)) #T<16
#print(paste("T=16", T4))
#hist(T4,nc=60,main="T=16",xlab="Inbreeding coefficient",xlim=c(0,1),col='tomato')



print("鉴定的HBD片段总数")
dim(results@hbdseg)[1]

print("HBD片段长度摘要")
summary(results@hbdseg$length)

print("HBD片段平均SNP数")
summary(results@hbdseg$number_snp)

# 关闭 PDF 文件
dev.off()
