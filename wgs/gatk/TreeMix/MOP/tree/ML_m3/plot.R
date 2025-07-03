setwd("/data2/projects/zwang/m.hk/Geneflow/TreeMix/MOP/tree/ML_m3")
source("/data2/projects/zwang/m.hk/Geneflow/TreeMix/MOP/plotting_funcs.R")

pdf("TreeMix.MOP_m3.pdf")
par(cex = 1.5)
plot_tree("migrate_3")
dev.off()
#plot_resid(stem = "migrate_1", pop_order = "poporder.txt" )
