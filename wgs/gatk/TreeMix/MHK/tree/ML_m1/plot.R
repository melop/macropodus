setwd("/data2/projects/zwang/m.hk/Geneflow/TreeMix/MHK/tree/ML_m1")
source("/data2/projects/zwang/m.hk/Geneflow/TreeMix/MHK/plotting_funcs.R")

pdf("TreeMix.MHK_m1.pdf")
par(cex = 1.5)
plot_tree("migrate_1")
dev.off()
#plot_resid(stem = "migrate_1", pop_order = "poporder.txt" )
