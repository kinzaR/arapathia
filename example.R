##### libraries #####
library("igraph")
#####################

##### Resources #####
source("main.R")
source("stats.R")
source("chart.R")
source("utils.R")
#####################

##### load graphs  #####
load("files/meta_graph_info_athPremium.RData")
#load("files/meta_graph_info_athKEGG.RData")
################################

##### load data #####
des <- read.experimental.design("examples/GSE77017_24__design.txt")
exp <- read.expression.matrix("examples/GSE77017_24__exp.txt")
#####################
exp <- cbind(exp[,7:9],exp[,10:12])
#####################

##### prepare data #####
exp <- normalize.data(exp, by.quantiles=F, by.gene=F, percentil=F)
########################

##### add missing genes #####
exp <- add.missing.genes(exp, genes=metaginfo$all.genes) 
#############################

##### araPathia #####
results <- prettyways(exp, metaginfo$pathigraphs, verbose=T)
####################

##### descriptive statistics #####
path_vals <- results$all$effector.path.vals
sample_group <- des[colnames(path_vals),"group"]
groups <- unique(sample_group)
n <- length(sample_group)
##################################

##### all paths #####
plot.heatmap(path_vals, sample_type=sample_group)
#####################

##### select most relevant using the T test #####
comp <- do.t.test(path_vals,sample_group,g1=groups[1],g2=groups[2], paired= F)
##### select most relevant using the wilcox test #####
#comp <- do.wilcox(path_vals, sample_group,g1=groups[1],g2=groups[2], paired= F)
##### select most relevant using the Limma test #####
#comp <- do.limma.test(path_vals, sample_group,g1=groups[1],g2=groups[2], paired= F)

ranked_path_vals <- path_vals[order(comp$p.value,decreasing = F),]
table(comp$FDRp.value<0.05)
plot.heatmap(ranked_path_vals[1:15,],sample_type=sample_group)
pca_model <- do.pca(ranked_path_vals[1:n,])
plot.multiple.pca(pca_model,colors = c("blue","red")[1+(sample_group==groups[1])])
#################################################

