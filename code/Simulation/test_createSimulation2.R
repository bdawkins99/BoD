# testing createSimulation2() function

setwd("C:/Users/bdawk/Documents/KNN_project_output/R_code_and_output/code")
source('createSimulation2.R')

num.samples <- 100
num.variables <- 10

dataset <- createSimulation2(num.samples=num.samples,
                             num.variables=num.variables,
                             pct.imbalance=0.5,
                             pct.signals=0.2,
                             main.bias=0.4,
                             interaction.bias=1,
                             hi.cor=0.95,
                             lo.cor=0.2,
                             mix.type="main-interactionScalefree",
                             label="class",
                             sim.type="mixed",
                             pct.mixed=0.5,
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             plot.graph=F,
                             verbose=T)
dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
dats <- dats[order(dats[,ncol(dats)]),]

mycor <- cor(dats[,-ncol(dats)])
mean(abs(mycor[upper.tri(mycor)]))

adj <- dataset$A.mat
colnames(adj) <- colnames(mycor)
row.names(adj) <- row.names(mycor)

d.mat <- as.matrix(dist(dats, method="euclidean", upper=T))
dist.vec <- d.mat[upper.tri(d.mat)]
hist(dist.vec, breaks=50, freq=F)

mean(d.mat[upper.tri(d.mat)])
sqrt(2*num.variables - 1)
var(d.mat[upper.tri(d.mat)])

########################################################################

dat.tmp <- dats[,-ncol(dats)]
dat.tmp <- dat.tmp[,sort(colnames(dat.tmp))]
adj.tmp <- adj[sort(colnames(adj)),sort(colnames(adj))]
mycor.ctrl <- cor(dat.tmp[1:50,])
mycor.case <- cor(dat.tmp[51:100,])

library(reshape2)
library(ggplot2)

blah <- melt(adj.tmp)
melted_cormat <- cbind(melt(mycor.ctrl), adj.value=blah[,3])
head(melted_cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_text(aes(Var1, Var2, label=adj.value), color="black", size=4)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+
  ggtitle("Control Data Feature Correlation Heatmap") +
  theme(axis.text.x = element_text(size=10,face="bold",angle=90,vjust=0.25),
        axis.text.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        legend.title.align = 0.5)

blah <- melt(adj.tmp)
melted_cormat <- cbind(melt(mycor.case), adj.value=blah[,3])
head(melted_cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_text(aes(Var1, Var2, label=adj.value), color="black", size=4)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+
  ggtitle("Case Data Feature Correlation Heatmap") +
  theme(axis.text.x = element_text(size=10,face="bold",angle=90,vjust=0.25),
        axis.text.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        legend.title.align = 0.5)