if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
library(devtools)

if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC") # build_vignettes = TRUE)
}
if (!("cncv" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/cncv", build_vignettes = TRUE)
}
library(privateEC)  # used to simulate data
library(cncv)
library(npdr) # for knnSURF()
library(caret)
library(nleqslv)
library(pracma)

dist_moments <- function(x){
  
  mu.d <- numeric()
  sd.d <- numeric()
  for(i in 1:nrow(x)){
    
    mu.d[i] <- mean(x[i,-i])
    sd.d[i] <- sd(x[i,-i])
    
  }
  
  list(mu=mu.d, sigma=sd.d)
  
}

letsSimulate <- T   # F to use previously simulated data
class.lab <- "class"
writeData <- F  # usually the same as letsSimulate
writeResults <- F

num.samp <- 100
num.attr <- 1000
pct.signals <- 0.1
bias <- 0.4
#sim.type <- "mainEffect"
sim.type <- "interactionErdos"

cv.acc.predk <- numeric()
cv.acc.naivek <- numeric()
validation.acc.predk <- numeric()
validation.acc.naivek <- numeric()

attr.recall.predk <- numeric()
attr.recall.naivek <- numeric()
attr.precision.predk <- numeric()
attr.precision.naivek <- numeric()

attr.tpr.predk <- numeric()
attr.tpr.naivek <- numeric()
attr.fpr.predk <- numeric()
attr.fpr.naivek <- numeric()

alphas <- numeric()
neighbor.probs <- numeric()
k.alphas <- numeric()

set.seed(1989)
for(iter in 1:30){
  cat("Begin Iteration: ", iter, "\n")
  
  if (letsSimulate == TRUE){
    sim.data <- createSimulation(num.samples = num.samp, num.variables = num.attr,
                                 pct.signals = pct.signals, pct.train = 1/2, pct.holdout = 1/2, 
                                 bias = bias, sim.type = sim.type, verbose = FALSE)
    dat <- rbind(sim.data$train, sim.data$holdout)
    predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
  } else { # optional: use provided data
    dat <- read.csv(pec_simFile)
    dat <- dat[,-1] # written file has first X column with subject names
    predictors.mat <- dat[, - which(colnames(dat) == class.lab)]
  }
  
  dat[, class.lab] <- as.factor(dat[, class.lab]) 
  pheno.class <- dat[, class.lab]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)
  
  outer_folds <- caret::createFolds(sim.data$train[,"class"], 5, list = FALSE)
  
  inner_folds <- caret::createFolds(sim.data$train[,"class"][outer_folds!=1], 5, list = TRUE)
  
  unlist(lapply(lapply(lapply(inner_folds, FUN=function(x){which(outer_folds != 1)[-x]}), FUN=function(x){table(sim.data$train[x,"class"])}), FUN=min))
  m.inner <- min(unlist(lapply(lapply(lapply(inner_folds, FUN=function(x){which(outer_folds != 1)[-x]}), FUN=function(x){table(sim.data$train[x,"class"])}), FUN=min)))
  m.inner <- m.inner - 1
  m.inner
  
  inner_idx <- which(outer_folds != 1)[-inner_folds[[1]]]
  
  knnSURF(m.inner, 0.5) # predicted multisurf k based on size of inner training folds
  floor((dim(sim.data$train[inner_idx, ])[1]-1)*0.154)
  
  d.mat <- as.matrix(dist(sim.data$train[,-ncol(sim.data$train)], upper=T, method="manhattan"))
  my.dist.moments <- dist_moments(d.mat)
  
  my.d <- density(d.mat[upper.tri(d.mat)], n=1000)
  #plot(my.d)
  
  myfn <- approxfun(x=my.d$x, y=my.d$y)
  prob.neighbor <- integrate(f=myfn, lower=min(my.d$x), upper=(mean(my.dist.moments$mu) - 0.7661796*mean(my.dist.moments$sigma)))
  
  prob.neighbor$value/2
  
  pnorm(-0.5)/2
  
  f <- function(k){
    A <- 2*exp(-x)-1
    B <- 2*exp(-y)-1
    sum((A*B)/(1+k*A*B))
  }
  
  set.seed(13)
  x <- runif(10)*10
  y <- runif(10)*5
  
  nleqslv(0,f)
  
  f <- function(alpha){
    
    val <- integrate(myfn, lower=min(my.d$x), upper=(mean(my.dist.moments$mu) - alpha*mean(my.dist.moments$sigma)))$value - 0.5*(1 - erf(0.5/sqrt(2)))
    
    val
  }
  
  myalpha <- nleqslv(0.7, f)
  
  myalpha$x
  knnSURF(m.inner, myalpha$x)
  
  prob.neighbor <- integrate(f=myfn, lower=min(my.d$x), upper=(mean(my.dist.moments$mu) - myalpha$x*mean(my.dist.moments$sigma)))
  
  alphas[iter] <- myalpha$x
  neighbor.probs[iter] <- prob.neighbor$value
  k.alphas[iter] <- knnSURF(m.inner, myalpha$x)
  
  my.rad <- mean(my.dist.moments$mu) - myalpha$x*mean(my.dist.moments$sigma)
  my.rad2 <- mean(my.dist.moments$mu) - 0.5*mean(my.dist.moments$sigma)
  
  par(mfrow=c(1,1), mar=c(4.1,4.5,1.3,0.8), xpd=F)
  plot(my.d, xlim=c(1.005*min(my.d$x), 0.99*mean(my.dist.moments$mu)),
       main="Estimated Distance Density")
  polygon(c(my.rad2, min(my.d$x), my.d$x[my.d$x <= my.rad2]), c(0, 0, my.d$y[my.d$x <= my.rad2]), col=rgb(1,0,0,0.3))
  abline(v=my.rad2, lty=1, lwd=2, col='red')
  
  #lines(my.d)
  polygon(c(my.rad, min(my.d$x), my.d$x[my.d$x <= my.rad]), c(0, 0, my.d$y[my.d$x <= my.rad]), col=rgb(0,0,1,0.6))
  abline(v=my.rad, lty=1, lwd=2, col='blue')
  legend('topleft', c("Adjusted alpha Prob", "alpha=0.5 Prob"), fill=c(rgb(0,0,1,0.6), rgb(1,0,0,0.3)), bg='white')
  legend('center', c('Adjusted alpha Radius','alpha=0.5 Radius'), lty=1, lwd=2, col=c('blue','red'), bg='white')

  
  # run consensus nested CV with predicted multisurf k
  cncv_result.predk <- consensus_nestedCV(train.ds = sim.data$train, 
                                  validation.ds =  sim.data$holdout, 
                                  label = sim.data$label,
                                  method.model = "classification",
                                  is.simulated = TRUE,
                                  ncv_folds = c(5, 5),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = "ReliefFequalK",
                                  relief.k.method = "k_half_sigma",                     # ReliefF knn
                                  wrapper = "relief",
                                  inner_selection_percent = 0.2,
                                  inner_selection_positivescores = T,
                                  tune.k = FALSE,
                                  tuneGrid = NULL,
                                  num_tree = 1000,
                                  verbose = T)


  # run consensus nested CV with naive k
  cncv_result.naivek <- consensus_nestedCV(train.ds = sim.data$train, 
                                  validation.ds =  sim.data$holdout, 
                                  label = sim.data$label,
                                  method.model = "classification",
                                  is.simulated = TRUE,
                                  ncv_folds = c(5, 5),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = "ReliefFequalK",
                                  #relief.k.method = knnSURF(m.inner, myalpha$x),             # ReliefF knn
                                  relief.k.method = 10,
                                  wrapper = "relief",
                                  inner_selection_percent = 0.2,
                                  inner_selection_positivescores = T,
                                  tune.k = FALSE,
                                  tuneGrid = NULL,
                                  num_tree = 1000,
                                  verbose = T)


  # naive k results
  cat("\n Naive K Nested Cross-Validation Accuracy [",cncv_result.naivek$cv.acc,"]\n")
  cat("\n Naive K Validation Accuracy [",cncv_result.naivek$Validation,"]\n")
  cat("\n Naive K Selected Features \n [",cncv_result.naivek$Features,"]\n")
  #cat("\n Elapsed Time [",cncv_result.naivek$Elapsed,"]\n")

  # predicted multisurf k results
  cat("\n Predicted MSURF K Nested Cross-Validation Accuracy [",cncv_result.predk$cv.acc,"]\n")
  cat("\n Predicted MSURF K Validation Accuracy [",cncv_result.predk$Validation,"]\n")
  cat("\n Predicted MSURF K Selected Features \n [",cncv_result.predk$Features,"]\n")
  #cat("\n Elapsed Time [",cncv_result.predk$Elapsed,"]\n")


  functional <- sim.data$signal.names

  cncv.predk.positives <- cncv_result.predk$Features
  TP <- sum(cncv.predk.positives %in% functional)
  FP <- sum((cncv.predk.positives %in% functional)==F)
  FN <- length(functional) - TP
  recall.predk <- TP/(TP+FN)
  cat("recall = ", recall.predk,"\n")
  precision.predk <- TP/(TP+FP)
  cat("precision = ", precision.predk,"\n")

  TPR.predk <- TP/length(cncv.predk.positives) #rate
  FPR.predk <- FP/length(cncv.predk.positives) #rate
  TPR.predk
  FPR.predk


  cncv.naivek.positives <- cncv_result.naivek$Features
  TP <- sum(cncv.naivek.positives %in% functional)
  FP <- sum((cncv.naivek.positives %in% functional)==F)
  FN <- length(functional) - TP
  recall.naivek <- TP/(TP+FN)
  cat("recall = ", recall.naivek,"\n")
  precision.naivek <- TP/(TP+FP)
  cat("precision = ", precision.naivek,"\n")

  TPR.naivek <- TP/length(cncv.naivek.positives) #rate
  FPR.naivek <- FP/length(cncv.naivek.positives) #rate
  TPR.naivek
  FPR.naivek

  cv.acc.predk[iter] <- cncv_result.predk$cv.acc
  cv.acc.naivek[iter] <- cncv_result.naivek$cv.acc
  validation.acc.predk[iter] <- cncv_result.predk$Validation
  validation.acc.naivek[iter] <- cncv_result.naivek$Validation

  attr.recall.predk[iter] <- recall.predk
  attr.recall.naivek[iter] <- recall.naivek
  attr.precision.predk[iter] <- precision.predk
  attr.precision.naivek[iter] <- precision.naivek

  attr.tpr.predk[iter] <- TPR.predk
  attr.tpr.naivek[iter] <- TPR.naivek
  attr.fpr.predk[iter] <- FPR.predk
  attr.fpr.naivek[iter] <- FPR.naivek
  
  cat("End Iteration: ", iter, "\n\n")

}

setwd("C:/Users/bdawk/Documents/KNN_project_output/R_code_and_output/code")

out.df <- data.frame(Iteration = rep(seq(1,30,by=1), length=60),
                     CV.Accuracy = c(cv.acc.predk, cv.acc.naivek),
                     Validation.Accuracy = c(validation.acc.predk, validation.acc.naivek),
                     Attribute.Recall = c(attr.recall.predk, attr.recall.naivek),
                     Attribute.Precision = c(attr.precision.predk, attr.precision.naivek),
                     Attribute.TPR = c(attr.tpr.predk, attr.tpr.naivek),
                     Attribute.FPR = c(attr.fpr.predk, attr.fpr.naivek),
                     Method=rep(c("PredK","NaiveK"), each=30))

CV.Accuracy <- c(cv.acc.predk, cv.acc.naivek)
Validation.Accuracy <- c(validation.acc.predk, validation.acc.naivek)
Attribute.Recall <- c(attr.recall.predk, attr.recall.naivek)
Attribute.Precision <- c(attr.precision.predk, attr.precision.naivek)
Attribute.TPR <- c(attr.tpr.predk, attr.tpr.naivek)
Attribute.FPR <- c(attr.fpr.predk, attr.fpr.naivek)
Method <- rep(c("PredK","NaiveK"), each=30)

out.df2 <- data.frame(Value=c(CV.Accuracy,
                               Validation.Accuracy,
                               Attribute.Recall,
                               Attribute.Precision),
                      Metric=c(rep('Training Accuracy',length=length(CV.Accuracy)),
                               rep('Validation Accuracy',length=length(Validation.Accuracy)),
                               rep('Attribute Recall',length=length(Attribute.Recall)),
                               rep('Attribute Precision',length=length(Attribute.Precision))),
                      Method=rep(Method,length=4*60))
out.df2[,2] <- as.factor(out.df2[,2])
out.df2[,3] <- as.factor(out.df2[,3])

ggplot(out.df2, aes(x=Metric, y=Value, fill=Method)) +
  geom_boxplot(outlier.shape=NA) +
  #geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
  #            aes(x=Metric, y=Value, color=Method), inherit.aes=F)+
  geom_point(shape=16, position=position_jitterdodge(),alpha=0.4, size=4.5,
             aes(x=Metric, y=Value, color=Method))+
  scale_y_continuous(name="Cross Validation Classification Accuracy", limits=c(0,1),breaks=seq(0,1,by=.25), 
                     labels=c("0%","25%","50%","75%","100%"))+
  #scale_x_discrete(breaks=c("PredK","NaiveK"),
  #                 labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values=c("magenta", "#0072B2"),
                     labels=c("Informed K", "Naive K"))+
  ggtitle("Performance comparison between distance-informed and distance-naive K")+
  scale_fill_manual(breaks=c("PredK","NaiveK"),
                    values = c("white", "white"))+
  guides(fill=FALSE)+
  theme_bw()+
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_blank(),
        plot.title=element_text(size=18,face='bold',hjust=0.5),
        legend.title=element_text(size=16,face='bold'),
        legend.text=element_text(size=16),
        legend.title.align=0.5,
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
ggsave("interaction_sim_m100-p1000_bias-pt4_pct-signals-pt1_predk_vs_naivek(grouped).pdf", height=7, width=11)


#write.csv(out.df, "consensus-nestedCV_iterates_mainEffect-m200_p1000_pct-signals-pt1_bias-pt6(naiveKmax).csv",row.names=F)
out.df[,1] <- as.factor(out.df[,1])
out.df[,2] <- as.numeric(out.df[,2])
out.df[,3] <- as.numeric(out.df[,3])
out.df[,4] <- as.numeric(out.df[,4])
out.df[,5] <- as.numeric(out.df[,5])
out.df[,6] <- as.numeric(out.df[,6])
out.df[,7] <- as.numeric(out.df[,7])
out.df[,8] <- as.factor(out.df[,8])

library(ggplot2)

# CV Accuracy
p1 <- ggplot(out.df, aes(x=Method, y=CV.Accuracy)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Cross Validation Classification Accuracy", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=17,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p1
#ggsave("cv-accuracy_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)


# Validation Accuracy
p2 <- ggplot(out.df, aes(x=Method, y=Validation.Accuracy)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Validation Data Classification Accuracy", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=17,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p2
#ggsave("validation-accuracy_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)

# Attribute Recall
p3 <- ggplot(out.df, aes(x=Method, y=Attribute.Recall)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Attribute Recall TP/(TP+FN)", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=17,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p3
#ggsave("attribute-recall_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)

# Attribute Precision
p4 <- ggplot(out.df, aes(x=Method, y=Attribute.Precision)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Attribute Precision TP/(TP+FP)", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=17,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p4
#ggsave("attribute-precision_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)

# Attribute TPR
ggplot(out.df, aes(x=Method, y=Attribute.TPR)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Attribute True Positive Rate (TP/Detected)", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=16,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
#ggsave("attribute-TPR_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)

# Attribute FPR
ggplot(out.df, aes(x=Method, y=Attribute.FPR)) +
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=length(table(out.df[,"Method"])))) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=Method), inherit.aes=T)+
  scale_y_continuous(name="Attribute False Positive Rate (FP/Detected)", limits=c(0,1),breaks=seq(0,1,by=.25))+
  scale_x_discrete(breaks=c("PredK","NaiveK"),
                   labels=c("Informed K","Naive K"))+
  scale_color_manual(breaks=c("PredK","NaiveK"),
                     values = c("magenta", "#0072B2"),
                     labels=c("Informed K","Naive K"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=16,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
#ggsave("attribute-FPR_cncv_naiveKmax-informedK(mainEffect_bias-pt6).pdf", height=5, width=6)

setwd("C:/Users/bdawk/Documents/KNN_project_output")
comb_plots <- cowplot::plot_grid(p1, p2, p3, p4, labels = 'AUTO',
                                 align = 'h', label_size=22, rel_widths=c(1,1), 
                                 rel_heights=c(1,1))
comb_plots

ggsave(filename="interaction_sim_m100-p1000_bias-pt4_pct-signals-pt1_predk_vs_naivek.pdf",
       height=12, width=12.1)
