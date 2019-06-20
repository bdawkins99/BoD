# nested cross-validation for tuning k in fixed-k STIR

# In inner folds, loop through k = 1, 2, ... , floor((m-1)/2) 
# to choose the value(s) of k that maximize test classification 
# accuracy, precision, and recall. 
#
# In outer folds, fit a rf-model with tuned k and chosen features. 
# Compute outer fold test classification accuracy, precision, and recall. 
#
# Record the value(s) of k for which the best outer fold test 
# classification accuracy, precision, and recall are maximized.

setwd("C:/Users/bdawk/Documents/KNN_project_output")
library(devtools)
#install_github("insilico/glmSTIR") 
library(privateEC)
library(ggplot2)
library(stir)
library(glmSTIR)
library(caret)
library(randomForest)

# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)  # helper function from STIR

attrDetectStats <- function(functional, positives){
  # Given a vector functional (true) attribute associations and a vector of positive
  # associations, returns detection statistics like true positive rate, recall and precision.
  # usage:
  # functional <- data.sets$signal.names
  # rncv.positives <- rncv_result$Features
  # rncv.detect <- attrDetectStats(functional, rncv.positives)
  TP <- sum(positives %in% functional)
  FP <- sum((positives %in% functional)==F)
  FN <- length(functional) - TP
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  num.positives <- length(positives)
  TPR <- TP/num.positives #rate
  FPR <- FP/num.positives #rate
  return(list(TP=TP, FP=FP, FN=FN, TPR=TPR, FPR=FPR, precision=precision, recall=recall))
}

mixed.ratio <- seq(0.1,0.9,by=0.1)
#mixed.ratio <- seq(0.3,0.9,by=0.1)
#mixed.ratio <- c(0.3)
#mixed.ratio <- c(0.1)
my.chosen.k <- numeric()
my.precision.list <- list()
my.recall.list <- list()
all.k.list <- list()
for(mixed.iter in 1:length(mixed.ratio)){
  cat("Mixed Percent = ",mixed.ratio[mixed.iter],"\n\n")
chosen.k <- numeric()
num.iter <- 100
k.accu.list <- list()
for(iter in 1:num.iter){
if((iter/10) == floor(iter/10)){
  cat("Mixed Percent = ",mixed.ratio[mixed.iter],"\n\n")
}
cat("Iteration: ",iter,"\n")
# create simulated data set

n.samples <- 100     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class"
type <- "mainEffect" # main effect simulations
#type <- "interactionErdos"
#type <- "interactionScalefree"
bias <- 0.4          # moderate effect size
pct.signals <- 0.1   # pct functional features
verbose <- FALSE
#data.sets <- createSimulation(num.samples = n.samples,
#                              num.variables = n.variables,
#                              pct.signals = pct.signals,
#                              label = label,
#                              bias = bias,
#                              pct.train = 1/3,
#                              pct.holdout = 1/3,
#                              pct.validation = 1/3,
#                              sim.type = type,
#                              save.file = NULL,
#                              verbose = verbose)

pct.mixed <- mixed.ratio[mixed.iter]
pct.imbalance <- 0.5
mixed.type <- c("mainEffect", "interactionScalefree")
data.sets <- createMixedSimulation(num.samples=n.samples,
                                   num.variables=n.variables,
                                   pct.signals=pct.signals,
                                   pct.imbalance=pct.imbalance,
                                   label = label,
                                   bias=bias,
                                   pct.train = 1/3,
                                   pct.holdout = 1/3,
                                   pct.validation = 1/3,
                                   pct.mixed=pct.mixed,
                                   mixed.type=mixed.type,
                                   save.file = NULL,
                                   verbose=FALSE)

# number of trees in rf-models
num_tree <- 500

# attributes that are functionally related to the phenotype
functional <- data.sets$signal.names

# number of folds in outer and inner loops
ncv_folds <- c(10,10)

# type of response
label <- "class"

# joing train and holdout into one training set
train.data <- rbind(data.sets$train,data.sets$holdout) # full training data set
train.ds <- train.data                                 # full training data set
validation.ds <- data.sets$validation                  # validation set to be predicted after tuning k and selecting features
predictors.mat <- apply(train.ds[,-which(colnames(train.data) == label)],2,as.numeric)

# number of samples
m <- n.samples

# values of k to consider inside inner folds
ks <- seq(1,floor((m-1)/4),by=1)

# relief method
RF.method <- "relieff"

# metric to use in distance calculations
metric <- "manhattan"

# create outer folds for nested CV
outer_folds <- createFolds(train.ds[,label],ncv_folds[1],list=FALSE)

relief_atts <- list() # list that will contain chosen features in each outer fold
best_k <- NULL        # will be assigned best value of k based on outer fold test classification accuracy
Train_accu <- NULL    # will be the full training data classification accuracy after k is chosen and an rf model is fit with features from nested loop
Test_accu <- NULL     # will be the validation set classification accuracy after k is chosen and an rf model is fit to with features from nested loop
tune_params <- NULL
k.outer <- numeric()
best_atts <- NULL
fold_atts_list <- list()

for(i in 1:ncv_folds[1]){
    cat("----------------------------------------------------------------------","\n")
    cat("Outer Fold: ",i,"\n")
    cat("----------------------------------------------------------------------","\n")
    
    outer_idx <- which(outer_folds != i)
    
    outer.trn <- predictors.mat[outer_idx,]
    outer.pheno <- as.factor(train.ds[,label][outer_idx])
    
    outer.rf.model <- randomForest(outer.trn, y=outer.pheno, mtry=max(floor(ncol(outer.trn)/3),1), ntree=num_tree,
                             importance=TRUE)
    var.importance <- importance(outer.rf.model,type=1)
    var.tmp <- as.data.frame(var.importance[order(var.importance,decreasing=T),])
    #relief_atts[[i]] <- row.names(var.importance)[which(var.importance > -1)]
    #relief_atts[[i]] <- row.names(var.tmp)[1:(nrow(var.tmp)-pct.signals*100)]
    relief_atts[[i]] <- row.names(var.tmp)

    trn.data.tmp <- as.matrix(train.ds[outer_idx,])
    trn.data <- as.matrix(train.ds[outer_idx, relief_atts[[i]]])
    tst.data <- as.matrix(train.ds[-outer_idx, relief_atts[[i]]])
    
    trn.pheno <- as.factor(train.ds[,label][outer_idx])
    tst.pheno <- as.factor(train.ds[,label][-outer_idx])
    
    inner_folds <- createFolds(trn.data.tmp[,label], ncv_folds[2], list=FALSE)
    
    k.accu.vec <- numeric()
    for(k.iter in 1:length(ks)){
         k <- ks[k.iter]
         
         #neighbor.idx.observed <- find.neighbors(trn.data, trn.pheno, k=k, method=RF.method)
         #results.list <- stir(trn.data, neighbor.idx.observed, k=k, metric = metric, method=RF.method)
         glmstir.mdd.rnaseq.results <- glmSTIR(trn.pheno, trn.data, regression.type="glm", attr.diff.type="numeric-abs", 
                                               nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                               fdr.method="bonferroni",k=k)
         #t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
         #t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
         #colnames(t_sorted_stir) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
         #t_sorted_stir$attribute <- rownames(t_sorted_stir) # adds a column for merge

         #stir.positives <- t_sorted_stir$attribute[1:30]
         stir.positives <- row.names(glmstir.mdd.rnaseq.results)[1:30]  # top 30, p-value sorted
       
         k.accu.vec.tmp <- numeric() # store classification accuracy for each of the inner folds
         for(j in 1:ncv_folds[2]){
             #cat("  Inner Fold: ",j,"\n")
             
             inner_idx <- which(inner_folds != j)
             
             trn.data.tmp.inner <- as.matrix(trn.data.tmp[inner_idx,])
             trn.data.inner <- as.matrix(trn.data[inner_idx, stir.positives])
             tst.data.inner <- as.matrix(trn.data[-inner_idx, stir.positives])
             
             trn.pheno.inner <- as.factor(trn.data.tmp[,label][inner_idx])
             tst.pheno.inner <- as.factor(trn.data.tmp[,label][-inner_idx])
             
             rf.model.inner <- randomForest(trn.data.inner, y=trn.pheno.inner, 
                                            mtry=max(floor(ncol(trn.data.inner)/3),1), ntree=num_tree)
             pred.class <- predict(rf.model.inner,tst.data.inner)
             k.accu.vec.tmp[j] <- 1 - mean(pred.class != tst.pheno.inner)
             
         }
         k.accu.vec[k.iter] <- mean(k.accu.vec.tmp) # test set CV classification accuracy for current k
         
    }
    k.accu.list[[iter]] <- k.accu.vec

    max.idx <- which(k.accu.vec == max(k.accu.vec))
    #k.outer[i] <- which.max(k.accu.vec) # value of k giving the highest classification accuracy from inner loop CV
    #k.outer[i] <- max(max.idx)
    #k.outer[i] <- min(max.idx)
    #k.outer[i] <- round(mean(max.idx))
    if(length(max.idx) > 1){
       k.outer[i] <- sample(max.idx,size=1,replace=F)
    }else{
       k.outer[i] <- max.idx
    }
    
    plot(ks,k.accu.vec,type='l',lwd=2,main=paste("Outer Fold: ",i,sep=""),
         xlab="k",ylab="Test Accuracy")
    points(ks,k.accu.vec,col='red',pch=16)
    max.k <- k.outer[i]
    points(max.k,k.accu.vec[max.k],pch=16,col='blue')
    
    
    #neighbor.idx.observed <- find.neighbors(trn.data, trn.pheno, k=k.outer[i], method=RF.method)
    #results.list <- stir(trn.data, neighbor.idx.observed, k=k.outer[i], metric = metric, method=RF.method)
    glmstir.mdd.rnaseq.results <- glmSTIR(trn.pheno, trn.data, regression.type="glm", attr.diff.type="numeric-abs", 
                                          nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                          fdr.method="bonferroni",k=k)
    #t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
    #t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
    #colnames(t_sorted_stir) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
    #t_sorted_stir$attribute <- rownames(t_sorted_stir) # adds a column for merge
    
    #stir.positives <- t_sorted_stir$attribute[1:30]
    stir.positives <- row.names(glmstir.mdd.rnaseq.results)[1:30]
    
    #n.positives <- length(t_sorted_stir$attribute[t_sorted_stir$t.pval.adj.stir < 1.1])
    #print(n.positives)
    #print(t_sorted_stir$attribute)
    #if(n.positives > 0){
    #  stir.positives <- t_sorted_stir$attribute[t_sorted_stir$t.pval.adj.stir < .05]
    #}else{
    #  stir.positives <- t_sorted_stir$attribute[1:30]
    #}
    fold_atts_list[[i]] <- stir.positives
    
    #print(ncol(trn.data[,stir.positives]))
    rf.model.tmp <- randomForest(trn.data[,stir.positives], y=trn.pheno, mtry=max(floor(ncol(trn.data[,stir.positives])/3),1),ntree=num_tree)
    pred.class <- predict(rf.model.tmp,tst.data[,stir.positives])
    accu <- 1 - mean(pred.class != tst.pheno)
    tune_params <- rbind(tune_params, data.frame(k.outer[i],accu))
}

max.idx <- which(tune_params$accu == max(tune_params$accu))
best_k <- max(tune_params[max.idx,1])
cat("Chosen k = ",best_k,"\n")

chosen.k[iter] <- best_k
###############################################################################################

len.list <- length(k.accu.list)
rep.ks <- rep(ks,length=(length(ks)*len.list))
if(iter == num.iter){
  mod <- lm(c(unlist(k.accu.list)) ~ rep.ks)
  s <- summary(mod)
  r.squ <- round(s$r.squared,3)
  F.stat <- round(c(s$fstatistic[1]),3)
  pval <- round(s$coefficients[,4][2],9)
  pdf(paste("glmSTIR_pct_mixed",mixed.ratio[mixed.iter],"_test_accu-vs-k.pdf",sep=""),height=7,width=7)
  plot(rep.ks,c(unlist(k.accu.list)),type='p',
       main=paste("All Iterations pct.mixed = ",mixed.ratio[mixed.iter],sep=""),
       xlab="k",ylab="Test Accuracy",ylim=c(0,1),font.lab=2)
  abline(mod,col='red',lwd=2)
  abline(v=round(mean(chosen.k)),lty=2,lwd=2)
  #points(rep.ks,c(unlist(k.accu.list)),col='red',pch=16)
  leg1 <- paste("Rsq. = ",r.squ,sep="")
  leg2 <- paste("Fstat = ",F.stat,sep="")
  leg3 <- paste("Pval = ",pval,sep="")
  legend("bottomright",
         c("LS-Fit",paste("Best k = ",round(mean(chosen.k)),sep=""),leg1,leg2,leg3),lwd=c(2,2,0,0,0),
         lty=c(1,2,0,0,0),col=c('red','black','white','white','white'),bg='white')
  dev.off()
}
################################################################################################

best_atts <- fold_atts_list[[max.idx[length(max.idx)]]]
best_atts

train.data <- as.matrix(train.ds[,best_atts])
test.data <- as.matrix(validation.ds[,best_atts])

train.pheno <- as.integer(train.ds[,label]) - 1
test.pheno <- as.integer(validation.ds[,label]) - 1

rf.model <- randomForest(train.data, y = as.factor(train.pheno), mtry=max(floor(ncol(train.data)/3),1),ntree=num_tree)

Train_accu <- 1 - mean(rf.model$confusion[,"class.error"])
cat("Train Accuracy = ",Train_accu,"\n")

test.pred <- predict(rf.model, newdata=test.data)

Test_accu <- 1 - mean(test.pred != test.pheno)
cat("Test Accuracy = ",Test_accu,"\n")

#neighbor.idx.observed <- find.neighbors(train.data, train.pheno, k=best_k, method=RF.method)
#results.list <- stir(train.data, neighbor.idx.observed, k=best_k, metric = metric, method=RF.method)
glmstir.mdd.rnaseq.results <- glmSTIR(train.pheno, train.data, regression.type="glm", attr.diff.type="numeric-abs", 
                                      nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=0.5,
                                      fdr.method="bonferroni",k=best_k)
#t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
#t_sorted_stir <- results.list$STIR_T[,-3]  # remove cohen-d
#colnames(t_sorted_stir) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
#t_sorted_stir$attribute <- rownames(t_sorted_stir) # adds a column for merge
#stir.positives <- t_sorted_stir$attribute[t_sorted_stir$t.pval.adj.stir < 5]
stir.positives <- row.names(glmstir.mdd.rnaseq.results)[1:30]

cat("\n","Stir Positives: ","\n",stir.positives,"\n")

positive.stats <- attrDetectStats(functional,stir.positives)
cat("\n","Precision = ",positive.stats$precision,"\n")
cat("\n","Fraction of Maximum Precision = ",(positive.stats$precision/(1/3)),"\n")
cat("\n","Recall = ",positive.stats$recall,"\n")

my.precision.list[[iter]] <- positive.stats$precision
my.recall.list[[iter]] <- positive.stats$recall
}
my.chosen.k[mixed.iter] <- round(mean(chosen.k))
all.k.list[[mixed.iter]] <- chosen.k
write.csv(data.frame(chosen.k=chosen.k),file=paste("glmSTIR_chosen_k_pct-mixed_",mixed.ratio[mixed.iter],".csv",sep=""),row.names=F)
write.csv(data.frame(precision=unlist(my.precision.list),recall=unlist(my.recall.list)),file=paste("glmSTIR_recall_precision_pct-mixed_",mixed.ratio[mixed.iter],".csv",sep=""),row.names=F)
pdf(paste("glmSTIR_chosen_k_hist_pct-mixed",mixed.ratio[mixed.iter],".pdf",sep=""),height=7,width=7)
hist(chosen.k,breaks=10,freq=F,
     main=paste("Histogram of Chosen k from Nested CV:","Mix Perc. = ",mixed.ratio[mixed.iter],sep=""),
     xlab="k")
dev.off()
}

out.df <- data.frame(pct.mixed=mixed.ratio,best.k=my.chosen.k)
write.csv(out.df,file="glmSTIR_nestedCV-k_mixedEffects_summary.csv",row.names=F)
round(mean(chosen.k))
pdf("glmSTIR_hist-chosen_k-mixedEffects.pdf",height=6,width=6)
hist(chosen.k,breaks=10,freq=F)
dev.off()

m <- 100
sig.frac <- 0.55
library(pracma)
floor(((m-1)/2)*((1 - erf(sig.frac/sqrt(2)))*0.5))
