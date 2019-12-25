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

library(privateEC)
library(ggplot2)
library(stir)

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
  num.positives <- length(rncv.positives)
  TPR <- TP/num.positives #rate
  FPR <- FP/num.positives #rate
  return(list(TP=TP, FP=FP, FN=FN, TPR=TPR, FPR=FPR, precision=precision, recall=recall))
}

chosen.k <- numeric()
atts.list <- list()
for(all.iter in 1:50){
    print(all.iter)

# create simulated data set

n.samples <- 100     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class"
type <- "mainEffect" # main effect simulations
bias <- 0.9          # moderate effect size
pct.signals <- 0.1   # pct functional features
verbose <- FALSE
data.sets <- createSimulation(num.samples = n.samples,
                              num.variables = n.variables,
                              pct.signals = pct.signals,
                              label = label,
                              bias = bias,
                              pct.train = 1/3,
                              pct.holdout = 1/3,
                              pct.validation = 1/3,
                              sim.type = type,
                              save.file = NULL,
                              verbose = verbose)

# number of trees in rf-models
num_tree <- 500

# attributes that are functionally related to the phenotype
functional <- data.sets$signal.names

# number of folds in outer and inner loops
ncv_folds <- c(10,10)

# type of response
label <- "class"

# joing train and holdout into one training set
train.data <- rbind(data.sets$train,data.sets$holdout)
train.ds <- train.data
validation.ds <- data.sets$validation

predictors.mat <- apply(train.data[,-which(colnames(train.data) == label)],2,as.numeric)
#get.distance(predictors.mat,metric="manhattan")

# create outer folds for nested CV
outer_folds <- createFolds(train.ds[,label],ncv_folds[1],list=FALSE)

# number of samples
m <- n.samples

# values of k to consider inside inner folds
ks <- seq(1,floor((m-1)/4),by=1)
#ks <- seq(1,10,by=1)
# numeric vector for test accuracy in outer folds for prevailing k
# from the inner loop

ks.test.accu <- numeric()

relief_atts <- list()
best.atts <- list()
out.accu <- numeric()

# nested cv loops

for(out.iter in 1:ncv_folds[1]){
  
    cat("Outer Fold: ",out.iter,"\n")
    
    RF.method <- "relieff"
    
    metric <- "manhattan"
    
    outer_idx <- which(outer_folds != out.iter)
    
    trn.data <- as.matrix(train.ds[outer_idx,])
    tst.data <- as.matrix(train.ds[-outer_idx,])
    
    predictors.mat <- apply(trn.data[,-which(colnames(trn.data)==label)],2,as.numeric)
    #print(class(predictors.mat))
    
    trn.pheno <- as.factor(train.ds[,label][outer_idx])
    tst.pheno <- as.factor(train.ds[,label][-outer_idx])
    
    rf.model <- randomForest(predictors.mat, y=trn.pheno, mtry=max(floor(ncol(trn.data)/3),1), ntree=num_tree,
                             importance=TRUE)
    var.importance <- importance(rf.model,type=1)
    #print(importance(rf.model))
    relief_atts[[out.iter]] <- row.names(var.importance)[which(var.importance > -0.5)]
    
    inner_folds <- createFolds(trn.data[,label], ncv_folds[2], list=FALSE)
    
    #############################################################################################
    
    ## start inner folds loop for each value of k
    #
    #############################################################################################
    ks.test.tmp <- numeric()
    atts.tmp <- list()
    for(k.iter in 1:length(ks)){
        #cat("k = ",k.iter,"\n")
        k <- ks[k.iter]
        predictors.mat <- predictors.mat[,relief_atts[[out.iter]]]
        neighbor.idx.observed <- find.neighbors(predictors.mat, trn.pheno, k=k, method=RF.method)
        #print(neighbor.idx.observed)
        results.list <- stir(predictors.mat, neighbor.idx.observed, k=k, metric = metric, method=RF.method)
        t_sorted_stir <- results.list$`STIR-t`[,-3]  # remove cohen-d
        t_sorted_stir <- results.list$`STIR-t`[,-3]  # remove cohen-d
        colnames(t_sorted_stir) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
        t_sorted_stir$attribute <- rownames(t_sorted_stir) # adds a column for merge
        n.positives <- length(t_sorted_stir$attribute[t_sorted_stir$t.pval.adj.stir < .5])
        if(n.positives > 0){
           stir.positives <- t_sorted_stir$attribute[t_sorted_stir$t.pval.adj.stir < .5]
        }else{
           stir.positives <- t_sorted_stir$attribute[1:(pct.signals*100+20)]
        }
        stir.positives <- t_sorted_stir$attribute[1:30]
        atts.tmp[[k.iter]] <- stir.positives
        #print(stir.positives)
        
        ks.test.inner <- numeric()
        for(in.iter in 1:ncv_folds[2]){
          inner_idx <- which(inner_folds != in.iter)
          trn.data.inner <- as.matrix(trn.data[inner_idx,c(stir.positives)])
          tst.data.inner <- as.matrix(trn.data[-inner_idx,c(stir.positives)])
          
          trn.pheno.inner <- as.factor(trn.data[,label][inner_idx])
          tst.pheno.inner <- as.factor(trn.data[,label][-inner_idx])
          
          rf.model.inner <- randomForest(trn.data.inner, y=trn.pheno.inner, 
                                         mtry=max(floor(ncol(trn.data.inner)/3),1), ntree=num_tree)
          pred.class <- predict(rf.model.inner,tst.data.inner)
          ks.test.inner[in.iter] <- 1 - mean(pred.class != tst.pheno.inner)
            
        }
        ks.test.tmp[k.iter] <- mean(ks.test.inner) # inner test set CV error for current value of k
    }
    plot(ks,ks.test.tmp,type='l',lwd=2,main=paste("Outer Fold: ",out.iter,sep=""),
         xlab="k",ylab="Test Accuracy")
    points(ks,ks.test.tmp,col='red',pch=16)
    
    max.idx <- which(ks.test.tmp == max(ks.test.tmp))
    #ks.test.accu[out.iter] <- which.max(ks.test.tmp) # value of k giving highest classification accuracy
    ks.test.accu[out.iter] <- max.idx[length(max.idx)]
    
    max.k <- ks.test.accu[out.iter]
    points(max.k,ks.test.tmp[max.k],pch=16,col='blue')
    
    best.atts[[out.iter]] <- atts.tmp[[max.k]]
    
    train.data <- as.matrix(train.ds[,best.atts[[out.iter]]])
    test.data <- as.matrix(validation.ds[,best.atts[[out.iter]]])
    
    train.pheno <- as.integer(train.ds[,label]) - 1
    test.pheno <- as.integer(validation.ds[,label]) - 1
    
    rf.model <- randomForest(train.data, y=as.factor(train.pheno), mtry=max(floor(ncol(train.data)/3),1),
                             ntree=num_tree)
    Train_accu <- 1 - mean(rf.model$confusion[,"class.error"])
    test.pred <- predict(rf.model, newdata=test.data)
    Test_accu <- 1 - mean(test.pred != test.pheno)
    out.accu[out.iter] <- Test_accu
}
best.fold <- which.max(out.accu)
best.k <- ks.test.accu[best.fold]
best.fold.accu <- out.accu[best.fold]
chosen.k[all.iter] <- best.k
best.fold.accu
best.atts[[best.fold]]

#for(i in 1:ncv_folds[1]){
#    cat("k = ",ks.test.accu[i],"\n")
  
#    cat("\n","Features: ","\n",best.atts[[i]],"\n")
    
#    cat("\n","Validation Accuracy: ",out.accu[i],"\n\n")
#}

functional <- data.sets$signal.names
rncv.positives <- best.atts[[best.fold]]
rncv.detect <- attrDetectStats(functional, rncv.positives)
cat("recall = ", rncv.detect$recall,"\n")
#cat("The recall is high (high TP) for regular nested: it gets most of the functional vars because it selects a lot of vars.\n") 
cat("precision = ", rncv.detect$precision,"\n")
#cat("But the precision is low (high FP) for regular nested: it selects so many vars that many are false.\n") 

#rncv.detect$TPR 
#rncv.detect$FPR 
}
round(mean(chosen.k))
consensus.atts <- unique(unlist(best.atts))
consensus.atts
rncv.positives <- consensus.atts
rncv.detect <- attrDetectStats(functional, rncv.positives)
cat("recall = ", rncv.detect$recall,"\n")
cat("precision = ", rncv.detect$precision,"\n")

