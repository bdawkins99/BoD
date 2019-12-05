setwd("C:/Users/bdawk/Documents/KNN_project_output")
library(devtools)
#install_github("insilico/npdr") 
library(privateEC)
library(ggplot2)
library(npdr)
library(caret)
library(randomForest)

# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
#check.packages(packages)  # helper function from STIR

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

nestedCV_best_k <- function(train.ds=NULL,
                            validation.ds=NULL,
                            ncv_folds=c(10,10),
                            num_tree=500,
                            verbose=FALSE,
                            label="class",
                            attr.diff.type="numeric-abs",
                            corr.attr.names=NULL,
                            dopar.nn=FALSE, dopar.reg=FALSE,
                            k.grid=NULL){
  
  if (sum(ncv_folds)>(dim(train.ds)[1])/3){
    stop("There are less than three observations in each fold")
  }
  
  m <- dim(train.ds)[1]
  if(attr.diff.type!="correlation-data"){
    p <- dim(train.ds[,-ncol(train.ds)])[2]
  }else{
    p <- dim(train.ds[,-1])
  }
  
  if(attr.diff.type=="correlation-data"){   # corrdata
    mynum <- dim(train.ds[,-1])[2]               # corrdata
    for(i in seq(1,mynum-1,by=1)){          # corrdata
      mydiv <- i                            # corrdata
      if((mydiv*(mydiv - 1)) == mynum){     # corrdata
        my.dimension <- mydiv               # corrdata
        break                               # corrdata
      }                                     # corrdata
    }                                       # corrdata
    num.attr <- my.dimension                # corrdata
  }else{
    num.attr <- ncol(train.ds[,-ncol(train.ds)])
  }
  num.samp <- nrow(train.ds)
  
  if(attr.diff.type=="correlation-data"){    # corrdata
    attr.idx.list <- list()                  # corrdata
    for(i in 1:num.attr){                    # corrdata
      lo.idx <- (i - 1)*(num.attr-1) + 1     # corrdata
      hi.idx <- i*(num.attr-1)               # corrdata
      attr.idx.list[[i]] <- c(lo.idx:hi.idx) # corrdata
    }                                        # corrdata
  }
  
  if(is.null(k.grid)){
    #ks <- seq(1,floor((m-1)/2),by=1)
    ks <- seq(1,(m-1),by=1)
  }else{
    ks <- k.grid
  }
  
  k.best.vec <- numeric()
  att.best.list <- list()
  pred.accu.vec <- numeric()
  
  outer_folds <- caret::createFolds(train.ds[,label], ncv_folds[1], list=FALSE)
  for(i in 1:ncv_folds[1]){
    cat("Outer Fold: ",i,"\n")
    outer.idx <- which(outer_folds != i)
    trn.data <- as.matrix(train.ds[outer.idx,-which(colnames(train.ds)==label)])
    tst.data <- as.matrix(train.ds[-outer.idx,-which(colnames(train.ds)==label)])
    trn.pheno <- as.factor(train.ds[,label][outer.idx])
    tst.pheno <- as.factor(train.ds[,label][-outer.idx])
    
    #print(dim(trn.data))
    
    inner_folds <- caret::createFolds(trn.pheno, ncv_folds[2], list=FALSE)
    
    k.inner.best.vec <- numeric()
    k.inner.best.accu <- numeric()
    att.best.inner.list <- list()
    for(j in 1:ncv_folds[2]){
      cat("      Inner Fold: ",j,"\n")
      inner.idx <- which(inner_folds != j)
      trn.data.inner <- trn.data[inner.idx,]
      tst.data.inner <- trn.data[-inner.idx,]
      trn.pheno.inner <- as.factor(trn.pheno[inner.idx])
      tst.pheno.inner <- as.factor(trn.pheno[-inner.idx])
      
      n.k <- length(trn.pheno.inner)
      #print(n.k)
      
      k.accu.vec <- numeric()
      att.best.tmp <- list()
      for(k.iter in 1:(n.k-1)){
        cat("      k = ",k.iter,"\n")
        k <- ks[k.iter]
        
        npdr.cc.results <- npdr(trn.pheno.inner, trn.data.inner, regression.type="binomial", attr.diff.type=attr.diff.type,
                                nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, knn=k,
                                neighbor.sampling="none",
                                padj.method="bonferroni", 
                                corr.attr.names=corr.attr.names,
                                verbose=verbose, dopar.nn=dopar.nn, dopar.reg=dopar.reg)
        
        npdr.positives <- as.character(npdr.cc.results$att)[1:30]
        
        if(attr.diff.type=="correlation-data"){
          idx.positive <- numeric()
          new.attr.idx.list <- list()
          for(pos in 1:length(npdr.positives)){
            idx.positive[pos] <- as.numeric(strsplit(npdr.positives[pos],split="ROI")[[1]][2])
            new.attr.idx.list[[pos]] <- attr.idx.list[[idx.positive[pos]]]
          }
        
          col.idx.npdr.positives <- unlist(new.attr.idx.list)
        
        
          att.best.tmp[[k.iter]] <- npdr.positives

          rf.model.inner <- randomForest::randomForest(trn.data.inner[,col.idx.npdr.positives], y=trn.pheno.inner, 
                                       mtry=max(floor(ncol(trn.data.inner[,col.idx.npdr.positives])/3),1), ntree=num_tree,
                                       na.action=na.omit)
          pred.class <- stats::predict(rf.model.inner,tst.data.inner[,col.idx.npdr.positives], na.action=na.omit)
        }else{
          att.best.tmp[[k.iter]] <- npdr.positives
          
          rf.model.inner <- randomForest::randomForest(trn.data.inner[,npdr.positives], y=trn.pheno.inner, 
                                                       mtry=max(floor(ncol(trn.data.inner[,npdr.positives])/3),1), ntree=num_tree,
                                                       na.action=na.omit)
          pred.class <- stats::predict(rf.model.inner,tst.data.inner[,npdr.positives], na.action=na.omit)
        }

        k.accu.vec[k.iter] <- 1 - mean(pred.class != tst.pheno.inner)
        #cat("Test Accu: ",k.accu.vec[k.iter],"\n")
        
      }
      
      max.accu <- max(k.accu.vec)
      cat("            Test Accuracy: ",max.accu,"\n")
      
      idx.max <- which(k.accu.vec==max.accu)
      #idx.max <- ifelse(length(idx.max)==1,idx.max,sample(idx.max,size=1,replace=F))
      if(length(idx.max)==1){
        idx.max <- idx.max
      }else{
        idx.max <- sample(idx.max,size=1,replace=F)
      }
      cat("                     k = ",idx.max,"\n")
      
      par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
      plot(c(1:length(k.accu.vec)),k.accu.vec,type='l',lwd=2,col='black')
      points(c(1:length(k.accu.vec)),k.accu.vec,pch=19,col='red',cex=1.3)
      points(idx.max,k.accu.vec[idx.max],pch=19,col='blue',cex=1.5)
      
      #cat("Inner Best k's: ",idx.max,"\n")
      k.inner.best.vec[j] <- idx.max
      k.inner.best.accu[j] <- k.accu.vec[which.max(k.accu.vec)]
      att.best.inner.list[[j]] <- att.best.tmp[[which.max(k.accu.vec)]]
      idx.att <- c(grep(att.best.inner.list[[j]],pattern="main",value=TRUE),
                   grep(att.best.inner.list[[j]],pattern="int",value=TRUE))
      cat("            Functional: ",idx.att,"\n")
    }
    max.accu <- max(k.inner.best.accu)
    idx.max <- which(k.inner.best.accu==max.accu)
    #idx.max <- ifelse(length(idx.max)==1,idx.max,sample(idx.max,size=1,replace=F))
    if(length(idx.max)==1){
      idx.max <- idx.max
    }else{
      idx.max <- sample(idx.max,size=1,replace=F)
    }
    
    k.best.vec[i] <- k.inner.best.vec[idx.max]
    cat("Outer Best k's: ",k.best.vec[i],"\n")
    
    npdr.positives <- att.best.inner.list[[idx.max]]
    
    if(attr.diff.type=="correlation-data"){
      idx.positive <- numeric()
      new.attr.idx.list <- list()
      for(pos in 1:length(npdr.positives)){
        idx.positive[pos] <- as.numeric(strsplit(npdr.positives[pos],split="ROI")[[1]][2])
        new.attr.idx.list[[pos]] <- attr.idx.list[[idx.positive[pos]]]
      }
    
      col.idx.npdr.positives <- unlist(new.attr.idx.list)
    
      rf.model.outer <- randomForest(trn.data[,col.idx.npdr.positives], y=trn.pheno,
                                   mtry=max(floor(ncol(trn.data[,col.idx.npdr.positives])/3),1), ntree=num_tree,
                                   na.action=na.omit)
      pred.class <- stats::predict(rf.model.outer,tst.data[,col.idx.npdr.positives],na.action=na.omit)
    }else{
      rf.model.outer <- randomForest(trn.data[,npdr.positives], y=trn.pheno,
                                     mtry=max(floor(ncol(trn.data[,npdr.positives])/3),1), ntree=num_tree,
                                     na.action=na.omit)
      pred.class <- stats::predict(rf.model.outer,tst.data[,npdr.positives],na.action=na.omit)
    }
    
    pred.accu.vec[i] <- 1 - mean(pred.class != tst.pheno)
    cat("Test Accuracy: ",pred.accu.vec[i],"\n")
    att.best.list[[i]] <- npdr.positives
    
    if(attr.diff.type != "correlation-data"){
      idx.att <- c(grep(npdr.positives,pattern="main",value=TRUE),
                 grep(npdr.positives,pattern="int",value=TRUE))
      cat("Functional: ",idx.att,"\n")
    }else{
      cat("ROIs: ",npdr.positives,"\n")
    }
  }
  max.accu <- max(pred.accu.vec)
  idx.max <- which(pred.accu.vec==max.accu)
  if(length(idx.max)==1){
    idx.max <- idx.max
  }else{
    idx.max <- sample(idx.max,size=1,replace=F)
  }
  
  idx.best.k <- idx.max
  best.k <- k.best.vec[idx.best.k]
  best.atts <- att.best.list[[idx.best.k]]

  if(attr.diff.type=="correlation-data"){
    idx.positive <- numeric()
    new.attr.idx.list <- list()
    for(pos in 1:length(best.atts)){
      idx.positive[pos] <- as.numeric(strsplit(best.atts[pos],split="ROI")[[1]][2])
      new.attr.idx.list[[pos]] <- attr.idx.list[[idx.positive[pos]]]
    }
  
    col.idx.npdr.positives <- unlist(new.attr.idx.list)
  
    train.data <- as.matrix(train.ds[,col.idx.npdr.positives])
    train.pheno <- as.factor(train.ds[,label])
  }else{
    train.data <- as.matrix(train.ds[,best.atts])
    train.pheno <- as.factor(train.ds[,label])
  }
  
  train.model <- randomForest::randomForest(train.data, y=train.pheno,
                              mtry=max(floor(ncol(train.data)/3),1), ntree=num_tree,
                              na.action=na.omit)
  train.accu <- 1 - mean(train.model$confusion[, "class.error"])
  
  if(!is.null(validation.ds)){
    
    if(attr.diff.type=="correlation-data"){
      test.data <- as.matrix(validation.ds[,col.idx.npdr.positives])
      
      npdr.results.tmp <- npdr(label, validation.ds, regression.type="binomial", attr.diff.type=attr.diff.type,
                               nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, knn=best.k,
                               neighbor.sampling="none",
                               padj.method="bonferroni", 
                               corr.attr.names=corr.attr.names,
                               verbose=verbose, dopar.nn=dopar.nn, dopar.reg=dopar.reg)
      #npdr.results.df <- data.frame(att=npdr.results$att,
      #                              betas=npdr.results$beta.Z.att,
      #                              pval.att=npdr.results$pval.adj)
    }else{
      test.data <- as.matrix(validation.ds[,best.atts])
      
      npdr.results.tmp <- npdr(label, validation.ds, regression.type="binomial", attr.diff.type=attr.diff.type,
                               nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, knn=best.k,
                               neighbor.sampling="none",
                               padj.method="bonferroni", 
                               corr.attr.names=corr.attr.names,
                               verbose=verbose, dopar.nn=dopar.nn, dopar.reg=dopar.reg)
      #npdr.results.df <- data.frame(att=npdr.results.tmp$att,
      #                              betas=npdr.results.tmp$beta.Z.att,
      #                              pval.att=npdr.results.tmp$pval.adj)
    }
    test.pheno <- as.factor(validation.ds[,label])

    test.pred <- stats::predict(train.model, newdata=test.data, na.action=na.omit)

    test.accu <- 1 - mean(test.pred != test.pheno)

    #test.accu <- confusionMatrix(as.factor(test.pred), as.factor(test.pheno))$byClass["Balanced Accuracy"]
  }

  list(cv.acc=train.accu, Validation=test.accu, Features=best.atts, Train.model=train.model, best.k=best.k, AllFeatures.NPDR=npdr.results.tmp)
  
}
