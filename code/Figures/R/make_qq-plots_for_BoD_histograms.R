# make QQ-plots for histograms Trang made for Manhattan/Euclidean with p=10,100,10000
set.seed(1618)

m <- 100  # number of subjects
ps <- c(10, 100, 10000) # numbers of attributes
dist.list.manhattan <- list()
dist.list.euclidean <- list()
for(i in 1:length(ps)){
  print(i)
  p <- ps[i] 
  
  # data matrix
  sub.mat <- matrix(rep(0, length = m*p), nrow = m, ncol = p)
  for(j in 1:m){
    sub.mat[j,] <- matrix(runif(p), nrow = 1, ncol = p)
  }
  
  # find range for each attribute
  maxs <- apply(sub.mat, 2, max) # attribute maximums
  mins <- apply(sub.mat, 2, min) # attribute minimums
  diffs <- maxs - mins           # attribute ranges
  diffs.inv <- 1/diffs           # inverse of attribute ranges
  sub.mat <- sub.mat %*% diag(diffs.inv) # scale the matrix before dist
  
  d.mat <- as.matrix(dist(sub.mat, method = 'manhattan', upper=T))
  d.mat2 <- as.matrix(dist(sub.mat, method = 'euclidean', upper=T))
  
  dist.list.manhattan[[i]] <- d.mat[upper.tri(d.mat)]
  
  dist.list.euclidean[[i]] <- d.mat2[upper.tri(d.mat2)]
}

Ws <- c("W = 0.9981", "W = 0.9991", "W = 0.9995")
Ps <- c("P = 1.12e-05", "P = 0.0128", "P = 0.262")

setwd("C:/Users/bdawk/Documents/KNN_project_output")
for(i in 1:length(dist.list.manhattan)){
  pdf(paste("manhattan_QQ-plot_p",ps[i],".pdf",sep=""),height=4.5,width=5)
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  qqnorm(dist.list.manhattan[[i]],main="",
         xlab="",ylab="",font=2,cex.lab=1.3,
         pch=21,bg='gray',col='black',cex=1.4,cex.axis=1.4,
         panel.first=grid(col = "lightgray", lty = "solid",
              lwd = 2, equilogs = TRUE))

  qqline(dist.list.manhattan[[i]],lwd=2,col='#cc79a7')
  my.point.W <- 0.9
  my.point.P <- 0.8
  y.W <- (1 - my.point.W)*min(dist.list.manhattan[[i]]) + my.point.W*max(dist.list.manhattan[[i]])
  y.P <- (1 - my.point.P)*min(dist.list.manhattan[[i]]) + my.point.P*max(dist.list.manhattan[[i]])
  text(x=-2,y=y.W,labels=Ws[i],cex=1.3)
  text(x=-2,y=y.P,labels=Ps[i],cex=1.3)
  dev.off()
}


Ws <- c("W = 0.9994", "W = 0.9996", "W = 0.9996")
Ps <- c("P = 0.145", "P = 0.421", "P = 0.483")

for(i in 1:length(dist.list.euclidean)){
  pdf(paste("euclidean_QQ-plot_p",ps[i],".pdf",sep=""),height=4.5,width=5)
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  qqnorm(dist.list.euclidean[[i]],main="",
         xlab="",ylab="",font=2,cex.lab=1.3,
         pch=21,bg='gray',col='black',cex=1.4,cex.axis=1.4,
         panel.first=grid(col = "lightgray", lty = "solid",
                          lwd = 2, equilogs = TRUE))
  qqline(dist.list.euclidean[[i]],lwd=2,col='#cc79a7')
  my.point.W <- 0.9
  my.point.P <- 0.8
  y.W <- (1 - my.point.W)*min(dist.list.euclidean[[i]]) + my.point.W*max(dist.list.euclidean[[i]])
  y.P <- (1 - my.point.P)*min(dist.list.euclidean[[i]]) + my.point.P*max(dist.list.euclidean[[i]])
  text(x=-2,y=y.W,labels=Ws[i],cex=1.3)
  text(x=-2,y=y.P,labels=Ps[i],cex=1.3)
  dev.off()
}

