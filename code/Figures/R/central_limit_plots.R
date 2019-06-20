# asymptotically normal distance distribution plots

ps <- c(10,50,100,1000,5000,10000)
dist.list.manhattan <- list()
dist.list.euclidean <- list()
for(i in 1:length(ps)){
  print(i)
p <- ps[i] # number of attributes

m <- 100  # number of subjects

# data matrix
sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
for(j in 1:m){
  sub.mat[j,] <- matrix(runif(p),nrow=1,ncol=p)
}

# find range for each attribute
maxs <- apply(sub.mat,2,max) # attribute maximums
mins <- apply(sub.mat,2,min) # attribute minimums
diffs <- maxs-mins           # attribute ranges
diffs.inv <- 1/diffs         # inverse of attribute ranges

# distance matrix
d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
d.mat2 <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
for(j in 1:m){
  for(k in j:m){
    
    d.mat2[j,k] <- sqrt(sum(((sub.mat[j,]-sub.mat[k,])^2)*diffs.inv^2))
    d.mat[j,k] <- sum(abs(sub.mat[j,]-sub.mat[k,])*diffs.inv)
    #d.mat[j,k] <- sum(abs(sub.mat[j,]-sub.mat[k,]))
    #d.mat2[j,k] <- sqrt(sum(((sub.mat[j,]-sub.mat[k,])^2)))
    
  }
  
}

tmp <- d.mat + t(d.mat)
d.mat <- tmp
dist.vec <- c(t(d.mat))
idx <- which(dist.vec==0)
dist.vec <- dist.vec[-idx] # distance vector
dist.list.manhattan[[i]] <- dist.vec

tmp <- d.mat2 + t(d.mat2)
d.mat2 <- tmp
dist.vec <- c(t(d.mat2))
idx <- which(dist.vec==0)
dist.vec <- dist.vec[-idx] # distance vector
dist.list.euclidean[[i]] <- dist.vec
}

setwd("C:/Users/bdawk/Documents/KNN_project_output")
# mar = c(bottom, left, top, right) 
pdf("central_limit_distances-Manhattan-Uniform-diff.pdf",height=10,width=8)
par(mfrow=c(3,2))
for(i in 1:length(dist.list.manhattan)){
  
  #if(i < 3){
  #  par(mar=c(2.5,4,1.4,0.8))
  #}else if(i >=3 && i < 5){
  #  par(mar=c(2.5,4,1.3,0.8))
  #}else{
  #  par(mar=c(4.5,4,1.3,0.8))
  #}
  
  if(i==1){
    par(mar=c(2.5,4,1.4,0.2))
  }else if(i==2){
    par(mar=c(2.5,2.6,1.4,1.3))
  }else if(i==3){
    par(mar=c(2.5,4,1.3,0.2))
  }else if(i==4){
    par(mar=c(2.5,2.6,1.3,1.3))
  }else if(i==5){
    par(mar=c(4.5,4,1.3,0.2))
  }else{
    par(mar=c(4.5,2.6,1.3,1.3))
  }

if(i==1){
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="",
       main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}else if(i==2){
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="",ylab="",
       main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}else if(i==3){
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="",main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}else if(i==4){
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="",ylab="",
       main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}else if(i==5){
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="Manhattan (normalized diff)",
       main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}else{
  hist(dist.list.manhattan[[i]],breaks=100,freq=F,
       xlab="Manhattan (normalized diff)",ylab="",
       main="",
       font.lab=2,cex.lab=1.2)
  title(paste("Histogram of Manhattan Distances"," (p = ",ps[i],")",sep=""),line=0.5)
}

x <- seq(min(dist.list.manhattan[[i]]),max(dist.list.manhattan[[i]]),by=0.01)
y <- dnorm(x,mean=mean(dist.list.manhattan[[i]]),sd=sd(dist.list.manhattan[[i]]))
lines(x,y,lty=1,col='red',lwd=2)
legend("topright",c("Normal Density",paste("mean = ",round(mean(dist.list.manhattan[[i]]),3),sep=""),paste("var = ",round(var(dist.list.manhattan[[i]]),3),sep="")),
       lty=1,col=c("red","white","white"),lwd=2,bg="white",cex=1)
box()
}
dev.off()

# mar = c(bottom, left, top, right) 
pdf("central_limit_distances-Euclidean-Uniform-diff.pdf",height=10,width=8)
par(mfrow=c(3,2))
for(i in 1:length(dist.list.euclidean)){
  
  #if(i < 3){
  #  par(mar=c(2.5,4,1.4,0.8))
  #}else if(i >=3 && i < 5){
  #  par(mar=c(2.5,4,1.3,0.8))
  #}else{
  #  par(mar=c(4.5,4,1.3,0.8))
  #}
  
  if(i==1){
    par(mar=c(2.5,4,1.4,0.2))
  }else if(i==2){
    par(mar=c(2.5,2.6,1.4,1.3))
  }else if(i==3){
    par(mar=c(2.5,4,1.3,0.2))
  }else if(i==4){
    par(mar=c(2.5,2.6,1.3,1.3))
  }else if(i==5){
    par(mar=c(4.5,4,1.3,0.2))
  }else{
    par(mar=c(4.5,2.6,1.3,1.3))
  }
  
  if(i==1){
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="",
         main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }else if(i==2){
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="",ylab="",
         main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }else if(i==3){
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="",main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }else if(i==4){
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="",ylab="",
         main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }else if(i==5){
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="Euclidean (normalized diff)",
         main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }else{
    hist(dist.list.euclidean[[i]],breaks=100,freq=F,
         xlab="Euclidean (normalized diff)",ylab="",
         main="",
         font.lab=2,cex.lab=1.2)
    title(paste("Histogram of Euclidean Distances"," (p = ",ps[i],")",sep=""),line=0.5)
  }
  
  x <- seq(min(dist.list.euclidean[[i]]),max(dist.list.euclidean[[i]]),by=0.01)
  y <- dnorm(x,mean=mean(dist.list.euclidean[[i]]),sd=sd(dist.list.euclidean[[i]]))
  lines(x,y,lty=1,col='red',lwd=2)
  legend("topright",c("Normal Density",paste("mean = ",round(mean(dist.list.euclidean[[i]]),3),sep=""),paste("var = ",round(var(dist.list.euclidean[[i]]),3),sep="")),
         lty=1,col=c("red","white","white"),lwd=2,bg="white",cex=1)
  box()
}
dev.off()

