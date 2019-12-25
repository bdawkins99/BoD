# make figures showing null vs correlated data distance distributions
setwd("C:/Users/bdawk/Documents/KNN_project_output")
library(npdr)

# standard normal
#####################################################################
save.fig <- T
lowers <- c(0.0,0.25,0.5,0.75)
uppers <- c(0.0,0.35,0.6,0.85)
avg.abs.cors <- c(0.0048,0.2622,0.4575,0.6124)

for(iter in 1:4){
  num.variables <- 100
  num.samples <- 100
  plot.graph <- F
  nbias <- round(0.1*num.variables)
  
  hi.cor <- 0
  lo.cor <- 0

  set.seed(1987)
  null.dats <- matrix(rnorm(num.variables*num.samples),nrow=num.samples,ncol=num.variables)
  
  e <- 1    # fudge factor to the number of nodes to avoid giant component
  prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component

  prob <- 0.95
  set.seed(1989)
  g <- erdos.renyi.game(num.variables, prob) # Erdos-Renyi network

  # generate correlation matrix from g
  set.seed(1991)
  network.atts <- npdr::generate_structured_corrmat(g=g,
                                              num.variables=num.variables, 
                                              hi.cor.tmp=hi.cor, 
                                              lo.cor.tmp=lo.cor, 
                                              hi.cor.fixed=hi.cor,
                                              graph.type="Erdos-Renyi",
                                              plot.graph=plot.graph,
                                              nbias=nbias)

  R <- as.matrix(network.atts$corrmat) # correlation matrix

  U <- t(chol(R))                             # upper tri cholesky
  #null.dats <- t(U %*% t(null.dats)) # correlated data

  hi.cor <- uppers[iter]
  lo.cor <- lowers[iter]

  # generate correlation matrix from g
  prob <- 0.2
  set.seed(1989)
  g <- erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
  
  set.seed(1991)
  network.atts <- npdr::generate_structured_corrmat(g=g,
                                              num.variables=num.variables, 
                                              hi.cor.tmp=hi.cor, 
                                              lo.cor.tmp=lo.cor,
                                              hi.cor.fixed=hi.cor,
                                              graph.type="Erdos-Renyi",
                                              plot.graph=plot.graph,
                                              nbias=nbias)

  R <- as.matrix(network.atts$corrmat) # correlation matrix

  set.seed(1987)
  tmp <- matrix(rnorm(num.variables*num.samples),nrow=num.samples,ncol=num.variables)
  U <- t(chol(R))                             # upper tri cholesky
  dats <- t(U %*% t(tmp)) # correlated data

  d.mat.null <- as.matrix(dist(null.dats,upper=T,method='euclidean'))
  d.mat.corr <- as.matrix(dist(dats,upper=T,method='euclidean'))

  dist.null <- d.mat.null[upper.tri(d.mat.null)]
  dist.corr <- d.mat.corr[upper.tri(d.mat.corr)]

  tmp <- c(dist.null,dist.corr)
  my.breaks <- seq(min(tmp)-1,max(tmp)+1,length.out=40)
  tmp.hist <- hist(tmp,breaks=my.breaks,plot=F)
  ymax <- 1.75*max(tmp.hist$density)
  ymax <- 0.5
  
  #print(mean(cor(null.dats)[upper.tri(cor(null.dats))]))
  print(mean(cor(dats)[upper.tri(cor(dats))]))
  
  mean.cor.null <- mean(abs(cor(null.dats)[upper.tri(cor(null.dats))]))
  mean.cor.corr <- mean(abs(cor(dats)[upper.tri(cor(dats))]))

  if(save.fig==T){
    pdf(paste("null_vs_correlated_density",iter,".pdf",sep=""),height=6,width=6)
  }
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  hist(dist.null,breaks=my.breaks,freq=F,ylim=c(0,ymax),col='white',border='white',
       ylab="",xlab="",main="",
       font=2,cex.lab=1.5,cex.main=1.7,xlim=c(0,35),cex.axis=1.4)
  hist(dist.corr,breaks=my.breaks,freq=F,add=T,col='white',border='white')
  x.null <- seq(min(my.breaks),max(my.breaks),by=0.01)
  y.null <- dnorm(x.null,mean=mean(dist.null),sd=sd(dist.null))
  x.tmp <- seq(min(my.breaks)-20,max(my.breaks)+20,by=0.01)
  y.tmp <- dnorm(x.tmp,mean=mean(dist.null),sd=1)
  lines(x=x.tmp,y=y.tmp,lwd=2,lty=1,col='orange')
  #lines(density(dist.null,adjust=3,from=min(my.breaks)-5,to=max(my.breaks)+5),lwd=2,lty=1,col='orange')
  lines(density(dist.corr,adjust=3,from=-10,to=50),lwd=2,lty=1,col='blue')
  abline(v=sqrt(2*num.variables - 1),lty=2,col='orange',lwd=2.5)
  abline(v=mean(dist.corr),lty=2,col='blue',lwd=2.5)
  
  abline(v=mean(dist.null)+1,lty=3,col='orange',lwd=3)
  abline(v=mean(dist.null)-1,lty=3,col='orange',lwd=3)
  
  #abline(v=mean(dist.null)+sd(dist.null),lty=3,col='orange',lwd=2)
  #abline(v=mean(dist.null)-sd(dist.null),lty=3,col='orange',lwd=2)
  abline(v=mean(dist.corr)+sd(dist.corr),lty=3,col='blue',lwd=3)
  abline(v=mean(dist.corr)-sd(dist.corr),lty=3,col='blue',lwd=3)
  legend("topright",c("Uncorrelated Data","Correlated Data"),lty=0,fill=c("orange","blue"),cex=1.5,bg='white')
  legend("right",c("Mean","SD"),lty=c(2,3),lwd=c(2,3),col='black',bg='white',cex=1.5)

  #text(x=20,y=0.1,labels=paste("|r| = ",round(mean.cor.corr,4),sep=""))
  box()
  if(save.fig==T){
    dev.off()
  }
}
#####################################################################

# GWAS - TiTv
#####################################################################
library(bindata)
library(npdr)

rmvBinomial2 <- function(n, size, probs, corrmat) {
  X <- replicate(n, {
    colSums(rmvbin(size, margprob=probs, sigma=corrmat))
  })
  t(X)
}

save.fig <- T
lowers <- c(0.0,0.25,0.5,0.75)
uppers <- c(0.0,0.35,0.6,0.85)
avg.abs.cors <- c(0.0846,0.1477,0.2333,0.3130)

for(iter in 1:4){

  num.variables <- 100
  plot.graph <- F
  nbias <- round(0.1*num.variables)
  hi.cor <- uppers[iter]
  lo.cor <- lowers[iter]

  e <- 1    # fudge factor to the number of nodes to avoid giant component
  prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
  prob <- 0.95

  set.seed(1989)
  g <- erdos.renyi.game(num.variables, prob) # Erdos-Renyi network

  # generate correlation matrix from g
  set.seed(1991)
  network.atts <- npdr::generate_structured_corrmat(g=g,
                                            num.variables=num.variables, 
                                            hi.cor.tmp=hi.cor, 
                                            lo.cor.tmp=lo.cor, 
                                            graph.type="Erdos-Renyi",
                                            plot.graph=plot.graph,
                                            nbias=nbias)

  R <- as.matrix(network.atts$corrmat) # correlation matrix

  set.seed(1987)
  f.a <- runif(num.variables,min=0.01,max=0.95)

  set.seed(1990)
  X <- rmvBinomial2(n=100,size=2,probs=f.a,corrmat=R)

  data.mat <- X
  corr.dats <- X

  m <- dim(X)[1]
  p <- dim(X)[2]

  my.eta <- 2
  g1 <- 1/(my.eta+1)

  set.seed(1964)
  g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
  g2 <- my.eta/(my.eta + 1) - g0

  set.seed(1962)
  PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0,g1,g2),replace=T)

  d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
  for(i in 1:m){
  
    for(j in i:m){
      diff <- abs(data.mat[i,]-data.mat[j,])
      diff.TiTv <- diff
      for(k in 1:p){
      
        if(diff[k]==1){
          if(PuPy.samp[k]==3){
            diff.TiTv[k] <- 0.25
          }else if(PuPy.samp[k]==5){
            diff.TiTv[k] <- 0.25
          }else{
            diff.TiTv[k] <- 0.5
          }
        }else if(diff[k]==2){
          if(PuPy.samp[k]==3){
            diff.TiTv[k] <- 0.75
          }else if(PuPy.samp[k]==5){
            diff.TiTv[k] <- 0.75
          }else{
            diff.TiTv[k] <- 1
          }
        }
      }
      d.mat[i,j] <- sum(diff.TiTv)
    
    }
  }
  tmp <- d.mat + t(d.mat)

  F.a <- ((1 - f.a)^3)*f.a + (f.a^3)*(1 - f.a)
  G.a <- ((1 - f.a)^2)*(f.a^2)
  mean.titv.corr <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
  var.titv.corr <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)

  dist.vec.corr <- tmp[upper.tri(tmp)]

  # null data
  data.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
  set.seed(1991)
  for(i in 1:m){
    samp <- matrix(as.double(rbinom(p,2,f.a)),nrow=1,ncol=p)
    data.mat[i,] <- samp
  }

  set.seed(1990)
  X <- rmvBinomial2(n=100,size=2,probs=f.a,corrmat=diag(length(f.a)))
  null.dats <- data.mat
  #data.mat <- X
  #null.dats <- X

  d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
  for(i in 1:m){
  
    for(j in i:m){
      diff <- abs(data.mat[i,]-data.mat[j,])
      diff.TiTv <- diff
      for(k in 1:p){
      
        if(diff[k]==1){
          if(PuPy.samp[k]==3){
            diff.TiTv[k] <- 0.25
          }else if(PuPy.samp[k]==5){
            diff.TiTv[k] <- 0.25
          }else{
            diff.TiTv[k] <- 0.5
          }
        }else if(diff[k]==2){
          if(PuPy.samp[k]==3){
            diff.TiTv[k] <- 0.75
          }else if(PuPy.samp[k]==5){
            diff.TiTv[k] <- 0.75
          }else{
            diff.TiTv[k] <- 1
          }
        }
      }
      d.mat[i,j] <- sum(diff.TiTv)
    
    }
  }
  tmp <- d.mat + t(d.mat)

  F.a <- ((1 - f.a)^3)*f.a + (f.a^3)*(1 - f.a)
  G.a <- ((1 - f.a)^2)*(f.a^2)
  mean.titv.null <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
  var.titv.null <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)

  dist.vec.null <- tmp[upper.tri(tmp)]

  # plots
  tmp <- c(dist.vec.null,dist.vec.corr)
  my.breaks <- seq(min(tmp)-1,max(tmp)+1,length.out=40)
  tmp.hist <- hist(tmp,breaks=my.breaks,plot=F)
  ymax <- 5*max(tmp.hist$density)
  ymax <- 0.3

  mean.cor.null <- mean(abs(cor(null.dats)[upper.tri(cor(null.dats))]))
  mean.cor.corr <- mean(abs(cor(corr.dats)[upper.tri(cor(corr.dats))]))

  #save.fig <- F
  if(save.fig==T){
    pdf(paste("null_vs_correlated_density-TiTv",iter,".pdf",sep=""),height=6,width=6)
  }
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  hist(dist.vec.null,breaks=my.breaks,freq=F,ylim=c(0,ymax),col='white',border='white',
      ylab="",xlab="",main="",
      font=2,cex.lab=1.5,cex.main=1.7,xlim=c(5,50),cex.axis=1.4)
  hist(dist.vec.corr,breaks=my.breaks,freq=F,add=T,col='white',border='white')
  x.null <- seq(min(my.breaks),max(my.breaks),by=0.01)
  y.null <- dnorm(x.null,mean=mean(dist.vec.null),sd=sd(dist.vec.null))
  x.tmp <- seq(min(my.breaks)-20,max(my.breaks)+20,by=0.01)
  y.tmp <- dnorm(x.tmp,mean=mean(dist.vec.null),sd=sd(dist.vec.null))
  lines(x=x.tmp,y=y.tmp,lwd=2,lty=1,col='orange')
  #lines(density(dist.null,adjust=3,from=min(my.breaks)-5,to=max(my.breaks)+5),lwd=2,lty=1,col='orange')
  lines(density(dist.vec.corr,adjust=3,from=-10,to=70),lwd=2,lty=1,col='blue')
  abline(v=mean(dist.vec.null),lty=2,col='orange',lwd=2.5)
  abline(v=mean(dist.vec.corr),lty=2,col='blue',lwd=2.5)

  abline(v=mean(dist.vec.null)+sd(dist.vec.null),lty=3,col='orange',lwd=3)
  abline(v=mean(dist.vec.null)-sd(dist.vec.null),lty=3,col='orange',lwd=3)

  #abline(v=mean(dist.null)+sd(dist.null),lty=3,col='orange',lwd=2)
  #abline(v=mean(dist.null)-sd(dist.null),lty=3,col='orange',lwd=2)
  abline(v=mean(dist.vec.corr)+sd(dist.vec.corr),lty=3,col='blue',lwd=3)
  abline(v=mean(dist.vec.corr)-sd(dist.vec.corr),lty=3,col='blue',lwd=3)
  legend("topright",c("Uncorrelated Data","Correlated Data"),lty=0,fill=c("orange","blue"),cex=1.5,bg='white')
  legend("right",c("Mean","SD"),lty=c(2,3),lwd=c(2,3),col='black',bg='white',cex=1.5)

  #text(x=20,y=0.1,labels=paste("|r| = ",round(mean.cor.corr,4),sep=""))
  box()
  if(save.fig==T){
    dev.off()
  }
  print(mean.cor.corr)

}

#####################################################################

# correlation data (rs-fMRI)
#####################################################################
r_to_z_fn <- function(x){
  r.z <- 0.5*log((1 + x)/(1 - x))
  r.z
}

stretch_mat <- function(M){
  
  mat <- numeric()
  for(k in 1:nrow(M)){
    mat <- c(mat,M[k,-k])
  }
  return(mat)
}

stretch_mat2 <- function(M){
  
  mat <- NULL
  for(k in 1:(nrow(M) - 1)){
    mat <- c(mat, M[k,(k + 1):ncol(M)])
  }
  return(mat)
}

save.fig <- T

num.variables <- 30
num.samples <- 100
plot.graph <- F
nbias <- round(0.1*((p*(p-1))/2))

lowers <- c(0.0,0.15,0.25,0.35)
uppers <- c(0.0,0.25,0.35,0.45)
avg.abs.cors <- c(0.0909,0.1157,0.1589,0.2136)

for(iter in 1:4){
  hi.cor <- uppers[iter]
  lo.cor <- lowers[iter]

  m <- num.samples
  p <- num.variables

  e <- 1    # fudge factor to the number of nodes to avoid giant component
  prob <- 1/(((p*(p-1))/2)+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
  #prob <- 0.9
  set.seed(1989)
  g <- erdos.renyi.game(((p*(p-1))/2), prob) # Erdos-Renyi network

  set.seed(1991)
  network.atts <- npdr::generate_structured_corrmat(g=g,
                                                  num.variables=((p*(p-1))/2), 
                                                  hi.cor.tmp=hi.cor, 
                                                  lo.cor.tmp=lo.cor, 
                                                  hi.cor.fixed=hi.cor,
                                                  graph.type="Erdos-Renyi",
                                                  plot.graph=plot.graph,
                                                  nbias=nbias)               

  A <- network.atts$A.mat

  R <- as.matrix(network.atts$corrmat) # correlation matrix
  #U <- matrix(rnorm(p^2),nrow=p,ncol=p)
  U <- t(chol(R))
  #Y <- matrix(rnorm(m*p*(p-1)),nrow=m,ncol=p*(p-1))
  Y <- matrix(0, nrow=m, ncol=((p*(p-1))/2))
  
  set.seed(1990)
  for(k in 1:m){
    
    null.mat <- matrix(rnorm(m*p),nrow=m,ncol=p)
    cor.tmp <- cor(null.mat) # correlation matrix
  
    zcorr <- apply(matrix(stretch_mat2(cor.tmp),ncol=1),1,r_to_z_fn)
    Y[k,] <- matrix(zcorr,ncol=((p*(p-1))/2),nrow=1)
  
  }
  Y <- scale(t(Y))
  Y <- t(Y)
  Y.new <- t(U %*% t(Y))

  Y.full <- matrix(0, nrow=m, ncol=(p*(p-1)))
  for(k in 1:m){
    mat.tmp <- Y.new[k,]
    zero.mat <- matrix(0,nrow=p,ncol=p)
  
    zero.mat[lower.tri(zero.mat,diag=F)] <- mat.tmp
    zero.mat <- t(zero.mat)
    zero.mat <- zero.mat + t(zero.mat)
  
    Y.full[k,] <- c(stretch_mat(zero.mat))
  
  }
  my.cor <- cor(Y.full)[upper.tri(cor(Y.full))]
  print(mean(abs(my.cor)))

  tmp <- as.matrix(dist(Y.full,upper=T,method="manhattan"))

  h <- hist(tmp[upper.tri(tmp)],breaks=30,plot=F)
  xmin <- 0.8*min(h$breaks)
  xmax <- 2.5*max(h$breaks)
  ymax <- 0.02

  if(save.fig==T){
    pdf(paste("null_vs_correlated_density-rsfMRI",iter,".pdf",sep=""),height=6,width=6)
  }
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  hist(tmp[upper.tri(tmp)],breaks=30,freq=F,ylim=c(0,ymax),xlim=c(700,1500),border='white',col='white',
       ylab="",xlab="",main="",
       font=2,cex.lab=1.5,cex.main=1.7,cex.axis=1.4)


  Y.null <- matrix(rnorm(m*p*(p-1)),nrow=m)
  d.mat2 <- as.matrix(dist(Y.null, upper=T, method="manhattan"))
  x <- seq(xmin,xmax,by=0.1)
  y <- dnorm(x, mean=mean(d.mat2[upper.tri(d.mat2)]), sd=sd(d.mat2[upper.tri(d.mat2)]))
  lines(x,y,lty=1,lwd=2,col='orange')

  lines(density(tmp[upper.tri(tmp)],adjust=3,from=xmin,to=xmax),lwd=2,col='blue')
  abline(v=mean(tmp[upper.tri(tmp)]),lty=2,lwd=2.5,col='blue')
  abline(v=mean(tmp[upper.tri(tmp)])+sd(tmp[upper.tri(tmp)]),lty=3,lwd=3,col='blue')
  abline(v=mean(tmp[upper.tri(tmp)])-sd(tmp[upper.tri(tmp)]),lty=3,lwd=3,col='blue')
  #abline(v=(2*p*(p-1))/sqrt(pi),lty=2,lwd=2,col='orange')
  abline(v=mean(d.mat2[upper.tri(d.mat2)]),lty=2,lwd=2.5,col='orange')
  abline(v=mean(d.mat2[upper.tri(d.mat2)])+sd(d.mat2[upper.tri(d.mat2)]),lty=3,lwd=3,col='orange')
  abline(v=mean(d.mat2[upper.tri(d.mat2)])-sd(d.mat2[upper.tri(d.mat2)]),lty=3,lwd=3,col='orange')
  
  legend("topright",c("Uncorrelated Data","Correlated Data"),lty=0,fill=c("orange","blue"),cex=1.5,bg='white')
  legend("right",c("Mean","SD"),lty=c(2,3),lwd=c(2,3),col='black',bg='white',cex=1.5)
  box()
  if(save.fig==T){
    dev.off()
  }
}


#####################################################################