# generate histograms for TiTv distances that explore the effect of:
# 
# (1) varying the Ti/Tv ratio and fixing the MAF probs
# (2) varying the MAF probs with fixed Ti/Tv ratio

partial_sum_func <- function(M,b){
  part.sum <- sum(M[1:(b-1)])
  return(part.sum)
}

# generate varied MAF histograms
#########################################################################################
lowers <- c(0.01,0.1,0.2,0.3)
uppers <- c(0.1,0.2,0.3,0.4)

save.fig <- T
for(maf in 1:4){
  samp.var <- numeric()
  theo.var1 <- numeric()
  theo.var2 <- numeric()
  theo.var3 <- numeric()
  theo.mean <- numeric()
  theo.mean2 <- numeric()
  theo.mean3 <- numeric()
  samp.mean <- numeric()
  samp.m2 <- numeric()
  theo.m2 <- numeric()
  long.dist.vec <- NULL

  m <- 100
  p <- 100

  set.seed(1989)
  for(myiter in 1:100){
    print(myiter)

    my.eta <- 2
    g1 <- 1/(my.eta+1)
    g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
    g2 <- my.eta/(my.eta + 1) - g0
  
    PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0,g1,g2),replace=T)
    table(PuPy.samp)/sum(table(PuPy.samp))
  
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    probs <- runif(p,lowers[maf],uppers[maf])
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    for(i in 1:m){
      samp <- matrix(as.double(rbinom(p,2,probs)),nrow=p,ncol=1)
      data.mat[,i] <- samp
    }
  
    w0 <- mean((1-probs)^2)
    w1 <- 2*mean(probs*(1-probs))
    w2 <- mean(probs^2)
  
    diff.mat <- matrix(rep(0,length=4*(m*(m-1)/2)),nrow=(m*(m-1)/2),ncol=5)
    iter <- 1
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    for(i in 1:m){

      for(j in i:m){
        diff <- abs(data.mat[,i]-data.mat[,j])
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
    long.dist.vec <- c(long.dist.vec,tmp[upper.tri(tmp)])
    dist.vec <- c(t(tmp))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
  
    E.diff0 <- mean((1 - probs)^4 + 4*(probs^2)*((1 - probs)^2) + probs^4)
    
    E.diff.25 <- 4*(g0 + g2)*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
  
    E.diff.5 <- 4*g1*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
  
    E.diff.75 <- 2*(g0 + g2)*mean(((1 - probs)^2)*(probs^2))
  
    E.diff1 <- 2*g1*mean(((1 - probs)^2)*(probs^2))
  
    E.Dij <- (0*E.diff0 + 0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + 1*E.diff1)*p
    exp.mean <- (g0 + g2 + 2*g1)*sum(probs*(1 - probs)^3 + (1 - probs)*probs^3) + (1.5*(g0 + g2) + 2*g1)*sum(((1 - probs)^2)*(probs^2))
    exp.mean
    E.Dij
    mean(dist.vec)
  
    ######################################################
    w.v <- 2*(g0*g1 + g1*g2 + g0*g2)
    w.i <- 1 - 2*w.v
    F.a <- ((1 - probs)^3)*probs + (probs^3)*(1 - probs)
    G.a <- ((1 - probs)^2)*(probs^2)
  
    new.mean <- (w.i+2*w.v)*sum(F.a) + (1.5*w.i + 2*w.v)*sum(G.a)
  
    new.mean2 <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
  
    ######################################################
  
    E.D2.1 <- (0.25*(g0+g2)+g1)*sum(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((9/8)*(g0+g2)+2*g1)*sum(((1-probs)^2)*(probs^2))
  
    v <- (g0+g2+2*g1)*(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((3/2)*(g0+g2)+2*g1)*((1-probs)^2)*((probs)^2)
    v1 <- (g0+g2+2*g1)*(((1-probs[-1])^3)*probs[-1]+(probs[-1]^3)*(1-probs[-1])) + ((3/2)*(g0+g2)+2*g1)*((1-probs[-1])^2)*((probs[-1])^2)
    v2 <- numeric()
    for(k in 2:p){
      v2[(k-1)] <- partial_sum_func(v,k)
    }
  
    E.D2.2 <- c(2*matrix(v1,nrow=1,ncol=length(v1)) %*% matrix(v2,nrow=length(v2),ncol=1))
    E.D2.2

    E.Dsqij.new <- E.D2.1 + E.D2.2
    E.Dsqij.new
  
    E.Dsqij <- (0*E.diff0 + (1/16)*E.diff.25 + 0.25*E.diff.5 + (9/16)*E.diff.75 + 1*E.diff1)*p + (p^2 - p)*(0*E.diff0 + 0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + 1*E.diff1)^2
  
    Var.Dij <- E.Dsqij - (E.Dij)^2
    Var.Dij
    Var.Dij.new <- E.Dsqij.new - (E.Dij)^2
    Var.Dij.new
    var(dist.vec)
  
    exp.var <- (0.25*(g0 + g2) + g1)*sum(probs*(1-probs)^3 + (1-probs)*probs^3) + ((9/8)*(g0 + g2) + 2*g1)*sum(((1-probs)^2)*(probs^2)) - (1/p)*(((g0 + g2 + 2*g1)*sum(probs*(1-probs)^3 + (1-probs)*probs^3) + (1.5*(g0 + g2) + 2*g1)*sum(((1-probs)^2)*(probs^2)))^2)
    exp.var
  
    new.var <- mean((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2)) - mean(((g0 + g2+ 2*g1)^2)*((probs*((1-probs)^3) + (1-probs)*(probs^3))^2) - 2*(g0 + g2 + 2*g1)*(1.5*(g0 + g2) + 2*g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3))*((1-probs)^2)*(probs^2) - ((1.5*(g0 + g2) + 2*g1)^2)*((1-probs)^4)*(probs^4))
    new.var2 <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)
  
    theo.var1[myiter] <- Var.Dij.new
    theo.var2[myiter] <- new.var
    theo.var3[myiter] <- new.var2
    samp.var[myiter] <- var(dist.vec)
    theo.mean[myiter] <- E.Dij
    theo.mean2[myiter] <- new.mean
    theo.mean3[myiter] <- new.mean2
    samp.mean[myiter] <- mean(dist.vec)
  
    samp.m2[myiter] <- mean(dist.vec^2)
    theo.m2[myiter] <- p*sum((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2))
  
  }

  my.breaks <- seq(0,40,length.out=30)
  setwd("C:/Users/bdawk/Documents/KNN_project_output")
  if(save.fig==T){
    pdf(paste("TiTv_distance_histogram_maf",maf,".pdf",sep=""),height=6,width=6)
  }
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  hist(long.dist.vec,breaks=my.breaks,main="",ylab="",xlab="",font.lab=2,
      cex.lab=1.5,cex.main=1.7,ylim=c(0,0.25),freq=F,cex.axis=1.4)
  abline(v=mean(long.dist.vec),lty=1,col='red',lwd=2)
  abline(v=mean(theo.mean3),lty=1,col='blue',lwd=2)

  abline(v=mean(long.dist.vec)+sd(long.dist.vec),lty=2,col='orange',lwd=2)
  abline(v=mean(long.dist.vec)-sd(long.dist.vec),lty=2,col='orange',lwd=2)

  abline(v=mean(theo.mean3)+sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
  abline(v=mean(theo.mean3)-sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
  legend("topright",c("sample mean","theoretical mean","sample SD","theoretical SD"),
        lty=c(1,1,2,2),col=c('red','blue','orange','purple'),lwd=2,bg='white',cex=1.5)
  box()
  if(save.fig==T){
    dev.off()
  }

}

#########################################################################################

# generate varied Ti/Tv histograms
#########################################################################################
lowers <- c(0.01,0.1,0.2,0.3)
uppers <- c(0.1,0.2,0.3,0.4)
etas <- c(2,1.5,1,0.5)

save.fig <- T
for(titv in 1:4){
  samp.var <- numeric()
  theo.var1 <- numeric()
  theo.var2 <- numeric()
  theo.var3 <- numeric()
  theo.mean <- numeric()
  theo.mean2 <- numeric()
  theo.mean3 <- numeric()
  samp.mean <- numeric()
  samp.m2 <- numeric()
  theo.m2 <- numeric()
  long.dist.vec <- NULL
  
  m <- 100
  p <- 100
  
  set.seed(1989)
  for(myiter in 1:100){
    print(myiter)
    
    my.eta <- etas[titv]
    g1 <- 1/(my.eta+1)
    g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
    g2 <- my.eta/(my.eta + 1) - g0
    
    PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0,g1,g2),replace=T)
    table(PuPy.samp)/sum(table(PuPy.samp))
    
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    probs <- runif(p,lowers[1],uppers[1])
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    for(i in 1:m){
      samp <- matrix(as.double(rbinom(p,2,probs)),nrow=p,ncol=1)
      data.mat[,i] <- samp
    }
    
    w0 <- mean((1-probs)^2)
    w1 <- 2*mean(probs*(1-probs))
    w2 <- mean(probs^2)
    
    diff.mat <- matrix(rep(0,length=4*(m*(m-1)/2)),nrow=(m*(m-1)/2),ncol=5)
    iter <- 1
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    for(i in 1:m){
      
      for(j in i:m){
        diff <- abs(data.mat[,i]-data.mat[,j])
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
    long.dist.vec <- c(long.dist.vec,tmp[upper.tri(tmp)])
    dist.vec <- c(t(tmp))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    E.diff0 <- mean((1 - probs)^4 + 4*(probs^2)*((1 - probs)^2) + probs^4)
    
    E.diff.25 <- 4*(g0 + g2)*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
    
    E.diff.5 <- 4*g1*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
    
    E.diff.75 <- 2*(g0 + g2)*mean(((1 - probs)^2)*(probs^2))
    
    E.diff1 <- 2*g1*mean(((1 - probs)^2)*(probs^2))
    
    E.Dij <- (0*E.diff0 + 0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + 1*E.diff1)*p
    exp.mean <- (g0 + g2 + 2*g1)*sum(probs*(1 - probs)^3 + (1 - probs)*probs^3) + (1.5*(g0 + g2) + 2*g1)*sum(((1 - probs)^2)*(probs^2))
    exp.mean
    E.Dij
    mean(dist.vec)
    
    ######################################################
    w.v <- 2*(g0*g1 + g1*g2 + g0*g2)
    w.i <- 1 - 2*w.v
    F.a <- ((1 - probs)^3)*probs + (probs^3)*(1 - probs)
    G.a <- ((1 - probs)^2)*(probs^2)
    
    new.mean <- (w.i+2*w.v)*sum(F.a) + (1.5*w.i + 2*w.v)*sum(G.a)
    
    new.mean2 <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
    
    ######################################################
    
    E.D2.1 <- (0.25*(g0+g2)+g1)*sum(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((9/8)*(g0+g2)+2*g1)*sum(((1-probs)^2)*(probs^2))
    
    v <- (g0+g2+2*g1)*(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((3/2)*(g0+g2)+2*g1)*((1-probs)^2)*((probs)^2)
    v1 <- (g0+g2+2*g1)*(((1-probs[-1])^3)*probs[-1]+(probs[-1]^3)*(1-probs[-1])) + ((3/2)*(g0+g2)+2*g1)*((1-probs[-1])^2)*((probs[-1])^2)
    v2 <- numeric()
    for(k in 2:p){
      v2[(k-1)] <- partial_sum_func(v,k)
    }
    
    E.D2.2 <- c(2*matrix(v1,nrow=1,ncol=length(v1)) %*% matrix(v2,nrow=length(v2),ncol=1))
    E.D2.2
    
    E.Dsqij.new <- E.D2.1 + E.D2.2
    E.Dsqij.new
    
    E.Dsqij <- (0*E.diff0 + (1/16)*E.diff.25 + 0.25*E.diff.5 + (9/16)*E.diff.75 + 1*E.diff1)*p + (p^2 - p)*(0*E.diff0 + 0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + 1*E.diff1)^2
    
    Var.Dij <- E.Dsqij - (E.Dij)^2
    Var.Dij
    Var.Dij.new <- E.Dsqij.new - (E.Dij)^2
    Var.Dij.new
    var(dist.vec)
    
    exp.var <- (0.25*(g0 + g2) + g1)*sum(probs*(1-probs)^3 + (1-probs)*probs^3) + ((9/8)*(g0 + g2) + 2*g1)*sum(((1-probs)^2)*(probs^2)) - (1/p)*(((g0 + g2 + 2*g1)*sum(probs*(1-probs)^3 + (1-probs)*probs^3) + (1.5*(g0 + g2) + 2*g1)*sum(((1-probs)^2)*(probs^2)))^2)
    exp.var
    
    new.var <- mean((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2)) - mean(((g0 + g2+ 2*g1)^2)*((probs*((1-probs)^3) + (1-probs)*(probs^3))^2) - 2*(g0 + g2 + 2*g1)*(1.5*(g0 + g2) + 2*g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3))*((1-probs)^2)*(probs^2) - ((1.5*(g0 + g2) + 2*g1)^2)*((1-probs)^4)*(probs^4))
    new.var2 <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)
    
    theo.var1[myiter] <- Var.Dij.new
    theo.var2[myiter] <- new.var
    theo.var3[myiter] <- new.var2
    samp.var[myiter] <- var(dist.vec)
    theo.mean[myiter] <- E.Dij
    theo.mean2[myiter] <- new.mean
    theo.mean3[myiter] <- new.mean2
    samp.mean[myiter] <- mean(dist.vec)
    
    samp.m2[myiter] <- mean(dist.vec^2)
    theo.m2[myiter] <- p*sum((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2))
    
  }
  
  my.breaks <- seq(0,20,length.out=30)
  setwd("C:/Users/bdawk/Documents/KNN_project_output")
  if(save.fig==T){
    pdf(paste("TiTv_distance_histogram_TiTv",titv,".pdf",sep=""),height=6,width=6)
  }
  par(mfrow=c(1,1),mar=c(4.5,4.1,1.1,0.8))
  hist(long.dist.vec,breaks=my.breaks,main="",ylab="",xlab="",font.lab=2,
       cex.lab=1.5,cex.main=1.7,ylim=c(0,0.35),freq=F,cex.axis=1.4)
  abline(v=mean(long.dist.vec),lty=1,col='red',lwd=2)
  abline(v=mean(theo.mean3),lty=1,col='blue',lwd=2)
  
  abline(v=mean(long.dist.vec)+sd(long.dist.vec),lty=2,col='orange',lwd=2)
  abline(v=mean(long.dist.vec)-sd(long.dist.vec),lty=2,col='orange',lwd=2)
  
  abline(v=mean(theo.mean3)+sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
  abline(v=mean(theo.mean3)-sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
  legend("topright",c("sample mean","theoretical mean","sample SD","theoretical SD"),
         lty=c(1,1,2,2),col=c('red','blue','orange','purple'),lwd=2,bg='white',cex=1.5)
  box()
  if(save.fig==T){
    dev.off()
  }
  
}

#########################################################################################