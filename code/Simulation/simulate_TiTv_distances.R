# TiTv Distances

PuPy_probs_func <- function(g0){
  g1a <- NA
  g2a <- NA
  g1b <- NA
  g2b <- NA
  
  while(is.na(g1a) || is.na(g1b) || is.na(g2a) || is.na(g2b)){
  if(g0 <= 0 || g0 >= 1) stop("Please enter a probability in (0,1).")
  
  probs <- list(5)
  g2a <- suppressWarnings((6*(1 - g0) + sqrt(36*((1 - g0)^2) - 24*(6*g0^2 - 6*g0 + 1)))/12)
  g2b <- suppressWarnings((6*(1 - g0) - sqrt(36*((1 - g0)^2) - 24*(6*g0^2 - 6*g0 + 1)))/12)
  g1a <- 1 - g2a - g0
  g1b <- 1 - g2b - g0
  
  checka <- g0 + g1a + g2a
  checkb <- g0 + g1b + g2b
  
  #if(abs(checka-1) > 0.0001) stop(paste("Sum of probabilities is ",checka,". The sum should be 1 dummy!"))
  #if(abs(checkb-1) > 0.0001) stop(paste("Sum of probabilities is ",checkb,". The sum should be 1 dummy!"))
  
  if(g1a < 0 || is.na(g1a < 0)){
    g1a <- NA
  }
  
  if(g2a < 0 || is.na(g2a < 0)){
    g2a <- NA
  }
  
  if(g1b < 0 || is.na(g1b < 0)){
    g1b <- NA
  }
  
  if(g2b < 0 || is.na(g2b < 0)){
    g2b <- NA
  }
  
  if(is.na(g1a) || is.na(g1b) || is.na(g2a) || is.na(g2b)){
     g0 <- runif(1)
  }
  
  }
  
  probs <- list(g0,g1a,g1b,g2a,g2b)
  names(probs) <- c("g0","g1a","g1b","g2a","g2b")
  probs
}

partial_sum_func <- function(M,b){
      part.sum <- sum(M[1:(b-1)])
      return(part.sum)
}

#set.seed(2)
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
#probs <- runif(p,0.1,0.9)
g1 <- 1/3
g0 <- runif(1,min=0.0001,max=(2/3 - 0.0001))
g2 <- 2/3 - g0
#g0 <- g1 <- g2 <- 1/3;
for(myiter in 1:100){
  print(myiter)
prob0 <- runif(1)
myprobs <- PuPy_probs_func(prob0)
#g0 <- myprobs$g0
#g0
#g1 <- myprobs$g1a
#g1
#g2 <- myprobs$g2a

my.eta <- 2
g1 <- 1/(my.eta+1)
g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
g2 <- my.eta/(my.eta + 1) - g0
#g2
#g0+g1+g2
#2*(g0*g1 + g1*g2 + g0*g2)

m <- 100
p <- 100

PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0,g1,g2),replace=T)
table(PuPy.samp)/sum(table(PuPy.samp))

data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
#probs <- runif(m*p,0.1,0.9)
bounds <- runif(2,min=0.001,max=0.99)
my.min <- bounds[which.min(bounds)]
my.max <- bounds[which.max(bounds)]
probs <- runif(p,0.3,0.4)
#data.mat <- matrix(as.double(rbinom(m*p,2,probs)),p,m)
data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
for(i in 1:m){
  #probs <- runif(p,0.001,0.5)
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
  #print(i)
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
    #idx0 <- which(diff.TiTv==0)
    #idx.25 <- which(diff.TiTv==0.25)
    #idx.5 <- which(diff.TiTv==0.5)
    #idx.75 <- which(diff.TiTv==0.75)
    #idx.1 <- which(diff.TiTv==1)
    #if(length(idx0)!=p){
    #  diff.mat[iter,] <- matrix(c((length(idx0)/p),(length(idx.25)/p),(length(idx.5)/p),(length(idx.75)/p),(length(idx.1)/p)),nrow=1,ncol=5)
    #  iter <- iter + 1
    #}
  }
}
tmp <- d.mat + t(d.mat)
long.dist.vec <- c(long.dist.vec,tmp[upper.tri(tmp)])
dist.vec <- c(t(tmp))
idx <- which(dist.vec==0)
dist.vec <- dist.vec[-idx]
#hist(dist.vec,breaks=100,freq=F)  

E.diff0 <- mean((1 - probs)^4 + 4*(probs^2)*((1 - probs)^2) + probs^4)
#E.diff0
E.diff.25 <- 4*(g0 + g2)*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
#E.diff.25
E.diff.5 <- 4*g1*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
#E.diff.5
E.diff.75 <- 2*(g0 + g2)*mean(((1 - probs)^2)*(probs^2))
#E.diff.75
E.diff1 <- 2*g1*mean(((1 - probs)^2)*(probs^2))
#E.diff1

#print(colMeans(diff.mat))
#print(c(E.diff0,E.diff.25,E.diff.5,E.diff.75,E.diff1))

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
#E.D2.2 <- 0
#for(k in 2:p){
    #print(k)
#    factor1 <- (g0+g2+2*g1)*(((1-probs[k])^3)*probs[k]+(probs[k]^3)*(1-probs[k])) + ((3/2)*(g0+g2)+2*g1)*((1-probs[k])^2)*((probs[k])^2)
#    for(j in 1:(k-1)){
#      factor2 <- (g0+g2+2*g1)*(((1-probs[j])^3)*probs[j]+(probs[j]^3)*(1-probs[j])) + ((3/2)*(g0+g2)+2*g1)*((1-probs[j])^2)*((probs[j])^2)
#      E.D2.2 <- E.D2.2 + 2*factor1*factor2
#    }
#}
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
mean(theo.var1)
#mean(theo.var2)
mean(samp.var)

mean(samp.mean)
mean(theo.mean)

hist(samp.var,breaks=30)
abline(v=mean(samp.var),lwd=2,lty=2,col='red')
abline(v=mean(theo.var1),lwd=2,lty=2,col='blue')

my.iter <- 2
hist(dist.vec,breaks=100)
abline(v=samp.mean[my.iter],lty=2,col='red',lwd=2)
abline(v=theo.mean3[my.iter],lty=2,col='blue',lwd=2)
abline(v=theo.mean3[my.iter]+sqrt(theo.var3[my.iter]),lty=2,col='green',lwd=2)
abline(v=theo.mean3[my.iter]-sqrt(theo.var3[my.iter]),lty=2,col='green',lwd=2)
abline(v=samp.mean[my.iter]+sqrt(samp.var[my.iter]),lty=2,col='purple',lwd=2)
abline(v=samp.mean[my.iter]-sqrt(samp.var[my.iter]),lty=2,col='purple',lwd=2)


######
my.breaks <- seq(0,40,length.out=30)
setwd("C:/Users/bdawk/Documents/KNN_project_output/paper_figs")
pdf("TiTv_distance_histogram_maf4.pdf",height=7,width=7)
par(mfrow=c(1,1),mar=c(4.5,4.1,1.8,0.8))
hist(long.dist.vec,breaks=my.breaks,main="Histogram of TiTv Distances (Avg. MAF = 0.35)",ylab="Frequency",xlab="TiTv Distance",font.lab=2,
     cex.lab=1.5,cex.main=1.7,ylim=c(0,0.25),freq=F)
abline(v=mean(long.dist.vec),lty=1,col='red',lwd=2)
abline(v=mean(theo.mean3),lty=1,col='blue',lwd=2)

abline(v=mean(long.dist.vec)+sd(long.dist.vec),lty=2,col='orange',lwd=2)
abline(v=mean(long.dist.vec)-sd(long.dist.vec),lty=2,col='orange',lwd=2)

abline(v=mean(theo.mean3)+sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
abline(v=mean(theo.mean3)-sqrt(mean(theo.var3)),lty=2,col='purple',lwd=2)
legend("topright",c("sample mean","theoretical mean","sample SD","theoretical SD"),
       lty=c(1,1,2,2),col=c('red','blue','orange','purple'),lwd=2,bg='white',cex=1.3)
box()
dev.off()

#############################################################################################



p <- 100
etas <- c(0.1,1,2)
means.list <- list()
sds.list <- list()
mafs.list <- list()

for(out.iter in 1:length(etas)){
  my.iterates <- seq(0.01,0.9,length.out=100)
  my.eta <- etas[out.iter]
  g1 <- 1/(my.eta+1)
  g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
  g2 <- my.eta/(my.eta + 1) - g0

  my.means <- numeric()
  my.sds <- numeric()
  my.mafs <- numeric()
  for(iter in c(1:(length(my.iterates)-1))){
  
    my.idx <- c(iter,iter+1)
  
    my.mafs[iter] <- mean(c(my.iterates[my.idx[1]],my.iterates[my.idx[2]]))
  
    probs <- runif(p,my.iterates[my.idx[1]],my.iterates[my.idx[2]])

    w.v <- 2*(g0*g1 + g1*g2 + g0*g2)
    w.i <- 1 - 2*w.v
    F.a <- ((1 - probs)^3)*probs + (probs^3)*(1 - probs)
    G.a <- ((1 - probs)^2)*(probs^2)

    my.means[iter] <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
    my.sds[iter] <- sqrt((0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2))
  
  }
  means.list[[out.iter]] <- my.means
  sds.list[[out.iter]] <- my.sds
  mafs.list[[out.iter]] <- my.mafs
}

setwd("C:/Users/bdawk/Documents/KNN_project_output/paper_figs")
pdf("predicted_TiTv_distance-vs-Average_maf.pdf",height=7,width=7)
par(mfrow=c(1,1),mar=c(4.5,4.1,1.8,0.8))
plot(mafs.list[[1]],means.list[[1]],type='l',main="",
     xlab="Average MAF",ylab="TiTv Distance",font.lab=2,cex.lab=1.5,cex.main=1.7,lwd=2,ylim=c(0,40))
lines(mafs.list[[2]],means.list[[2]],type='l',col='orange',lwd=2,lty=2)
lines(mafs.list[[3]],means.list[[3]],type='l',col='blue',lwd=3,lty=3)
abline(v=0.5,lty=2,lwd=2,col='gray')
ex1 <- expression(paste(eta," = ",0.1,sep=""))
ex2 <- expression(paste(eta," = ",1,sep=""))
ex3 <- expression(paste(eta," = ",2,sep=""))
#text(x=0.56,y=39,labels=c("Max Distance"),font=2,cex=1.2)
#text(x=0.56,y=37.5,labels=c(expression(paste(f[a]," = 0.5",sep=""))),font=2,cex=1.2)
legend("bottomright",c(ex1,ex2,ex3),lty=c(1,2,3),lwd=c(2,2,3),col=c('black','orange','blue'),bg='white',cex=1.5)
dev.off()
#lines(my.mafs,my.means,type='l',col=4,lwd=2,lty=4)
#lines(my.mafs,my.means,type='l',col=5,lwd=2,lty=5)
#lines(my.mafs,my.means,type='l',col=6,lwd=2,lty=6)
#lines(my.mafs,my.means,type='l',col=7,lwd=2,lty=7)
#lines(my.mafs,my.means,type='l',col=8,lwd=2,lty=8)
#lines(my.mafs,(my.means + my.sds),type='l',col='red',lwd=2)
#lines(my.mafs,(my.means - my.sds),type='l',col='red',lwd=2)
