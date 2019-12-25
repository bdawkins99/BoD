############################################################
#
# simulations for all distance distributions with several different
# pairs for m and p
#
############################################################
#
# (1) Manhattan Distances (non-normalized diff)
#    (a) Standard Normal Data
#    (b) Uniform (0,1) Data
#
# (2) Manhattan Distances (min-max normalized diff)
#    (a) Standard Normal Data
#    (b) Uniform (0,1) Data
#
# (3) Euclidean Distances (non-normalized diff)
#    (a) Standard Normal Data
#    (b) Uniform (0,1) Data
#
# (4) Euclidean Distances (min-max normalized diff)
#    (a) Standard Normal Data
#    (b) Uniform (0,1) Data
#
############################################################

# (1a) Manhattan Distances (non-normalized diff) - Normal Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
samp.var <- numeric()
samp.mean <- numeric()
for(iter in 1:100){
  #print(iter)

p <- ps[outiter] # number of attributes

m <- ms[outiter]  # number of subjects

# data matrix
sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
for(i in 1:m){
  sub.mat[i,] <- matrix(rnorm(p,sd=1),nrow=1,ncol=p)
}

# find range for each attribute
#maxs <- apply(sub.mat,2,max) # attribute maximums
#mins <- apply(sub.mat,2,min) # attribute minimums
#diffs <- maxs-mins           # attribute ranges
#diffs.inv <- 1/diffs         # inverse of attribute ranges

# distance matrix
d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
for(i in 1:m){
  for(k in i:m){
    
    #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
    #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
    d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
    
  }
  
}

tmp <- d.mat + t(d.mat)
d.mat <- tmp
dist.vec <- c(t(d.mat))
idx <- which(dist.vec==0)
dist.vec <- dist.vec[-idx] # distance vector

# theoretical mean
theo.mean <- 2*p/sqrt(pi)
theo.mean

# sample mean
samp.mean[iter] <- mean(dist.vec)
samp.mean

# theoretical variance
theo.var <- (2*p*(pi - 2))/pi
theo.var

# sample variance
samp.var[iter] <- var(dist.vec)
samp.var
}

samp.mean.vec[outiter] <- mean(samp.mean)
samp.var.vec[outiter] <- mean(samp.var)
}

mean(samp.mean)
mean(samp.var)

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
    print(i)
    p <- ps[i]
    m <- ms[i]
    theo.mean.vec[i] <- 2*p/sqrt(pi)
    theo.var.vec[i] <- (2*p*(pi - 2))/pi
}
theo.mean <- 2*p/sqrt(pi)
theo.mean

theo.var <- (2*p*(pi - 2))/pi
theo.var

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# (1b) Manhattan Distances (non-normalized diff) - Uniform Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)

p <- ps[outiter] # number of attributes

m <- ms[outiter]  # number of subjects

# data matrix
sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
for(i in 1:m){
  sub.mat[i,] <- matrix(runif(p),nrow=1,ncol=p)
}

# find range for each attribute
#maxs <- apply(sub.mat,2,max) # attribute maximums
#mins <- apply(sub.mat,2,min) # attribute minimums
#diffs <- maxs-mins           # attribute ranges
#diffs.inv <- 1/diffs         # inverse of attribute ranges

# distance matrix
d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
for(i in 1:m){
  for(k in i:m){
    
    #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
    #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
    d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
    
  }
  
}

tmp <- d.mat + t(d.mat)
d.mat <- tmp
dist.vec <- c(t(d.mat))
idx <- which(dist.vec==0)
dist.vec <- dist.vec[-idx] # distance vector

# theoretical mean
theo.mean[iter] <- p/3
theo.mean

# sample mean
samp.mean[iter] <- mean(dist.vec)
samp.mean

# theoretical variance
theo.var <- p/18
theo.var

# sample variance
samp.var[iter] <- var(dist.vec)
samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  theo.mean.vec[i] <- p/3
  theo.var.vec[i] <- p/18
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# (2a) Manhattan Distances (min-max normalized diff) - Normal Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:10){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(rnorm(p,sd=1),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    maxs <- apply(sub.mat,2,max) # attribute maximums
    mins <- apply(sub.mat,2,min) # attribute minimums
    diffs <- maxs-mins           # attribute ranges
    diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    # expected range of attribute
    mu.m <- log(log(2))/qnorm(1/m) - qnorm(1/m)
    
    # theoretical mean
    theo.mean <- p/(sqrt(pi)*mu.m)
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- (p*(pi - 2))/(2*pi*mu.m^2)
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec1 <- numeric()
theo.mean.vec2 <- numeric()
theo.var.vec1 <- numeric()
theo.var.vec2 <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  mu.m <- log(log(2))/qnorm(1/m) - qnorm(1/m)
  mu.m2 <- qnorm(1 - (1/m)*exp(digamma(1)))
  EX2 <- (pi^2 + 24*(mu.m^2)*log(m))/(6*log(m))
  theo.mean.vec1[i] <- p/(sqrt(pi)*mu.m)
  theo.mean.vec2[i] <- p/(sqrt(pi)*mu.m2)
  theo.var.vec1[i] <- (p*(pi - 2))/(2*pi*mu.m^2)
  theo.var.vec2[i] <- (2*p*(pi-2)/pi)/EX2
  #theo.var.vec2[i] <- (p*(pi - 2))/(2*pi*mu.m2^2)
  
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean1=theo.mean.vec1,and3=ands,
                     theo.mean2=theo.mean.vec2,and4=ands,samp.var=samp.var.vec,and5=ands,theo.var1=theo.var.vec1,
                     and6=ands,theo.var2=theo.var.vec2,dash=dashes)
out.df

# (2b) Manhattan Distances (min-max normalized diff) - Uniform Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(runif(p),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    maxs <- apply(sub.mat,2,max) # attribute maximums
    mins <- apply(sub.mat,2,min) # attribute minimums
    diffs <- maxs-mins           # attribute ranges
    diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    # theoretical mean
    theo.mean <- ((m+1)*p)/(3*(m-1))
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- ((m^5 - 18*m^2 - 5*m + 2)*p)/(18*(m^3 + m^2 + 2)*((m-1)^2))
    theo.var <- (((m+1)^2)*(m^3 - 7*m + 2)*p)/(18*(m^3 - m + 2)*((m-1)^2))
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec1 <- numeric()
theo.var.vec2 <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  theo.mean.vec[i] <- ((m+1)*p)/(3*(m-1))
  theo.var.vec1[i] <- ((m^5 - 18*m^2 - 5*m + 2)*p)/(18*(m^3 + m^2 + 2)*((m-1)^2))
  theo.var.vec2[i] <- (((m+1)^2)*(m^3 - 7*m + 2)*p)/(18*(m^3 - m + 2)*((m-1)^2))
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var1=theo.var.vec1,and5=ands,
                     theo.var2=theo.var.vec2,dash=dashes)

# (3a) Euclidean Distances (non-normalized diff) - Normal Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(rnorm(p,sd=1),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    #maxs <- apply(sub.mat,2,max) # attribute maximums
    #mins <- apply(sub.mat,2,min) # attribute minimums
    #diffs <- maxs-mins           # attribute ranges
    #diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)))
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    # theoretical mean
    theo.mean <- sqrt(2*p - 1)
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- 1
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  theo.mean.vec[i] <- sqrt(2*p-1)
  theo.var.vec[i] <- 1
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# (3b) Euclidean Distances (non-normalized diff) - Uniform Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(runif(p),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    #maxs <- apply(sub.mat,2,max) # attribute maximums
    #mins <- apply(sub.mat,2,min) # attribute minimums
    #diffs <- maxs-mins           # attribute ranges
    #diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)))
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    # theoretical mean
    theo.mean <- sqrt((p/6) - (7/120))
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- 7/120
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  theo.mean.vec[i] <- sqrt((p/6) - (7/120))
  theo.var.vec[i] <- 7/120
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# (4a) Euclidean Distances (min-max normalized diff) - Normal Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(rnorm(p,sd=1),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    maxs <- apply(sub.mat,2,max) # attribute maximums
    mins <- apply(sub.mat,2,min) # attribute minimums
    diffs <- maxs-mins           # attribute ranges
    diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)))
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    mu.m <- log(log(2))/qnorm(1/m) - qnorm(1/m)
    
    # theoretical mean
    theo.mean <- sqrt(2*p - 1)/(2*mu.m)
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- (3*log(m))/(pi^2 + 12*(mu.m^2)*log(m))
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  mu.m <- log(log(2))/qnorm(1/m) - qnorm(1/m)
  
  theo.mean.vec[i] <- sqrt(2*p - 1)/(2*mu.m)
  theo.var.vec[i] <- (3*log(m))/(pi^2 + 12*(mu.m^2)*log(m))
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# (4b) Euclidean Distances (min-max normalized diff) - uniform Data

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:100){
    #print(iter)
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    # data matrix
    sub.mat <- matrix(rep(0,length=m*p),nrow=m,ncol=p)
    for(i in 1:m){
      sub.mat[i,] <- matrix(runif(p),nrow=1,ncol=p)
    }
    
    # find range for each attribute
    maxs <- apply(sub.mat,2,max) # attribute maximums
    mins <- apply(sub.mat,2,min) # attribute minimums
    diffs <- maxs-mins           # attribute ranges
    diffs.inv <- 1/diffs         # inverse of attribute ranges
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m,byrow=TRUE)
    for(i in 1:m){
      for(k in i:m){
        
        d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)*diffs.inv^2))
        #d.mat[i,k] <- sqrt(sum(((sub.mat[i,]-sub.mat[k,])^2)))
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,])*diffs.inv)
        #d.mat[i,k] <- sum(abs(sub.mat[i,]-sub.mat[k,]))
        
      }
      
    }
    
    tmp <- d.mat + t(d.mat)
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx] # distance vector
    
    # theoretical mean
    theo.mean <- (sqrt((p/6) - (7/120))*(m+1))/(m-1)
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- (7*((m+1)^2)*(m+2))/(120*(m^3 + m^2 + 2))
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]

  theo.mean.vec[i] <- (sqrt((p/6) - (7/120))*(m+1))/(m-1)
  theo.var.vec[i] <- (7*((m+1)^2)*(m+2))/(120*(m^3 + m^2 + 2))
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)

# GWAS data

# diffAM

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
w0 <- 0.3
w1 <- 0.3
w2 <- 0.4
N <- 10
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
probs.list <- list()
for(outiter in 1:length(ms)){
  print(outiter)
  
  p <- ps[outiter] # number of attributes
  
  m <- ms[outiter]  # number of subjects
  
  probs <- runif(p,0.1,0.9)
  
  samp.var <- numeric()
  samp.mean <- numeric()
  probs.mat.tmp <- matrix(rep(0,length=N*p),nrow=p,ncol=N)
  for(iter in 1:N){
    cat("Iteration: ",iter,"\n")
    
    #p <- ps[outiter] # number of attributes
    
    #m <- ms[outiter]  # number of subjects
    
    #data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    #for(i in 1:m){
    #  samp <- sample(c(0,1,2),size=p,prob=c(w0,w1,w2),replace=T)
    #  data.mat[,i] <- matrix(samp,nrow=p,ncol=1)
    #}
    
    #probs <- runif(p,0.1,0.9)
    probs.mat.tmp[,iter] <- probs
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    for(i in 1:m){
      samp <- matrix(as.double(rbinom(p,2,probs)),nrow=p,ncol=1)
      data.mat[,i] <- samp
    }
    
    #diff.func <- function(x){
    #  if(x==0){
    #    f <- 0
    #  }else if(x==1){
    #    f <- 0.5
    #  }else{
    #    f <- 1
    #  }
    #}
    
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    place <- 1
    for(i in 1:m){
      for(k in i:m){
        #my.diff <- apply(matrix(abs(data.mat[,i] - data.mat[,k]),ncol=1),1,diff.func)
        my.diff <- abs(data.mat[,i] - data.mat[,k])
        place <- place + 1 
        d.mat[i,k] <- sum(0.5*my.diff)
      }
      
    }
    
    tmp <- t(d.mat) + d.mat
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
  probs.list[[outiter]] <- rowMeans(probs.mat.tmp)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  q <- probs.list[[i]]
  #p.diff1 <- 2*mean(2*((1 - q)^3)*q + 2*(q^3)*(1 - q) + ((1 - q)^2)*(q^2))
  term1 <- sum(((1-q)^3)*q + (q^3)*(1-q) + 2*((1-q)^2)*(q^2))
  #print(term1)
  
  term2 <- 0
  for(k in 2:p){
    
    myfactor <- sum(((1-q[1:(k-1)])^3)*q[1:(k-1)] + (q[1:(k-1)]^3)*(1-q[1:(k-1)]) + ((1-q[1:(k-1)])^2)*(q[1:(k-1)]^2))
    term2 <- term2 + 8*(((1-q[k])^3)*q[k] + (q[k]^3)*(1-q[k]) + ((1-q[k])^2)*(q[k]^2))*myfactor
    
  }
  #print(term2)
  
  theo.mean.vec[i] <- 2*sum(((1-q)^3)*q + (q^3)*(1-q) + ((1-q)^2)*(q^2))
  theo.var.vec[i] <- term1 + term2 - theo.mean.vec[i]^2
  
  t1 <- ((1-q)^3)*q + (q^3)*(1-q) + 2*((1-q)^2)*(q^2)
  t2 <- (((1-q)^3)*q + (q^3)*(1-q) + ((1-q)^2)*(q^2))^2
  print(sum(t1 - 4*t2))
  #theo.mean.vec[i] <- p*(w1*(w0 + w2) + 2*w0*w2)
  #theo.var.vec[i] <- (0.5*w1*(w0 + w2) + 2*w0*w2)*p - p*(w1*(w0 + w2) + 2*w0*w2)^2
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)
out.df

# diffGM

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
w0 <- 0.3
w1 <- 0.3
w2 <- 0.4
N <- 10
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
probs.list <- list()
for(outiter in 1:length(ms)){
  print(outiter)
  
  p <- ps[outiter] # number of attributes
  
  m <- ms[outiter]  # number of subjects
  
  probs <- runif(p,0.1,0.9)
  
  samp.var <- numeric()
  samp.mean <- numeric()
  probs.mat.tmp <- matrix(rep(0,length=N*p),nrow=p,ncol=N)
  for(iter in 1:N){
    #print(iter)
    cat("Iteration: ",iter,"\n")
    
    #p <- ps[outiter] # number of attributes
    
    #m <- ms[outiter]  # number of subjects
    
    #data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    #for(i in 1:m){
    #  samp <- sample(c(0,1,2),size=p,prob=c(w0,w1,w2),replace=T)
    #  data.mat[,i] <- matrix(samp,nrow=p,ncol=1)
    #}
    
    #probs <- runif(p,0.1,0.9)
    probs.mat.tmp[,iter] <- probs
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    for(i in 1:m){
      samp <- matrix(as.double(rbinom(p,2,probs)),nrow=p,ncol=1)
      data.mat[,i] <- samp
    }
    
    #diff.func <- function(x){
    #  if(x==0){
    #    f <- 0
    #  }else{
    #    f <- 1
    #  }
    #}
    
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    place <- 1
    for(i in 1:m){
      for(k in i:m){
        #my.diff <- apply(matrix(abs(data.mat[,i] - data.mat[,k]),ncol=1),1,diff.func)
        my.diff <- abs(data.mat[,i] - data.mat[,k])
        idx <- which(my.diff==2)
        my.diff[idx] <- 1
        place <- place + 1 
        d.mat[i,k] <- sum(my.diff)
      }
      
    }
    
    tmp <- t(d.mat) + d.mat
    d.mat <- tmp
    dist.vec <- c(t(d.mat))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
  probs.list[[outiter]] <- rowMeans(probs.mat.tmp)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(probs.list)){
#for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  q <- probs.list[[i]]
  p.diff1 <- 2*mean(2*((1 - q)^3)*q + 2*(q^3)*(1 - q) + ((1 - q)^2)*(q^2))
  term1 <- 2*sum(2*((1-q)^3)*q + 2*(q^3)*(1-q) + ((1-q)^2)*(q^2))
  
  term2 <- 0
  for(k in 2:p){
      
        myfactor <- sum(2*((1-q[1:(k-1)])^3)*q[1:(k-1)] + 2*(q[1:(k-1)]^3)*(1-q[1:(k-1)]) + ((1-q[1:(k-1)])^2)*(q[1:(k-1)]^2))
        term2 <- term2 + 8*(2*((1-q[k])^3)*q[k] + 2*(q[k]^3)*(1-q[k]) + ((1-q[k])^2)*(q[k]^2))*myfactor
  
  }
  w0 <- mean(((1-q)^2))
  w1 <- 2*mean(q*(1-q))
  w2 <- mean(q^2)
  
  #theo.mean.vec[i] <- 2*sum((2*(1 - q)^3)*q + 2*(q^3)*(1 - q) + 0.33*((1 - q)^2)*q^2)
  theo.mean.vec[i] <- p.diff1*p
  #theo.var.vec[i] <- theo.mean.vec[i] - (4/p)*(0.5*theo.mean.vec[i])^2
  theo.var.vec[i] <- term1 + term2 - theo.mean.vec[i]^2
  #theo.mean.vec[i] <- 2*p*(w0*w1 + w1*w2 + w0*w2)
  #theo.var.vec[i] <- 2*p*(w0*w1 + w1*w2 + w0*w2)*(1 - 2*(w0*w1 + w1*w2 + w0*w2))
  fact1 <- 2*((1-q)^3)*q + 2*(q^3)*(1-q) + ((1-q)^2)*(q^2)
  fact2 <- 1 - 2*(2*((1-q)^3)*q + 2*(q^3)*(1-q) + ((1-q)^2)*(q^2))
  print(2*sum(fact1*fact2))
}

ands <- data.frame(and=rep("&",length=length(ms[1:5])))
dashes <- data.frame(dash=rep("\\",length=length(ms[1:5])))

dims <- paste("$(",ms[1:5],",",ps[1:5],")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)
out.df

#############################################################################
#############################################################################

w0 <- 0.3
w1 <- 0.3
w2 <- 0.4

p <- 1000 # number of attributes

m <- 100  # number of subjects

data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
for(i in 1:m){
  samp <- sample(c(0,1,2),size=p,prob=c(w0,w1,w2),replace=T)
  data.mat[,i] <- matrix(samp,nrow=p,ncol=1)
}

data.vec <- c(data.mat)
idx0 <- which(data.vec==0)
idx1 <- which(data.vec==1)
idx2 <- which(data.vec==2)

w0bar <- length(idx0)/(m*p)
w0bar

w1bar <- length(idx1)/(m*p)
w1bar

w2bar <- length(idx2)/(m*p)
w2bar

mean(data.vec)

p <- 1000
m <- 100

PuPy_probs_func <- function(g0){
    
    if(g0 <= 0 || g0 >= 1) stop("Please enter a probability in (0,1).")
    
    probs <- list(5)
    g2a <- (6*(1 - g0) + sqrt(36*((1 - g0)^2) - 24*(6*g0^2 - 6*g0 + 1)))/12
    g2b <- (6*(1 - g0) - sqrt(36*((1 - g0)^2) - 24*(6*g0^2 - 6*g0 + 1)))/12
    g1a <- 1 - g2a - g0
    g1b <- 1 - g2b - g0
    
    checka <- g0 + g1a + g2a
    checkb <- g0 + g1b + g2b
     
    if(abs(checka-1) > 0.0001) stop(paste("Sum of probabilities is ",checka,". The sum should be 1 dummy!"))
    if(abs(checkb-1) > 0.0001) stop(paste("Sum of probabilities is ",checkb,". The sum should be 1 dummy!"))
    
    probs <- list(g0,g1a,g1b,g2a,g2b)
    names(probs) <- c("g0","g1a","g1b","g2a","g2b")
    probs
}

myprobs <- PuPy_probs_func(0.2)

P.Tva <- 2*(myprobs$g0*myprobs$g1a + myprobs$g1a*myprobs$g2a + myprobs$g0*myprobs$g2a)
P.Tvb <- 2*(myprobs$g0*myprobs$g1b + myprobs$g1b*myprobs$g2b + myprobs$g0*myprobs$g2b)
P.Tia <- 1 - P.Tva
P.Tib <- 1 - P.Tvb

#probs <- runif(m*p,0.1,0.9)
probs <- runif(p,0.1,0.9)
#data.mat <- matrix(as.double(rbinom(m*p,2,probs)),p,m)
data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
for(i in 1:m){
    samp <- matrix(as.double(rbinom(p,2,runif(p,0.1,0.9))),nrow=p,ncol=1)
    data.mat[,i] <- samp
}

MAs <- sample(c("A","T","C","G"),size=p,replace=T)
SAs <- numeric()
for(i in 1:length(MAs)){
    samp.idx <- c("A","T","C","G") %in% MAs[i]
    samp.idx <- !samp.idx
    samp <- c("A","T","C","G")[samp.idx]
    SAs[i] <- sample(samp,size=1,replace=T)
}
probsSAs <- 1 - probs

PuPy <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
for(i in 1:p){
    print(i)
    MA <- MAs[i]
    SA <- SAs[i]
    if(MA=="A" && SA=="G"){
       PuPy[i,] <- matrix(rep(3,length=m),nrow=1,ncol=m)
    }else if(MA=="A" && SA=="T"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 5
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 3
          }
       }
    }else if(MA=="A" && SA=="C"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 5
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 3
          }
       }
    }else if(MA=="T" && SA=="C"){
       PuPy[i,] <- matrix(rep(5,length=m),nrow=1,ncol=m)
    }else if(MA=="T" && SA=="A"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 3
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 5
          }
       }
    }else if(MA=="T" && SA=="G"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 3
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 5
          }
       }
    }else if(MA=="C" && SA=="T"){
       PuPy[i,] <- matrix(rep(5,length=m),nrow=1,ncol=m)
    }else if(MA=="C" && SA=="A"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 3
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 5
          }
       }
    }else if(MA=="C" && SA=="G"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 3
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 5
          }
       }
    }else if(MA=="G" && SA=="A"){
       PuPy[i,] <- matrix(rep(3,length=m),nrow=1,ncol=m)
    }else if(MA=="G" && SA=="T"){
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 5
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 3
          }
       }
    }else{
       for(j in 1:m){
          if(data.mat[i,j]==0){
            PuPy[i,j] <- 5
          }else if(data.mat[i,j]==1){
            PuPy[i,j] <- 4
          }else{
            PuPy[i,j] <- 3
          }
       }
    }
}

diff.mat <- matrix(rep(0,length=4*(m*(m-1)/2)),nrow=(m*(m-1)/2),ncol=5)
d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
iter <- 1
for(i in 1:m){
  print(i)
    for(j in i:m){
        diff <- abs(data.mat[,i]-data.mat[,j])
        diff.TiTv <- diff
        gens <- as.matrix(cbind(PuPy[,i],PuPy[,j]))
        
           for(k in 1:p){
             if(diff[k]==1){
                if(gens[k,1]==3 && gens[k,2]==3){
                   diff.TiTv[k] <- 0.25
                }else if(gens[k,1]==5 && gens[k,2]==5){
                   diff.TiTv[k] <- 0.25
                }else if(gens[k,1]==3 && gens[k,2]==4){
                   diff.TiTv[k] <- 0.5
                }else if(gens[k,1]==4 && gens[k,2]==3){
                   diff.TiTv[k] <- 0.5
                }else if(gens[k,1]==4 && gens[k,2]==5){
                   diff.TiTv[k] <- 0.5
                }else if(gens[k,1]==5 && gens[k,2]==4){
                   diff.TiTv[k] <- 0.5
                }
             }else if(diff[k]==2){
                if(gens[k,1]==3 && gens[k,2]==3){
                   diff.TiTv[k] <- 0.75
                }else if(gens[k,1]==5 && gens[k,2]==5){
                   diff.TiTv[k] <- 0.75
                }else if(gens[k,1]==3 && gens[k,2]==5){
                   diff.TiTv[k] <- 1
                }else if(gens[k,1]==5 && gens[k,2]==3){
                   diff.TiTv[k] <- 1
                }
             }
           }
        d.mat[i,j] <- sum(diff.TiTv)
        idx0 <- which(diff.TiTv==0)
        idx.25 <- which(diff.TiTv==0.25)
        idx.5 <- which(diff.TiTv==0.5)
        idx.75 <- which(diff.TiTv==0.75)
        idx.1 <- which(diff.TiTv==1)
        if(length(idx0)!=p){
        diff.mat[iter,] <- matrix(c((length(idx0)/p),(length(idx.25)/p),(length(idx.5)/p),(length(idx.75)/p),(length(idx.1)/p)),nrow=1,ncol=5)
        iter <- iter + 1
        }
    }
}

tmp <- d.mat
tmp <- tmp + t(tmp)

d.vec <- c(t(tmp))
idx <- which(d.vec==0)
d.vec <- d.vec[-idx]
hist(d.vec,breaks=100,freq=F)

PuPy.vec <- c(PuPy)
idx3 <- which(PuPy.vec==3)
idx4 <- which(PuPy.vec==4)
idx5 <- which(PuPy.vec==5)

g0bar <- length(idx3)/(m*p)
g0bar

g1bar <- length(idx4)/(m*p)
g1bar

g2bar <- length(idx5)/(m*p)
g2bar

#data.mat <- matrix(as.double(rbinom(m*p,2,probs)),p,m)
#data.mat <- matrix(as.double(rbinom(m*p,2,rep(0.5,length=m*p))),p,m)
#probs <- runif(m*p,0.1,0.9)
mean(probs)

data.vec <- c(data.mat)
idx0 <- which(data.vec==0)
idx1 <- which(data.vec==1)
idx2 <- which(data.vec==2)

w0bar <- length(idx0)/(m*p)
w0bar
mean((1-probs)^2)

w1bar <- length(idx1)/(m*p)
w1bar
2*mean(probs*(1-probs))

w2bar <- length(idx2)/(m*p)
w2bar
mean(probs^2)

mean(data.vec)

P.diff0 <- w0bar^2 + w1bar^2 + w2bar^2                                           # P(diff=0)
P.diff.25 <- 2*(g0bar^2+g2bar^2)*(w0bar*w1bar + w1bar*w2bar + w0bar*w2bar)      # P(diff=0.25)
P.diff.5 <- 4*(w0bar*w1bar + w1bar*w2bar)*(g0bar*g1bar + g1bar*g2bar + g0bar*g2bar) # P(diff=0.5)
P.diff.75 <- 2*w0bar*w2bar*(g0bar^2+g2bar^2)                                     # P(diff=0.75)
P.diff1 <- 4*w0bar*w2bar*(g0bar*g1bar+g1bar*g2bar+g0bar*g2bar)                  # P(diff=1)
P.diff0 + P.diff.25 + P.diff.5 + P.diff.75 + P.diff1

colMeans(diff.mat)
matrix(c(P.diff0,P.diff.25,P.diff.5,P.diff.75,P.diff1),nrow=1)

mean(d.vec)
#(0.5*(w0bar*w1bar + w1bar*w2bar + 3*w0bar*w2bar)*(g0bar^2+g2bar^2) + 2*(w0bar*w1bar + w1bar*w2bar+2*w0bar*w2bar)*(g0bar*g1bar+g1bar*g2bar+g0bar*g2bar))*p
E.Dij <- (0*P.diff0+0.25*P.diff.25+0.5*P.diff.5+0.75*P.diff.75+P.diff1)*p
E.Dij
abline(v=c(mean(d.vec),E.Dij),lty=2,col=c('red','blue'),lwd=2)

exp <- (0.5*(w0bar*w1bar + w1bar*w2bar + 4*w0bar*w2bar)*(g0bar^2 + g2bar^2) + 2*(w0bar*w1bar + w1bar*w2bar + 2*w0bar*w2bar)*(g0bar*g1bar + g1bar*g2bar + g0bar*g2bar))*p
exp
V.Dij <- ((1/8)*(w0bar*w1bar + w1bar*w2bar + 10*w0bar*w2bar)*(g0bar^2 + g2bar^2) + (w0bar*w1bar + w1bar*w2bar + 4*w0bar*w2bar)*(g0bar*g1bar + g1bar*g2bar + g0bar*g2bar))*p - ((0.5*(w0bar*w1bar + w1bar*w2bar + 4*w0bar*w2bar)*(g0bar^2 + g2bar^2) + 2*(w0bar*w1bar + w1bar*w2bar + 2*w0bar*w2bar)*(g0bar*g1bar + g1bar*g2bar + g0bar*g2bar))^2)*p
V.Dij

#############################################################################
#############################################################################

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

ps <- rep(c(seq(1000,5000,by=1000)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
#ps <- rep(c(seq(1000,5000,by=1000)),length=5)
#ms <- c(rep(100,length=5))

samp.mean.vec <- numeric()
samp.var.vec <- numeric()
w0 <- numeric()
w1 <- numeric()
w2 <- numeric()
g0 <- numeric()
g1 <- numeric()
g2 <- numeric()
E.diff0 <- numeric()
E.diff.25 <- numeric()
E.diff.5 <- numeric()
E.diff.75 <- numeric()
E.diff1 <- numeric()
E.Dsqij <- numeric()
theo.var.vec <- numeric()
N <- 10
probs.list <- list()
gs.list <- list()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  w0.vec <- numeric()
  w1.vec <- numeric()
  w2.vec <- numeric()
  g0.vec <- numeric()
  g1.vec <- numeric()
  g2.vec <- numeric()
  E.diff0.vec <- numeric()
  E.diff.25.vec <- numeric()
  E.diff.5.vec <- numeric()
  E.diff.75.vec <- numeric()
  E.diff1.vec <- numeric()
  E.Dsqij.vec <- numeric()
  theo.var <- numeric()
  
  p <- ps[outiter] # number of attributes
  
  m <- ms[outiter]  # number of subjects
  
  probs <- runif(p,0.1,0.9)
  
  prob0 <- runif(1)
  myprobs <- PuPy_probs_func(prob0)
  g0.tmp <- myprobs$g0
  g1.tmp <- myprobs$g1a
  g2.tmp <- myprobs$g2a
  g0.tmp
  g1.tmp
  g2.tmp
  #gs.list[[outiter]] <- c(g0.tmp,g1.tmp,g2.tmp)
  #PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0.tmp,g1.tmp,g2.tmp),replace=T)
  
  probs.mat.tmp <- matrix(rep(0,length=N*p),nrow=p,ncol=N)
  for(iter in 1:N){
    #print(iter)
    cat("Iteration: ",iter,"\n")
    
    #p <- ps[outiter] # number of attributes
    
    #m <- ms[outiter]  # number of subjects
    
    #prob0 <- runif(1)
    #myprobs <- PuPy_probs_func(prob0)
    #g0.vec[iter] <- myprobs$g0
    #g1.vec[iter] <- myprobs$g1a
    #g2.vec[iter] <- myprobs$g2a
    g0.vec[iter] <- g0.tmp
    g1.vec[iter] <- g1.tmp
    g2.vec[iter] <- g2.tmp
    
    PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0.vec[iter],g1.vec[iter],g2.vec[iter]),replace=T)
    
    #probs <- runif(p,0.1,0.9)
    probs.mat.tmp[,iter] <- probs
    data.mat <- matrix(rep(0,length=m*p),nrow=p,ncol=m)
    for(i in 1:m){
      samp <- matrix(as.double(rbinom(p,2,probs)),nrow=p,ncol=1)
      data.mat[,i] <- samp
    }
    
    w0.vec[iter] <- mean((1-probs)^2)
    w1.vec[iter] <- 2*mean(probs*(1-probs))
    w2.vec[iter] <- mean(probs^2)
    
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
      }
    }
    
    tmp <- d.mat + t(d.mat)
    dist.vec <- c(t(tmp))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    
    #E.diff0.vec[iter] <- mean((1 - probs)^4 + 4*(probs^2)*((1 - probs)^2) + probs^4)
    #E.diff.25.vec[iter] <- 4*(g0.vec[iter] + g2.vec[iter])*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
    #E.diff.5.vec[iter] <- 4*g1.vec[iter]*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
    #E.diff.75.vec[iter] <- 2*(g0.vec[iter] + g2.vec[iter])*mean(((1 - probs)^2)*(probs^2))
    #E.diff1.vec[iter] <- 2*g1.vec[iter]*mean(((1 - probs)^2)*(probs^2))
    
    #E.Dij <- (0*E.diff0.vec[iter] + 0.25*E.diff.25.vec[iter] + 0.5*E.diff.5.vec[iter] + 0.75*E.diff.75.vec[iter] + 1*E.diff1.vec[iter])*p
    
    #E.D2.1 <- (0.25*(g0.vec[iter]+g2.vec[iter])+g1.vec[iter])*sum(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((9/8)*(g0.vec[iter]+g2.vec[iter])+2*g1.vec[iter])*sum(((1-probs)^2)*(probs^2))
    
    #v <- (g0.vec[iter]+g2.vec[iter]+2*g1.vec[iter])*(((1-probs)^3)*probs+(probs^3)*(1-probs)) + ((3/2)*(g0.vec[iter]+g2.vec[iter])+2*g1.vec[iter])*((1-probs)^2)*((probs)^2)
    #v1 <- (g0.vec[iter]+g2.vec[iter]+2*g1.vec[iter])*(((1-probs[-1])^3)*probs[-1]+(probs[-1]^3)*(1-probs[-1])) + ((3/2)*(g0.vec[iter]+g2.vec[iter])+2*g1.vec[iter])*((1-probs[-1])^2)*((probs[-1])^2)
    #v2 <- numeric()
    #for(k in 2:p){
    #  v2[(k-1)] <- partial_sum_func(v,k)
    #}
    
    #E.D2.2 <- c(2*matrix(v1,nrow=1,ncol=length(v1)) %*% matrix(v2,nrow=length(v2),ncol=1))
    
    #E.Dsqij.vec[iter] <- E.D2.1 + E.D2.2
    #theo.var[iter] <- E.Dsqij.vec[iter] - (E.Dij)^2
  }
  #w0[outiter] <- mean(w0.vec)
  #w1[outiter] <- mean(w1.vec)
  #w2[outiter] <- mean(w2.vec)
  g0[outiter] <- mean(g0.vec)
  g1[outiter] <- mean(g1.vec)
  g2[outiter] <- mean(g2.vec)
  
  #E.diff0[outiter] <- mean(E.diff0.vec)
  #E.diff.25[outiter] <- mean(E.diff.25.vec)
  #E.diff.5[outiter] <- mean(E.diff.5.vec)
  #E.diff.75[outiter] <- mean(E.diff.75.vec)
  #E.diff1[outiter] <- mean(E.diff1.vec)
  #E.Dsqij[outiter] <- mean(E.Dsqij.vec)
  #theo.var.vec[outiter] <- mean(theo.var)
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
  probs.list[[outiter]] <- rowMeans(probs.mat.tmp)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  q <- probs.list[[i]]
  g0.tmp <- g0[i]
  g1.tmp <- g1[i]
  g2.tmp <- g2[i]
  
  #E.diff.25 <- 4*(g0.tmp + g2.tmp)*mean(((1-q)^3)*q + (q^3)*(1-q))
  #E.diff.5 <- 4*g1.tmp*mean(((1-q)^3)*q + (q^3)*(1-q))
  #E.diff.75 <- 2*(g0.tmp + g2.tmp)*mean(((1-q)^2)*(q^2))
  #E.diff1 <- 2*g1.tmp*mean(((1-q)^2)*(q^2))
  #print((0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + E.diff1)*p)
  
  term1 <- (g0.tmp + g2.tmp + 2*g1.tmp)*sum(((1-q)^3)*q + (q^3)*(1-q))
  term2 <- (1.5*(g0.tmp + g2.tmp) + 2*g1.tmp)*sum(((1-q)^2)*(q^2))
  theo.mean.vec[i] <- term1 + term2
  
  term1 <- (0.25*(g0.tmp + g2.tmp) + g1.tmp)*sum(((1-q)^3)*q + (q^3)*(1-q))
  term2 <- ((9/8)*(g0.tmp + g2.tmp) + 2*g1.tmp)*sum(((1-q)^2)*(q^2))
  term3 <- sum(((g0.tmp + g2.tmp + 2*g1.tmp)*(((1-q)^3)*q + (q^3)*(1-q)) + (1.5*(g0.tmp + g2.tmp) + 2*g1.tmp)*((1-q)^2)*(q^2))^2)
  
  theo.var.vec[i] <- term1 + term2 - term3
  #theo.mean.vec[i] <- (0*E.diff0[i] + 0.25*E.diff.25[i] + 0.5*E.diff.5[i] + 0.75*E.diff.75[i] + 1*E.diff1[i])*p
  #E.Dsqij <- (0*E.diff0[i] + (1/16)*E.diff.25[i] + 0.25*E.diff.5[i] + (9/16)*E.diff.75[i] + 1*E.diff1[i])*p + (p^2 - p)*(0*E.diff0[i] + 0.25*E.diff.25[i] + 0.5*E.diff.5[i] + 0.75*E.diff.75[i] + 1*E.diff1[i])^2
  
  #theo.var.vec[i] <- E.Dsqij[i] - theo.mean.vec[i]^2
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)
out.df

#######################################################################
#######################################################################

# rs-fMRI correlation metric distances

myfisher <- function(rho){
  
  z <- 0.5*log((1+rho)/(1-rho))
  return(z)
}

stretch_mat <- function(M){
  
  mat <- numeric()
  for(k in 1:nrow(M)){
    mat <- c(mat,M[k,-k])
  }
  return(mat)
}

library(Matrix) 

ps <- rep(c(seq(100,300,by=50)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:5){
    cat("Iteration: ",iter,"\n")
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    data.mat <- matrix(rep(0,length=m*(p*(p-1))),nrow=(p*(p-1)),ncol=m)
    for(k in 1:m){
      #print(k)
      
      p <- 1000
      X <- matrix(runif(p^2,min=-1,max=1), ncol=p) 
      
      cov <- X %*% t(X) 
      
      mycorr <- cov2cor(cov) 
      
      #zcorr <- apply(matrix(mycorr[upper.tri(mycorr)],ncol=1),1,myfisher)
      zcorr <- apply(matrix(stretch_mat(mycorr),ncol=1),1,myfisher)
      #hist(zcorr1,breaks=100,freq=F)
      #data.mat[,k] <- matrix(zcorr,nrow=(p*(p-1)/2),ncol=1)
      data.mat[,k] <- matrix(zcorr,nrow=(p*(p-1)),ncol=1)
    }
    
    #maxs <- numeric()
    #mins <- numeric()
    #diffs <- numeric()
    #diffs.inv <- numeric()
    #for(i in 1:p){
      #print(i)
    #  idx1 <- (i-1)*(p-1)+1
    #  idx2 <- i*(p-1)
    #  dats <- c(data.mat[c(idx1:idx2),])
    #  maxs[i] <- max(dats)
    #  mins[i] <- min(dats)
    #  diffs[i] <- max(dats)-min(dats)
    #  diffs.inv[i] <- 1/diffs[i]
    #}
    
    #diffs.inv.vec <- numeric()
    #for(i in 1:length(diffs.inv)){
    #  #print(i)
    #  idx1 <- (i-1)*(p-1)+1
    #  idx2 <- i*(p-1)
    #  diffs.inv.vec[c(idx1:idx2)] <- rep(diffs.inv[i],length=(p-1))
    #}
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    for(i in 1:m){
      for(j in i:m){
        diff <- abs(data.mat[,i]-data.mat[,j])
        #diff <- abs(data.mat[,i]-data.mat[,j])*diffs.inv.vec
        d.mat[i,j] <- sum(diff)
        }
    }
    
    tmp <- d.mat + t(d.mat)
    dist.vec <- c(t(tmp))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    #mu.m <- (-qnorm(1/(m*(p-1))) + digamma(1)/qnorm(1/(m*(p-1))))*(1/sqrt(m-3))
    
    # theoretical mean
    #theo.mean <- p*(p-1)/(mu.m*sqrt(pi*(m-3)))
    theo.mean <- 2*p*(p-1)*(1/sqrt(pi*(p-3)))
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    theo.var <- 4*((pi-2)/pi)*(p*(p-1)/(p-3))
    #theo.var <- 2*(6*(m-3)*(mu.m^2)*log(m*(p-1))*(pi - 2) - pi^2)*p/(pi*(m-3)*(mu.m^2)*(pi^2 + 12*(m-3)*(mu.m^2)*log(m*(p-1))))*(p-1)
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  
  theo.mean.vec[i] <- 2*p*(p-1)*(1/sqrt(pi*(p-3)))
  #theo.mean.vec[i] <- p*(p-1)/(mu.m*sqrt(pi*(m-3)))
  theo.var.vec[i] <- 4*((pi-2)/pi)*(p*(p-1)/(p-3))
  #theo.var.vec[i] <- 2*(6*(m-3)*(mu.m^2)*log(m*(p-1))*(pi - 2) - pi^2)*p/(pi*(m-3)*(mu.m^2)*(pi^2 + 12*(m-3)*(mu.m^2)*log(m*(p-1))))*(p-1)
  
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)
out.df

# min-max normalized diff

start.time <- Sys.time()

ps <- rep(c(seq(100,300,by=50)),length=25)
ms <- c(rep(100,length=5),rep(200,length=5),rep(300,length=5),
        rep(400,length=5),rep(500,length=5))
ps <- ps[1:5]
ms <- ms[1:5]
samp.mean.vec <- numeric()
samp.var.vec <- numeric()
for(outiter in 1:length(ms)){
  print(outiter)
  samp.var <- numeric()
  samp.mean <- numeric()
  for(iter in 1:50){
    cat("Simulation: ",iter,"\n")
    
    p <- ps[outiter] # number of attributes
    
    m <- ms[outiter]  # number of subjects
    
    data.mat <- matrix(rep(0,length=m*(p*(p-1))),nrow=(p*(p-1)),ncol=m)
    for(k in 1:m){
      #print(k)
      
      X <- matrix(runif(p^2,min=-1,max=1), ncol=p) 
      
      cov <- X %*% t(X) 
      
      mycorr <- cov2cor(cov) 
      
      #zcorr <- apply(matrix(mycorr[upper.tri(mycorr)],ncol=1),1,myfisher)
      zcorr <- apply(matrix(stretch_mat(mycorr),ncol=1),1,myfisher)
      #hist(zcorr1,breaks=100,freq=F)
      #data.mat[,k] <- matrix(zcorr,nrow=(p*(p-1)/2),ncol=1)
      data.mat[,k] <- matrix(zcorr,nrow=(p*(p-1)),ncol=1)
    }
    
    maxs <- numeric()
    mins <- numeric()
    diffs <- numeric()
    diffs.inv <- numeric()
    for(i in 1:p){
      #print(i)
      idx1 <- (i-1)*(p-1)+1
      idx2 <- i*(p-1)
      dats <- c(data.mat[c(idx1:idx2),])
      maxs[i] <- max(dats)
      mins[i] <- min(dats)
      diffs[i] <- max(dats)-min(dats)
      diffs.inv[i] <- 1/diffs[i]
    }
    
    diffs.inv.vec <- numeric()
    for(i in 1:length(diffs.inv)){
      #print(i)
      idx1 <- (i-1)*(p-1)+1
      idx2 <- i*(p-1)
      diffs.inv.vec[c(idx1:idx2)] <- rep(diffs.inv[i],length=(p-1))
    }
    
    # distance matrix
    d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
    for(i in 1:m){
      for(j in i:m){
        #diff <- abs(data.mat[,i]-data.mat[,j])
        diff <- abs(data.mat[,i]-data.mat[,j])*diffs.inv.vec
        d.mat[i,j] <- sum(diff)
      }
    }
    
    tmp <- d.mat + t(d.mat)
    dist.vec <- c(t(tmp))
    idx <- which(dist.vec==0)
    dist.vec <- dist.vec[-idx]
    
    mu.m <- (-qnorm(1/(m*(p-1))) + digamma(1)/qnorm(1/(m*(p-1))))*(1/sqrt(p-3)) # + (1/sqrt(p-3))*sqrt(((pi^2)/(12*(p-3)*log(m*(p-1)))))

    
    # theoretical mean
    theo.mean <- p*(p-1)/(mu.m*sqrt(pi*(p-3)))
    #theo.mean <- 2*p*(p-1)*(1/sqrt(pi*(p-3)))
    theo.mean
    
    # sample mean
    samp.mean[iter] <- mean(dist.vec)
    samp.mean
    
    # theoretical variance
    #theo.var <- 4*((pi-2)/pi)*(p*(p-1)/(m-3))
    #theo.var <- 2*(6*(p-3)*(mu.m^2)*log(m*(p-1))*(pi - 2) - pi^2)*p/(pi*(p-3)*(mu.m^2)*(pi^2 + 12*(p-3)*(mu.m^2)*log(m*(p-1))))*(p-1) #+ (p-1)*(1/sqrt(p-3))
    theo.var <- (p*(p-1)*(pi - 2))/(pi*(mu.m^2)*(p-3))
    theo.var
    
    # sample variance
    samp.var[iter] <- var(dist.vec)
    samp.var
  }
  
  samp.mean.vec[outiter] <- mean(samp.mean)
  samp.var.vec[outiter] <- mean(samp.var)
}

theo.mean.vec <- numeric()
theo.var.vec <- numeric()
for(i in 1:length(ms)){
  print(i)
  p <- ps[i]
  m <- ms[i]
  
  mu.m <- (-qnorm(1/(m*(p-1))) + digamma(1)/qnorm(1/(m*(p-1))))*(1/sqrt(p-3)) # + (1/sqrt(p-3))*sqrt(((pi^2)/(12*(p-3)*log(m*(p-1)))))
  
  #theo.mean.vec[i] <- 2*p*(p-1)*(1/sqrt(pi*(p-3)))
  theo.mean.vec[i] <- p*(p-1)/(mu.m*sqrt(pi*(p-3)))
  #theo.var.vec[i] <- 4*((pi-2)/pi)*(p*(p-1)/(p-3))
  #theo.var.vec[i] <- 2*(6*(p-3)*(mu.m^2)*log(m*(p-1))*(pi - 2) - pi^2)*p/(pi*(p-3)*(mu.m^2)*(pi^2 + 12*(p-3)*(mu.m^2)*log(m*(p-1))))*(p-1) #+ (p-1)*(1/sqrt(p-3))
  theo.var.vec[i] <- (p*(p-1)*(pi - 2))/(pi*(mu.m^2)*(p-3))
}

ands <- data.frame(and=rep("&",length=length(ms)))
dashes <- data.frame(dash=rep("\\",length=length(ms)))

dims <- paste("$(",ms,",",ps,")$",sep="")
out.df <- data.frame(dims=dims,and1=ands,samp.mean=samp.mean.vec,and2=ands,theo.mean=theo.mean.vec,
                     and3=ands,samp.var=samp.var.vec,and4=ands,theo.var=theo.var.vec,dash=dashes)
out.df

end.time <- Sys.time()
end.time - start.time
