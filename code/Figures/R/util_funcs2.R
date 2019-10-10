# Functions to make figures 6-8 and compute means and sds for fig 9, 10

draw_intervals <- function(intervals_df, x_limits){
  intervals_df %>%
    ggplot(aes(x = Type, color = Type)) +
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = .1) +
    geom_point(aes(y = Means), size = 2) +
    coord_flip() +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(legend.position = 'None') +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(limits = x_limits)
}

draw_densities <- function(sim_df, x_limits){
  sim_df %>%
    ggplot(aes(Value, fill = Type)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    scale_fill_viridis_d() +
    labs(x = NULL, y = NULL) +
    theme(legend.position = 'None',
          legend.title = element_blank()) +
    scale_x_continuous(label = NULL, limits = x_limits)
}

compute_dmat <- function(m, data.mat, PuPy.samp) {
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
  d.mat <- tmp
  d.mat
}

get_means_sds <- function(my.eta, m, p, probs, gammas, PuPy) {
  
  g0 <- gammas[1]
  g1 <- gammas[2]
  g2 <- gammas[3]
  
  PuPy.samp <- PuPy
  
  data.mat <- matrix(rep(0, length = m * p), nrow = p, ncol = m)
  for (i in 1:m) {
    samp <- matrix(as.double(rbinom(p, 2, probs)), nrow = p, ncol = 1)
    data.mat[, i] <- samp
  }
  
  d.mat <- compute_dmat(m, data.mat, PuPy.samp)

  long.dist.vec <- c(long.dist.vec,d.mat[upper.tri(d.mat)])

  samp.mean <- mean(d.mat[upper.tri(d.mat)])
  samp.var <- var(d.mat[upper.tri(d.mat)])

  F.a <- ((1 - probs)^3)*probs + (probs^3)*(1 - probs)
  G.a <- ((1 - probs)^2)*(probs^2)
  
  theo.mean <- (g0+g2+2*g1)*sum(F.a) + (1.5*(g0 + g2) + 2*g1)*sum(G.a)
  
  theo.var <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)

  list(
    theo.var = theo.var,
    samp.var = samp.var,
    theo.mean = theo.mean,
    samp.mean = samp.mean,
    long.dist.vec = long.dist.vec
  )
}

# Function for correlated rs-fMRI:

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

# functions for correlate Ti/Tv distances
rmvBinomial2 <- function(n, size, probs, corrmat) {
  X <- replicate(n, {
    colSums(rmvbin(size, margprob=probs, sigma=corrmat))
  })
  t(X)
}
