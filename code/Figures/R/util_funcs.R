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
  d.mat <- matrix(rep(0, length = m * m), nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in i:m) {
      diff <- abs(data.mat[i,] - data.mat[j,])
      diff.TiTv <- diff
      for (k in 1:p) {
        if (diff[k] == 1) {
          if (PuPy.samp[k] == 3) {
            diff.TiTv[k] <- 0.25
          } else if (PuPy.samp[k] == 5) {
            diff.TiTv[k] <- 0.25
          } else{
            diff.TiTv[k] <- 0.5
          }
        } else if (diff[k] == 2) {
          if (PuPy.samp[k] == 3) {
            diff.TiTv[k] <- 0.75
          } else if (PuPy.samp[k] == 5) {
            diff.TiTv[k] <- 0.75
          } else{
            diff.TiTv[k] <- 1
          }
        }
      }
      d.mat[i, j] <- sum(diff.TiTv)
    }
  }
  d.mat
}

get_means_sds <- function(my.eta, p, lowers, uppers, maf_idx) {
  g0 <- runif(1, min = 0.01, max = (my.eta / (my.eta + 1) - 0.01))
  g1 <- 1/(my.eta+1)
  g2 <- my.eta / (my.eta + 1) - g0
  
  PuPy.samp <- sample(
    c(3, 4, 5),
    size = p,
    prob = c(g0, g1, g2),
    replace = T
  )
  
  probs <- runif(p, lowers[maf_idx], uppers[maf_idx])
  data.mat <- matrix(rep(0, length = m * p), nrow = p, ncol = m)
  for (i in 1:m) {
    samp <- matrix(as.double(rbinom(p, 2, probs)), nrow = p, ncol = 1)
    data.mat[, i] <- samp
  }
  
  diff.mat <-
    matrix(rep(0, length = 4 * (m * (m - 1) / 2)), nrow = (m * (m - 1) / 2), ncol = 5)
  d.mat <- compute_dmat(m, data.mat, PuPy.samp)
  # d.mat <- matrix(rep(0, length = m * m), nrow = m, ncol = m)
  # for (i in 1:m) {
  #   for (j in i:m) {
  #     diff <- abs(data.mat[, i] - data.mat[, j])
  #     diff.TiTv <- diff
  #     for (k in 1:p) {
  #       if (diff[k] == 1) {
  #         if (PuPy.samp[k] == 3) {
  #           diff.TiTv[k] <- 0.25
  #         } else if (PuPy.samp[k] == 5) {
  #           diff.TiTv[k] <- 0.25
  #         } else{
  #           diff.TiTv[k] <- 0.5
  #         }
  #       } else if (diff[k] == 2) {
  #         if (PuPy.samp[k] == 3) {
  #           diff.TiTv[k] <- 0.75
  #         } else if (PuPy.samp[k] == 5) {
  #           diff.TiTv[k] <- 0.75
  #         } else{
  #           diff.TiTv[k] <- 1
  #         }
  #       }
  #     }
  #     d.mat[i, j] <- sum(diff.TiTv)
  #   }
  # }
  tmp <- d.mat + t(d.mat)
  long.dist.vec <- c(long.dist.vec,tmp[upper.tri(tmp)])
  dist.vec <- c(t(tmp))
  dist.vec <- dist.vec[dist.vec!=0]
  
  E.diff0 <- mean((1 - probs)^4 + 4*(probs^2)*((1 - probs)^2) + probs^4)
  E.diff.25 <- 4*(g0 + g2)*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
  E.diff.5 <- 4*g1*mean(probs*(1 - probs)^3 + (1 - probs)*probs^3)
  E.diff.75 <- 2*(g0 + g2)*mean(((1 - probs)^2)*(probs^2))
  E.diff1 <- 2*g1*mean(((1 - probs)^2)*(probs^2))
  
  E.Dij <- (0*E.diff0 + 0.25*E.diff.25 + 0.5*E.diff.5 + 0.75*E.diff.75 + 1*E.diff1)*p
  
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
  v2 <- cumsum(v)[-1]
  
  E.D2.2 <- c(2 * matrix(v1, nrow = 1, ncol = length(v1)) %*% 
                matrix(v2, nrow = length(v2), ncol = 1))
  E.Dsqij.new <- E.D2.1 + E.D2.2
  Var.Dij.new <- E.Dsqij.new - (E.Dij) ^ 2
  
  new.var <- mean((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2)) - mean(((g0 + g2+ 2*g1)^2)*((probs*((1-probs)^3) + (1-probs)*(probs^3))^2) - 2*(g0 + g2 + 2*g1)*(1.5*(g0 + g2) + 2*g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3))*((1-probs)^2)*(probs^2) - ((1.5*(g0 + g2) + 2*g1)^2)*((1-probs)^4)*(probs^4))
  new.var2 <- (0.25*(g0+g2)+g1)*sum(F.a)+((9/8)*(g0+g2)+2*g1)*sum(G.a)-sum(((g0+g2+2*g1)*F.a+(1.5*(g0+g2)+2*g1)*G.a)^2)
  theo.m2 <- p*sum((0.25*(g0 + g2) + g1)*(probs*((1-probs)^3) + (1-probs)*(probs^3)) + ((9/8)*(g0 + g2) + 2*g1)*((1-probs)^2)*(probs^2))
  
  list(
    theo.var1 = Var.Dij.new,
    theo.var2 = new.var,
    theo.var3 = new.var2,
    samp.var = var(dist.vec),
    theo.mean = E.Dij,
    theo.mean2 = new.mean,
    theo.mean3 = new.mean2,
    samp.mean = mean(dist.vec),
    samp.m2 = mean(dist.vec ^ 2),
    theo.m2 = theo.m2,
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
