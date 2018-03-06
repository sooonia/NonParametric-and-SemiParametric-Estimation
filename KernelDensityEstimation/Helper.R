# Silverman's Rule of Thumb Functions
Silverman <- function(obs, n=1000, ker = 'epanechnikov'){
  s <- summary(obs)
  iqr <- as.numeric(s[5] - s[2])
  use <- min(c(iqr/1.34, sd(obs)))
  h.rot <- 1.06*use*(n^(-1/5))
  
  return(density(obs, bw = h.rot, kernel = ker))
}


# Cross Validation Functions
ep.kernel <- function(u){
  val <- 0
  if(abs(u) <= 1){
    val <- .75*(1-(u*u))
  }
  return(val)
}
k_h <- function(h,u){
  return(ep.kernel(u/h)/h)
}
#K*K(u) = 3/160*(32-40*u^2+20*|u|^3-|u|^5) * I(|u|<=2)
conv <- function(u){
  val <- 0
  if(abs(u) <= 2){
    val <- (3/160)*(32-(40*u^2)+(20*abs(u^3))-abs(u^5))
  }
  return(val)
}
CV <- function(h,X){
  n = length(X)
  sum1 = 0
  for(obs.i in X){
    for(obs.j in X){
      sum1 = sum1 + conv((obs.j - obs.i)/h)
    }
  }
  term1 = sum1 / (n^2*h)
  sum2 = 0
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        sum2 = sum2 + k_h(h, X[i]-X[j])
      }
    }
  }
  term2 = 2*sum2 / (n*(n-1))
  return(term1-term2)
}


# Evaluation functions
get.f.hat <- function(n, h, data, x){
  sum = 0
  for(obs in data){
    sum = sum + ep.kernel((x-obs)/h)
  }
  return(sum/(n*h))
}