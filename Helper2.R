# Silverman's Rule of Thumb Functions
Silverman <- function(obs, n=10000, ker = 'epanechnikov'){
  n = length(obs)
  s <- summary(obs)
  iqr <- as.numeric(s[5] - s[2])
  use <- min(c(iqr/1.34, sd(obs)))
  h.rot <- 1.06*use*(n^(-1/5))
  
  return(density(obs, bw = h.rot, kernel = ker))
}


SSE <- function(preds, obs){
  SSE = 0
  for(i in seq(1,length(preds))){
    SSE = SSE + (preds[i] - obs[i])^2
  }
  
  return(SSE)
}