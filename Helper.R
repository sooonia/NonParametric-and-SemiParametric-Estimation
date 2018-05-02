# Silverman's Rule of Thumb Functions
Silverman <- function(obs, n=10000, ker = 'epanechnikov'){
  n = length(obs)
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
# speed up loops
CV <- function(h,X){
  n = length(X)
  sum1 = 0
  
  sum2 = 0
  for(i in 1:n){
    obs.i = X[i]
    for(j in 1:n){
      obs.j = X[j]
      sum1 = sum1 + conv((obs.j - obs.i)/h)
      if(i != j){
        sum2 = sum2 + k_h(h, X[i]-X[j])
      }
    }
  }
  term1 = sum1 / (n^2*h)
  term2 = 2*sum2 / (n*(n-1))
  return(term1-term2)
}


# Evaluation functions

# epanechnikov
get.f.hat.ep <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + ep.kernel((x-obs)/bw)
  }
  return(sum/(n*bw))
}

get.f.hat.ep.c <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + ep.kernel((x-obs)/(bw *(1.7188/0.7764)))
  }
  return(sum/(n*bw*(1.7188/0.7764)))
}

#Uniform
unif.kernel <- function(u){
  if(abs(u) <= 1){
    return(1/2)
  }  
  else{return(0)}
}

get.f.hat.unif.c <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + unif.kernel((x-obs)/bw * (1.3510/0.7764))
  }
  return(sum/(n*bw*1.3510/0.7764))
}

# Triweight
triw.kernel <- function(u){
  if(abs(u) <= 1){
    return(35/32*((1-u^2)^3))
  }  
  else{return(0)}
}

get.f.hat.triw.c <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + triw.kernel((x-obs)/(bw*2.3122/0.7764))
  }
  return(sum/(n*bw*2.3122/0.7764))
}

# Quartic
quart.kernel <- function(u){
  if(abs(u) <= 1){
    return(15/16*((1-u^2)^2))
  }  
  else{return(0)}
}

get.f.hat.quart.c <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + quart.kernel((x-obs)/(bw*2.0362/0.7764))
  }
  return(sum/(n*bw*2.0362/0.7764))
}


# Gaussian
gauss.kernel <- function(u){
  return(1/sqrt(2*pi)*exp(-.5*(u^2)))
}

get.f.hat.gauss.c <- function(bw, data, x){
  n = length(data)
  sum = 0
  for(obs in data){
    sum = sum + gauss.kernel((x-obs)/bw)
  }
  return(sum/(n*bw))
}

get.true.f <- function(d_name, x, para){
  if(d_name == 'ChiSq'){
    df <- para[1]
    return(dchisq(x, df))
  }
  else if(d_name == 'Normal'){
    mu = para[1]
    sig = para[2]
    return(dnorm(x, mu, sig))
  }
  else if(d_name == 'Beta'){
    alpha = para[1]
    beta = para[2]
    return(dbeta(x, alpha, beta))
  }
  else{return('error: d_name not valid')}
  
}




# ASE Loop