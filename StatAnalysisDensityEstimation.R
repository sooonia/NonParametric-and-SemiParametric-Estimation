library(KernSmooth)
library(ggplot2)
source('Helper.R')

num_runs <- 1000
num_observations <- 100

results <- data.frame(density = factor(), bandwidth = double(), 
                      method = factor(), sample_num = integer(), ASE = double())


# NORMAL
norm.params <- c(10, 2)
for(samp in 1:num_runs){
  d <- rnorm(num_observations, norm.params[1], norm.params[2])
  CV.results <- data.frame(bw=double(), CV = double())
  # sped up runtime --> O(n^2)
  bws <- c(Silverman(d)$bw, dpik(d, kernel = 'epanech'))
  lower <- max(min(bws)-sd(bws), .01)
  upper <- max(bws)+sd(bws)
  for(h in seq(from=lower,to=upper, by=sd(bws)/3)){
    CV.results = rbind(CV.results, data.frame(bw = h, CV = CV(h, d)))
  }
  bws <- c(bws, CV.results[which(CV.results$CV == min(CV.results$CV)) , 1])
  
  sum.silverman = 0
  sum.plug = 0
  sum.cv = 0
  
  for(obs in d){
    sum.silverman = sum.silverman + ((get.f.hat.ep(bws[1], d, obs) - dnorm(obs, norm.params))^2)
    sum.plug = sum.plug + ((get.f.hat.ep(bws[2], d, obs)- dnorm(obs, norm.params))^2)
    sum.cv = sum.cv + ((get.f.hat.ep(bws[3], d, obs) - dnorm(obs, norm.params))^2)
  }
  ase.silv <- sum.silverman/num_observations
  results <- rbind(results, data.frame(density = 'Normal', bandwidth= bws[1], 
                                       method = 'Silverman', sample_num=samp, ASE = ase.silv))
  
  ase.plug <- sum.plug/num_observations
  results <- rbind(results, data.frame(density = 'Normal', bandwidth= bws[2], 
                                       method = 'Plug In', sample_num=samp, ASE = ase.plug))
  
  ase.cv <- sum.cv/num_observations
  results <- rbind(results, data.frame(density = 'Normal', bandwidth= bws[3], 
                                       method = 'CV', sample_num=samp, ASE = ase.cv))
  
}



# Chi Square
chisq.params <- c(2)
for(samp in 1:num_runs){
  d <- rchisq(num_observations, chisq.params[1])
  CV.results <- data.frame(bw=double(), CV = double())
  # sped up runtime --> O(n^2)
  bws <- c(Silverman(d)$bw, dpik(d, kernel = 'epanech'))
  lower <- max(min(bws)-sd(bws), .01)
  upper <- max(bws)+sd(bws)
  for(h in seq(from=lower,to=upper, by=sd(bws)/3)){
    CV.results = rbind(CV.results, data.frame(bw = h, CV = CV(h, d)))
  }
  bws <- c(bws, CV.results[which(CV.results$CV == min(CV.results$CV)) , 1])
  
  sum.silverman = 0
  sum.plug = 0
  sum.cv = 0
  
  for(obs in d){
    sum.silverman = sum.silverman + ((get.f.hat.ep(bws[1], d, obs) - dchisq(obs, chisq.params))^2)
    sum.plug = sum.plug + ((get.f.hat.ep(bws[2], d, obs) - dchisq(obs, chisq.params))^2)
    sum.cv = sum.cv + ((get.f.hat.ep(bws[3], d, obs) - dchisq(obs, chisq.params))^2)
  }
  ase.silv <- sum.silverman/num_observations
  results <- rbind(results, data.frame(density = 'ChiSq', bandwidth= bws[1], 
                                       method = 'Silverman', sample_num=samp, ASE = ase.silv))
  
  ase.plug <- sum.plug/num_observations
  results <- rbind(results, data.frame(density = 'ChiSq', bandwidth= bws[2], 
                                       method = 'Plug In', sample_num=samp, ASE = ase.plug))
  
  ase.cv <- sum.cv/num_observations
  results <- rbind(results, data.frame(density = 'ChiSq', bandwidth= bws[3], 
                                       method = 'CV', sample_num=samp, ASE = ase.cv))
  
}



# BETA
beta.params <- c(2, 5)
for(samp in 1:num_runs){
  d <- rbeta(num_observations, beta.params[1], beta.params[2])
  CV.results <- data.frame(bw=double(), CV = double())
  # sped up runtime --> O(n^2)
  bws <- c(Silverman(d)$bw, dpik(d, kernel = 'epanech'))
  lower <- max(min(bws)-sd(bws), .01)
  upper <- max(bws)+sd(bws)
  for(h in seq(from=lower,to=upper, by=sd(bws)/3)){
    CV.results = rbind(CV.results, data.frame(bw = h, CV = CV(h, d)))
  }
  bws <- c(bws, CV.results[which(CV.results$CV == min(CV.results$CV)) , 1])
  
  sum.silverman = 0
  sum.plug = 0
  sum.cv = 0
  
  for(obs in d){
    sum.silverman = sum.silverman + ((get.f.hat.ep(bws[1], d, obs) - dbeta(obs, beta.params[1], beta.params[2]))^2)
    sum.plug = sum.plug + ((get.f.hat.ep(bws[2], d, obs) - dbeta(obs, beta.params[1], beta.params[2]))^2)
    sum.cv = sum.cv + ((get.f.hat.ep(bws[3], d, obs) - dbeta(obs, beta.params[1], beta.params[2]))^2)
  }
  ase.silv <- sum.silverman/num_observations
  results <- rbind(results, data.frame(density = 'Beta', bandwidth= bws[1], 
                                       method = 'Silverman', sample_num=samp, ASE = ase.silv))
  
  ase.plug <- sum.plug/num_observations
  results <- rbind(results, data.frame(density = 'Beta', bandwidth= bws[2], 
                                       method = 'Plug In', sample_num=samp, ASE = ase.plug))
  
  ase.cv <- sum.cv/num_observations
  results <- rbind(results, data.frame(density = 'Beta', bandwidth= bws[3], 
                                       method = 'CV', sample_num=samp, ASE = ase.cv))
  
}

save(results, file = 'densityResults.rda')

load('densityResults.rda')

n <- results[results$density == 'Normal',]
b <- results[results$density == 'Beta',]
c <- results[results$density == 'Beta',]

for(k in levels(n$method)){
  print(median(c[c$method == k, 'ASE']))
  print(paste(k, ":", sep = ''))
  print(paste("Mean:", mean(c[c$method == k, 'ASE'])))
  print(paste("Standard deviation:", sd(c[c$method == k, 'ASE'])))
  print(paste("Min:", min(c[c$method == k, 'ASE']), 
              "Max:", max(c[c$method == k, 'ASE'])))
  print(' ')
}

ggplot(n, aes(ASE, color = method))+
  geom_density() +
  ggtitle('Normal Distrubtion ASE Density')

ggplot(b, aes(ASE, color = method))+
  geom_density() +
  ggtitle('Beta Distrubtion ASE Density')

ggplot(c, aes(ASE, color = method))+
  geom_density() +
  ggtitle('Chi Square Distrubtion ASE Density')




# Comparing kernels
results.c <- data.frame(kernel = factor(), sample_num = integer(), ASE = double())
num_runs <- 1000
num_observations <- 100
beta.params <- c(2,5)

for(i in 1:num_runs){
  d <- rbeta(num_observations, beta.params[1], beta.params[2])
  bw <- dpik(d, kernel = 'normal')
  sum.unif = 0
  sum.ep = 0
  sum.quart = 0
  sum.triw = 0
  sum.gauss = 0
  for(obs in d){
    sum.unif = sum.unif + (get.f.hat.unif.c(bw, d, obs)- dbeta(obs, beta.params[1], beta.params[2]))^2
    sum.ep = sum.ep + (get.f.hat.ep.c(bw, d, obs)- dbeta(obs, beta.params[1], beta.params[2]))^2
    sum.quart = sum.quart + (get.f.hat.quart.c(bw, d, obs)- dbeta(obs, beta.params[1], beta.params[2]))^2
    sum.triw = sum.triw + (get.f.hat.triw.c(bw, d, obs)- dbeta(obs, beta.params[1], beta.params[2]))^2
    sum.gauss = sum.gauss + (get.f.hat.gauss.c(bw, d, obs)- dbeta(obs, beta.params[1], beta.params[2]))^2
  }
  
  results.c <- rbind(results.c, data.frame(kernel = 'Unif', sample_num=i, ASE = sum.unif/length(d)))
  results.c <- rbind(results.c, data.frame(kernel = 'Epanech.', sample_num=i, ASE = sum.ep/length(d)))
  results.c <- rbind(results.c, data.frame(kernel = 'Quartic', sample_num=i, ASE = sum.quart/length(d)))
  results.c <- rbind(results.c, data.frame(kernel = 'Triweight', sample_num=i, ASE = sum.triw/length(d)))
  results.c <- rbind(results.c, data.frame(kernel = 'Gaussian', sample_num=i, ASE = sum.gauss/length(d)))
}


#save(results.c, file = 'ComparingKernelsResults.rda')
#load('ComparingKernelsResults.rda')

for(k in levels(results.c$kernel)){
  print(median(results.c[results.c$kernel == k, 'ASE']))
  #print(paste(k, ":", sep = ''))
  #print(paste("Mean:", mean(results.c[results.c$kernel == k, 'ASE'])))
  #print(paste("Standard deviation:", sd(results.c[results.c$kernel == k, 'ASE'])))
  #print(paste("Min:", min(results.c[results.c$kernel == k, 'ASE']), 
   #     "Max:", max(results.c[results.c$kernel == k, 'ASE'])))
  print(' ')
}

ggplot(results.c, aes(ASE, color = kernel))+
  geom_density() +
  ggtitle('ChiSq(2) Distrubtion ASE Canonical Kernel Comparison')

ggplot(results.c, aes(ASE, color = kernel))+
  geom_density() +
  scale_x_continuous(limits = c(0, 500)) +
  ggtitle('Normal Distrubtion ASE Canonical Kernel Comparison \nZoomed')

#Strange:
d <- rchisq(num_observations, 2)
plot(density(d, bw=dpik(d, kernel = 'normal')
             *1.7188/0.7764, kernel="epanech"), main = 'Epanechnikov Kernel')
plot(density(d, bw=dpik(d, kernel = 'normal')
             *1.3510/0.7764, kernel="rect"), main = 'Uniform Kernel')
plot(density(d, bw=dpik(d, kernel = 'normal')* 2, kernel="gaussian"), main='Gaussian Kernel')
plot(density(d, bw=dpik(d, kernel = 'normal')
             *1.3510/0.7764, kernel="rect"), main = 'Uniform Kernel')
x <- seq(-5, 20, length=500)
hx <- dchisq(x,2)
plot(x, hx, type="l", lty=1, ylab="Density", main="True PDF of ChiSq(2)")
