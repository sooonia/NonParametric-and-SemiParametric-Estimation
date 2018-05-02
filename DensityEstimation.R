# generate data
norm1 <- rnorm(500, 10, 2)
norm2 <- rnorm(500, 3, 2)
nonnorm <- rbeta(1000, .5, .5)


# estimate density with histogram using different bin sizes
for (i in c(1,5, 10, 20, 50, 100, 200, 500)){
  hist(nonnorm, breaks = i, main = paste('Num Bins=', i))
}
h.opt <- 3.5*1000^(-1/3)
range <- max(nonnorm) - min(nonnorm)
num.bins <- round(range/h.opt)
hist(nonnorm, breaks = num.bins, main = paste('Num Bins=', num.bins))


# plot density
#plot(density(nonnorm, bw= NUM, kernel = c("gaussian", "epanechnikov", "rectangular",
#                                        "triangular", "biweight",
#                                        "cosine", "optcosine")))
#
plot(density(nonnorm))
for(b in c(.1, .5, 1, 5, 10)){
  plot(density(nonnorm, bw = b, kernel = 'gaussian'), main=paste(paste('bandwidth=',b),'\nKernel= gaussian'))
}

for(b in c(.1, .5, 1, 5, 10)){
  plot(density(nonnorm, bw = b, kernel = 'triangular'), main=paste(paste('bandwidth=',b),'\nKernel= triangular'))
}

for(b in c(.01, .025, .05, .1, .5, 1)){
  plot(density(nonnorm, bw = b, kernel = 'epanechnikov'), main=paste(paste('bandwidth=',b),'\nKernel= epanechnikov'))
}

for(b in c(.1, .5, 1, 5, 10)){
  plot(density(nonnorm, bw = b, kernel = 'rectangular'), main=paste(paste('bandwidth=',b),'\nKernel= rectangular'))
}
# using optimal bandwidth from book
s <- summary(nonnorm)
iqr <- as.numeric(s[5] - s[2])
use <- min(c(iqr/1.34, sd(nonnorm)))
h.rot <- 1.06*use*(1000^(-1/5))
plot(density(nonnorm, bw = h.rot, kernel = 'gaussian'), main=paste(paste('bandwidth=',h.rot),'\nKernel= gaussian'))
plot(density(nonnorm, bw = h.rot, kernel = 'triangular'), main=paste(paste('bandwidth=',h.rot),'\nKernel= triangular'))
plot(density(nonnorm, bw = h.rot, kernel = 'epanechnikov'), main=paste(paste('bandwidth=',h.rot),'\nKernel= epanechnikov'))
plot(density(nonnorm, bw = h.rot, kernel = 'rectangular'), main=paste(paste('bandwidth=',h.rot),'\nKernel= rectangular'))
plot(density(nonnorm))

