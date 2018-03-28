library(KernSmooth)
source('Helper.R')

# generate test distributions
# need to figure out how to get true pdf
#d1 <- samplepdf(1000, pdf1)
n= 100
d1 <- rchisq(n, 2)
d2 <- rnorm(n, 10, 2)
d3 <- rbeta(n, .5, .5)

plot(density(d1))
plot(density(d2))
plot(density(d3))

results <- data.frame(density = factor(), bandwidth = double(), method = factor())

# Silverman's Rule of Thumb
plot(Silverman(d1), main=paste('Silverman\'s Better Rule of Thumb', '\nDensity 1'))
plot(Silverman(d2), main=paste('Silverman\'s Better Rule of Thumb', '\nDensity 2'))
plot(Silverman(d3), main=paste('Silverman\'s Better Rule of Thumb', '\nDensity 3'))

results <- rbind(results, data.frame(density = 'd1', bandwidth= Silverman(d1)$bw, 
                                     method = 'Silverman'))
results <- rbind(results, data.frame(density = 'd2', bandwidth= Silverman(d2)$bw, 
                                     method = 'Silverman'))
results <- rbind(results, data.frame(density = 'd3', bandwidth= Silverman(d3)$bw, 
                                     method = 'Silverman'))

# Plug in
plot(density(d1, bw = dpik(d1, kernel = 'epanech'), kernel = 'epanechnikov'), 
     main=paste('Refined Plug in Method', '\nDensity 1'))
plot(density(d2, bw = dpik(d2, kernel = 'epanech'), kernel = 'epanechnikov'), 
     main=paste('Refined Plug in Method', '\nDensity 2'))
plot(density(d3, bw = dpik(d3, kernel = 'epanech'), kernel = 'epanechnikov'), 
     main=paste('Refined Plug in Method', '\nDensity 3'))
results <- rbind(results, data.frame(density = 'd1', bandwidth= dpik(d1, kernel = 'epanech'),
                                     method = 'PlugIn'))
results <- rbind(results, data.frame(density = 'd2', bandwidth= dpik(d2, kernel = 'epanech'), 
                                     method = 'PlugIn'))
results <- rbind(results, data.frame(density = 'd3', bandwidth= dpik(d3, kernel = 'epanech'), 
                                     method = 'PlugIn'))




# Cross Validation
CV.results <- data.frame(bw=double(), CV = double(), density = factor())
# sped up runtime --> O(n^2)
for(h in seq(from=.05,to=2.5, by=.3)){
  CV.results = rbind(CV.results, data.frame(bw = h, CV = CV(h, d1), density = d1 ))
}
plot(CV.results[,c(1,2)])

# Choosing Kernel
# change to canonical
plot(density(d1, bw = dpik(d1, kernel = 'epanech', canonical = TRUE), kernel = 'epanechnikov'), 
     main=paste('Density 1' ,'\nKernel= epanechnikov'))

plot(density(d1, bw = dpik(d1, kernel = 'box', canonical = TRUE), kernel = 'rectangular'),
     main=paste('Density 1' ,'\nKernel= uniform'))

plot(density(d1, bw = dpik(d1, kernel = 'normal', canonical = TRUE), kernel = 'gaussian'),
     main=paste('Density 1' ,'\nKernel= gaussian'))
plot(density(d1, bw = dpik(d1, kernel = 'biweight', canonical = TRUE), kernel = 'biweight'), 
     main=paste('Density 1' ,'\nKernel= Quartic'))

plot(density(d2, bw = dpik(d2, kernel = 'epanech', canonical = TRUE), kernel = 'epanechnikov'), 
     main=paste('Density 2' ,'\nKernel= epanechnikov'))
plot(density(d2, bw = dpik(d2, kernel = 'box', canonical = TRUE), kernel = 'rectangular'),
     main=paste('Density 2' ,'\nKernel= uniform'))
plot(density(d2, bw = dpik(d2, kernel = 'triweight', canonical = TRUE), kernel = 'triangular'), 
     main=paste('Density 2' ,'\nKernel= triweight'))
plot(density(d2, bw = dpik(d2, kernel = 'normal', canonical = TRUE), kernel = 'gaussian'),
     main=paste('Density 2' ,'\nKernel= gaussian'))
plot(density(d2, bw = dpik(d2, kernel = 'biweight', canonical = TRUE), kernel = 'biweight'), 
     main=paste('Density 2' ,'\nKernel= Quartic'))

plot(density(d3, bw = dpik(d3, kernel = 'epanech', canonical = TRUE), kernel = 'epanechnikov'), 
     main=paste('Density 3' ,'\nKernel= epanechnikov'))
plot(density(d3, bw = dpik(d3, kernel = 'box', canonical = TRUE), kernel = 'rectangular'),
     main=paste('Density 3' ,'\nKernel= uniform'))
plot(density(d3, bw = dpik(d3, kernel = 'triweight', canonical = TRUE), kernel = 'triangular'), 
     main=paste('Density 3' ,'\nKernel= triweight'))
plot(density(d3, bw = dpik(d3, kernel = 'normal', canonical = TRUE), kernel = 'gaussian'),
     main=paste('Density 3' ,'\nKernel= gaussian'))
plot(density(d3, bw = dpik(d3, kernel = 'biweight', canonical = TRUE), kernel = 'biweight'), 
     main=paste('Density 3' ,'\nKernel= Quartic'))

