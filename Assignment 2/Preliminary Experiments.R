source("Density Estimators.R")

# Generate data for the experiment
set.seed(360)
# The sample one
mu = 1; sigma = 2
sample1 = rnorm(5000, mean = mu, sd = sigma)
#The sample two
alpha = 2; beta = 4
sample2 = rbeta(5000, shape1 = alpha, shape2 = beta)

# The density functions for the distributions of the simulated data
# The density function for the distribution of sample one
density_sam1 = function(x){
  exp(-(x - mu)^2/(2*sigma^2))/(sqrt(2*pi)*sigma)
}
# The density function for the distribution of sample two
density_sam2 = function(x){
  x^(alpha - 1)*(1 - x)^(beta - 1)/beta(alpha, beta)
}

# Create the visualizations to illustrate the performances of different density estimators
# The visualization based on sample one
X1 = seq(-6, 9, 0.02)
pdf("Part1_fig1.pdf")
plot(X1, density_sam1(X1), type = "l", col = "blue", lwd = 1.6,
     main = "Performances of different density estimators\n on data from Normal distribution",
     xlab = "X", ylab = "Density", xlim = c(-6, 9), ylim = c(0, 0.22), xaxt = "n")
axis(side = 1, at = seq(-6, 9, 3))
lines(density(sample1), col = "green", lwd = 1.6)
lines(X1, sapply(X1, Histogram, h = 0.16, data = sample1), col = "purple", lwd = 1.6)
lines(X1, sapply(X1, OS_Hermite, data = sample1), col = "red", lwd = 1.6)
legend(x = "topright", legend = c("Accurate density", "Kernel density estimator",
                                 "Histogram estimator", "Orthogonal series estimator"),
       lty = c(1, 1, 1, 1), col = c("blue", "green", "purple", "red"), lwd = 1.8, bty = "n")
dev.off()

# The visualization based on sample two
X2 = seq(0, 1, 0.01)
pdf("Part1_fig2.pdf")
plot(X2, density_sam2(X2), type = "l", col = "blue", lwd = 1.6,
     main = "Performances of different density estimators\n on data from Beta distribution",
     xlab = "X", ylab = "Density", xlim = c(0, 1), ylim = c(0, 2.2))
lines(density(sample2), col = "green", lwd = 1.6)
lines(X2, sapply(X2, Histogram, h = 0.02, data = sample2), col = "purple", lwd = 1.6)
lines(X2, sapply(X2, OS_trigo, data = sample2), col = "red", lwd = 1.6)
legend(x = "topright", legend = c("Accurate density", "Kernel density estimator",
                                  "Histogram estimator", "Orthogonal series estimator"),
       lty = c(1, 1, 1, 1), col = c("blue", "green", "purple", "red"), lwd = 1.8, bty = "n")
dev.off()
