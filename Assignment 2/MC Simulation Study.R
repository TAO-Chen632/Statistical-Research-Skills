source("Density Estimators.R")
# Initialize some variables
# R is the times that we repeat the one-shot experiment
R = 1000
n = c(250, 500, 1000)
df = 5
ISE = list()
ISE$n3 = ISE$n2 = ISE$n1 = data.frame("kernel" = rep(0, R), "hist" = rep(0, R), "OS" = rep(0, R))
set.seed(600)

# The density function of the distribution of the sample
density_sam = function(x, df){
  (1 + x^2/df)^(-(df + 1)/2)*gamma((df + 1)/2)/gamma(df/2)/sqrt(df*pi)
}

# The Monte Carlo Simulation study for sample size n=250
for(i in 1:R){
  sample = rt(n[1], df = df)
  # The calculation of the ISE for the kernel density estimator
  ker_den = density(sample)
  interval = ker_den$x[-1] - ker_den$x[-length(ker_den$x)]
  x_pts = (ker_den$x[-1] + ker_den$x[-length(ker_den$x)])/2
  y_pts = (ker_den$y[-1] + ker_den$y[-length(ker_den$y)])/2
  ISE$n1$kernel[i] = sum((y_pts - density_sam(x_pts, df = df))^2*interval)
  # The calculation of the ISE for the histogram estimator
  ISE$n1$hist[i] = integrate(function(x){
    (Histogram(x, h = 0.1, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
  # The calculation of the ISE for the orthogonal series estimator
  ISE$n1$OS[i] = integrate(function(x){
    (OS_Hermite(x, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
}

# The Monte Carlo Simulation study for sample size n=500
for(i in 1:R){
  sample = rt(n[2], df = df)
  # The calculation of the ISE for the kernel density estimator
  ker_den = density(sample)
  interval = ker_den$x[-1] - ker_den$x[-length(ker_den$x)]
  x_pts = (ker_den$x[-1] + ker_den$x[-length(ker_den$x)])/2
  y_pts = (ker_den$y[-1] + ker_den$y[-length(ker_den$y)])/2
  ISE$n2$kernel[i] = sum((y_pts - density_sam(x_pts, df = df))^2*interval)
  # The calculation of the ISE for the histogram estimator
  ISE$n2$hist[i] = integrate(function(x){
    (Histogram(x, h = 0.1, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
  # The calculation of the ISE for the orthogonal series estimator
  ISE$n2$OS[i] = integrate(function(x){
    (OS_Hermite(x, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
}

# The Monte Carlo Simulation study for sample size n=1000
for(i in 1:R){
  sample = rt(n[3], df = df)
  # The calculation of the ISE for the kernel density estimator
  ker_den = density(sample)
  interval = ker_den$x[-1] - ker_den$x[-length(ker_den$x)]
  x_pts = (ker_den$x[-1] + ker_den$x[-length(ker_den$x)])/2
  y_pts = (ker_den$y[-1] + ker_den$y[-length(ker_den$y)])/2
  ISE$n3$kernel[i] = sum((y_pts - density_sam(x_pts, df = df))^2*interval)
  # The calculation of the ISE for the histogram estimator
  ISE$n3$hist[i] = integrate(function(x){
    (Histogram(x, h = 0.1, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
  # The calculation of the ISE for the orthogonal series estimator
  ISE$n3$OS[i] = integrate(function(x){
    (OS_Hermite(x, data = sample) - density_sam(x, df = df))^2},
    lower = -Inf, upper = Inf, subdivisions = 10^6)$value
}

# Create the visualization of the integrated squared errors (ISE) of different estimators
# on samples with sizes n = 250, 500 and 1000 respectively
library(dplyr)
library(tidyr)
library(ggplot2)

# The visualization of the ISE of different estimators on samples with size n = 250
pdf("Part2_fig1.pdf")
ISE$n1 %>% rename("Kernel density estimator" = kernel, "Histogram estimator" = hist,
                  "Orthogonal series estimator" = OS) %>%
  pivot_longer(cols = 1:3, names_to = "Estimator", values_to = "ISE") %>%
  ggplot(
    aes(
      x = ISE,
      fill = Estimator
    )
  ) +
  geom_histogram(bins = 25) +
  facet_grid(Estimator ~ .) +
  guides(fill = FALSE) +
  labs(
    title = "Integrated squared errors of different estimators on samples with size n = 250",
    x = "Integrated squared error (ISE)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("purple", "green3", "red2")) +
  scale_x_continuous(breaks = seq(0, 0.1, 0.01))
dev.off()

# The visualization of the ISE of different estimators on samples with size n = 500
pdf("Part2_fig2.pdf")
ISE$n2 %>% rename("Kernel density estimator" = kernel, "Histogram estimator" = hist,
                  "Orthogonal series estimator" = OS) %>%
  pivot_longer(cols = 1:3, names_to = "Estimator", values_to = "ISE") %>%
  ggplot(
    aes(
      x = ISE,
      fill = Estimator
    )
  ) +
  geom_histogram(bins = 25) +
  facet_grid(Estimator ~ .) +
  guides(fill = FALSE) +
  labs(
    title = "Integrated squared errors of different estimators on samples with size n = 500",
    x = "Integrated squared error (ISE)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("purple", "green3", "red2")) +
  scale_x_continuous(breaks = seq(0, 0.04, 0.005))
dev.off()

# The visualization of the ISE of different estimators on samples with size n = 1000
pdf("Part2_fig3.pdf")
ISE$n3 %>% rename("Kernel density estimator" = kernel, "Histogram estimator" = hist,
                  "Orthogonal series estimator" = OS) %>%
  pivot_longer(cols = 1:3, names_to = "Estimator", values_to = "ISE") %>%
  ggplot(
    aes(
      x = ISE,
      fill = Estimator
    )
  ) +
  geom_histogram(bins = 25) +
  facet_grid(Estimator ~ .) +
  guides(fill = FALSE) +
  labs(
    title = "Integrated squared errors of different estimators on samples with size n = 1000",
    x = "Integrated squared error (ISE)",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("purple", "green3", "red2")) +
  scale_x_continuous(breaks = seq(0, 0.02, 0.0025))
dev.off()
