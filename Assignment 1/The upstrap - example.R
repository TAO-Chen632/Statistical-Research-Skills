library(dplyr)

# Prepare the data set
load("SHHS.Rdata")
df = SHHS %>%
  mutate(slp_apnea = (rdi4p >= 15)) %>%
  select(pptid, slp_apnea, gender, age_s1, bmi_s1, HTNDerv_s1) %>%
  rename(ID = pptid, age = age_s1, BMI = bmi_s1, HTN = HTNDerv_s1)

# Perform the upstrap algorithm
U = 10000
n = dim(df)[1]
r_seq = seq(1, 5, 0.2)
frequ = c()
alpha = 0.05
for (r in r_seq){
  count = 0
  for (u in 1:U){
    data = df[sample(1:n, size = n*r, replace = TRUE), ]
    sleep_logi = glm(slp_apnea ~ gender + age + BMI + HTN + HTN*age, data = data, family = "binomial")
    summary_glm = summary(sleep_logi)
    if(summary_glm$coefficients[5, 4] < alpha){
      count = count + 1
    }
  }
  frequ = c(frequ, count/U)
}

# Create the visualization of the power to detect the HTN effect
pdf("Figure.pdf")
plot(r_seq, frequ, type = "l", col = "blue", lwd = 2,
     xlab = "Multiplier of the sample size, r",
     ylab = "Power to detect the HTN effect",
     bty = "l", cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0.8,  col = "red", lwd = 2, lty = 1)
abline(h = 0.9, col = "purple", lwd = 2, lty = 1)
abline(v = 3.06, col = "green3", lwd = 2, lty = 2)
dev.off()
