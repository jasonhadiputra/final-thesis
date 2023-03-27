library(tidyverse)
library(readxl)
library(statmod)
library(tweedie)
library(fitdistrplus)

# ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd']
linewidth <- 0.75
theme <- theme_light() +
  theme(
    line=element_line(linewidth = 0.75),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_line(colour='lightgray', linewidth=0.6),
    panel.border=element_rect(colour='lightgray', linewidth=1.2),
    axis.ticks=element_blank(),
    axis.text=element_text(colour='black'),
    text=element_text(colour='black', size=14)
  )

base_dir <- r"{G:\.shortcut-targets-by-id\1_UrMWiRHhBVY6BmkXpUmghuXVDupcUm8\Jason\}"
df <- read_excel(paste(base_dir, r"{Datasets\Medical Cost Personal.xlsx}", sep=""))


y <- df$charges

df <- df %>%
  mutate(bmi3 = case_when (
    bmi < 18.5 ~ "underweight",
    bmi >= 18.5 & bmi < 24.9 ~ "healthy",
    bmi >= 24.9 & bmi < 29.9 ~ "overweight",
    bmi >= 29.9 & bmi < 34.9 ~ "class I obese",
    bmi >= 34.9 & bmi < 39.9 ~ "class II obese",
    TRUE ~ "class III obese"
  ))


### ESTIMASI PARAMETER DISTRIBUSI GAMMA
(est_gamma <- fitdist(y / 1000, "gamma", method='mle')$estimate)
est_gamma['rate'] <- est_gamma['rate'] / 1000

(est_shape <- est_gamma['shape'])
(est_rate <- est_gamma['rate'])
est_shape/est_rate
est_shape/est_rate^2

# uji Kolmogorov-Smirnov
library(goftest)
ks.test(y, "pgamma", shape=est_shape, rate=est_rate)

# uji Cramer-von Mises
library(goftest)
cvm.test(y, "pgamma", shape=est_shape, rate=est_rate)


### ESTIMASI PARAMETER DISTRIBUSI INVERSE GAUSSIAN
(est_ig <- fitdist(y / 1e6, "invgauss", method='mle', start = list(mean = 5, shape = 1))$estimate)
est_ig['mean'] <- est_ig['mean'] * 1e6
est_ig['shape'] <- est_ig['shape'] * 1e6
(est_mean <- est_ig['mean'])
(est_shape <- est_ig['shape'])
est_mean
1/est_shape

# uji Kolmogorov-Smirnov
library(goftest)
ks.test(y, "pinvgauss", mean=est_mean, shape=est_shape)

# uji Cramer-von Mises
library(goftest)
cvm.test(y, "pinvgauss", mean=est_mean, shape=est_shape)


### ESTIMASI PARAMETER DISTRIBUSI TWEEDIE

library(statmod)
library(tweedie)
out <- tweedie.profile(y~1,
                       p.vec=seq(2, 3, by=0.1),
                       do.smooth=TRUE)

est_p <- out$p.max
Lmax <- out$L.max
est_phi <- out$phi.max
est_mu <- mean(y)
(est_tweedie <- list(xi=est_p, phi=est_phi, mu=est_mu))

# Plot profiling
(plot <- ggplot() +
    geom_line(aes(x=out$x, y=out$y), linewidth=linewidth, colour='#4c72b0') +
    geom_point(aes(x=out$p, y=out$L), size=3, colour='#4c72b0') +
    geom_hline(yintercept=out$L.max, linetype='longdash', colour='#dd8452') +
    geom_vline(xintercept=out$p.max, linetype='longdash', colour='#dd8452') +
    geom_point(aes(x=out$p.max, y=out$L.max), size=3, colour='#dd8452') +
    labs(x='p index (95% confidence interval)', y='Log-likelihood') +
    theme)

plot + ggsave(paste(base_dir, r"{Pictures\profile log likelihood.png}", sep=""), width=15, height=9, units='cm', dpi=600)

# uji Kolmogorov-Smirnov
library(goftest)
ks.test(y, "ptweedie", xi=out$p.max, mu = mean(y),
        phi=out$phi.max, power=out$p.max)

# uji Cramer-von Mises
library(goftest)
cvm.test(y, "ptweedie", xi=out$p.max, mu = mean(y),
         phi=out$phi.max, power=out$p.max)

### PLOTTING DISTRIBUSI


# PDF No Tweedie
(plot <- ggplot() +
    geom_histogram(aes(x=y, y=..density..), bins=30, alpha=0.8, fill="#8c8c8c", colour="white") +
    geom_function(aes(colour='Gamma'), fun=dgamma, args=est_gamma, linewidth=linewidth, n=1000) +
    geom_function(aes(colour='IG'), fun=dinvgauss, args=est_ig, linewidth=linewidth, n=1000) +
    labs(x='Besar premi (USD)', y='Fungsi kepadatan peluang') +
    scale_colour_manual(
      name='Distribusi',
      breaks=c('Gamma', 'IG', 'Tweedie', 'Empiris'),
      values=c("#4c72b0", "#dd8452", "#55a868", "#8c8c8c")
    ) +
    coord_cartesian(ylim = c(0, 9e-5), expand=FALSE) +
    theme
    )

plot + ggsave(paste(base_dir, r"{Pictures\pdf gamma, ig, no tweedie.png}", sep=""), width=15, height=9, units='cm', dpi=600)

# CDF No Tweedie
(plot <- ggplot() +
    geom_function(aes(colour='Gamma'), fun=pgamma, args=est_gamma, linewidth=linewidth) +
    geom_function(aes(colour='IG'), fun=pinvgauss, args=est_ig, linewidth=linewidth) +
#    geom_function(aes(colour='Tweedie'), fun=ptweedie, args=est_tweedie, linewidth=linewidth) +
    stat_ecdf(aes(x=y, colour="Empiris"), linewidth=linewidth) +
    labs(x='Besar premi (USD)', y='Fungsi distribusi kumulatif') +
    scale_colour_manual(
      name='Distribusi',
      breaks=c('Gamma', 'IG', 'Tweedie', 'Empiris'),
      values=c("#4c72b0", "#dd8452", "#55a868", "#8c8c8c")
    ) +
    coord_cartesian(ylim = c(0, NA), expand=FALSE) +
    theme)

plot + ggsave(paste(base_dir, r"{Pictures\cdf gamma, ig, no tweedie.png}", sep=""), width=15, height=9, units='cm', dpi=600)



# PDF Tweedie
(plot <- ggplot() +
    geom_histogram(aes(x=y, y=..density..), bins=30, alpha=0.8, fill="#8c8c8c", colour="white") +
    geom_function(aes(colour='Gamma'), fun=dgamma, args=est_gamma, linewidth=linewidth, n=1000) +
    geom_function(aes(colour='IG'), fun=dinvgauss, args=est_ig, linewidth=linewidth, n=1000) +
    geom_function(aes(colour='Tweedie'), fun=dtweedie, args=est_tweedie, linewidth=linewidth, n=1000) +
    labs(x='Besar premi (USD)', y='Fungsi kepadatan peluang') +
    scale_colour_manual(
      name='Distribusi',
      breaks=c('Gamma', 'IG', 'Tweedie', 'Empiris'),
      values=c("#4c72b0", "#dd8452", "#55a868", "#8c8c8c")
    ) +
    coord_cartesian(ylim = c(0, 9e-5), expand=FALSE) +
    theme
)

plot + ggsave(paste(base_dir, r"{Pictures\pdf gamma, ig, tweedie.png}", sep=""), width=15, height=9, units='cm', dpi=600)


# Tweedie
(plot <- ggplot() +
    geom_function(aes(colour='Gamma'), fun=pgamma, args=est_gamma, linewidth=linewidth) +
    geom_function(aes(colour='IG'), fun=pinvgauss, args=est_ig, linewidth=linewidth) +
    geom_function(aes(colour='Tweedie'), fun=ptweedie, args=est_tweedie, linewidth=linewidth) +
    stat_ecdf(aes(x=y, colour="Empiris"), linewidth=linewidth) +
    labs(x='Besar premi (USD)', y='Fungsi distribusi kumulatif') +
    scale_colour_manual(
      name='Distribusi',
      breaks=c('Gamma', 'IG', 'Tweedie', 'Empiris'),
      values=c("#4c72b0", "#dd8452", "#55a868", "#8c8c8c")
    ) +
    coord_cartesian(ylim = c(0, NA), expand=FALSE) +
    theme)

plot + ggsave(paste(base_dir, r"{Pictures\cdf gamma, ig, tweedie.png}", sep=""), width=15, height=9, units='cm', dpi=600)


### PARTISI DATA

set.seed(17)
data_shuff <- df[sample(nrow(df)), ]
# write_csv(data_shuff, r"{G:\.shortcut-targets-by-id\1_UrMWiRHhBVY6BmkXpUmghuXVDupcUm8\Jason\Datasets\Shuffled Medical Cost Personal.csv}")
training_index <- seq(1, nrow(df)*0.8, 1)
train_data <- data_shuff[training_index, ]
test_data <- data_shuff[-training_index, ]



### PENENTUAN BASE LEVEL

library(dplyr)
train_data <- train_data %>%
  mutate(sex = relevel(factor(sex), ref = "male"))
train_data <- train_data %>%
  mutate(smoker = relevel(factor(smoker), ref = "no"))
train_data <- train_data %>%
  mutate(region = relevel(factor(region), ref = "southeast"))
train_data <- train_data %>%
  mutate(bmi2 = relevel(factor(bmi3), ref = "class I obese"))

### PEMODELAN
# Full
mod1 <- glm(charges ~ age + sex + bmi3 + children + smoker + region,
            data = train_data,
            family = tweedie(var.power = est_p, link.power = 0))
summary(mod1, dispersion=est_phi)
AICtweedie(mod1)
# Remove sex
mod2 <- glm(charges ~ age + bmi3 + children + smoker + region,
            data = train_data,
            family = tweedie(var.power = est_p, link.power = 0))
summary(mod2, dispersion=est_phi)
AICtweedie(mod2)
# Remove region
mod3 <- glm(charges ~ age + sex + bmi3 + children + smoker,
            data = train_data,
            family = tweedie(var.power = est_p, link.power = 0))
summary(mod3, dispersion=est_phi)
AICtweedie(mod3)

# Remove both
mod4 <- glm(charges ~ age + bmi3 + children + smoker,
            data = train_data,
            family = tweedie(var.power = est_p, link.power = 0))
AICtweedie(mod4)
summary(mod4, dispersion = est_phi) # mod2 adalah model terbaik

y <- train_data$charges
y_pred <- exp(predict(mod2, dispersion = est_phi, train_data))

sum(tweedie.dev(y, y_pred, est_p))

### PENGHITUNGAN NILAI RMSE PADA TEST SET

actual <- test_data$charges
pred <- exp(predict(mod2, dispersion = est_phi, test_data))
mse <- mean((actual-pred)^2)
(rmse <- sqrt(mse))
