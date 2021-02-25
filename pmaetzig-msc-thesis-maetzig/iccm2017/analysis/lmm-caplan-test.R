library(dplyr)
library(ggplot2)
library(lme4)
#library(plotly)
#library(plot3D)         ## for 3D scatterplots
#library(scatterplot3d)  ## for 3D scatterplots
#library(ggdendro)
library(factoextra)

## Caplan et al. (2015) data for subject vs. object relative sentences
caplan_data <- read.table('../data/CaplanEtAl_2015_data.txt', header=TRUE)
#xtabs(~subj+seg, caplan_data[caplan_data$grp=="EC",])
caplan_data <- subset(caplan_data, caplan_data$seg == "b")
caplan_data$condition <- as.factor(ifelse(caplan_data$sent == "SS", -1, 1))
#head(caplan_data)

# simple glmers, this doesn't give us the uncertainty of praameter estimates
m_acc_controls <- glmer(acc ~ condition + (1|subj) + (1|item), family = "binomial", subset(caplan_data, caplan_data$grp=="EC"))
m_acc_iwa <- glmer(acc ~ condition + (1|subj) + (1|item), family = "binomial", subset(caplan_data, caplan_data$grp!="EC"))

### LOGISTIC REGRESSION

library(arm)
library(MASS)
library(car)

beetle <- read.table("./beetle.txt", header = T) 
deads <- c()
doses <- c()
for (j in 1:length(beetle$killed)) {
  deads <- append(deads, rep(1, beetle$killed[j]))
  deads <- append(deads, rep(0, beetle$number[j] - beetle$killed[j]))
  doses <- append(doses, rep(beetle$dose[j], beetle$number[j]))
}
beetle_zero_one <- cbind.data.frame(deads, doses)
colnames(beetle_zero_one) <- c('dead', 'dose')

beetle_m1 <- glm(dead ~ dose, binomial(logit), data = beetle_zero_one)

# preprocessing for Bayesian modeling
caplan_acc_data <- list(N_obs = length(caplan_data$acc),
                        N_subj = length(unique(caplan_data$subj)),
                        acc = caplan_data$acc,
                        x = caplan_data$condition)

library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

fit_caplan_wmc <- stan(file='./caplan-wmc.stan',
                       data = caplan_acc_data,
                       chains = 4,
                       iter = 2000)

alphas <- rstan::extract(fit_caplan_wmc)$alpha
betas <- rstan::extract(fit_caplan_wmc)$beta

