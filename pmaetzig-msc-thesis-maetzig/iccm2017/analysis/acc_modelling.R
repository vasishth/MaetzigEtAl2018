library(dplyr)
library(ggplot2)
library(plot3D)         ## for 3D scatterplots
library(scatterplot3d)  ## for 3D scatterplots
 
caplan_data <- read.table('../data/CaplanEtAl_2015_data.txt', header=TRUE)
model_data_SR <- read.csv('../data/response_accuracies_GG_SR_2017016_1736.csv', header=TRUE)
model_data_OR <- read.csv('../data/response_accuracies_GG_OR_2017018_1827.csv', header=TRUE)
caplan_data$sent <- as.character(caplan_data$sent)

caplan_data_SR <- filter(caplan_data, sent=='SS')
caplan_data_OR <- filter(caplan_data, sent=='SO')

# data are not balanced, see:
xtabs(~subj+item, caplan_data_SR)
xtabs(~subj+item, caplan_data_OR)

caplan_data_SR_acc <- summarise(group_by(caplan_data_SR, subj,item), sum(acc) / n())
#head(caplan_data_SR_acc)
#summary(caplan_data_SR_acc)
#xtabs(~subj+item, caplan_data_SR_acc)
colnames(caplan_data_SR_acc) <- c('subj', 'item', 'acc')
caplan_data_SR_acc <- summarise(group_by(caplan_data_SR_acc, subj), sum(acc) / n())
colnames(caplan_data_SR_acc)<-c('subj', 'acc')
summary(caplan_data_SR_acc)

caplan_data_OR_acc <- summarise(group_by(caplan_data_OR, subj,item), sum(acc) / n())
colnames(caplan_data_OR_acc) <- c('subj', 'item', 'acc')
#head(caplan_data_OR_acc)
#summary(caplan_data_OR_acc)
#xtabs(~subj+item, caplan_data_OR_acc)
caplan_data_OR_acc <- summarise(group_by(caplan_data_OR_acc, subj), sum(acc) / n())
colnames(caplan_data_OR_acc)<-c('subj', 'acc')
summary(caplan_data_OR_acc)

slots_SR <- matrix(rep(NA, dim(model_data_SR)[1]*dim(caplan_data_SR_acc)[1]), ncol=dim(caplan_data_SR_acc)[1])
slots_OR <- matrix(rep(NA, dim(model_data_OR)[1]*dim(caplan_data_OR_acc)[1]), ncol=dim(caplan_data_OR_acc)[1])

for (i in 1:length(caplan_data_SR_acc$subj)) {
  tmp_acc <- caplan_data_SR_acc$acc[i]
  slots_SR[,i] <- abs(model_data_SR$acc - tmp_acc)
}

for (i in 1:length(caplan_data_OR_acc$subj)) {
  tmp_acc <- caplan_data_OR_acc$acc[i]
  slots_OR[,i] <- abs(model_data_OR$acc - tmp_acc)
}

colnames(slots_SR) <- paste("s", caplan_data_SR_acc$subj, sep = "")
colnames(slots_OR) <- paste("s", caplan_data_OR_acc$subj, sep = "")

model_data_SR <- cbind(model_data_SR, slots_SR)
model_data_OR <- cbind(model_data_OR, slots_OR)

get_min_params <- function(data, subject) {
  ref <- which(data[c(subject)]==min(data[c(subject)]))
  params <- data[ref, c("GA", "DAT", "ANS")]
  return(params)
}

subjects_SR <- colnames(model_data_SR)[9:dim(model_data_SR)[2]]
subjects_OR <- colnames(model_data_OR)[9:dim(model_data_OR)[2]]

results_SR <- data.frame()
for (i in 1:length(subjects_SR)) {
  tmp <- get_min_params(model_data_SR, subjects_SR[i])
  tmp$subj <- rep(subjects_SR[i], dim(tmp)[1])
  results_SR <- rbind(results_SR, tmp, make.row.names=FALSE)
}

#subjects<-unique(results_SR$subj)
#for(i in 1:length(subjects)){
#print(dim(subset(results_SR,subj==subjects[i]))[1])
#}
  
results_OR <- data.frame()
for (i in 1:length(subjects_OR)) {
  tmp <- get_min_params(model_data_OR, subjects_OR[i])
  tmp$subj <- rep(subjects_OR[i], dim(tmp)[1])
  results_OR <- rbind(results_OR, tmp, make.row.names=FALSE)
}

## splitting the data by subject/object RC and controls/IWA
results_SR_controls <- results_SR[1:(min(which(results_SR$subj=="s54001"))-1), ]
results_SR_iwa <- results_SR[min(which(results_SR$subj=="s54001")):dim(results_SR)[1], ]
results_OR_controls <- results_OR[1:(min(which(results_OR$subj=="s54001"))-1), ]
results_OR_iwa <- results_OR[min(which(results_OR$subj=="s54001")):dim(results_OR)[1], ]

## averaging values for each subject:
control_subj<-unique(results_SR_controls$subj)
results_SR_c_mn <- data.frame()
for(i in control_subj){
      tmp<-subset(results_SR_controls,subj==i)
      temp<-data.frame(GA=mean(tmp$GA),DAT=mean(tmp$DAT),ANS=mean(tmp$ANS),subj=i)
      results_SR_c_mn<-rbind(results_SR_c_mn,
                        temp,make.row.names=FALSE)
}

control_subj<-unique(results_OR_controls$subj)
results_OR_c_mn <- data.frame()
for(i in control_subj){
      tmp<-subset(results_OR_controls,subj==i)
      temp<-data.frame(GA=mean(tmp$GA),DAT=mean(tmp$DAT),ANS=mean(tmp$ANS),subj=i)
      results_OR_c_mn<-rbind(results_OR_c_mn,
                        temp,make.row.names=FALSE)
}

iwa_subj<-unique(results_SR_iwa$subj)
results_SR_i_mn <- data.frame()
for(i in iwa_subj){
      tmp<-subset(results_SR_iwa,subj==i)
      temp<-data.frame(GA=mean(tmp$GA),DAT=mean(tmp$DAT),ANS=mean(tmp$ANS),subj=i)
      results_SR_i_mn<-rbind(results_SR_i_mn,
                        temp,make.row.names=FALSE)
}

iwa_subj<-unique(results_OR_iwa$subj)
results_OR_i_mn <- data.frame()
for(i in iwa_subj){
      tmp<-subset(results_OR_iwa,subj==i)
      temp<-data.frame(GA=mean(tmp$GA),DAT=mean(tmp$DAT),ANS=mean(tmp$ANS),subj=i)
      results_OR_i_mn<-rbind(results_OR_i_mn,
                        temp,make.row.names=FALSE)
}

results_SR_controls_avg <- results_SR_c_mn
results_OR_controls_avg <- results_OR_c_mn
results_SR_iwa_avg <- results_SR_i_mn
results_OR_iwa_avg <- results_OR_i_mn

# 
# # plot of model accuracies
# plot(1:dim(model_data_SR)[1], jitter(model_data_SR$acc, 10), pch=19, xlab="number of parameter vector", ylab="accuracy", main="Accuracy simulations, SR")
# plot(1:dim(model_data_OR)[1], jitter(model_data_OR$acc, 10), pch=19, xlab="number of parameter vector", ylab="accuracy", main="Accuracy simulations, OR")
# 
# # 3D scatterplots, using plot3D package
# results_SR_controls$subj <- as.integer(factor(results_SR_controls$subj))
# results_OR_controls$subj <- as.integer(factor(results_OR_controls$subj))
# 
# results_SR_iwa$subj <- as.integer(factor(results_SR_iwa$subj))
# results_OR_iwa$subj <- as.integer(factor(results_OR_iwa$subj))
# 
# ## Trying out how plot with subj number (as integers) looks like
# ## However, not very readable and too much visual information, and clustering
# ## information is lost
# #tmpSRcontrol <- cbind(jitter(results_SR_controls[,1], factor=10),
# #                      jitter(results_SR_controls[,2], factor=4),
# #                      jitter(results_SR_controls[,3], factor=4) )
# #s3d <- scatterplot3d(tmpSRcontrol[, 1:3], pch = "")
# #text(s3d$xyz.convert(tmpSRcontrol[, 1:3]), labels = results_SR_controls$subj,
# #     col = "steelblue",cex=1.5)
# 
# pdf(file="../manuscript/figures/SR_controls.pdf", paper="a4")
# #par(mfrow=c(2,2))
# scatter3D(jitter(results_SR_controls$GA,factor=4),
#           jitter(results_SR_controls$DAT,factor=4),
#           jitter(results_SR_controls$ANS), colvar=NULL,
#           col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")
# 
# # base R plots without colours and jitter (can be added at will)
# # plot(results_SR_controls$GA, results_SR_controls$DAT, pch=19, xlab="GA", ylab="DAT")
# # plot(results_SR_controls$GA, results_SR_controls$ANS, pch=19, xlab="GA", ylab="ANS")
# # plot(results_SR_controls$DAT, results_SR_controls$ANS, pch=19, xlab="DAT", ylab="ANS")
# print(p_SR_controls_GA_DAT <- ggplot(data=results_SR_controls, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, GA and DAT"))
# 
# print(p_SR_controls_GA_ANS <- ggplot(data=results_SR_controls, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, GA and ANS"))
# print(p_SR_controls_DAT_ANS <- ggplot(data=results_SR_controls, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, DAT and ANS"))
# dev.off()
# 
# pdf(file="../manuscript/figures/SR_iwa.pdf", paper="a4")
# #par(mfrow=c(2,2))
# scatter3D(jitter(results_SR_iwa$GA),
#           jitter(results_SR_iwa$DAT),
#           jitter(results_SR_iwa$ANS),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")
# 
# # base R plots without colours and jitter (can be added at will)
# # plot(results_SR_iwa$GA, results_SR_iwa$DAT, pch=19, xlab="GA", ylab="DAT")
# # plot(results_SR_iwa$GA, results_SR_iwa$ANS, pch=19, xlab="GA", ylab="ANS")
# # plot(results_SR_iwa$DAT, results_SR_iwa$ANS, pch=19, xlab="DAT", ylab="ANS")
# print(p_SR_control_GA_DAT <- ggplot(data=results_SR_iwa, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, GA and DAT"))
# print(p_SR_control_GA_ANS <- ggplot(data=results_SR_iwa, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, GA and ANS"))
# print(p_SR_control_DAT_ANS <- ggplot(data=results_SR_iwa, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, DAT and ANS"))
# dev.off()
# 
# pdf(file="../manuscript/figures/OR_controls.pdf", paper="a4")
# #par(mfrow=c(2,2))
# scatter3D(jitter(results_OR_controls$GA),
#           jitter(results_OR_controls$DAT),
#           jitter(results_OR_controls$ANS),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")
# 
# # base R plots without colours and jitter (can be added at will)
# # plot(results_OR_controls$GA, results_OR_controls$DAT, pch=19, xlab="GA", ylab="DAT")
# # plot(results_OR_controls$GA, results_OR_controls$ANS, pch=19, xlab="GA", ylab="ANS")
# # plot(results_OR_controls$DAT, results_OR_controls$ANS, pch=19, xlab="DAT", ylab="ANS")
# print(p_OR_controls_GA_DAT <- ggplot(data=results_OR_controls, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, GA and DAT"))
# print(p_OR_controls_GA_ANS <- ggplot(data=results_OR_controls, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, GA and ANS"))
# print(p_OR_controls_DAT_ANS <- ggplot(data=results_OR_controls, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, DAT and ANS"))
# dev.off()
# 
# pdf(file="../manuscript/figures/OR_iwa.pdf", paper="a4")
# #par(mfrow=c(2,2))
# scatter3D(jitter(results_OR_iwa$GA),
#           jitter(results_OR_iwa$DAT),
#           jitter(results_OR_iwa$ANS),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")
# 
# # base R plots without colours and jitter (can be added at will)
# # plot(results_OR_iwa$GA, results_OR_iwa$DAT, pch=19, xlab="GA", ylab="DAT")
# # plot(results_OR_iwa$GA, results_OR_iwa$ANS, pch=19, xlab="GA", ylab="ANS")
# # plot(results_OR_iwa$DAT, results_OR_iwa$ANS, pch=19, xlab="DAT", ylab="ANS")
# print(p_OR_iwa_GA_DAT <- ggplot(data=results_OR_iwa, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, GA and DAT"))
# print(p_OR_iwa_GA_ANS <- ggplot(data=results_OR_iwa, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, GA and ANS"))
# print(p_OR_iwa_DAT_ANS <- ggplot(data=results_OR_iwa, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, DAT and ANS"))
# dev.off()


## comparing distributions:
op<-par(mfrow=c(2,3),pty="s")
plot(density(results_SR_controls$GA),xlab="GA",main="controls, SR")
plot(density(results_SR_controls$DAT),xlab="DAT",main="controls, SR")
plot(density(results_SR_controls$ANS),xlab="ANS",main="controls, SR")

plot(density(results_SR_iwa$GA),xlab="GA",main="IWA, SR")
plot(density(results_SR_iwa$DAT),xlab="DAT",main="IWA, SR")
plot(density(results_SR_iwa$ANS),xlab="ANS",main="IWA, SR")

op<-par(mfrow=c(2,3),pty="s")
plot(density(results_OR_controls$GA),xlab="GA",main="controls, OR")
plot(density(results_OR_controls$DAT),xlab="DAT",main="controls, OR")
plot(density(results_OR_controls$ANS),xlab="ANS",main="controls, OR")

plot(density(results_OR_iwa$GA),xlab="GA",main="IWA, OR")
plot(density(results_OR_iwa$DAT),xlab="DAT",main="IWA, OR")
plot(density(results_OR_iwa$ANS),xlab="ANS",main="IWA, OR")

## distributions with truncation of parameters
## comparing distributions:
op<-par(mfrow=c(2,3),pty="s")
plot(density(results_SR_controls$GA, from=0.2, to=1.1),xlab="GA",main="controls, SR")
plot(density(results_SR_controls$DAT, from=0.05, to=0.1),xlab="DAT",main="controls, SR")
plot(density(results_SR_controls$ANS, from=0.15, to=0.45),xlab="ANS",main="controls, SR")

plot(density(results_SR_iwa$GA, from=0.2, to=1.1),xlab="GA",main="IWA, SR")
plot(density(results_SR_iwa$DAT, from=0.05, to=0.1),xlab="DAT",main="IWA, SR")
plot(density(results_SR_iwa$ANS, from=0.15, to=0.45),xlab="ANS",main="IWA, SR")

op<-par(mfrow=c(2,3),pty="s")
plot(density(results_OR_controls$GA, from=0.2, to=1.1),xlab="GA",main="controls, OR")
plot(density(results_OR_controls$DAT, from=0.05, to=0.1),xlab="DAT",main="controls, OR")
plot(density(results_OR_controls$ANS, from=0.15, to=0.45),xlab="ANS",main="controls, OR")

plot(density(results_OR_iwa$GA, from=0.2, to=1.1),xlab="GA",main="IWA, OR")
plot(density(results_OR_iwa$DAT, from=0.05, to=0.1),xlab="DAT",main="IWA, OR")
plot(density(results_OR_iwa$ANS, from=0.15, to=0.45),xlab="ANS",main="IWA, OR")

# mixtools
library(mixtools)

## In SRs, higher proportion of low GA in iwa than controls:
mixmodSRcGA <- normalmixEM(results_SR_controls$GA, 
                           lambda = c(1/2), ## inits 
                           mu = c(.3,.9), ## inits
                           sigma = 1)
summary(mixmodSRcGA)

mixmodSRiGA <- normalmixEM(results_SR_iwa$GA, 
                           lambda = c(1/2), ## inits 
                           mu = c(.2,.9), ## inits
                           sigma = 1)
summary(mixmodSRiGA)

## In ORs iwa have lower goal activations and more variance than controls:
mixmodORcGA <- normalmixEM(results_OR_controls$GA, 
                      lambda = c(1/2), ## inits 
                      mu = c(.04,.08), ## inits
                      sigma = 1)
summary(mixmodORcGA)

mixmodORiGA <- normalmixEM(results_OR_iwa$GA, 
                           lambda = c(1/2), ## inits 
                           mu = c(.04,.1), ## inits
                           sigma = 1)
summary(mixmodORiGA)

## stan

## truncations/bounds:
#GA in {0.2, 0.3, ..., 1.1}
#DAT in {0.05, 0.06, ..., 0.1}
#ANS in {0.15, 0.2, ..., 0.45}

## Note: None of the Bayesian stuff works yet
if(0){
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat<-list(y=results_SR_controls$GA,
          N=length(results_SR_controls$GA),
          K=2)

## test code:
option<-stanc("mix1.stan")
## run code:
fit<-stan(file="mix1.stan",data=dat,iter=2000,
          control = list(adapt_delta = 0.999))

## print results:
print(fit,probs=c(0.025,0.975))

dat<-list(y=results_OR_controls$GA,
          N=length(results_OR_controls$GA))

library(rjags)
mix1jags <- jags.model( file = "mix1.jag",
                        data = dat,
                        n.chains = 4,
                        n.adapt =20000 , quiet=T)

track.variables<-c("mu","s","p")

res <- coda.samples( mix1jags,
                     var = track.variables,
                     n.iter = 40000,
                     thin = 1 )

gelman.diag(res)
summary(res)
}

## hierarchical clustering
## shows greater clustering in controls vs aphasics

clustersSRc<-hclust(dist(results_SR_controls[,1:3],
                         method="euclidean"),method="single")

clustersSRi<-hclust(dist(results_SR_iwa[,1:3],
                         method="euclidean"),method="single")

plot(clustersSRi)

clustersORc<-hclust(dist(results_OR_controls[,1:3],
                         method="euclidean"),method="single")

plot(clustersORc)

clustersORi<-hclust(dist(results_OR_iwa[,1:3],
                         method="euclidean"),method="single")

plot(clustersORi)

## PCA
pcSRc<-prcomp(results_SR_controls[,1:3],scale=TRUE)
summary(pcSRc)
pcSRc$sdev
plot(pcSRc)
round(head(pcSRc$rotation, 5), 2)
pc<-predict(pcSRc)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcSRc)

pcSRi<-prcomp(results_SR_iwa[,1:3],scale=TRUE)
summary(pcSRi)
pcSRi$sdev
plot(pcSRi)
round(head(pcSRi$rotation, 5), 2)
pc<-predict(pcSRi)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcSRi)

pcORc<-prcomp(results_SR_controls[,1:3],scale=TRUE)
summary(pcORc)
pcORc$sdev
plot(pcORc)
round(head(pcORc$rotation, 5), 2)
pc<-predict(pcORc)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcORc)

pcORi<-prcomp(results_SR_iwa[,1:3],scale=TRUE)
summary(pcORi)
pcORi$sdev
plot(pcORi)
round(head(pcORi$rotation, 5), 2)
pc<-predict(pcORi)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcORi)

## combining controls and patient data:
results_SR_controls$type<-"control"
results_SR_iwa$type<-"iwa"

results_SR<-rbind(results_SR_controls,results_SR_iwa)

pcSR<-prcomp(results_SR[,1:3],scale=TRUE)
summary(pcSR)
pcSR$sdev
plot(pcSR)
round(head(pcSR$rotation, 5), 2)
pc<-predict(pcSR)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcSR)

clusters <- hclust(dist(results_SR[, 1:3]),method="mcquitty")
plot(clusters)

clusterCut <- cutree(clusters, 2)
## not able to identify aphasics well
table(clusterCut, results_SR$type)

results_OR_controls$type<-"control"
results_OR_iwa$type<-"iwa"

results_OR<-rbind(results_OR_controls,results_OR_iwa)

pcOR<-prcomp(results_OR[,1:3],scale=TRUE)
summary(pcOR)
pcOR$sdev
plot(pcOR)
round(head(pcOR$rotation, 5), 2)
pc<-predict(pcOR)
plot(pc[,1:2])
plot(pc[,2:3])
biplot(pcOR)

clusters <- hclust(dist(results_OR[, 1:3]),method="complete")
plot(clusters)

clusterCut <- cutree(clusters, 2)
## so so discrimination ability for aphasics:
table(clusterCut, results_OR$type)


## ggplot2
ggplot(results_SR_controls,aes(x=GA)) + geom_density(alpha=0.25)
ggplot(results_SR_iwa,aes(x=GA)) + geom_density(alpha=0.25)

ggplot(results_OR_controls,aes(x=GA)) + geom_density(alpha=0.25)
ggplot(results_OR_iwa,aes(x=GA)) + geom_density(alpha=0.25)

## FactoMineR
library(FactoMineR)
pca_SR<-PCA(results_SR[,1:3])
pca_SR$eig

barplot(pca_SR$eig[,"eigenvalue"], border = NA, col = "gray80", names.arg = rownames(pca_SR$eig))

## these are the *correlations* between the variables and the PC
round(pca_SR$var$coord[,1:2], 4)
## circle of correlations:
plot(pca_SR,choix="var")

## contributions of variables in each PC:
print(rbind(pca_SR$var$contrib,
            TOTAL = colSums(pca_SR$var$contrib)), print.gap = 3)

library(RColorBrewer)
# color palette
colpal = brewer.pal(n = 4, name = "Blues")[4:1]
barplot(t(pca_SR$var$contrib), beside = TRUE,
        border = NA, ylim = c(0, 90), col = colpal, legend.text = colnames(pca_ir$var$contrib), args.legend = list(x = "top", ncol = 4, bty = 'n'))
abline(h = 25, col = "#ff572255", lwd = 2)
