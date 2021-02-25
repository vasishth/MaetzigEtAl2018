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

## using get_min_params.R, helper function found in "../R/" dir
## returns all the parameter vector that minimise the optimisation criterion
source('../R/get_min_params.R')

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

## trying to visualise param combinations
results_SR_controls <- results_SR[1:(min(which(results_SR$subj=="s54001"))-1), ]
results_SR_iwa <- results_SR[min(which(results_SR$subj=="s54001")):dim(results_SR)[1], ]
results_OR_controls <- results_OR[1:(min(which(results_OR$subj=="s54001"))-1), ]
results_OR_iwa <- results_OR[min(which(results_OR$subj=="s54001")):dim(results_OR)[1], ]

## plot of model accuracies
#plot(1:dim(model_data_SR)[1], jitter(model_data_SR$acc, 10), pch=19, xlab="number of parameter vector", ylab="accuracy", main="Accuracy simulations, SR")
#plot(1:dim(model_data_OR)[1], jitter(model_data_OR$acc, 10), pch=19, xlab="number of parameter vector", ylab="accuracy", main="Accuracy simulations, OR")

## 3D scatterplots, using plot3D package
source("../R/jitter_scatter3D.R")
#results_SR_controls$subj <- as.integer(factor(results_SR_controls$subj))
#results_OR_controls$subj <- as.integer(factor(results_OR_controls$subj))

#results_SR_iwa$subj <- as.integer(factor(results_SR_iwa$subj))
#results_OR_iwa$subj <- as.integer(factor(results_OR_iwa$subj))

### Trying out how plot with subj number (as integers) looks like
### However, not very readable and too much visual information, and clustering
### information is lost
##tmpSRcontrol <- cbind(jitter(results_SR_controls[,1], factor=10),
##                      jitter(results_SR_controls[,2], factor=4),
##                      jitter(results_SR_controls[,3], factor=4) )
##s3d <- scatterplot3d(tmpSRcontrol[, 1:3], pch = "")
##text(s3d$xyz.convert(tmpSRcontrol[, 1:3]), labels = results_SR_controls$subj,
##     col = "steelblue",cex=1.5)

#
##scatter3D_to_2D_proj <- function(data, x, y, main_title) {
#  
#pdf(file="../manuscript/figures/SR_controls.pdf", paper="a4")
##par(mfrow=c(2,2))
#scatter3D(jitter(results_SR_controls$GA,factor=4),
#          jitter(results_SR_controls$DAT,factor=4),
#          jitter(results_SR_controls$ANS), colvar=NULL,
#          col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")

## base R plots without colours and jitter (can be added at will)
## plot(results_SR_controls$GA, results_SR_controls$DAT, pch=19, xlab="GA", ylab="DAT")
## plot(results_SR_controls$GA, results_SR_controls$ANS, pch=19, xlab="GA", ylab="ANS")
## plot(results_SR_controls$DAT, results_SR_controls$ANS, pch=19, xlab="DAT", ylab="ANS")
#print(p_SR_controls_GA_DAT <- ggplot(data=results_SR_controls, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, GA and DAT"))
#print(p_SR_controls_GA_ANS <- ggplot(data=results_SR_controls, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, GA and ANS"))
#print(p_SR_controls_DAT_ANS <- ggplot(data=results_SR_controls, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR controls, DAT and ANS"))
#dev.off()

#pdf(file="../manuscript/figures/SR_iwa.pdf", paper="a4")
##par(mfrow=c(2,2))
#scatter3D(jitter(results_SR_iwa$GA),
#          jitter(results_SR_iwa$DAT),
#          jitter(results_SR_iwa$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")

## base R plots without colours and jitter (can be added at will)
## plot(results_SR_iwa$GA, results_SR_iwa$DAT, pch=19, xlab="GA", ylab="DAT")
## plot(results_SR_iwa$GA, results_SR_iwa$ANS, pch=19, xlab="GA", ylab="ANS")
## plot(results_SR_iwa$DAT, results_SR_iwa$ANS, pch=19, xlab="DAT", ylab="ANS")
#print(p_SR_control_GA_DAT <- ggplot(data=results_SR_iwa, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, GA and DAT"))
#print(p_SR_control_GA_ANS <- ggplot(data=results_SR_iwa, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, GA and ANS"))
#print(p_SR_control_DAT_ANS <- ggplot(data=results_SR_iwa, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("SR iwa, DAT and ANS"))
#dev.off()

#pdf(file="../manuscript/figures/OR_controls.pdf", paper="a4")
##par(mfrow=c(2,2))
#scatter3D(jitter(results_OR_controls$GA),
#          jitter(results_OR_controls$DAT),
#          jitter(results_OR_controls$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")

## base R plots without colours and jitter (can be added at will)
## plot(results_OR_controls$GA, results_OR_controls$DAT, pch=19, xlab="GA", ylab="DAT")
## plot(results_OR_controls$GA, results_OR_controls$ANS, pch=19, xlab="GA", ylab="ANS")
## plot(results_OR_controls$DAT, results_OR_controls$ANS, pch=19, xlab="DAT", ylab="ANS")
#print(p_OR_controls_GA_DAT <- ggplot(data=results_OR_controls, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, GA and DAT"))
#print(p_OR_controls_GA_ANS <- ggplot(data=results_OR_controls, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, GA and ANS"))
#print(p_OR_controls_DAT_ANS <- ggplot(data=results_OR_controls, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR controls, DAT and ANS"))
#dev.off()

#pdf(file="../manuscript/figures/OR_iwa.pdf", paper="a4")
##par(mfrow=c(2,2))
#scatter3D(jitter(results_OR_iwa$GA),
#          jitter(results_OR_iwa$DAT),
#          jitter(results_OR_iwa$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g")

## base R plots without colours and jitter (can be added at will)
## plot(results_OR_iwa$GA, results_OR_iwa$DAT, pch=19, xlab="GA", ylab="DAT")
## plot(results_OR_iwa$GA, results_OR_iwa$ANS, pch=19, xlab="GA", ylab="ANS")
## plot(results_OR_iwa$DAT, results_OR_iwa$ANS, pch=19, xlab="DAT", ylab="ANS")
#print(p_OR_iwa_GA_DAT <- ggplot(data=results_OR_iwa, aes(x=GA, y=DAT)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, GA and DAT"))
#print(p_OR_iwa_GA_ANS <- ggplot(data=results_OR_iwa, aes(x=GA, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, GA and ANS"))
#print(p_OR_iwa_DAT_ANS <- ggplot(data=results_OR_iwa, aes(x=DAT, y=ANS)) + geom_jitter(aes(colour=subj)) + ggtitle("OR iwa, DAT and ANS"))
#dev.off()

# ## all 4 scatterplots, in a grid, for paper
# pdf(file="../manuscript/figures/3Dscatterplots.pdf")
# par(mfrow=c(2,2))
# scatter3D(jitter(results_SR_controls$GA,factor=4),
#           jitter(results_SR_controls$DAT,factor=4),
#           jitter(results_SR_controls$ANS), colvar=NULL,
#           col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS",main="Subject relatives, controls",
#           ticktype="detailed", bty="g")
# scatter3D(jitter(results_OR_controls$GA, factor=4),
#           jitter(results_OR_controls$DAT, factor=4),
#           jitter(results_OR_controls$ANS),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", main="Object relatives, controls",
#           ticktype="detailed", bty="g")
# scatter3D(jitter(results_SR_iwa$GA, factor=4),
#           jitter(results_SR_iwa$DAT, factor=4),
#           jitter(results_SR_iwa$ANS),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", main="Subject relatives, IWA",
#           ticktype="detailed", bty="g")
# scatter3D(jitter(results_OR_iwa$GA, factor=4),
#           jitter(results_OR_iwa$DAT, factor=4),
#           jitter(results_OR_iwa$ANS, factor=4),
#           colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#           xlab="GA", ylab="DAT", zlab="ANS", main="Object relatives, IWA",
#           ticktype="detailed", bty="g")
# dev.off()

## Test: taking averages of multiple parameter estimates, across subjects
## taking averages over all parameter vectors by subject yields:

results_SR_controls_avg <- summarise_each(group_by(results_SR_controls, subj), c("mean"))
results_SR_iwa_avg <- summarise_each(group_by(results_SR_iwa, subj), c("mean"))
results_OR_controls_avg <- summarise_each(group_by(results_OR_controls, subj), c("mean"))
results_OR_iwa_avg <- summarise_each(group_by(results_OR_iwa, subj), c("mean"))

library(plotly)
sr_controls_plotly <- plot_ly(results_SR_controls_avg, x=~GA, y=~DAT, z=~ANS) %>%
  add_markers()
print(sr_controls_plotly)
sr_iwa_plotly <- plot_ly(results_SR_iwa_avg, x=~GA, y=~DAT, z=~ANS) %>%
  add_markers()
print(sr_iwa_plotly)
or_controls_plotly <- plot_ly(results_OR_controls_avg, x=~GA, y=~DAT, z=~ANS) %>%
  add_markers()
print(or_controls_plotly)
or_iwa_plotly <- plot_ly(results_OR_iwa_avg, x=~GA, y=~DAT, z=~ANS) %>%
  add_markers()
print(or_iwa_plotly)

#par(mfrow=c(2,2))
#scatter3D(jitter(results_SR_controls_avg$GA,factor=4),
#          jitter(results_SR_controls_avg$DAT,factor=4),
#          jitter(results_SR_controls_avg$ANS), colvar=NULL,
#          col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", main="Subject relatives, controls_avg",
#          ticktype="detailed", bty="g")

#scatter3D(jitter(results_OR_controls_avg$GA),
#          jitter(results_OR_controls_avg$DAT),
#          jitter(results_OR_controls_avg$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", main="Object relatives, controls_avg",
#          ticktype="detailed", bty="g")

#scatter3D(jitter(results_SR_iwa_avg$GA),
#          jitter(results_SR_iwa_avg$DAT),
#          jitter(results_SR_iwa_avg$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", main="Subject relatives, iwa_avg",
#          ticktype="detailed", bty="g")

#scatter3D(jitter(results_OR_iwa_avg$GA),
#          jitter(results_OR_iwa_avg$DAT),
#          jitter(results_OR_iwa_avg$ANS),
#          colvar=NULL, col='black', pch = 19, cex = 1, phi=0,
#          xlab="GA", ylab="DAT", zlab="ANS", main="Object relatives, iwa_avg",
#          ticktype="detailed", bty="g")

## CLUSTER ANALYSIS

## K-Means clustering
## 
## using my_kmeans.R, a helper function found in ../R/ directory
source('../R/my_kmeans.R')

k_max <- 15
kmeans_SR_controls <- my_kmeans(results_SR_controls[,1:3], k_max=k_max, iter.max=100000, maintitle = "SR controls")
kmeans_OR_controls <- my_kmeans(results_SR_controls[,1:3], k_max=k_max, iter.max=100000, maintitle = "OR controls")
kmeans_SR_iwa <- my_kmeans(results_SR_iwa[,1:3], k_max=k_max, iter.max=100000, maintitle = "SR IWA")
kmeans_OR_iwa <- my_kmeans(results_SR_iwa[,1:3], k_max=k_max, iter.max=100000, maintitle = "OR IWA")

plot(1:k_max, kmeans_SR_controls$kmeans_totss)
kmeans_SR_controls$kmeans_results

plot(1:k_max, kmeans_SR_iwa$kmeans_totss)
plot(1:k_max, kmeans_OR_controls$kmeans_totss)
plot(1:k_max, kmeans_OR_iwa$kmeans_totss)

# kmeans() documentation recommends using multiple starting points, this is tested here
for (i in 100:110) {
  my_kmeans(results_SR_iwa[,1:3], k_max=k_max, iter.max=100000, nstart=i, maintitle = paste("SR iwa", ",", i))
}

par(mfrow=c(1,1))
## Hierarchical clustering
hclust_SR_controls <- hclust(dist(results_SR_controls[,1:3]), method='ward.D')
plot(hclust_SR_controls)
hclust_OR_controls <- hclust(dist(results_OR_controls[,1:3]), method='ward.D')
plot(hclust_OR_controls)
hclust_SR_iwa <- hclust(dist(results_SR_iwa[,1:3]), method='ward.D')
plot(hclust_SR_iwa)
hclust_OR_iwa <- hclust(dist(results_OR_iwa[,1:3]), method='ward.D')
plot(hclust_OR_iwa)

results_SR_controls$type<-"control"
results_SR_iwa$type<-"iwa"
results_SR<-rbind(results_SR_controls,results_SR_iwa)
hclust_SR <- hclust(dist(results_SR[,1:3]), method="ward.D")
plot(hclust_SR)
clusterCut <- cutree(hclust_SR, 2)
table(clusterCut, results_SR$type)
e

par(mfrow=c(2,1))
plot(hclust_OR_controls)
plot(hclust_OR_iwa)
par(mfrow=c(1,1))
plot(hclust_OR_controls)
plot(hclust_OR_iwa)
