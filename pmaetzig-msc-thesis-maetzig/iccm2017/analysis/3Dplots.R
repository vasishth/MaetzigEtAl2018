library(dplyr)
library(ggplot2)
library(plotly)
library(plot3D)         ## for 3D scatterplots
library(scatterplot3d)  ## for 3D scatterplots
library(ggplot2)
library(ggdendro)
library(factoextra)

caplan_data <- read.table('../data/CaplanEtAl_2015_data.txt', header=TRUE)
model_data_SR <- read.csv('../data/response_accuracies_GG_SR_2017016_1736.csv', header=TRUE)
model_data_OR <- read.csv('../data/response_accuracies_GG_OR_2017018_1827.csv', header=TRUE)
caplan_data$sent <- as.character(caplan_data$sent)

caplan_data_SR <- filter(caplan_data, sent=='SS')
caplan_data_OR <- filter(caplan_data, sent=='SO')

# data are not balanced, see:
#xtabs(~subj+item, caplan_data_SR)
#xtabs(~subj+item, caplan_data_OR)

caplan_data_SR_acc <- summarise(group_by(caplan_data_SR, subj,item), sum(acc) / n())
#head(caplan_data_SR_acc)
#summary(caplan_data_SR_acc)
#xtabs(~subj+item, caplan_data_SR_acc)
colnames(caplan_data_SR_acc) <- c('subj', 'item', 'acc')
caplan_data_SR_acc <- summarise(group_by(caplan_data_SR_acc, subj), sum(acc) / n())
colnames(caplan_data_SR_acc)<-c('subj', 'acc')
#summary(caplan_data_SR_acc)

caplan_data_OR_acc <- summarise(group_by(caplan_data_OR, subj,item), sum(acc) / n())
colnames(caplan_data_OR_acc) <- c('subj', 'item', 'acc')
#head(caplan_data_OR_acc)
#summary(caplan_data_OR_acc)
#xtabs(~subj+item, caplan_data_OR_acc)
caplan_data_OR_acc <- summarise(group_by(caplan_data_OR_acc, subj), sum(acc) / n())
colnames(caplan_data_OR_acc)<-c('subj', 'acc')
#summary(caplan_data_OR_acc)

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

## splitting data by SR/OR condition and by control/IWA group
results_SR_controls <- results_SR[1:(min(which(results_SR$subj=="s54001"))-1), ]
results_SR_iwa <- results_SR[min(which(results_SR$subj=="s54001")):dim(results_SR)[1], ]
results_OR_controls <- results_OR[1:(min(which(results_OR$subj=="s54001"))-1), ]
results_OR_iwa <- results_OR[min(which(results_OR$subj=="s54001")):dim(results_OR)[1], ]

## averaging the values for each subject, renaming
if(1){
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

results_SR_controls<-results_SR_c_mn
results_OR_controls<-results_OR_c_mn
results_SR_iwa<-results_SR_i_mn
results_OR_iwa<-results_OR_i_mn
}

results_SR_controls$subj <- as.integer(factor(results_SR_controls$subj))
results_OR_controls$subj <- as.integer(factor(results_OR_controls$subj))

results_SR_iwa$subj <- as.integer(factor(results_SR_iwa$subj))
results_OR_iwa$subj <- as.integer(factor(results_OR_iwa$subj))

# combining the data into SR / OR dataframes
results_SR_controls$type <- 'control'
results_SR_iwa$type <- 'iwa'
results_SR <- rbind(results_SR_controls, results_SR_iwa)
results_SR$type <- as.factor(results_SR$type)

results_OR_controls$type <- "control"
results_OR_iwa$type <- "iwa"
results_OR <- rbind(results_OR_controls,results_OR_iwa)
results_OR$type <- as.factor(results_OR$type)

## 3D plots with colour coding for control/iwa in one plot

library(plot3Drgl)

#pdf(file="../manuscript/figures/SRaverage3D.pdf")
par(mar=c(5.1, 4.1, 4.1, 0))
scatter3D(jitter(results_SR$GA, factor=4),
          jitter(results_SR$DAT, factor=4),
          jitter(results_SR$ANS, factor=4), colvar=as.integer(results_SR$type),
          #col=ramp.col(col = c('red', 'black'), 2),
          col='black',
          colkey = FALSE,
          pch = c(rep(16, 46), rep(4, 56)), cex = 1, phi=0,
          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g",main="Subject relative clauses",
          xlim=c(0.2, 1.1), nticks=3)
legend('right', legend=c('control', 'IWA'), pch=c(16,4), box.lty=0, inset=c(0,-10), cex=1.2)
#dev.off()
plotrgl() ## run this to get an interactive view of the 3Dplot

#pdf(file="../manuscript/figures/ORaverage3D.pdf")
scatter3D(jitter(results_OR$GA, factor=4),
          jitter(results_OR$DAT, factor=4),
          jitter(results_OR$ANS, factor=4), #colvar=as.integer(results_OR$type),
          #col=ramp.col(col = c('red', 'black'), 2),
          col='black',
          colkey = FALSE,
          pch = c(rep(16, 46), rep(4, 56)), cex = 1, phi=0,
          xlab="GA", ylab="DAT", zlab="ANS", ticktype="detailed", bty="g",main="Object relative clauses",
          xlim=c(0.2, 1.1), nticks=3)
#dev.off()
#plotrgl() ## run this to get an interactive view of the 3Dplot


### I copied this from https://plot.ly/r/axes/#modifying-axes-for-3d-plots
## Create lists for axis properties
#f1 <- list(
#  family = "Arial, sans-serif",
#  size = 12)

#f2 <- list(
#  family = "Arial, sans-serif",
#  size = 11)

#axis_GA <- list(
#  title='GA',
#  range = c(0, 1.1),
#  titlefont = f1,
#  tickfont = f2,
#  tick0 = 0.2,
#  dtick = 0.4,
#  tickangle = 30,
#  showgrid = TRUE
#)

#axis_DAT <- list(
#  title='DAT',
#  range = c(0.04, 0.1),
#  titlefont = f1,
#  tickfont = f2,
#  tick0 = 0.05,
#  dtick = 0.02,
#  tickangle = -30,
#  showgrid = TRUE
#)

#axis_ANS <- list(
#  title='ANS',
#  range = c(0.15, 0.55),
#  titlefont = f1,
#  tickfont = f2,
#  tick0 = 0.15,
#  dtick = 0.1,
#  showgrid = TRUE
#)

#my_marker <- list(size = 4,
#              color = 'black')

#p1_SR_controls <- plot_ly(results_SR_controls, 
#                          x=jitter(results_SR_controls$GA, factor=4),
#                          y=jitter(results_SR_controls$DAT, factor=4),
#                          z=jitter(results_SR_controls$ANS, factor=4),
#                          marker = my_marker,
#                          showlegend=FALSE,
#                          scene = "scene1")
#p2_OR_controls <- plot_ly(results_OR_controls, 
#                          x=jitter(results_OR_controls$GA, factor=4),
#                          y=jitter(results_OR_controls$DAT, factor=4),
#                          z=jitter(results_OR_controls$ANS, factor=4),
#                          marker = my_marker,
#                          showlegend=FALSE,
#                          scene = "scene2")
#p3_SR_iwa <- plot_ly(results_SR_iwa, 
#                     x=jitter(results_SR_iwa$GA, factor=4),
#                     y=jitter(results_SR_iwa$DAT, factor=4),
#                     z=jitter(results_SR_iwa$ANS, factor=4),
#                     marker = my_marker,
#                     showlegend=FALSE,
#                     scene = "scene3")
#p4_OR_iwa <- plot_ly(results_OR_iwa, 
#                     x=jitter(results_OR_iwa$GA, factor=4),
#                     y=jitter(results_OR_iwa$DAT, factor=4),
#                     z=jitter(results_OR_iwa$ANS, factor=4),
#                     marker = my_marker,
#                     showlegend=FALSE,
#                     scene = "scene4")

#avg_3D <- subplot(p1_SR_controls, p2_OR_controls, p3_SR_iwa, p4_OR_iwa) %>% 
#  layout(scene1 = list(domain=list(x=c(0,0.5),y=c(0.5,1)), 
#                      xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                      camera = list(eye = list(x = 2, y = 1, z = 0.2)),
#                      aspectmode='cube'),
#         scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)), 
#                       xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                       camera = list(eye = list(x = 2, y = 1, z = 0.2)),
#                       aspectmode='cube'),
#         scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)), 
#                       xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                       camera = list(eye = list(x = 2, y = 1, z = 0.2)),
#                       aspectmode='cube'),
#         scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)), 
#                       xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                       camera = list(eye = list(x = 2, y = 1, z = 0.2)),
#                       aspectmode='cube'),
#         annotations = list(
#           list(x = 0.15, y = 0.45, showarrow = FALSE,
#                text = "SR, iwa",
#                font = list(family="Arial", size=18)),
#           list(x = 0.75, y = 0.45, showarrow = FALSE,
#                text = "OR, iwa",
#                font = list(family="Arial", size=18)),
#           list(x = 0.15, y = 0.97, showarrow = FALSE,
#                text = "SR, controls",
#                font = list(family="Arial", size=18)),
#           list(x = 0.8, y = 0.97, showarrow = FALSE,
#                text = "OR, controls",
#                font = list(family="Arial", size=18))))
#print(avg_3D)
##
##         annotations = list(
##           list(x = 0.15, y = 0.45, showarrow = FALSE,
##                text = "SR, iwa",
##                font = list(family="Arial", size=18)),
##           list(x = 0.75, y = 0.45, showarrow = FALSE,
##                text = "OR, iwa",
##                font = list(family="Arial", size=18)),
##           list(x = 0.15, y = 0.97, showarrow = FALSE,
##                text = "SR, controls",
##                font = list(family="Arial", size=18)),
##           list(x = 0.8, y = 0.97, showarrow = FALSE,
##                text = "OR, controls",
##                font = list(family="Arial", size=18)),
##           ###########################################
##           list(x = 0.445, y = 0.57, showarrow = FALSE,
##                text = "GA",
##                font = list(family="Arial", size=12)),
##           list(x = 0.15, y = 0.52, showarrow=FALSE,
##                text = "DAT",
##                font = list(family="Arial", size=12)),
##           list(x = -0.03, y = 0.75, showarrow=FALSE,
##                text = "ANS", textangle = -90,
##                font = list(family="Arial", size=12)),
##           ###########################################
##           list(x = 0.445, y = 0.05, showarrow = FALSE,
##                text = "GA",
##                font = list(family="Arial", size=12)),
##           list(x = 0.15, y = 0, showarrow=FALSE,
##                text = "DAT",
##                font = list(family="Arial", size=12)),
##           list(x = 0.47, y = 0.75, showarrow=FALSE,
##                text = "ANS", textangle = -90,
##                font = list(family="Arial", size=12)),
##           ###########################################
##           list(x = 0.95, y = 0.05, showarrow = FALSE,
##                text = "GA",
##                font = list(family="Arial", size=12)),
##           list(x = 0.65, y = 0, showarrow=FALSE,
##                text = "DAT",
##                font = list(family="Arial", size=12)),
##           list(x = -0.03, y = 0.2, showarrow=FALSE,
##                text = "ANS", textangle = -90,
##                font = list(family="Arial", size=12)),
##           ###########################################
##           list(x = 0.95, y = 0.57, showarrow = FALSE,
##                text = "GA",
##                font = list(family="Arial", size=12)),
##           list(x = 0.65, y = 0.52, showarrow=FALSE,
##                text = "DAT",
##                font = list(family="Arial", size=12)),
##           list(x = 0.47, y = 0.2, showarrow=FALSE,
##                text = "ANS", textangle = -90,
##                font = list(family="Arial", size=12))
##           )

### single plots
#p1_SR_controls <- plot_ly(results_SR_controls,
#                          x=jitter(results_SR_controls$GA, factor=4),
#                          y=jitter(results_SR_controls$DAT, factor=4),
#                          z=jitter(results_SR_controls$ANS, factor=4),
#                          marker = my_marker) %>%
#  layout(title='SR, controls', scene=list(xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#         camera = list(eye = list(x = 2, y = 1.5, z = 0.2))))
#print(p1_SR_controls)

#p2_OR_controls <- plot_ly(results_OR_controls,
#                          x=jitter(results_OR_controls$GA, factor=4),
#                          y=jitter(results_OR_controls$DAT, factor=4),
#                          z=jitter(results_OR_controls$ANS, factor=4),
#                          marker = my_marker) %>%
#  layout(scene=list(xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#         camera = list(eye = list(x = 2, y = 1.5, z = 0.2))))
#print(p2_OR_controls)

#p3_SR_iwa <- plot_ly(results_SR_iwa,
#                          x=jitter(results_SR_iwa$GA, factor=4),
#                          y=jitter(results_SR_iwa$DAT, factor=4),
#                          z=jitter(results_SR_iwa$ANS, factor=4),
#                          marker = my_marker) %>%
#  layout(scene=list(xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                    camera = list(eye = list(x = 2, y = 1.5, z = 0.2))))
#print(p3_SR_iwa)

#p4_OR_iwa <- plot_ly(results_OR_iwa,
#                          x=jitter(results_OR_iwa$GA, factor=4),
#                          y=jitter(results_OR_iwa$DAT, factor=4),
#                          z=jitter(results_OR_iwa$ANS, factor=4),
#                          marker = my_marker) %>%
#  layout(scene=list(xaxis=axis_GA, yaxis=axis_DAT, zaxis=axis_ANS,
#                    camera = list(eye = list(x = 2, y = 1.5, z = 0.2))))
#print(p4_OR_iwa)
