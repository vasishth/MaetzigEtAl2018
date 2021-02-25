library(ACTRPredictions2016)
library(ggplot2)

results <- preprocess_data("data/GG-EXP1-20161212-0343-GA-ANS/")
results <- do.call("rbind", results)

results_EmbV_OR <- subset(results, results$roi=="EmbV" & results$condition=="OR")
results_EmbV_SR <- subset(results, results$roi=="EmbV" & results$condition=="SR")
TFT_EmbV <- results_EmbV_OR[, 5:length(colnames(results_EmbV_OR))]
TFT_EmbV$mu <- results_EmbV_OR$`mean(TFT)` - results_EmbV_SR$`mean(TFT)`
TFT_EmbV$ans <- round(TFT_EmbV$ans, digits=2)

#TFT_EmbV <- reshape2::melt(TFT_EmbV, c("ga", "ans"), c("mu"), variable.name="measure")
#TFT_EmbV <- reshape2::dcast(TFT_EmbV, ga ~ ans, mean)

myHeatMap <- ggplot(TFT_EmbV, aes(ga, ans))
myHeatMap <- myHeatMap + geom_tile(aes(fill=round(mu)), colour="white")
myHeatMap <- myHeatMap + scale_fill_gradient(low="white", high="brown")
(myHeatMap)

results_ga_ans_repl <- preprocess_data("data/GG-EXP1-20161219-0910-GA-ANS-replication/")
results_ga_ans_repl <- do.call("rbind", results_ga_ans_repl)
results_ga_ans_EmbV_OR_repl <- subset(results_ga_ans_repl, results_ga_ans_repl$roi=="EmbV" & results_ga_ans_repl$condition=="OR")
results_ga_ans_EmbV_SR_repl <- subset(results_ga_ans_repl, results_ga_ans_repl$roi=="EmbV" & results_ga_ans_repl$condition=="SR")
ga_ans_repl <- results_ga_ans_EmbV_OR_repl[, 5:length(colnames(results_ga_ans_EmbV_OR_repl))]
ga_ans_repl$mu <- results_ga_ans_EmbV_OR_repl$`mean(TFT)` - results_ga_ans_EmbV_SR_repl$`mean(TFT)`
ga_ans_repl$ans <- round(ga_ans_repl$ans, digits=2)

myHeatMap_ga_ans_repl <- ggplot(ga_ans_repl, aes(ga, ans))
myHeatMap_ga_ans_repl <- myHeatMap_ga_ans_repl + geom_tile(aes(fill=round(mu)), colour="white")
myHeatMap_ga_ans_repl <- myHeatMap_ga_ans_repl + scale_fill_gradient(low="white", high="brown")
(myHeatMap_ga_ans_repl)

results_ga_ans_fine <- preprocess_data("data/GG-EXP1-20161219-0303-GA-ANS-finer/")
results_ga_ans_fine <- do.call("rbind", results_ga_ans_fine)
results_ga_ans_EmbV_OR_fine <- subset(results_ga_ans_fine, results_ga_ans_fine$roi=="EmbV" & results_ga_ans_fine$condition=="OR")
results_ga_ans_EmbV_SR_fine <- subset(results_ga_ans_fine, results_ga_ans_fine$roi=="EmbV" & results_ga_ans_fine$condition=="SR")
ga_ans_fine <- results_ga_ans_EmbV_OR_fine[, 5:length(colnames(results_ga_ans_EmbV_OR_fine))]
ga_ans_fine$mu <- results_ga_ans_EmbV_OR_fine$`mean(TFT)` - results_ga_ans_EmbV_SR_fine$`mean(TFT)`
ga_ans_fine$ans <- round(ga_ans_fine$ans, digits=2)

myHeatMap_ga_ans_fine <- ggplot(ga_ans_fine, aes(ga, ans))
myHeatMap_ga_ans_fine <- myHeatMap_ga_ans_fine + geom_tile(aes(fill=round(mu)), colour="white")
myHeatMap_ga_ans_fine <- myHeatMap_ga_ans_fine + scale_fill_gradient(low="white", high="brown")
(myHeatMap_ga_ans_fine)
