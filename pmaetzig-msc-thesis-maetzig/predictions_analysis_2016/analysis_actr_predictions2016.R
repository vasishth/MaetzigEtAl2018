library(ACTRPredictions2016)
library(dplyr)
library(ggplot2)
library(cowplot)


## different method: combine all dfs into a single large one, introduce additional
## id variable for the param value combination set, and plot via ggplot and
## facet_grid(.~sim_set)
#plot_by_region(results_ga_complete, maintitle="ga-0.2") +
#  facet_grid(.~sim_set)
#ggplot(results_ga_complete, aes(x=roi, y=results_ga_complete$`mean(TFT)`, shape=condition, group=condition)) 

## put this function into the ACTRPredictions package
## maybe add functionality to calc both, embV and MV with argument
calculate_effects_embV <- function(results_list, indices, param_combinations) {
  effects_list <- rep(0.0, length(indices))
  tmp <- rep(0, 4)
  
  for (i in 1:length(indices)) {
    for (j in 1:length(indices[[i]])) {
      k <- indices[[i]][j]
      tmp[j] <- results_list[[k]][1,]$`mean(TFT)` - results_list[[k]][2,]$`mean(TFT)`
    }
    effects_list[i] <- mean(tmp)
  }
  return(effects_list)
}

###################################################################################################
##############  MONDAY, 2016-12-12 ################################################################
###################################################################################################

## only GA (:ga 0.5 1.5 0.25)
results_ga_rep2 <- preprocess_data(path="data/GG-EXP1-20161212-0219-GA/")
indices_ga_rep2 <- get_reps_same_params(results_ga_rep2)[[1]]
param_combinations_ga_rep2 <- get_reps_same_params(results_ga_rep2)[[2]]
effect_embV_ga_rep2 <- calculate_effect_ORSR(results_ga_rep2, indices_ga_rep2, param_combinations_ga_rep2)
effects_ga_rep2 <- dplyr::data_frame(params=param_combinations_ga_rep2, `OR-SR`=effect_embV_ga_rep2)


## only DAT (:dat 0.05 0.1 0.025)
results_dat_rep2 <- preprocess_data(path="data/GG-EXP1-20161212-0243-DAT/")
indices_dat_rep2 <- get_reps_same_params(results_dat_rep2)[[1]]
param_combinations_dat_rep2 <- get_reps_same_params(results_dat_rep2)[[2]]
effect_embV_dat_rep2 <- calculate_effect_ORSR(results_dat_rep2, indices_dat_rep2, param_combinations_dat_rep2)
effects_dat_rep2 <- dplyr::data_frame(params=param_combinations_dat_rep2, `OR-SR`=effect_embV_dat_rep2)


## GA and DAT 
results_ga_dat_rep2 <- preprocess_data("data/GG-EXP1-20161212-0250-GA-DAT/")
indices_ga_dat_rep2 <- get_reps_same_params(results_ga_dat_rep2)[[1]]
param_combinations_ga_dat_rep2 <- get_reps_same_params(results_ga_dat_rep2)[[2]]
effect_embV_ga_dat_rep2 <- calculate_effect_ORSR(results_ga_dat_rep2, indices_ga_dat_rep2, param_combinations_ga_dat_rep2)
effects_ga_dat_rep2 <- dplyr::data_frame(params=param_combinations_ga_dat_rep2, `OR-SR`=effect_embV_ga_dat_rep2)



## GA and ANS
results_ga_ans_rep2 <- preprocess_data("data/GG-EXP1-20161212-0343-GA-ANS/")
indices_ga_ans_rep2 <- get_reps_same_params(results_ga_ans_rep2)[[1]]
param_combinations_ga_ans_rep2 <- get_reps_same_params(results_ga_ans_rep2)[[2]]
effect_embV_ga_ans_rep2 <- calculate_effect_ORSR(results_ga_ans_rep2, indices_ga_ans_rep2, param_combinations_ga_ans_rep2)
effects_ga_ans_rep2 <- dplyr::data_frame(params=param_combinations_ga_ans_rep2, `OR-SR`=effect_embV_ga_ans_rep2)

plot_list_ga_ans_rep2 <- vector(mode="list", length=length(indices_ga_ans_rep2))

for (i in 1:length(indices_ga_ans_rep2)) {
  tmp <- multiplot_to_list(results_ga_ans_rep2, indices_ga_ans_rep2[[i]])
  plot_list_ga_ans_rep2[[i]] <- assign(param_combinations_ga_ans_rep2[[i]], tmp)
}

cowplot::plot_grid(plotlist=plot_list_ga_ans_rep2[[1]])
cowplot::plot_grid(plotlist=plot_list_ga_ans_rep2[[26]])

## 

## DAT and ANS
results_dat_ans_rep2 <- preprocess_data("data/GG-EXP1-20161212-0525-DAT-ANS/")
indices_dat_ans_rep2 <- get_reps_same_params(results_dat_ans_rep2)[[1]]
param_combinations_dat_ans_rep2 <- get_reps_same_params(results_dat_ans_rep2)[[2]]
effect_embV_dat_ans_rep2 <- calculate_effect_ORSR(results_dat_ans_rep2, indices_dat_ans_rep2, param_combinations_dat_ans_rep2)
effects_dat_ans_rep2 <- dplyr::data_frame(params=param_combinations_dat_ans_rep2, `OR-SR`=effect_embV_dat_ans_rep2)


## GA, DAT and ANS
results_ga_dat_ans_rep2 <- preprocess_data("data/GG-EXP1-20161212-0626-GA-DAT-ANS/")
indices_ga_dat_ans_rep2 <- get_reps_same_params(results_ga_dat_ans_rep2)[[1]]
param_combinations_ga_dat_ans_rep2 <- get_reps_same_params(results_ga_dat_ans_rep2)[[2]]
effect_embV_ga_dat_ans_rep2 <- calculate_effect_ORSR(results_ga_dat_ans_rep2, indices_ga_dat_ans_rep2, param_combinations_ga_dat_ans_rep2)
effects_ga_dat_ans_rep2 <- dplyr::data_frame(params=param_combinations_ga_dat_ans_rep2, `OR-SR`=effect_embV_ga_dat_ans_rep2)

###################################################################################################
###################################################################################################


## 4 repetitions of (:ga 0.2 1.0 0.2)
## 2016-12-11: this time ensured that 
## MAS = 3.5
## MP = 1.5 (see below for MP = NIL)
results_ga_MP15 <- preprocess_data(path="data/GG-EXP1-20161211-1727-GA-MAS35-MP15/")
indices_ga_MP15 <- get_reps_same_params(results_ga_MP15)[[1]]
param_combinations_ga_MP15 <- get_reps_same_params(results_ga_MP15)[[2]]
plot_list_ga_MP15 <- vector(mode="list", length=length(indices_ga_MP15))

for (i in 1:length(indices_ga_MP15)) {
  tmp <- multiplot_to_list(results_ga_MP15, indices_ga_MP15[[i]])
  plot_list_ga_MP15[[i]] <- assign(param_combinations_ga_MP15[[i]], tmp)
}
cowplot::plot_grid(plotlist=plot_list_ga_MP15[[1]])
cowplot::plot_grid(plotlist=plot_list_ga_MP15[[5]])

## The effects, i.e. RT for OR - SR, are as follows:
effect_EmbV_ga <- calculate_effect_ORSR(results_ga_MP15, indices_ga_MP15, param_combinations_ga_MP15)
effects_GA <- dplyr::data_frame(parameter=param_combinations_ga_MP15, `OR - SR`=effect_EmbV_ga)

## 4 repetitions of (:ga 0.2 1.0 0.2)
## 2016-12-11:
## MAS = 3.5
## MP = NIL
results_ga_MPNIL <- preprocess_data(path="data/GG-EXP1-20161211-1750-GA-MAS35-MPNIL")
indices_ga_MPNIL <- get_reps_same_params(results_ga_MPNIL)[[1]]
param_combinations_ga_MPNIL <- get_reps_same_params(results_ga_MPNIL)[[2]]
plot_list_ga_MPNIL <- vector(mode="list", length=length(indices_ga_MPNIL))

for (i in 1:length(indices_ga_MPNIL)) {
  tmp <- multiplot_to_list(results_ga_MPNIL, indices_ga_MPNIL[[i]])
  plot_list_ga_MPNIL[[i]] <- assign(param_combinations_ga_MPNIL[[i]], tmp)
}
cowplot::plot_grid(plotlist=plot_list_ga_MPNIL[[1]])
cowplot::plot_grid(plotlist=plot_list_ga_MPNIL[[5]])

## Again, repetitions of ga (:ga 0.5 1.5 0.25), more resembling Felix PhD Thesis
## MAS = 3.5, MP = 1.5 -- DISCUSS ENORMOUS FAILE RATE FOR MP = NIL
results_ga_rep1 <- preprocess_data(path="data/GG-EXP1-20161211-2310-GA-rep1")
indices_ga_rep1 <- get_reps_same_params(results_ga_rep1)[[1]]
param_combinations_ga_rep1 <- get_reps_same_params(results_ga_rep1)[[2]]


## 4 repetitions of (:dat 0.05 0.1 0.01)
results_dat_4 <- preprocess_data(path="data/GG-EXP1-2016113-1820-DAT-4/")
indices_dat_4 <- get_reps_same_params(results_dat_4)[[1]]
param_combinations_dat_4 <- get_reps_same_params(results_dat_4)[[2]]

plot_list_dat_4 <- vector(mode="list", length=length(indices_dat_4))

for (i in 1:length(indices_dat_4)) {
  tmp <- multiplot_to_list(results_dat_4, indices_dat_4[[i]])
  plot_list_dat_4[[i]] <- assign(param_combinations_dat_4[[i]], tmp)
}

cowplot::plot_grid(plotlist=plot_list_dat_4[[6]])


# 4 repetitions of (:ans 0.15 0.3 0.01)
results_ans_4 <- preprocess_data(path="data/GG-EXP1-2016119-1700-ANS-4")
indices_ans_4 <- get_reps_same_params(results_ans_4)[[1]]
param_combinations_ans_4 <- get_reps_same_params(results_ans_4)[[2]]

# plot_list_ans_ <- multiplot_to_list(results_ans_4, indices_ans_4[[1]])

cowplot::plot_grid(plotlist=plot_list_ans4)

## 4 repetitions of ga-ans
results_ga_ans_4 <- preprocess_data(path="data/GG-EXP1-20161115-1817-GA-ANS-4/")
indices_ga_ans_4 <- get_reps_same_params(results_ga_ans_4)[[1]]
param_combinations_ga_ans_4 <- get_reps_same_params(results_ga_ans_4)[[2]]

## 4 repetitions of ga-dat 
results_ga_dat_4 <- preprocess_data(path="data/GG-EXP1-20161115-1716-GA-DAT-4/")
indices_ga_dat_4 <- get_reps_same_params(results_ga_dat_4)[[1]]
param_combinations_ga_dat_4 <- get_reps_same_params(results_ga_dat_4)[[2]]
print(param_combinations_ga_dat_4)

## 4 repetitions of ga-ans-dat
results_ga_ans_dat_4 <- preprocess_data(path="data/GG-EXP1-20161115-2106-GA-DAT-ANS-4/")
indices_ga_ans_dat_4 <- get_reps_same_params(results_ga_ans_dat_4)[[1]]
param_combinations_ga_ans_dat_4 <- get_reps_same_params(results_ga_ans_dat_4)[[2]]

###############################################################################
###############################################################################
### 2016-12-


