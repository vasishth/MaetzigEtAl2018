regionmeans_new <- function(dataset) {

  dataset.molten <- melt(dataset, id=c("subj", "sent", "roi"),
                         measure=c("rt"),
                         na.rm=TRUE)

  # exploit plyr::ddply to rearrange df with 'mean rt' and 'n' columns
  dataset.molten.plied <- ddply(dataset.molten, .(subj, sent, roi), summarise, rt=mean(value), N=length(value))

  # melt again to unify the variable
  dataset.doubleMolten <- melt(dataset.molten.plied, id=c("subj", "sent", "roi"),
                               measure=c("rt", "N"),
                               na.rm=TRUE)

  # casting, finally
  dataset.id <- dcast(dataset.doubleMolten, subj + sent + roi ~ variable)

  # change 'sent' to 'condition' (for regionplot function)
  colnames(dataset.id) <- c("subj", "condition", "roi", "rt", "N")

  # get grand mean
  dataset.id$GM <- mean(dataset.id$rt)

  # deviation from the grandmean after considering intra-subject variability
  # removing between subject variance
  dataset.id <- ddply(dataset.id, .(subj),
                      transform, rt.w = rt - mean(rt) + GM)

  temp <- melt(dataset.id, id.var=c("subj","condition","roi"),
               measure.var="rt.w")

  M.id.w <- ddply(temp, .(condition, roi), summarise, M=mean(value), SE=sd(value)/sqrt(length(value)), N=length(value))
  M.id.w$M <- round(M.id.w$M)
  M.id.w$SE <- round(M.id.w$SE)
  M.id.w
}
