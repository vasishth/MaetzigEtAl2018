library(reshape2)

dat <- read.table('rcdata.txt', header = TRUE, dec = ',')

# ste coding of 'seg' to region type
# This has to be done in two steps, once for the SR and once for the OR clauses
# better way?
dat$roi <- ifelse(dat$seg == "a", "NP1", 
                  ifelse(dat$seg == "b", "who",
                         ifelse(dat$seg == "c", ifelse(dat$senttype == "SO", "NP2", "V1"), 
                                ifelse(dat$seg == "d", ifelse(dat$senttype == "SO", "V1", "NP2"),
                                       ifelse(dat$seg == "e", "V2", 
                                              ifelse(dat$seg == "f", "NP3", NA))))))
dat$roi <- as.factor(dat$roi)
dat$rt <- as.numeric(dat$rt)

# sanity check (no NA should appear as factor level)
levels(dat$roi)

# playing around with reshape2 (melt+cast)
dat.molten <- melt(dat, id=c("subj", "sent", "roi"),
                   measure=c("rt"),
                   na.rm=TRUE)

# this doesn't work because of some dimension error
# maybe it's due to the float - integer diff.ce between RT and N?
dat.id <- dcast(dat.molten, subj+sent+roi ~ variable,
                function(x) c(rt = mean(x), N = length(x)) )

# this also doesn't work
dat.id <- dcast(dat.molten, subj+sent+roi ~ .,
                function(x) c(rt = mean(x), N = length(x)) )

# this works but it doesn't recognise the list as two columns
dat.id <- dcast(dat.molten, subj+sent+roi ~ variable,
               function(x) list(c(rt=mean(x), N=length(x))) )
head(dat.id)

# ... but ONLY list() doesn't work
dat.id <- dcast(dat.molten, subj+sent+roi ~ variable,
                function(x) list(rt=mean(x), N=length(x)) )