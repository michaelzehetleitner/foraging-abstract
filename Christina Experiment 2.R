## Packages installieren bzw. laden

rm(list=ls()) # Löschen des gesamten Workspace

# Packages beim ersten Mal installieren, kann später weggelassen werden
#install.packages("ggplot2")
#install.packages("plyr")
#install.packages("lsr")
#update.packages("ggplot2")

library(ggplot2)
library(plyr)
library(lsr)
library(BayesFactor)

#########################################################################################################

## Vorbereitung der später folgenden deskriptive Analysen

## Summarizes data, handling within-subjects variables by removing inter-subject variability. 
## It will still work if there are no within-S variables. 
## Gives count, un-normed mean, normed mean (with same between-group mean), 
##   standard deviation, standard error of the mean, and confidence interval. 
## If there are within-subject variables, calculate adjusted values using method from Morey (2008). 
##   data: a data frame. 
##   measurevar: the name of a column that contains the variable to be summariezed 
##   betweenvars: a vector containing names of columns that are between-subjects variables 
##   withinvars: a vector containing names of columns that are within-subjects variables 
##   idvar: the name of a column that identifies each subject (or matched subjects) 
##   na.rm: a boolean that indicates whether to ignore NA's 
##   conf.interval: the percent range of the confidence interval (default is 95%) 


summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

########################################################################################

# Current Working Directory festlegen
setwd("/Users/Breitenladner/Desktop/Experimente/Foraging_Intensive_Work/Auswertung/Data")

# Daten aus csv einlesen
data <- read.csv(file="Data_For_Intensive_Work.csv", header = T)

#############################################################################################

## Daten aufbereiten
#### Berechnung der relevanten Kenngrößen und Ausreißerbereinigung

# Benötigte Variablen generieren
data <- subset(data, data$mouse_3.rightButton ==1) # Betrachtung der Wechsel, anderes ist irrelevant

# first for each trial & condition get the number of foraging actions
agg.all <- aggregate(data$counter, by = list(participant = data$participant,condition = data$condition,trial = as.factor(data$trials_gesamt.thisRepN)), FUN = mean)
agg.all$participant <- as.factor(agg.all$participant)
head(agg.all)


## Ausreißerbereinigung, zweistufig: 5 % extreme Werte pro Proband, dann "extreme" Probanden
# calculate mean termination time per participant and condition, exclude 5% extreme data 
agg <- ddply(agg.all, .(participant, condition), function(x){
  # rank order observations
  obs <- sort(x$x) # x ist die ursprüngliche Variable Counter aus agg.all
  # remove 5% smallest and largest
  n.remove <- floor(.05 * length(obs))
  obs <- obs[c(n.remove : (length(obs) - n.remove))] # Anzahl der Trials, die nach Bereinigung übrig bleiben
  # calculate mean
  mean.ori <- mean(obs)
  n.ori <- length(obs)
  
  data.frame(mean.ori = mean.ori, n.ori = n.ori)
})

# calculate outliers via boxplot per condition: 
op <- par(mfrow = c(2,2))
agg$remove <- FALSE

agg <- ddply(agg, .(condition), function(x){
  # do the plot
  boxplot(x$mean.ori, range = 2.0, ylim = c(0, 10), main = paste("N =", NROW(x)), do.out = TRUE);
  # remove outliers
  out <- boxplot.stats(x$mean.ori, coef = 2.0)$out
  if(length(out)>0){
    for (i in 1:length(out)){
      x$remove[x$mean.ori == out[i]] <- TRUE
    }
  }
  # return only the number of outliers
  x
})


# remove all conditions of that participant
out <- subset(agg, remove == TRUE)
for(i in 1:NROW(out)){
  agg$remove[agg$participant == out$participant[i]] <- TRUE
}

data.frame(agg$participant, agg$condition, agg$remove)
par(op)
agg$participant <- as.factor(agg$participant)


#############################################################################################

## Vorbereitung für die Analysen

# Erstellung der Variable travel = Travel Cost
agg$travel[agg$condition == "baseline"] <- "short"
agg$travel[agg$condition == "quality"] <- "short"
agg$travel[agg$condition == "travel"] <- "long"
agg$travel <- as.factor(agg$travel)

# Erstellung der Variable quality = Gain Curve
agg$quality[agg$condition == "baseline"] <- "high"
agg$quality[agg$condition == "quality"] <- "low"
agg$quality[agg$condition == "travel"] <- "high"
agg$quality <- as.factor(agg$quality)

## Boxplot gesamt
par(mfrow = c(1,1))
agg_boxplot <- boxplot(agg$mean.ori ~ as.factor(agg$condition) ,
                      xlab="Condition",
                      ylab="# Foraging Actions", ylim = c(0,8), range = 2.0)


###################################################################################

###  Hypothesen prüfen

## Prüfen von Hypothese 1 (Effekt Travel Time) und 2 (Effekt Patch Quality) & Optimalität

########## H1: Effekt Travel Time
par(mfrow=c(1,1))
e1 <- subset(agg, quality == "high"  & !remove) 

e1_bild <- ggplot(data=e1, aes(x=travel, y=mean.ori)) +
    ylim(0,8) +
    geom_point(size = 3, alpha = 0.25)
e1_bild

e1_boxplot <- boxplot(e1$mean.ori ~ e1$travel ,
                      xlab="Travel Time",
                      ylab="# Foraging Actions")

# Deskriptive Statistik
descr <- summarySEwithin(data=e1, measurevar = "mean.ori",  withinvars="travel",
                         idvar="participant", na.rm=FALSE, conf.interval=.95, .drop=TRUE) 
descr ###results

# 95% confidence intervals of the means
ggplot(descr, aes(x=factor(travel), y=mean.ori, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=mean.ori-ci, ymax=mean.ori+ci)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(0,8) +
  ylab("# Foraging Actions") +
  xlab("Travel Time")

# test within one-sided BF ttest
bf <- ttestBF(x = e1$mean.ori[e1$travel == "long"], y = e1$mean.ori[e1$travel == "short"], paired = TRUE, nullInterval=c(0, +Inf), rscale=1)
bf  ###results
chains = posterior(bf[1], iterations = 10000)
summary(chains)
diff <- e1$mean.ori[e1$travel == "long"] -  e1$mean.ori[e1$travel == "short"]
mean(diff) ###results
sd(diff) ###results

# Nochmal 2-sided gegen 0 für richtige Effektstärke
bf <- ttestBF(diff, mu = 0,rscale=1)
bf
chains = posterior(bf[1], iterations = 10000)
summary(chains) ###results


### Prüfen auf Optimalität
e1short <- subset(e1, travel == "short"  & !remove) 
e1long <- subset(e1, travel == "long"  & !remove) 

# e1short: check against optimum: 3 Mal sammeln
e1short.opt <- ddply(e1short, .(participant), function(x){
  data.frame(x = mean(x$mean.ori))
})
bf <- ttestBF(e1short.opt$x, mu = 3.0, rscale=1)
bf ###results
chains = posterior(bf, iterations = 10000)
summary(chains)

# Nochmal Differenz gegen 0 für richtige Effektstärke
diff_short <- e1short$mean.ori- rep(3,47)
bf <- ttestBF(diff_short, mu = 0, rscale=1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains) ###results

# e1long: check against optimum: 4 Mal sammeln
e1long.opt <- ddply(e1long, .(participant), function(x){
  data.frame(x = mean(x$mean.ori))
})
bf <- ttestBF(e1long.opt$x, mu = 4.0, rscale=1)
bf ###results
chains = posterior(bf, iterations = 10000)
summary(chains)

# Nochmal Differenz gegen 0 für richtige Effektstärke
diff_long <- e1long$mean.ori- rep(4,47)
bf <- ttestBF(diff_long, mu = 0, rscale=1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains) ###results


########## H2: Effekt Patch Quality
e2 <- subset(agg, travel == "short"  & !remove)

e2_boxplot <- boxplot(e2$mean.ori ~ e2$quality ,
                      xlab="Quality",
                      ylab="# Foraging Actions")

# Deskriptive Statistik
descr <- summarySEwithin(data=e2, measurevar = "mean.ori",  withinvars="quality",
                         idvar="participant", na.rm=FALSE, conf.interval=.95, .drop=TRUE) 
descr ###results

# 95% confidence intervals of the means
ggplot(descr, aes(x=factor(quality), y=mean.ori, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=mean.ori-ci, ymax=mean.ori+ci)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(0,8) +
  ylab("# Foraging Actions") +
  xlab("Patch Quality")

# test within one-sided BF ttest
bf <- ttestBF(x = e2$mean.ori[e2$quality == "high"], y = e2$mean.ori[e2$quality == "low"], paired = TRUE, nullInterval=c(0, +Inf), rscale = 1)
bf ###results
chains = posterior(bf[1], iterations = 10000)
summary(chains)
diff <- e2$mean.ori[e2$quality == "high"] -  e2$mean.ori[e2$quality == "low"]
mean(diff) ###results
sd(diff) ###results

# Nochmal 2-sided gegen 0 für richtige Effektstärke
bf <- ttestBF(diff, mu = 0, rscale=1)
bf
chains = posterior(bf[1], iterations = 10000)
summary(chains) ###results


# Prüfen auf Optimalität
e2low <- subset(e2, quality == "low"  & !remove) 
e2high <- subset(e2, quality == "high"  & !remove) 

# e2low: check against optimum: 2 Mal sammeln
e2low.opt <- ddply(e2low, .(participant), function(x){
  data.frame(x = mean(x$mean.ori))
})
bf <- ttestBF(e2low.opt$x, mu = 2.0, rscale=1)
bf ###results
chains = posterior(bf, iterations = 10000)
summary(chains)

# Nochmal Differenz gegen 0 für richtige Effektstärke
diff_low <- e2low$mean.ori- rep(2,47)
bf <- ttestBF(diff_low, mu = 0, rscale=1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains) ###results


# e2high: check against optimum: 3 Mal Sammeln
e2high.opt <- ddply(e2high, .(participant), function(x){
  data.frame(x = mean(x$mean.ori))
})
bf <- ttestBF(e2high.opt$x, mu = 3.0, rscale=1)
bf ###results
chains = posterior(bf, iterations = 10000)
summary(chains)

# Nochmal Differenz gegen 0 für richtige Effektstärke
diff_high <- e2high$mean.ori- rep(3,47)
bf <- ttestBF(diff_high, mu = 0, rscale=1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains) ###results

#######################################################################


# Between-subject confidence intervals
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Deskriptive Statistik
descr <- summarySE(data=e1, measurevar="mean.ori", groupvars="travel", na.rm=FALSE, conf.interval=.95, .drop=TRUE)
descr ###results
descr <- summarySE(data=e2, measurevar="mean.ori", groupvars="quality", na.rm=FALSE, conf.interval=.95, .drop=TRUE)
descr ###results

