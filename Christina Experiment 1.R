
#install.packages("BayesFactor")
#install.packages("ggplot2")
#install.packages("plyr")

rm(list = ls())

library(BayesFactor)
library(ggplot2)
library(plyr)

## Vorbereitung der sp?ter folgenden deskriptive Analysen inkl. Bargraphen (unter Experimente 1 & 2)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

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
########################################################################################
# Current Working Directory festlegen
#setwd("/Users/Breitenladner/Desktop/Experimente/Foraging_I/Auswertung/Data")

#Daten aus csv einlesen
data <- read.csv(file="Data_For_I.csv", header = T)



#############################################################################################


## Prüfen von Hypothese 1 (Effekt Travel Time) und 2 (Effekt Patch Quality) & Optimalität

#### Berechnung der relevanten Kenngr??en und Ausrei?erbereinigung

# MEAN NUMBER OF FORAGING-ACTIONS PER TRIAL (= pro K?stchen) PER CONDITION 

# first for each trial & condition get the number of foraging actions
agg.all <- evalq(aggregate(state, list(pbn = pbn, cond = cond,
                                       trialNumber = trialNumber, scoreReduction = scoreReduction, 
                                       travelTime = travelTime), length), subset(data, state=="foraging"))
agg.all$pbn <- factor(agg.all$pbn)

		
# calculate mean termination time per participant and condition, exclude 5% extreme data
agg <- ddply(agg.all, .(pbn, cond, scoreReduction = scoreReduction, travelTime = travelTime), function(x){
	# rank order observations
	obs <- sort(x$x)
	# remove 5% smallest and largest
	n.remove <- floor(.05 * length(obs))
	obs <- obs[c(n.remove : (length(obs) - n.remove))]
	# calculate mean
	mean.ori <- mean(obs)
	n.ori <- length(obs)
	data.frame(x = mean(obs), mean.ori = mean.ori, n.ori = n.ori, n = length(obs))
	})
	
	
## Ausreißerbereinigung	
	
	# calculate outliers via boxplot per condition: 
	op <- par(mfrow = c(2,2))
	agg$remove <- FALSE
	agg <- ddply(agg, .(cond), function(x){
			# do the plot
			boxplot(x$x, range = 2.0, ylim = c(0, 15), main = paste("N =", NROW(x)));
		# remove outliers
		out <- boxplot.stats(x$x, coef = 2.0)$out
		if(length(out)>0){
			for (i in 1:length(out)){
				x$remove[x$x == out[i]] <- TRUE
				
			}
		}
				# return only the number of outliers
			
			x
		})
		# remove all conditions of that participant
		out <- subset(agg, remove == TRUE)
		for(i in 1:NROW(out)){
			agg$remove[agg$pbn == out$pbn[i]] <- TRUE
		}
		table(agg$pbn, agg$remove)
		
	par(op)
	
	agg$pbn <- factor(agg$pbn)

	
### H1: Effekt Travel Time
	e1 <- subset(agg, scoreReduction == 0.15  & !remove) 
	

 # descriptive

  descr <- summarySE(data=e1, measurevar="mean.ori", groupvars="travelTime", na.rm=FALSE, conf.interval=.95, .drop=TRUE)
  descr ### results

  # 95% confidence intervals of the mean
  ggplot(descr, aes(x=factor(travelTime), y=mean.ori, group=1)) + 
    geom_errorbar(aes(ymin=mean.ori-ci, ymax=mean.ori+ci), width=.1) +
    geom_line() +
    geom_point() +
  ylab("Number of foraging actions") +
  xlab("Travel time")

		
	# test between one-sided BF ttest
	bf <- ttestBF(x = e1$mean.ori[e1$travelTime == 1089], y = e1$mean.ori[e1$travelTime == 4905], paired = FALSE, r = 1, nullInterval=c(-Inf,0))
	bf ###results
	# BF nochmal nicht one-sided, um richtige Effektstärke zu erhalten
  bf <- ttestBF(x = e1$mean.ori[e1$travelTime == 1089], y = e1$mean.ori[e1$travelTime == 4905], paired = FALSE, r = 1)
  bf 
  chains = posterior(bf, iterations = 10000)
  summary(chains) ### results
	diff <- e1$mean.ori[e1$travelTime == 1089] -  e1$mean.ori[e1$travelTime == 4905]
	mean(diff)### results
	sd(diff)### results
	
# Prüfen auf Optimalität
e1short <- subset(e1, travelTime == 1089  & !remove) 
e1long <- subset(e1, travelTime == 4905  & !remove) 

	# e1short: check against optimum: 2
	e1short.opt <- ddply(e1short, .(pbn), function(x){
		
		data.frame(x = mean(x$x))
	})
	
	dev <- e1short.opt$x -2
	
	bf <- ttestBF(dev, mu = 0, r = 1) 
	bf ### results
	chains = posterior(bf, iterations = 10000)
	summary(chains) ### results
	
  # e1long: check against optimum: 5
  e1long.opt <- ddply(e1long, .(pbn), function(x){
    
    data.frame(x = mean(x$x))
  })
  
  dev <- e1long.opt$x -5
  
  bf <- ttestBF(dev, mu = 0, r = 1)
  bf ### results
  chains = posterior(bf, iterations = 10000)
  summary(chains) ### results
	
	
###### H2: Effekt Patch Quality

e2 <- subset(agg, travelTime == 4905  & !remove) 


# descriptive

descr <- summarySE(data=e2, measurevar="mean.ori", groupvars="scoreReduction", na.rm=FALSE, conf.interval=.95, .drop=TRUE)
descr ###results

# 95% confidence intervals of the mean
ggplot(descr, aes(x=factor(scoreReduction), y=mean.ori, group=1)) + 
  geom_errorbar(aes(ymin=mean.ori-ci, ymax=mean.ori+ci), width=.1) +
  geom_line() +
  geom_point() +
  ylab("Number of foraging actions") +
  xlab("Gain per foraging action in %")


# test between one-sided BF ttest
bf <- ttestBF(x = e2$mean.ori[e2$scoreReduction ==0.238], y = e2$mean.ori[e2$scoreReduction ==0.15], paired = FALSE, r = 1, nullInterval=c(-Inf,0))
bf ### results
# BF nochmal nicht one-sided, um richtige Effektstärke zu erhalten
bf <- ttestBF(x = e2$mean.ori[e2$scoreReduction == 0.238], y = e2$mean.ori[e2$scoreReduction == 0.15], paired = FALSE, r = 1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains) ### results
diff <- e2$mean.ori[e2$scoreReduction == 0.15] -  e2$mean.ori[e2$scoreReduction == 0.238]
mean(diff) ### results
sd(diff)### results


# Prüfen auf Optimalität
e2lowqual <- subset(e2, scoreReduction == 0.15  & !remove) 
e2highqual <- subset(e2, scoreReduction == 0.238  & !remove) 

# e2lowqual: check against optimum: 5
e2lowqual.opt <- ddply(e2lowqual, .(pbn), function(x){
  
  data.frame(x = mean(x$x))
})

dev <- e2lowqual.opt$x -5

bf <- ttestBF(dev, mu = 0, r = 1)
bf
chains = posterior(bf, iterations = 10000)
summary(chains)

# e2highqual: check against optimum: 4
e2highqual.opt <- ddply(e2highqual, .(pbn), function(x){
  
  data.frame(x = mean(x$x))
})

dev <- e2highqual.opt$x -4

bf <- ttestBF(dev, mu = 0, r = 1)
bf ### results
chains = posterior(bf, iterations = 10000)
summary(chains) ### results

