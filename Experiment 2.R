#=========================================================#
######    Analysis Script for Experiment 2          #######
#=========================================================#

## Preamble    ####
# Check if packages are installed; if not, install them
if (!require(Rmisc)) install.packages("Rmisc")
if (!require(BayesFactor)) install.packages("BayesFactor")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(tidyverse)) install.packages("tidyverse")

rm(list = ls())

# Define formatting function for BFs
format_bfs <- function(bf, digits) {
  if (!is.numeric(bf)) bf <- extractBF(bf)[1,1]
  if (bf > 1000) {
    BF_power = floor(log(bf, base=10))
    BF_coeff = round(bf/10^BF_power, digits = digits)
    BF_label=paste("BF[10]==",format(BF_coeff, nsmall=1), "*x*10^", BF_power, sep="") #"BF==",
  } else {
    BF_label <- paste("BF[10]==",format(round(bf, digits), nsmall=1), sep="")
  }
  return(BF_label)
}

# Set right working directory
#setwd("/Users/Breitenladner/Desktop/Experimente/Foraging_I/Auswertung/Data")

#=========================================================#
#####       Data processing and preparation           #####
#=========================================================#

##* Read data   ####
dataE2 <- read.csv(file="Data_For_Intensive_Work.csv", header = T)

###*  Data cleaning  ####
dataE2 <- dataE2 %>%
  filter(mouse_3.rightButton == 1) %>% # select only those lines in the data, where the there is no foraging-action but a change-patch decision
  mutate(travelTime = case_when( # recode conditions
    condition %in% c("baseline", "quality") ~ "short",
    condition == "travel" ~ "long"
  ),
  patchQuality = case_when(
    condition == "baseline" ~ "high",
    condition == "quality" ~ "low",
    condition == "travel" ~ "high"
  ),
  condition = factor(condition, levels=c("travel", "baseline", "quality"),
                     labels=c("Long travel", "Baseline", "Low quality"))) %>%
  mutate(across(c(handedness = v3_handedness, expName, sex = v2_sex, vision = v4_vision, condition, travelTime, patchQuality, participant, trial = trials_gesamt.thisRepN), as.factor)) %>% # change names and factor
  select(-where(~all(is.na(.)))) %>% # remove empty cols
   rename(age = v1_age, patchLeavingTime = counter)   %>% # rename 
select(-c(X.2, X.1, trials_gesamt.thisTrialN, trials_gesamt.thisN, 
           trials_gesamt.thisIndex, starts_with("mouse_"), 
           frameRate, date, trials_gesamt.thisTrial, v3_handedness, v2_sex, v4_vision, trials_gesamt.thisRepN))  # remove unused cols
  
###* Filtering on trial level  ####
aggE2 <-dataE2 %>%
group_by(participant, condition, patchQuality, travelTime) %>%
  mutate(n.original = n()) %>% # count the original number of trials per group
  arrange(patchLeavingTime) %>% # sort patchLeavingTimes
  mutate(lower = floor(n() * 0.05),
         upper = n() - floor(n()*0.05),
         length = length(floor(n() * 0.05)  : (n() - floor(n()*0.05) ))) %>%
  slice(floor(n() * 0.05)  : (n() - floor(n()*0.05) )) %>% # leave out the first and last 5% of trials
  mutate(n = n()) %>% # count the new number of trials
  summarize(patchLeavingTime = mean(patchLeavingTime),
            n = mean(n),
            n.original = mean(n.original),
            lower = mean(lower), 
            upper = mean(upper),
            length = mean(length)) %>%
  ungroup()

###* Filtering on participant level   ####
aggE2 <- aggE2 %>%
	  mutate(remove = FALSE) %>%
	  group_by(condition, patchQuality, travelTime) %>%
	  mutate(remove = patchLeavingTime %in% boxplot.stats(patchLeavingTime, coef = 2.0)$out) %>%
  filter(!remove) %>%
	  ungroup()

###* Sample characteristics  ####
sample_characteristics <- list(Gender = table((dataE2 %>% filter(participant %in% unique(aggE2$participant)) %>%
                                                 select(participant, sex) %>% distinct())$sex), 
                               Age = dataE2 %>% filter(participant %in% unique(aggE2$participant)) %>%
                                 select(participant, age) %>% distinct() %>%
                                 summarise(AgeMean = mean(age), 
                                           SD = sd(age),
                                           Min = min(age),
                                           Max = max(age)),
                               Handedness = table((dataE2 %>% filter(participant %in% unique(aggE2$participant)) %>%
                                                     select(participant, handedness) %>% distinct())$handedness),
                               Vision = table((dataE2 %>% filter(participant %in% unique(aggE2$participant)) %>%
                                                 select(participant, vision) %>% distinct())$vision))
### Descriptive results  ####
descriptives <- Rmisc::summarySEwithin(aggE2,
                                       measurevar = "patchLeavingTime", 
                                       idvar="participant", 
                                       withinvars = "condition") %>%
  mutate(ci.low=patchLeavingTime-ci, ci.high=patchLeavingTime+ci)


#=========================================================#
#####          Hypothesis Testing                     #####
#=========================================================#

###* H1: Effect of travel time		 ####
# Subset the data frame to get the relevant number of foraging actions
	x_long <- aggE2$patchLeavingTime[aggE2$condition == "Long travel"]
	x_short <- aggE2$patchLeavingTime[aggE2$condition == "Baseline"]
	
	# Calculate the Bayes factor for a one-sided t-test
	bfH1 <- ttestBF(x = x_long, y = x_short, paired = TRUE, r = 1, nullInterval = c(0, Inf))
	# Posterior distribution of the effect size
	posteriorH1 <- summary(posterior(ttestBF(x = x_long, y = x_short, paired = FALSE, r = 1), iterations = 100000))
	posterior_deltaH1 <- format(round(c(Mean=posteriorH1$statistics["delta","Mean"], posteriorH1$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
	
	
	
	###* H2: Effect of patch quality   ####
	# Subset the data frame to get the relevant number of foraging actions
	x_high <- aggE2$patchLeavingTime[aggE2$condition == "Baseline"]
	x_low <- aggE2$patchLeavingTime[aggE2$condition == "Low quality"]
	
	# Calculate the Bayes factor for a one-sided t-test
	bfH2 <- ttestBF(x = x_high, y =  x_low, paired = TRUE, r = 1, nullInterval = c( 0, Inf))
	# Posterior distribution of the effect size
	posteriorH2 <- summary(posterior(ttestBF(x = x_high, y = x_low, paired = FALSE, r = 1), iterations = 100000))
	posterior_deltaH2 <- format(round(c(Mean=posteriorH2$statistics["delta","Mean"], posteriorH2$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
	

	###* H3: Comparison to optimality ####
	## Check against optimal foraging actions
	# Subtract optimal number of foragin actions: 4 for long travel, 3 for baseline, 2 for low quality conditions
	aggE2 <- mutate(aggE2, 
	                 opt_diff = patchLeavingTime - case_when(condition=="Long travel"~4, 
	                                                     condition=="Baseline"~3, 
	                                                     condition=="Low quality"~2))
	# Calculate the Bayes factor for the deviation of the mean values from the optimal value using `ttestBF()`
	bfBaseline <- ttestBF(aggE2$opt_diff[aggE2$condition == "Baseline"], mu = 0, rscale=1)
	bfQuality <- ttestBF(aggE2$opt_diff[aggE2$condition == "Low quality"], mu = 0, rscale=1)
	bfTravel <- ttestBF(aggE2$opt_diff[aggE2$condition == "Long travel"], mu = 0, rscale=1)
	
	# Posterior distribution of effect sizes
	posteriorH3Travel <- summary(posterior(bfTravel, iterations = 100000))
	posteriorH3Baseline <- summary(posterior(bfBaseline, iterations = 100000))
	posteriorH3Quality <- summary(posterior(bfQuality , iterations = 100000))
	posterior_deltaH3Travel <- format(round(c(Mean=posteriorH3Travel$statistics["delta","Mean"], 
	                                               posteriorH3Travel$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
	posterior_deltaH3Baseline <- format(round(c(Mean=posteriorH3Baseline$statistics["delta","Mean"], 
	                                            posteriorH3Baseline$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
	posterior_deltaH3Quality <- format(round(c(Mean=posteriorH3Quality$statistics["delta","Mean"], 
	                                               posteriorH3Quality$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
	
	
	#=========================================================#
	#####      Visualization of results in barplot        #####
	#=========================================================#
	text_size=12
	###* Plot for observed mean foraging actions    ####
brackets_hyp = data.frame(x=rep(c(1, 1.97, 2.03, 3),each=2),
                      y=c(descriptives$ci.high+0.2,5.6)[c(1,4,4,2,2,4,4,3)],
                      group=rep(c("travel", "quality"), each=4))
BF_labels <- lapply(list(bfH1, bfH2), FUN=format_bfs, digits=2) %>% 
  unlist()
labels_df <- data.frame(x=c(1.5, 2.5), y=c(5.7, 5.7), 
                        label=BF_labels)
opt_df <- data.frame(x=c(0.5, 1.5, 1.5, 2.5, 2.5, 3.5), y=rep(c(4, 3, 2), each=2),
                     group=rep(c("travel", "baseline", "quality"), each=2))

pl_cond <- ggplot(descriptives)+
  geom_bar(aes(x=condition, y=patchLeavingTime),
           stat="identity", fill="gray", color="black")+
  geom_errorbar(aes(x=condition,ymin=ci.low, ymax = ci.high), width=0.1)+
  geom_line(data=brackets_hyp, aes(x=x, y=y, group=group))+
  geom_line(data=opt_df, aes(x=x, y=y, group=group), color="red")+
  geom_text(data=labels_df, aes(x=x, y=y, label=label), parse = TRUE, vjust=0,
             size=text_size/(1.5*.pt))+
  xlab("Condition")+ylab("# Foraging Actions")+ylim(c(0, 6.1))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))

###* Plot for deviation to optimal number of foraging actions ####

descriptives$diffopt <- descriptives$patchLeavingTime - c(4,3,2)
optBF_labels <- lapply(list(bfTravel, bfBaseline, bfQuality), 
                       FUN= format_bfs, digits=2) %>%
  unlist()
labelsopt_df <- data.frame(x=1:3, y=rep(1.4,3), 
                        label=optBF_labels)
pl_opt <- ggplot(descriptives)+
  geom_bar(aes(x=condition, y=diffopt),
           stat="identity", fill="gray", color="black")+
  geom_errorbar(aes(x=condition,ymin=diffopt-ci, 
                 ymax = diffopt+ci), width=0.1)+
  geom_text(data=labelsopt_df, aes(x=x, y=y, label=label), parse = TRUE, vjust=0,
             size=text_size/(1.5*.pt))+
  xlab("Condition")+ylab("Difference to optimal\nforaging actions")+
  geom_hline(yintercept = 0, color="red")+
  ylim(c(0, 1.5))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))
ggpubr::ggarrange(pl_cond, pl_opt, nrow=1)
ggsave(file="fig_exp2.jpg", width=18, height=6, unit="cm", dpi=600)





#=========================================================#
#####      Print all numbers to console        #####
#=========================================================#

print(sample_characteristics)
descriptives
paste("Hypothesis 1 (Travel time): ",format_bfs(bfH1, 2), "; delta: ",
      paste(names(posterior_deltaH1), posterior_deltaH1,sep="=",collapse=", "), sep="")
paste("Hypothesis 2 (Quality): ",format_bfs(bfH2, 2), "; delta: ",
      paste(names(posterior_deltaH2), posterior_deltaH2,sep="=",collapse=", "), sep="")

paste("Hypothesis 3 (Long travel condition): ",format_bfs(bfTravel, 2), "; delta: ",
      paste(names(posterior_deltaH3Travel), posterior_deltaH3Travel,sep="=",collapse=", "), sep="")
paste("Hypothesis 3 (Baseline condition): ",format_bfs(bfBaseline, 2), "; delta: ",
      paste(names(posterior_deltaH3Baseline), posterior_deltaH3Baseline,sep="=",collapse=", "), sep="")
paste("Hypothesis 3 (Low Quality condition): ", format_bfs(bfQuality, 2), "; delta: ",
      paste(names(posterior_deltaH3Quality), posterior_deltaH3Quality,sep="=",collapse=", "), sep="")


