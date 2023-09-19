#=========================================================#
######    Analysis Script for Experiment 3          #######
######     (Honey pot /continuous foraging)         #######
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
#setwd("/Users/Breitenladner/Desktop/Experimente/1_1_Foraging_Continuous/")

#=========================================================#
#####       Data processing and preparation           #####
#=========================================================#

##* Read data  ####
dataE3 <- read.csv("Data_Honig.csv", header = TRUE)[,-1]

##* Re-code variables                         ####
dataE3 <- dataE3 %>% 
  # define conditions
  mutate(cond = factor(paste(parameterQuality, parameterTravel, parameterMain, sep=""), levels=c("001", "010", "100"),
                       labels=c("main", "travel", "quality"))) %>%
  # Combine and clean columns for dependent variable foraging time (patchLeavingTime)
  mutate(patchLeavingTime = as.numeric(gsub("NA", "", paste(mouse_3.time, mouse_4.time, mouse_5.time, sep="")))) %>%
  select(-starts_with("mouse")) %>%
  # recode condition variables
  mutate( 
    travelTime = ifelse(cond=="travel", "long", "short"),
    patchQuality = ifelse(cond == "quality", "high", "low"),
    condition = factor(cond, levels = c("travel", "main", "quality"),
      labels = c("Long travel", "Baseline", "High quality"))) %>%
  rename(age = v1_age, sex = v2_sex, handedness = v3_handedness, vision = v4_vision)


###* Filtering on trial level  ####
aggE3 <-dataE3 %>%
  group_by(participant, condition, patchQuality, travelTime) %>%
  mutate(n.original = n()) %>% # count the original number of trials per group
  arrange(patchLeavingTime) %>% # sort patchLeavingTimes
  mutate(lower = ceiling(n() * 0.05),
         upper = n() - floor(n()*0.05),
         length = length(ceiling(n() * 0.05)  : (n() - floor(n()*0.05) ))) %>%
  slice(ceiling(n() * 0.05)  : (n() - floor(n()*0.05) )) %>% # leave out the first and last 5% of trials
  mutate(n = n()) %>% # count the new number of trials
  summarize(patchLeavingTime = mean(patchLeavingTime),
            n = mean(n),
            n.original = mean(n.original),
            lower = mean(lower), 
            upper = mean(upper),
            length = mean(length), .groups = "drop") 

###* Filtering on participant level   ####
aggE3 <- aggE3 %>%
  mutate(remove = FALSE) %>%
  group_by(condition, patchQuality, travelTime) %>%
  mutate(low_ex_thr = boxplot.stats(patchLeavingTime, coef = 2.0)$stats[1],
         upp_ex_thr = boxplot.stats(patchLeavingTime, coef = 2.0)$stats[5],
    remove = patchLeavingTime %in% boxplot.stats(patchLeavingTime, coef = 2.0)$out) %>%
  group_by(participant) %>% 
  mutate(remove = any(remove)) # remove all conditions from outliers
filter(aggE3, remove)[,c("participant", "condition", "patchLeavingTime","low_ex_thr", "upp_ex_thr")]
aggE3 <- aggE3 %>% filter(!remove) %>%
  ungroup()
# sanity check
table((aggE3 %>% group_by(participant) %>% summarise(N=n()))$N)
# --> Remaining participants: 29 out of 30


###* Sample characteristics  ####
sample_characteristics <- (dataE3 %>% 
                             filter(participant %in% unique(aggE3$participant)) %>%
                             select(participant, age, sex, handedness, vision) %>% 
                             distinct() %>%
                             list(Gender=table(.$sex), 
                                  Age=summarise(., AgeMean=mean(age), SD=sd(age),
                                                Min = min(age), Max=max(age)),
                                  Handedness = table(.$handedness),
                                  Vision = table(.$vision)))[-1]
sample_characteristics

### Descriptive results  ####
descriptives <- Rmisc::summarySEwithin(aggE3,
                                       measurevar = "patchLeavingTime", 
                                       idvar="participant", 
                                       withinvars = "condition") %>%
  mutate(ci.low=patchLeavingTime-ci, ci.high=patchLeavingTime+ci)


#=========================================================#
#####          Hypothesis Testing                     #####
#=========================================================#

###* H1: Effect of travel time		 ####
# Subset the data frame to get the relevant number of foraging actions
x_long <- aggE3$patchLeavingTime[aggE3$condition == "Long travel"]
x_short <- aggE3$patchLeavingTime[aggE3$condition == "Baseline"]

# Calculate the Bayes factor for a one-sided t-test
bfH1 <- ttestBF(x = x_long, y = x_short, paired = TRUE, r = 1, nullInterval = c(0, Inf))
# Posterior distribution of the effect size
posteriorH1 <- summary(posterior(ttestBF(x = x_long, y = x_short, paired = FALSE, r = 1), iterations = 100000))
posterior_deltaH1 <- format(round(c(Mean=posteriorH1$statistics["delta","Mean"], posteriorH1$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)


###* H2: Effect of patch quality   ####
# Subset the data frame to get the relevant number of foraging actions
x_low <- aggE3$patchLeavingTime[aggE3$condition == "Baseline"]
x_high <- aggE3$patchLeavingTime[aggE3$condition == "High quality"]

# Calculate the Bayes factor for a one-sided t-test
bfH2 <- ttestBF(x = x_low, y =  x_high, paired = TRUE, r = 1, nullInterval = c( 0, Inf))
# Posterior distribution of the effect size
posteriorH2 <- summary(posterior(ttestBF(x = x_low, y = x_high, paired = FALSE, r = 1), iterations = 100000))
posterior_deltaH2 <- format(round(c(Mean=posteriorH2$statistics["delta","Mean"], posteriorH2$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)



###* H3: Comparison to optimality ####
## Check against optimal foraging time
# Subtract optimal number of foraging time (seconds): 
# 4 for long travel, 3 for baseline, 2 for high quality conditions
aggE3 <- mutate(aggE3, 
                opt_diff = patchLeavingTime - case_when(condition=="Long travel"~4, 
                                                        condition=="Baseline"~3, 
                                                        condition=="High quality"~2))
# Calculate the Bayes factor for the deviation of the mean values from the optimal value using `ttestBF()`
bfBaseline <- ttestBF(aggE3$opt_diff[aggE3$condition == "Baseline"], mu = 0, rscale=1)
bfQuality <- ttestBF(aggE3$opt_diff[aggE3$condition == "High quality"], mu = 0, rscale=1)
bfTravel <- ttestBF(aggE3$opt_diff[aggE3$condition == "Long travel"], mu = 0, rscale=1)

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
                          y=c(descriptives$ci.high+0.15,6.42)[c(1,4,4,2,2,4,4,3)],
                          group=rep(c("travel", "quality"), each=4))
BF_labels <- lapply(list(bfH1, bfH2), FUN=format_bfs, digits=2) %>% 
  unlist()
labels_df <- data.frame(x=c(1.5, 2.5), y=c(6.5, 6.5), 
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
  xlab("Condition")+ylab("# Foraging Time")+ylim(c(0, 6.74))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))

###* Plot for deviation to optimal number of foraging actions ####
descriptives$opt <- c(4,3,2)
descriptives$diffopt <- descriptives$patchLeavingTime - descriptives$opt
optBF_labels <- lapply(list(bfTravel, bfBaseline, bfQuality), 
                       FUN= format_bfs, digits=2) %>%
  unlist()
labelsopt_df <- data.frame(x=1:3, y=rep(2.7,3), 
                           label=optBF_labels)
pl_opt <- ggplot(descriptives)+
  geom_bar(aes(x=condition, y=diffopt),
           stat="identity", fill="gray", color="black")+
  geom_errorbar(aes(x=condition,ymin=diffopt-ci, 
                    ymax = diffopt+ci), width=0.1)+
  geom_text(data=labelsopt_df, aes(x=x, y=y, label=label), parse = TRUE, vjust=0,
            size=text_size/(1.5*.pt))+
  xlab("Condition")+ylab("Difference to optimal\nforaging time")+
  geom_hline(yintercept = 0, color="red")+
  ylim(c(0, 2.8))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))
ggpubr::ggarrange(pl_cond, pl_opt, nrow=1)
ggsave(file="fig_exp3.jpg", width=18, height=6, unit="cm", dpi=600)





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
paste("Hypothesis 3 (High Quality condition): ", format_bfs(bfQuality, 2), "; delta: ",
      paste(names(posterior_deltaH3Quality), posterior_deltaH3Quality,sep="=",collapse=", "), sep="")




