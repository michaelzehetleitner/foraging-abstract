#=========================================================#
######    Analysis Script for Experiment 1          #######
#=========================================================#

## Preamble    ####
# Check if packages are installed; if not, install them
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
data <- read.csv(file="Data_For_I.csv", header = T)

###*  Data cleaning  ####
agg.allNEW <- data %>%
  filter(state == "foraging") %>%
  group_by(participant = pbn, cond, trialNumber, scoreReduction, travelTime) %>%
  summarise(x = n()) %>%
  ungroup()
agg.allNEW$participant <- factor(agg.allNEW$participant)
agg.allNEW$cond <- factor(agg.allNEW$cond)

###* Filtering on trial level  ####
aggNEW <-agg.allNEW %>%
  group_by(participant, cond, scoreReduction, travelTime) %>%
  mutate(n.original = n()) %>% # count the original number of trials per group
  arrange(x) %>% # sort patchLeavingTimes
  mutate(lower = floor(n() * 0.05),
         upper = n() - floor(n()*0.05),
         length = length(floor(n() * 0.05)  : (n() - floor(n()*0.05) ))) %>%
  slice(floor(n() * 0.05)  : (n() - floor(n()*0.05) )) %>% # leave out the first and last 5% of trials
  mutate(n = n()) %>% # count the new number of trials
  summarize(mean.ori = mean(x),
            n = mean(n),
            n.original = mean(n.original),
            lower = mean(lower), 
            upper = mean(upper),
            length = mean(length)) %>%
  ungroup()

###* Filtering on participant level   ####
aggNEW <- aggNEW %>%
  mutate(remove = FALSE) %>%
  group_by(cond, scoreReduction, travelTime) %>%
  mutate(remove = mean.ori %in% boxplot.stats(mean.ori, coef = 2.0)$out) %>%
  filter(!remove) %>%
  mutate(condition = factor(ifelse(cond=="0.238-4905", "quality",  
                                    ifelse(cond=="0.15-1089", "travel", "baseline")), 
                            levels=c("travel", "baseline", "quality"),
                            labels=c("Short travel", "Baseline", "High quality"))) %>%
  rename(patchLeaving = mean.ori) %>%
  ungroup()
# number of participants included in analyses
NROW(aggNEW)

###* Sample characteristics  ####
sample_characteristics <- list(Gender = table((data %>% filter(pbn %in% unique(aggNEW$participant)) %>%
                                                 select(pbn, sex) %>% distinct())$sex), 
                               Age = data %>% filter(pbn %in% unique(aggNEW$participant)) %>%
                                 select(pbn, age) %>% distinct() %>%
                                 summarise(AgeMean = mean(age), 
                                           SD = sd(age),
                                           Min = min(age),
                                           Max = max(age)),
                               Handedness = table((data %>% filter(pbn %in% unique(aggNEW$participant)) %>%
                                                     select(pbn, handedness) %>% distinct())$handedness),
                               Vision = table((data %>% filter(pbn %in% unique(aggNEW$participant)) %>%
                                                 select(pbn, opticalaids) %>% distinct())$opticalaids))
### Descriptive results  ####
descriptives <- aggNEW %>%
  group_by(travelTime, scoreReduction, condition) %>%
  summarise(MLeaving = mean(patchLeaving),
            se = sd(patchLeaving)/sqrt(n()), 
            ci = qt(0.975, n()-1)*se, ci.low=MLeaving-ci, ci.high=MLeaving+ci,.groups = "drop") %>%
  arrange(travelTime, scoreReduction) 


#=========================================================#
#####          Hypothesis Testing                     #####
#=========================================================#

###* H1: Effect of travel time		 ####
# Subset the data frame to get the relevant number of foraging actions
x_long <- aggNEW$patchLeaving[aggNEW$condition=="Baseline"]
x_short <- aggNEW$patchLeaving[aggNEW$condition=="Short travel"]

# Calculate the Bayes factor for a one-sided t-test
bfH1 <- ttestBF(x = x_long, y = x_short, paired = FALSE, r = 1, nullInterval = c(0, Inf))
# Posterior distribution of the effect size
posteriorH1 <- summary(posterior(ttestBF(x = x_long, y = x_short, paired = FALSE, r = 1), iterations = 100000))
posterior_deltaH1 <- format(round(c(Mean=posteriorH1$statistics["delta","Mean"], posteriorH1$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)


###* H2: Effect of patch quality   ####
# Subset the data frame to get the relevant number of foraging actions
x_high <- aggNEW$patchLeaving[aggNEW$condition=="High quality"]
x_low <- aggNEW$patchLeaving[aggNEW$condition=="Baseline"]

# Calculate the Bayes factor for a one-sided t-test
bfH2 <- ttestBF(x = x_low, y = x_high, paired = FALSE, r = 1, nullInterval = c( 0, Inf))
# Posterior distribution of the effect size
posteriorH2 <- summary(posterior(ttestBF(x = x_low, y = x_high, paired = FALSE, r = 1), iterations = 100000))
posterior_deltaH2 <- format(round(c(Mean=posteriorH2$statistics["delta","Mean"], posteriorH2$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)

###* H3: Comparison to optimality ####
## Check against optimal foraging actions
# Subtract optimal number of foragin actions: 2 for short travel, 5 for baseline, 4 for high quality conditions
aggNEW <- mutate(aggNEW, 
                 opt_diff = patchLeaving - case_when(condition=="Short travel"~2, 
                                                     condition=="Baseline"~5, 
                                                     condition=="High quality"~4))

# Calculate the Bayes factor for the deviation of the mean values from the optimal value using `ttestBF()`
bfShortTravel <- ttestBF(aggNEW$opt_diff[aggNEW$condition=="Short travel"], mu = 0, r = 1)
bfLongTravelBaseline <- ttestBF(aggNEW$opt_diff[aggNEW$condition=="Baseline"], mu = 0, r = 1)
bfHighQuality <- ttestBF(aggNEW$opt_diff[aggNEW$condition=="High quality"], mu = 0, r = 1)

# Posterior distribution of effect sizes
posteriorH3ShortTravel <- summary(posterior(bfShortTravel, iterations = 100000))
posteriorH3Baseline <- summary(posterior(bfLongTravelBaseline, iterations = 100000))
posteriorH3HighQuality <- summary(posterior(bfHighQuality , iterations = 100000))
posterior_deltaH3ShortTravel <- format(round(c(Mean=posteriorH3ShortTravel$statistics["delta","Mean"], 
                                               posteriorH3ShortTravel$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
posterior_deltaH3Baseline <- format(round(c(Mean=posteriorH3Baseline$statistics["delta","Mean"], 
                                           posteriorH3Baseline$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)
posterior_deltaH3HighQuality <- format(round(c(Mean=posteriorH3HighQuality$statistics["delta","Mean"], 
                                               posteriorH3HighQuality$quantiles["delta",c("2.5%", "97.5%")]), 2), nsmall=2)


#=========================================================#
#####      Visualization of results in barplot        #####
#=========================================================#
text_size=12
###* Plot for observed mean foraging actions    ####
brackets_hyp = data.frame(x=rep(c(1, 1.97, 2.03, 3),each=2),
                          y=c(descriptives$ci.high+0.2,6.85)[c(1,4,4,2,2,4,4,3)],
                          group=rep(c("travel", "quality"), each=4))
BF_labels <- lapply(list(bfH1, bfH2), FUN=format_bfs, digits=2) %>% 
  unlist()
labels_df <- data.frame(x=c(1.5, 2.5), y=c(6.95, 6.95), 
                        label=BF_labels)
opt_df <- data.frame(x=c(0.5, 1.5, 1.5, 2.5, 2.5, 3.5), y=rep(c(2,5,4), each=2),
                     group=rep(c("travel", "baseline", "quality"), each=2))
pl_cond <- ggplot(descriptives)+
  geom_bar(aes(x=condition, y=MLeaving),
           stat="identity", fill="gray", color="black")+
  geom_errorbar(aes(x=condition,ymin=ci.low, ymax=ci.high), width=0.1)+
  geom_line(data=brackets_hyp, aes(x=x, y=y, group=group))+
  geom_line(data=opt_df, aes(x=x, y=y, group=group), color="red")+
  geom_text(data=labels_df, aes(x=x, y=y, label=label), parse = TRUE, vjust=0,
            size=text_size/(1.5*.pt))+
  xlab("Condition")+ylab("# Foraging Actions")+ylim(c(0, 7.5))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))

###* Plot for deviation to optimal number of foraging actions ####
descriptives$diffopt <- descriptives$MLeaving - c(2,5,4)
optBF_labels <- lapply(list(bfShortTravel, bfLongTravelBaseline, bfHighQuality), 
                       FUN=format_bfs, digits=2) %>%  unlist()
labelsopt_df <- data.frame(x=1:3, y=rep(2.1,3), 
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
  ylim(c(-0.7, 2.3216))+
  theme_minimal()+
  theme(text = element_text(size=text_size),
        axis.text = element_text(size=text_size))
ggpubr::ggarrange(pl_cond, pl_opt, nrow=1)
ggsave(file="fig_exp1.jpg", width=17, height=6, unit="cm", dpi=600)


#=========================================================#
#####      Print all numbers to console        #####
#=========================================================#

print(sample_characteristics)
descriptives
paste("Hypothesis 1 (Travel time): ",format_bfs(bfH1, 2), "; delta: ",
      paste(names(posterior_deltaH1), posterior_deltaH1,sep="=",collapse=", "), sep="")
paste("Hypothesis 2 (Quality): ",format_bfs(bfH2, 2), "; delta: ",
      paste(names(posterior_deltaH2), posterior_deltaH2,sep="=",collapse=", "), sep="")

paste("Hypothesis 3 (Short travel condition): ",format_bfs(bfShortTravel, 2), "; delta: ",
      paste(names(posterior_deltaH3ShortTravel), posterior_deltaH3ShortTravel,sep="=",collapse=", "), sep="")
paste("Hypothesis 3 (Baseline condition): ",format_bfs(bfLongTravelBaseline, 2), "; delta: ",
      paste(names(posterior_deltaH3Baseline), posterior_deltaH3Baseline,sep="=",collapse=", "), sep="")
paste("Hypothesis 3 (High Quality condition): ",format_bfs(bfHighQuality, 2), "; delta: ",
      paste(names(posterior_deltaH3HighQuality), posterior_deltaH3HighQuality,sep="=",collapse=", "), sep="")

