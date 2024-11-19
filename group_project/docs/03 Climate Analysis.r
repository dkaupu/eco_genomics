library(pcadapt)
library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)


meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal.csv", row.names = "X")
metaclim <- meta[,c("subpop","region","TempM","TempR","Prec")]
metaclim %>% group_by(subpop) 
metaclim <- plyr::arrange(metaclim, region, subpop)
unique(metaclim$subpop)

ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = TempM, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Mean Temperature by Region",
       x = "Subpopulation",
       y = "Mean Temperature") +
  theme_minimal()

ggplot(metaclim, aes(x = subpop, y = Prec, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Annual Precipitation by Region",
       x = "Subpopulation",
       y = "Annual Precipitation") +
  theme_minimal()
