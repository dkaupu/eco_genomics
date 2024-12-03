library(pcadapt)
library(tidyverse)
library(ggplot2)
library(gridExtra)

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")
metaclim <- meta[,c("subpop","region","TempM","MDiurnalR","TempR","Prec","PrecDM","PrecWmQ")]
metaclim %>% group_by(subpop) 
metaclim <- plyr::arrange(metaclim, region, subpop)
unique(metaclim$subpop)

## even though I've arranged by region then subpop, when we use ggplot they are not visually
## grouping in that order, so the level variable is manually sorting them

#### Var = TempM
a1 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = TempM, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "A)",
       x = "Subpopulation",
       y = "Mean Temperature") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 

#### Var = MDiurnalR
a2 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = MDiurnalR, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "B)",
       x = "Subpopulation",
       y = "Mean Diurnal Range") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5))

#### Var = TempR
a3 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = TempR, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "C)",
       x = "Subpopulation",
       y = "Temperature Range") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 

#### Var = Prec
a4 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = Prec, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "D)",
       x = "Subpopulation",
       y = "Annual Precipitation") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 

#### Var = PrecDM
a5 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = PrecDM, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "E)",
       x = "Subpopulation",
       y = "Mean Precipiation during Driest Month") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 

#### Var = PrecWmQ
a6 <- ggplot(metaclim, aes(x = factor(subpop, level= c("GER1","GER2","GER3","GER4","GER5","GER6","GER7","GER8","LUX1","CB","CH","CO","DF","FC","FL",  
                                                 "FP","HV","JF","JV","LH","MC","MP","PR","RM","SH","SP","TH","WR","WV","NOR1","NOR2","NOR3","NOR4",
                                                 "NOR5","NOR6","NOR7","NOR8","NOR9","AP","BJ","CG","DR","HL","HO","HR","IC","IL","JC","KR","OE","SE",
                                                 "TR","WC","WE","WF","WL","SP1","SP2","SP3","SP4","SP5","SP6","FR1","FR2","FR3","FR4","FR5","FR6","SW1")), y = PrecWmQ, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "F)",
       x = "Subpopulation",
       y = "Mean Precipiation during Warmest Quarter") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=5)) 

##### combining into one mega graph!!! #####
## this mega plot allows the reader to try and look for any similar/different patterns
## between our climactic variables; just helps visualize the differences

combined_clims <- grid.arrange(a1, a2, a3, a4, a5, a6, nrow = 2, ncol = 3)
ggsave("~/projects/eco_genomics/group_project/figures/climate1.png",combined_clims, width=18, height=8, units="in")
