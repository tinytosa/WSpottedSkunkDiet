#############
#summary stats

######
#load packages
require(ggplot2)

#####
#load data
data.all <- read.table("Data_Raw/dataall_meta_man_2022-04-18.txt", sep=",") #data with metabarcoding and manual sorting; only has 128 samples
sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193, only clupea
sk.final <- sk.final[!sk.final$ScatID %in% c("MT904","F37-193"),] #remove MT904 and F37-193 from sk.final, only clupea

#
classes <- sk.final[,c("Amphibia","Aves","Mammalia","Reptilia","Gastropoda","Arachnida","myriapoda","Insecta","Streptophyta")]
classes[classes > 0] <- 1
colSums(classes)

# Amphibia         Aves     Mammalia     Reptilia   Gastropoda    Arachnida    myriapoda      Insecta Streptophyta 
#       17           18           60            3            9            3           45           92           37

table(sk.final$scatseason)
colSums(classes[sk.final$scatseason == "wet",]) #71 samples
colSums(classes[sk.final$scatseason == "dry",]) #57 samples

table(sk.final$logged)
colSums(classes[sk.final$logged == 1,]) #57 samples
colSums(classes[sk.final$logged == 0,]) #71 samples

#number of vertebrate skunk scats
nrow(sk.final[sk.final$Aves > 0 | sk.final$Mammalia > 0 | sk.final$Amphibia > 0 | sk.final$Reptilia > 0,]) #75 samples
nrow(sk.final[sk.final$Aves > 0 | sk.final$Mammalia > 0 | sk.final$Amphibia > 0 | sk.final$Reptilia > 0,])/nrow(sk.final) # 0.5859 proportion of scats that were vertebrates

nrow(sk.final[sk.final$Mammalia > 0,])/nrow(sk.final) #0.46875 proportion of scats that were mammals
60/75
17/75
18/75

#number of invertebrate skunk scats
nrow(sk.final[sk.final$Gastropoda > 0 | sk.final$Arachnida > 0 | sk.final$myriapoda > 0 | sk.final$Insecta > 0,]) #109 samples
nrow(sk.final[sk.final$Gastropoda > 0 | sk.final$Arachnida > 0 | sk.final$myriapoda > 0 | sk.final$Insecta > 0,])/nrow(sk.final) #0.8516 proportion of scats that were invertebrates

nrow(data.all[grep(data.all$lineage, pattern="Vespula"),]) #73
73/128
73/109
nrow(data.all[grep(data.all$lineage, pattern="Diplopoda"),]) #44
44/109

#number of plant skunk scats
nrow(sk.final[sk.final$Streptophyta > 0,]) #37 samples
nrow(sk.final[sk.final$Streptophyta > 0,])/nrow(sk.final) #0.2891 proportion of scats that were plants

stat1 <- merge(sk.final, data.all[grep(data.all$lineage, pattern="Vespula"),], all.y=T, by.x="ScatID", by.y="sample")
stat1[stat1$locus == "manual",]
table(stat1[stat1$locus == "manual",]$scatseason)
# dry wet 
#  42  29 

stat2 <- merge(sk.final, data.all[grep(data.all$lineage, pattern="Diplopoda"),], all.y=T, by.x="ScatID", by.y="sample")
stat2[stat2$locus == "manual",]
table(stat2[stat2$locus == "manual",]$scatseason)
# dry wet 
#  21  22 

##############
#statistical comparison of seasons, logging
#all together
##############

#load packages
require(mvabund)
require(dplyr)
require(reshape2)

#load data
data.all <- read.table("Data_Raw/dataall_meta_man_2022-01-26.txt", sep=",") #data with metabarcoding and manual sorting; only has 128 samples
sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193
sk.final <- sk.final[!sk.final$ScatID %in% c("MT904","F37-193"),] #remove MT904 and F37-193 from sk.final

#######
#summarize dates of collection
sk.final$CollectionDate <- as.POSIXct(as.character(sk.final$CollectionDate), format="%Y-%m-%d %H:%M:%S")

sk.final$vert <- sk.final$Amphibia + sk.final$Aves + sk.final$Mammalia + sk.final$Reptilia
sk.final$invert <- sk.final$Gastropoda + sk.final$Arachnida + sk.final$myriapoda + sk.final$Insecta
sk.final$plant <- sk.final$Streptophyta
sk.final$season <- as.numeric(as.factor(sk.final$scatseason))-1 #summer = 0, fall = 1

#######
#run manyglm for vert, plant, inverts
data <- mvabund(sk.final[,c("vert","invert","plant")]) #taxonomic data
data[data > 0] <- 1

vars <- sk.final[,c("season","logged","hja","class")] #covariates

family="binomial" #if values are 1s and 0s
mod.s <- manyglm(formula = data ~ season, data=vars, family=family) 
mod.l <- manyglm(formula = data ~ logged, data=vars, family=family)

#additive
mod.sl <- manyglm(formula = data ~ season + logged, data=vars, family=family)
#interaction
# mod.sxl <- manyglm(formula = data ~ season*logged, data=vars, family="binomial")

#check assumptions, is the family correct? check to see if points are spread randomly
plot(mod.sl)

#interpret the results
anova(mod.sl)
#scatseason     126       1 10.101    0.023 *

#examine which taxa is different
anova(mod.sl, p.uni="adjusted")
#scatseason invert: 8.24    0.018

######
table(factor(data.all[data.all$season == "S" & data.all$k == "Viridiplantae",]$g))
table(factor(data.all[data.all$season == "F" & data.all$k == "Viridiplantae",]$g))

table(factor(data.all[data.all$season == "S" & data.all$p == "Chordata",]$g))
table(factor(data.all[data.all$season == "F" & data.all$p == "Chordata",]$g))

###############
#split data by class
vars <- sk.final[,c("season","logged","hja","class")] #covariates

#summary stats
sk <- sk.final[,c("Amphibia","Aves","Mammalia","Reptilia","Gastropoda","Arachnida","myriapoda","Insecta","Streptophyta")]
sk[sk > 0] = 1 #make presence/absence
sk[is.na(sk)] = 0 #replace NA with 0
colSums(sk)

data <- mvabund(sk)
mod.sl <- manyglm(formula = data ~ season + logged, data=vars, family="binomial")
plot(mod.sl)
anova(mod.sl)
#season:         126       1 32.04    0.002 **
#logged:         125       1 10.70    0.356   

# Model: data ~ scatseason + logged
# Multivariate test:
#               Res.Df Df.diff   Dev Pr(>Dev)    
#   (Intercept)    127                           
#   scatseason     126       1 32.04    0.001 ***
#   logged         125       1 10.70    0.364    

anova(mod.sl, p.uni="adjusted")
#Insecta by scatseason: 25.395    0.001
#Mammalia by scatseason: 5.792    0.126

# Multivariate test:
#               Res.Df Df.diff   Dev Pr(>Dev)    
#   (Intercept)    127                           
#   scatseason     126       1 32.04    0.001 ***
#   logged         125       1 10.70    0.378    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Univariate Tests:
#                      Amphibia           Aves          Mammalia          Reptilia          Gastropoda          Arachnida          myriapoda          Insecta          Streptophyta         
#                  Dev Pr(>Dev)   Dev Pr(>Dev)      Dev Pr(>Dev)      Dev Pr(>Dev)        Dev Pr(>Dev)       Dev Pr(>Dev)       Dev Pr(>Dev)     Dev Pr(>Dev)          Dev Pr(>Dev)
# (Intercept)                                                                                                                                                                      
# scatseason     0.051    0.998     0    0.999    5.792    0.119     0.16    0.998          0    0.999      0.16    0.998     0.128    0.998  25.395    0.001        0.356    0.998
# logged         0.512    0.966 2.582    0.639    0.743    0.966    5.451    0.200      0.524    0.966     0.105    0.987     0.732    0.966   0.048    0.987        0.007    0.987
# Arguments:
#   Test statistics calculated assuming uncorrelated response (for faster computation) 
# P-value calculated using 999 iterations via PIT-trap resampling.