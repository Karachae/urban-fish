# R version 3.6.3 (server)
# Bayesian generalized linear mixed effects models

library(brms)
library(performance)
library(bayestestR)
library(tidyr)

### Read in data ###
fishydata<- read.csv("ReProj_FinalUrbanizationdata.csv", header = TRUE)

#Heterozygosity aka Genetic Diversity
fishy.he <- drop_na(fishydata, gene_diversity)

#FST
fishy.fst <- drop_na(fishydata, global_fst)

#Allelic Richness
fishy.alri <- drop_na(fishydata, allelic_richness)

#Ne
fishy.ne <- drop_na(fishydata, Ne)
## Log transform Ne ##
fishy.ne$Ne<- log(fishy.ne$Ne)

### BRM Function ###
brmfun <- function(model_spec, data_group)
{brm(model_spec, cores = 4, chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta = 0.999, max_treedepth = 15),
     data = data_group)}

### Models Without Priors ###

#### Heterozygosity ####

# Human Population Density Models
m_het_humanpopdens200 <- bf(scale(gene_diversity) ~ scale(hupopdens200km) + 
                                         (scale(hupopdens200km)|species))
m_het_humanpopdens100 <- bf(scale(gene_diversity) ~ scale(hupopdens100km) + 
                                         (scale(hupopdens100km)|species))
m_het_humanpopdens50 <- bf(scale(gene_diversity) ~ scale(hupopdens50km) + 
                                         (scale(hupopdens50km)|species))
m_het_humanpopdens25 <- bf(scale(gene_diversity) ~ scale(hupopdens25km) + 
                                         (scale(hupopdens25km)|species))

# Cumulative Human Impacts Models
m_het_impacts200 <- bf(scale(gene_diversity) ~ scale(impacts200km) + 
                                    (scale(impacts200km)|species))
m_het_impacts100 <- bf(scale(gene_diversity) ~ scale(impacts100km) + 
                                    (scale(impacts100km)|species))
m_het_impacts50 <- bf(scale(gene_diversity) ~ scale(impacts50km) + 
                                    (scale(impacts50km)|species))
m_het_impacts25 <- bf(scale(gene_diversity) ~ scale(impacts25km) + 
                                    (scale(impacts25km)|species))

# Run Human Population Density Models
MH_humanpopdens200 <-brmfun(m_het_humanpopdens200, fishy.he) 
print(summary(MH_humanpopdens200)) 
plot(MH_humanpopdens200)
p_direction(MH_humanpopdens200) # Probability of Direction
bayes_R2(MH_humanpopdens200) # Bayesian R2
waic(MH_humanpopdens200) # waic
pp_check(MH_humanpopdens200 , ndraws = 500) 
saveRDS(MH_humanpopdens200, "brms_models_nopriors/MH_humanpopdens200.rds") 
MH_humanpopdens200 <- readRDS("brms_models_nopriors/MH_humanpopdens200.rds")

MH_humanpopdens100 <-brmfun(m_het_humanpopdens100, fishy.he) 
print(summary(MH_humanpopdens100)) 
plot(MH_humanpopdens100)
p_direction(MH_humanpopdens100) # Probability of Direction
bayes_R2(MH_humanpopdens100) # Bayesian R2
waic(MH_humanpopdens100) # waic
pp_check(MH_humanpopdens100 , ndraws = 500) 
saveRDS(MH_humanpopdens100, "brms_models_nopriors/MH_humanpopdens100.rds") 
MH_humanpopdens100 <- readRDS("brms_models_nopriors/MH_humanpopdens100.rds")

MH_humanpopdens50 <-brmfun(m_het_humanpopdens50, fishy.he) 
print(summary(MH_humanpopdens50)) 
plot(MH_humanpopdens50)
p_direction(MH_humanpopdens50) # Probability of Direction
bayes_R2(MH_humanpopdens50) # Bayesian R2 
waic(MH_humanpopdens50) # waic
pp_check(MH_humanpopdens50 , ndraws = 500) 
saveRDS(MH_humanpopdens50, "brms_models_nopriors/MH_humanpopdens50.rds") 
MH_humanpopdens50 <- readRDS("brms_models_nopriors/MH_humanpopdens50.rds")

MH_humanpopdens25 <-brmfun(m_het_humanpopdens25, fishy.he) 
print(summary(MH_humanpopdens25)) 
plot(MH_humanpopdens25)
p_direction(MH_humanpopdens25) # Probability of Direction 
pp_check(MH_humanpopdens25 , ndraws = 500) 
bayes_R2(MH_humanpopdens25) # Bayesian R2  
waic(MH_humanpopdens25) # waic
saveRDS(MH_humanpopdens25, "brms_models_nopriors/MH_humanpopdens25.rds") 
MH_humanpopdens25 <- readRDS("brms_models_nopriors/MH_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
MH_huimpacts200 <-brmfun(m_het_impacts200, fishy.he) 
print(summary(MH_huimpacts200))
plot(MH_huimpacts200)
p_direction(MH_huimpacts200) # Probability of Direction
bayes_R2(MH_huimpacts200) # Bayesian R2
waic(MH_huimpacts200) # waic
pp_check(MH_huimpacts200, ndraws = 500) 
saveRDS(MH_huimpacts200, "brms_models_nopriors/MH_huimpacts200.rds") 
MH_huimpacts200 <- readRDS("brms_models_nopriors/MH_huimpacts200.rds")

MH_huimpacts100 <-brmfun(m_het_impacts100, fishy.he) 
print(summary(MH_huimpacts100))
plot(MH_huimpacts100)
p_direction(MH_huimpacts100) # Probability of Direction
bayes_R2(MH_huimpacts100) # Bayesian R2
waic(MH_huimpacts100) # waic
pp_check(MH_huimpacts100, ndraws = 500) 
saveRDS(MH_huimpacts100, "brms_models_nopriors/MH_huimpacts100.rds") 
MH_huimpacts100 <- readRDS("brms_models_nopriors/MH_huimpacts100.rds")

MH_huimpacts50 <-brmfun(m_het_impacts50, fishy.he) 
print(summary(MH_huimpacts50))
plot(MH_huimpacts50)
p_direction(MH_huimpacts50) # Probability of Direction
bayes_R2(MH_huimpacts50) # Bayesian R2
waic(MH_huimpacts50) # waic
pp_check(MH_huimpacts50, ndraws = 500) 
saveRDS(MH_huimpacts50, "brms_models_nopriors/MH_huimpacts50.rds") 
MH_huimpacts50 <- readRDS("brms_models_nopriors/MH_huimpacts50.rds")

MH_huimpacts25 <-brmfun(m_het_impacts25, fishy.he) 
print(summary(MH_huimpacts25))
plot(MH_huimpacts25) 
p_direction(MH_huimpacts25) # Probability of Direction
bayes_R2(MH_huimpacts25) # Bayesian R2
waic(MH_huimpacts25) # waic
pp_check(MH_huimpacts25, ndraws = 500) 
saveRDS(MH_huimpacts25, "brms_models_nopriors/MH_huimpacts25.rds") 
MH_huimpacts25 <- readRDS("brms_models_nopriors/MH_huimpacts25.rds")

#### Allelic Richness ####

# Human Population Density Models
m_alri_humanpopdens200 <- bf(scale(allelic_richness) ~ scale(hupopdens200km) 
                                          + (scale(hupopdens200km)|species))
m_alri_humanpopdens100 <- bf(scale(allelic_richness) ~ scale(hupopdens100km)
                                          + (scale(hupopdens100km)|species))
m_alri_humanpopdens50 <- bf(scale(allelic_richness) ~ scale(hupopdens50km) 
                                          + (scale(hupopdens50km)|species))
m_alri_humanpopdens25 <- bf(scale(allelic_richness) ~ scale(hupopdens25km)
                                          + (scale(hupopdens25km)|species))

# Cumulative Human Impacts Models
m_alri_impacts200 <- bf(scale(allelic_richness) ~ scale(impacts200km) 
                                    + (scale(impacts200km)|species))
m_alri_impacts100 <- bf(scale(allelic_richness) ~ scale(impacts100km)
                                    + (scale(impacts100km)|species))
m_alri_impacts50 <- bf(scale(allelic_richness) ~ scale(impacts50km)  
                                    + (scale(impacts50km)|species))
m_alri_impacts25 <- bf(scale(allelic_richness) ~ scale(impacts25km)
                                    + (scale(impacts25km)|species))

# Run Human Population Density Models
MAR_humanpopdens200 <-brmfun(m_alri_humanpopdens200, fishy.alri) 
print(summary(MAR_humanpopdens200)) 
plot(MAR_humanpopdens200)
p_direction(MAR_humanpopdens200) # Probability of Direction
bayes_R2(MAR_humanpopdens200) # Bayesian R2
waic(MAR_humanpopdens200) # waic
pp_check(MAR_humanpopdens200, ndraws = 500)
saveRDS(MAR_humanpopdens200, "brms_models_nopriors/MAR_humanpopdens200.rds") 
MAR_humanpopdens200<- readRDS("brms_models_nopriors/MAR_humanpopdens200.rds")

MAR_humanpopdens100 <-brmfun(m_alri_humanpopdens100, fishy.alri) 
print(summary(MAR_humanpopdens100)) 
plot(MAR_humanpopdens100)
p_direction(MAR_humanpopdens100) # Probability of Direction
bayes_R2(MAR_humanpopdens100) # Bayesian R2
waic(MAR_humanpopdens100) # waic
pp_check(MAR_humanpopdens100, ndraws = 500)
saveRDS(MAR_humanpopdens100, "brms_models_nopriors/MAR_humanpopdens100.rds") 
MAR_humanpopdens100<- readRDS("brms_models_nopriors/MAR_humanpopdens100.rds")

MAR_humanpopdens50 <-brmfun(m_alri_humanpopdens50, fishy.alri) 
print(summary(MAR_humanpopdens50)) 
plot(MAR_humanpopdens50)
p_direction(MAR_humanpopdens50) # Probability of Direction
bayes_R2(MAR_humanpopdens50) # Bayesian R2
waic(MAR_humanpopdens50) # waic
pp_check(MAR_humanpopdens50, ndraws = 500)
saveRDS(MAR_humanpopdens50, "brms_models_nopriors/MAR_humanpopdens50.rds")
MAR_humanpopdens50<- readRDS("brms_models_nopriors/MAR_humanpopdens50.rds")

MAR_humanpopdens25 <-brmfun(m_alri_humanpopdens25, fishy.alri) 
print(summary(MAR_humanpopdens25)) 
plot(MAR_humanpopdens25)
p_direction(MAR_humanpopdens25) # Probability of Direction
bayes_R2(MAR_humanpopdens25) # Bayesian R2
waic(MAR_humanpopdens25) # waic
pp_check(MAR_humanpopdens25, ndraws = 500)
saveRDS(MAR_humanpopdens25, "brms_models_nopriors/MAR_humanpopdens25.rds") 
MAR_humanpopdens25<- readRDS("brms_models_nopriors/MAR_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
MAR_huimpacts200 <-brmfun(m_alri_impacts200, fishy.alri) 
print(summary(MAR_huimpacts200))
plot(MAR_huimpacts200)
p_direction(MAR_huimpacts200) # Probability of Direction
bayes_R2(MAR_huimpacts200) # Bayesian R2
waic(MAR_huimpacts200) # waic
pp_check(MAR_huimpacts200, ndraws = 500)
saveRDS(MAR_huimpacts200, "brms_models_nopriors/MAR_huimpacts200.rds") 
MAR_huimpacts200<- readRDS("brms_models_nopriors/MAR_huimpacts200.rds")

MAR_huimpacts100 <-brmfun(m_alri_impacts100, fishy.alri) 
print(summary(MAR_huimpacts100))
plot(MAR_huimpacts100)
p_direction(MAR_huimpacts100) # Probability of Direction
bayes_R2(MAR_huimpacts100) # Bayesian R2
waic(MAR_huimpacts100) # waic
pp_check(MAR_huimpacts100, ndraws = 500)
saveRDS(MAR_huimpacts100, "brms_models_nopriors/MAR_huimpacts100.rds") 
MAR_huimpacts100<- readRDS("brms_models_nopriors/MAR_huimpacts100.rds")

MAR_huimpacts50 <-brmfun(m_alri_impacts50, fishy.alri)  
print(summary(MAR_huimpacts50))
plot(MAR_huimpacts50)
p_direction(MAR_huimpacts50) # Probability of Direction
bayes_R2(MAR_huimpacts50) # Bayesian R2
waic(MAR_huimpacts50) # waic
pp_check(MAR_huimpacts50, ndraws = 500)
saveRDS(MAR_huimpacts50, "brms_models_nopriors/MAR_huimpacts50.rds") 
MAR_huimpacts50<- readRDS("brms_models_nopriors/MAR_huimpacts50.rds")

MAR_huimpacts25 <-brmfun(m_alri_impacts25, fishy.alri) 
print(summary(MAR_huimpacts25))
plot(MAR_huimpacts25)
p_direction(MAR_huimpacts25) # Probability of Direction
bayes_R2(MAR_huimpacts25) # Bayesian R2
waic(MAR_huimpacts25) # waic
pp_check(MAR_huimpacts25, ndraws = 500)
saveRDS(MAR_huimpacts25, "brms_models_nopriors/MAR_huimpacts25.rds") 
MAR_huimpacts25<- readRDS("brms_models_nopriors/MAR_huimpacts25.rds")

#### FST ####

# Human Population Density Models
m_fst_humanpopdens200 <- bf(scale(global_fst) ~ scale(hupopdens200km) + (scale(hupopdens200km)|species))
m_fst_humanpopdens100 <- bf(scale(global_fst) ~ scale(hupopdens100km) + (scale(hupopdens100km)|species))
m_fst_humanpopdens50 <- bf(scale(global_fst) ~ scale(hupopdens50km) + (scale(hupopdens50km)|species))
m_fst_humanpopdens25 <- bf(scale(global_fst) ~ scale(hupopdens25km) + (scale(hupopdens25km)|species))

# Cumulative Human Impacts Models
m_fst_impacts200 <- bf(scale(global_fst) ~ scale(impacts200km) + (scale(impacts200km)|species))
m_fst_impacts100 <- bf(scale(global_fst) ~ scale(impacts100km) + (scale(impacts100km)|species))
m_fst_impacts50 <- bf(scale(global_fst) ~ scale(impacts50km) + (scale(impacts50km)|species))
m_fst_impacts25 <- bf(scale(global_fst) ~ scale(impacts25km) + (scale(impacts25km)|species))

# Run Human Population Density Models
MFST_humanpopdens200 <-brmfun(m_fst_humanpopdens200, fishy.fst)  
print(summary(MFST_humanpopdens200)) 
plot(MFST_humanpopdens200)
p_direction(MFST_humanpopdens200) # Probability of Direction
bayes_R2(MFST_humanpopdens200) # Bayesian R2
waic(MFST_humanpopdens200) # waic
pp_check(MFST_humanpopdens200, ndraws = 500)
saveRDS(MFST_humanpopdens200, "brms_models_nopriors/MFST_humanpopdens200.rds") 
MFST_humanpopdens200<- readRDS("brms_models_nopriors/MFST_humanpopdens200.rds")

MFST_humanpopdens100 <-brmfun(m_fst_humanpopdens100, fishy.fst) 
print(summary(MFST_humanpopdens100)) 
plot(MFST_humanpopdens100)
p_direction(MFST_humanpopdens100) # Probability of Direction
bayes_R2(MFST_humanpopdens100) # Bayesian R2
waic(MFST_humanpopdens100) # waic
pp_check(MFST_humanpopdens100, ndraws = 500)
saveRDS(MFST_humanpopdens100, "brms_models_nopriors/MFST_humanpopdens100.rds") 
MFST_humanpopdens100<- readRDS("brms_models_nopriors/MFST_humanpopdens100.rds")

MFST_humanpopdens50 <-brmfun(m_fst_humanpopdens50, fishy.fst) 
print(summary(MFST_humanpopdens50)) 
plot(MFST_humanpopdens50)
p_direction(MFST_humanpopdens50) # Probability of Direction
bayes_R2(MFST_humanpopdens50) # Bayesian R2
waic(MFST_humanpopdens50) # waic
pp_check(MFST_humanpopdens50, ndraws = 500)
saveRDS(MFST_humanpopdens50, "brms_models_nopriors/MFST_humanpopdens50.rds") 
MFST_humanpopdens50<- readRDS("brms_models_nopriors/MFST_humanpopdens50.rds")

MFST_humanpopdens25 <-brmfun(m_fst_humanpopdens25, fishy.fst) 
print(summary(MFST_humanpopdens25)) 
plot(MFST_humanpopdens25)
p_direction(MFST_humanpopdens25) # Probability of Direction
bayes_R2(MFST_humanpopdens25) # Bayesian R2
waic(MFST_humanpopdens25) # waic
pp_check(MFST_humanpopdens25, ndraws = 500)
saveRDS(MFST_humanpopdens25, "brms_models_nopriors/MFST_humanpopdens25.rds") 
MFST_humanpopdens25<- readRDS("brms_models_nopriors/MFST_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
MFST_huimpacts200 <-brmfun(m_fst_impacts200, fishy.fst) 
print(summary(MFST_huimpacts200))
plot(MFST_huimpacts200)
p_direction(MFST_huimpacts200) # Probability of Direction
bayes_R2(MFST_huimpacts200) # Bayesian R2
waic(MFST_huimpacts200) # waic
pp_check(MFST_huimpacts200, ndraws = 500)
saveRDS(MFST_huimpacts200, "brms_models_nopriors/MFST_huimpacts200.rds") 
MFST_huimpacts200<- readRDS("brms_models_nopriors/MFST_huimpacts200.rds")

MFST_huimpacts100 <-brmfun(m_fst_impacts100, fishy.fst) 
print(summary(MFST_huimpacts100))
plot(MFST_huimpacts100)
p_direction(MFST_huimpacts100) # Probability of Direction
bayes_R2(MFST_huimpacts100) # Bayesian R2  
waic(MFST_huimpacts100) # waic
pp_check(MFST_huimpacts100, ndraws = 500)
saveRDS(MFST_huimpacts100, "brms_models_nopriors/MFST_huimpacts100.rds")
MFST_huimpacts100<- readRDS("brms_models_nopriors/MFST_huimpacts100.rds")

MFST_huimpacts50 <-brmfun(m_fst_impacts50, fishy.fst) 
print(summary(MFST_huimpacts50))
plot(MFST_huimpacts50)
p_direction(MFST_huimpacts50) # Probability of Direction
bayes_R2(MFST_huimpacts50) # Bayesian R2 
waic(MFST_huimpacts50) # waic
pp_check(MFST_huimpacts50, ndraws = 500)
saveRDS(MFST_huimpacts50, "brms_models_nopriors/MFST_huimpacts50.rds") 
MFST_huimpacts50<- readRDS("brms_models_nopriors/MFST_huimpacts50.rds")

MFST_huimpacts25 <-brmfun(m_fst_impacts25, fishy.fst) 
print(summary(MFST_huimpacts25))
plot(MFST_huimpacts25)
p_direction(MFST_huimpacts25) # Probability of Direction
bayes_R2(MFST_huimpacts25) # Bayesian R2
waic(MFST_huimpacts25) # waic
pp_check(MFST_huimpacts25, ndraws = 500)
saveRDS(MFST_huimpacts25, "brms_models_nopriors/MFST_huimpacts25.rds") 
MFST_huimpacts25<- readRDS("brms_models_nopriors/MFST_huimpacts25.rds")

#### Ne ####

# Human Population Density Models
m_ne_humanpopdens200 <- bf(scale(Ne) ~ scale(hupopdens200km) + (scale(hupopdens200km)|species))
m_ne_humanpopdens100 <- bf(scale(Ne) ~ scale(hupopdens100km) + (scale(hupopdens100km)|species))
m_ne_humanpopdens50 <- bf(scale(Ne) ~ scale(hupopdens50km) + (scale(hupopdens50km)|species))
m_ne_humanpopdens25 <- bf(scale(Ne) ~ scale(hupopdens25km) + (scale(hupopdens25km)|species))

# Cumulative Human Impacts Models
m_ne_impacts200 <- bf(scale(Ne) ~ scale(impacts200km) + (scale(impacts200km)|species))
m_ne_impacts100 <- bf(scale(Ne) ~ scale(impacts100km) + (scale(impacts100km)|species))
m_ne_impacts50 <- bf(scale(Ne) ~ scale(impacts50km) + (scale(impacts50km)|species))
m_ne_impacts25 <- bf(scale(Ne) ~ scale(impacts25km) + (scale(impacts25km)|species))

# Run Human Population density Models
MNE_humanpopdens200 <-brmfun(m_ne_humanpopdens200, fishy.ne) 
print(summary(MNE_humanpopdens200)) 
plot(MNE_humanpopdens200)
p_direction(MNE_humanpopdens200) # Probability of Direction
bayes_R2(MNE_humanpopdens200) # Bayesian R2
waic(MNE_humanpopdens200) # waic
pp_check(MNE_humanpopdens200, ndraws = 500)
saveRDS(MNE_humanpopdens200, "brms_models_nopriors/MNE_humanpopdens200.rds") 
MNE_humanpopdens200 <- readRDS("brms_models_nopriors/MNE_humanpopdens200.rds")

MNE_humanpopdens100 <-brmfun(m_ne_humanpopdens100, fishy.ne) 
print(summary(MNE_humanpopdens100)) 
plot(MNE_humanpopdens100)
p_direction(MNE_humanpopdens100) # Probability of Direction
bayes_R2(MNE_humanpopdens100) # Bayesian R2
waic(MNE_humanpopdens100) # waic
pp_check(MNE_humanpopdens100, ndraws = 500)
saveRDS(MNE_humanpopdens100, "brms_models_nopriors/MNE_humanpopdens100.rds") 
MNE_humanpopdens100 <- readRDS("brms_models_nopriors/MNE_humanpopdens100.rds")

MNE_humanpopdens50 <-brmfun(m_ne_humanpopdens50, fishy.ne) 
print(summary(MNE_humanpopdens50)) 
plot(MNE_humanpopdens50)
p_direction(MNE_humanpopdens50) # Probability of Direction
bayes_R2(MNE_humanpopdens50) # Bayesian R2
waic(MNE_humanpopdens50) # waic
pp_check(MNE_humanpopdens50, ndraws = 500)
saveRDS(MNE_humanpopdens50, "brms_models_nopriors/MNE_humanpopdens50.rds") 
MNE_humanpopdens50 <- readRDS("brms_models_nopriors/MNE_humanpopdens50.rds")

MNE_humanpopdens25 <-brmfun(m_ne_humanpopdens25, fishy.ne) 
print(summary(MNE_humanpopdens25)) 
plot(MNE_humanpopdens25)
p_direction(MNE_humanpopdens25) # Probability of Direction
bayes_R2(MNE_humanpopdens25) # Bayesian R2
waic(MNE_humanpopdens25) # waic
pp_check(MNE_humanpopdens25, ndraws = 500)
saveRDS(MNE_humanpopdens25, "brms_models_nopriors/MNE_humanpopdens25.rds") 
MNE_humanpopdens25 <- readRDS("brms_models_nopriors/MNE_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
MNE_huimpacts200 <-brmfun(m_ne_impacts200, fishy.ne)  
print(summary(MNE_huimpacts200))
plot(MNE_huimpacts200)
p_direction(MNE_huimpacts200) # Probability of Direction
bayes_R2(MNE_huimpacts200) # Bayesian R2
waic(MNE_huimpacts200) # waic
pp_check(MNE_huimpacts200, ndraws = 500)
saveRDS(MNE_huimpacts200, "brms_models_nopriors/MNE_huimpacts200.rds") 
MNE_huimpacts200 <- readRDS("brms_models_nopriors/MNE_huimpacts200.rds")

MNE_huimpacts100 <-brmfun(m_ne_impacts100, fishy.ne) 
print(summary(MNE_huimpacts100))
plot(MNE_huimpacts100)
p_direction(MNE_huimpacts100) # Probability of Direction
bayes_R2(MNE_huimpacts100) # Bayesian R2
waic(MNE_huimpacts100) # waic
pp_check(MNE_huimpacts100, ndraws = 500)
saveRDS(MNE_huimpacts100, "brms_models_nopriors/MNE_huimpacts100.rds") 
MNE_huimpacts100 <- readRDS("brms_models_nopriors/MNE_huimpacts100.rds")

MNE_huimpacts50 <-brmfun(m_ne_impacts50, fishy.ne)  
print(summary(MNE_huimpacts50))
plot(MNE_huimpacts50)
p_direction(MNE_huimpacts50) # Probability of Direction
bayes_R2(MNE_huimpacts50) # Bayesian R2
waic(MNE_huimpacts50) # waic 
pp_check(MNE_huimpacts50, ndraws = 500)
saveRDS(MNE_huimpacts50, "brms_models_nopriors/MNE_huimpacts50.rds") 
MNE_huimpacts50 <- readRDS("brms_models_nopriors/MNE_huimpacts50.rds")

MNE_huimpacts25 <-brmfun(m_ne_impacts25, fishy.ne) 
print(summary(MNE_huimpacts25))
plot(MNE_huimpacts25)
p_direction(MNE_huimpacts25) # Probability of Direction
bayes_R2(MNE_huimpacts25) # Bayesian R2
waic(MNE_huimpacts25) # waic
pp_check(MNE_huimpacts25, ndraws = 500)
saveRDS(MNE_huimpacts25, "brms_models_nopriors/MNE_huimpacts25.rds") 
MNE_huimpacts25 <- readRDS("brms_models_nopriors/MNE_huimpacts25.rds")

### Models with Priors ###

#### Heterozygosity with priors ####

# Run Human Population density Models
p.MH_humanpopdens200 <-brm(m_het_humanpopdens200, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens200km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data =fishy.he) 
print(summary(p.MH_humanpopdens200)) 
plot(p.MH_humanpopdens200)
p_direction(p.MH_humanpopdens200) # Probability of Direction
pp_check(p.MH_humanpopdens200, ndraws = 500)
bayes_R2(p.MH_humanpopdens200) # Bayesian R2
waic(p.MH_humanpopdens200) # waic
saveRDS(p.MH_humanpopdens200, "brms_models/p.MH_humanpopdens200.rds")
p.MH_humanpopdens200 <- readRDS("brms_models/p.MH_humanpopdens200.rds")

p.MH_humanpopdens100 <-brm(m_het_humanpopdens100, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens100km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data =fishy.he) 
print(summary(p.MH_humanpopdens100)) 
plot(p.MH_humanpopdens100)
p_direction(p.MH_humanpopdens100) # Probability of Direction
pp_check(p.MH_humanpopdens100, ndraws = 500)
bayes_R2(p.MH_humanpopdens100) # Bayesian R2
waic(p.MH_humanpopdens100) # waic
saveRDS(p.MH_humanpopdens100, "brms_models/p.MH_humanpopdens100.rds")
p.MH_humanpopdens100 <- readRDS("brms_models/p.MH_humanpopdens100.rds")

p.MH_humanpopdens50 <-brm(m_het_humanpopdens50, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                          prior = prior(normal(0,1), class = b, coef = scalehupopdens50km),
                          control = list(adapt_delta = 0.95, max_treedepth = 10), 
                          data =fishy.he)  
print(summary(p.MH_humanpopdens50)) 
plot(p.MH_humanpopdens50)
p_direction(p.MH_humanpopdens50) # Probability of Direction
pp_check(p.MH_humanpopdens50, ndraws = 500)
bayes_R2(p.MH_humanpopdens50) # Bayesian R2
waic(p.MH_humanpopdens50) # waic
saveRDS(p.MH_humanpopdens50, "brms_models/p.MH_humanpopdens50.rds")
p.MH_humanpopdens50 <- readRDS("brms_models/p.MH_humanpopdens50.rds")

p.MH_humanpopdens25 <-brm(m_het_humanpopdens25, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                          prior = prior(normal(0,1), class = b, coef = scalehupopdens25km),
                          control = list(adapt_delta = 0.95, max_treedepth = 10), 
                          data =fishy.he) 
print(summary(p.MH_humanpopdens25)) 
plot(p.MH_humanpopdens25)
p_direction(p.MH_humanpopdens25) # Probability of Direction
pp_check(p.MH_humanpopdens25, ndraws = 500)
bayes_R2(p.MH_humanpopdens25) # Bayesian R2
waic(p.MH_humanpopdens25) # waic
saveRDS(p.MH_humanpopdens25, "brms_models/p.MH_humanpopdens25.rds")
p.MH_humanpopdens25 <- readRDS("brms_models/p.MH_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
p.MH_huimpacts200 <-brm(m_het_impacts200, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                        prior = prior(normal(0,1), class = b, coef = scaleimpacts200km),
                        control = list(adapt_delta = 0.95, max_treedepth = 10), 
                        data =fishy.he) 
print(summary(p.MH_huimpacts200))
plot(p.MH_huimpacts200)
p_direction(p.MH_huimpacts200) # Probability of Direction
pp_check(p.MH_huimpacts200, ndraws = 500)
bayes_R2(p.MH_huimpacts200) # Bayesian R2
waic(p.MH_huimpacts200) # waic
saveRDS(p.MH_huimpacts200, "brms_models/p.MH_huimpacts200.rds")
p.MH_huimpacts200 <- readRDS("brms_models/p.MH_huimpacts200.rds")

p.MH_huimpacts100 <-brm(m_het_impacts100, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                        prior = prior(normal(0,1), class = b, coef = scaleimpacts100km),
                        control = list(adapt_delta = 0.99, max_treedepth = 10), 
                        data =fishy.he) 
print(summary(p.MH_huimpacts100))
plot(p.MH_huimpacts100)
p_direction(p.MH_huimpacts100) # Probability of Direction
pp_check(p.MH_huimpacts100, ndraws = 500)
bayes_R2(p.MH_huimpacts100) # Bayesian R2
waic(p.MH_huimpacts100) # waic
saveRDS(p.MH_huimpacts100, "brms_models/p.MH_huimpacts100.rds")
p.MH_huimpacts100 <- readRDS("brms_models/p.MH_huimpacts100.rds")

p.MH_huimpacts50 <-brm(m_het_impacts50, cores = 4, chains = 4, iter = 6000, warmup = 3000,
                       prior = prior(normal(0,1), class = b, coef = scaleimpacts50km),
                       control = list(adapt_delta = 0.95, max_treedepth = 10), 
                       data =fishy.he) 
print(summary(p.MH_huimpacts50))
plot(p.MH_huimpacts50)
p_direction(p.MH_huimpacts50) # Probability of Direction
pp_check(p.MH_huimpacts50, ndraws = 500)
bayes_R2(p.MH_huimpacts50) # Bayesian R2
waic(p.MH_huimpacts50) # waic
saveRDS(p.MH_huimpacts50, "brms_models/p.MH_huimpacts50.rds")
p.MH_huimpacts50 <- readRDS("brms_models/p.MH_huimpacts50.rds")

p.MH_huimpacts25 <-brm(m_het_impacts25, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                       prior = prior(normal(0,1), class = b, coef = scaleimpacts25km),
                       control = list(adapt_delta = 0.95, max_treedepth = 10), 
                       data =fishy.he) 
print(summary(p.MH_huimpacts25))
plot(p.MH_huimpacts25)
p_direction(p.MH_huimpacts25) # Probability of Direction
pp_check(p.MH_huimpacts25, ndraws = 500)
bayes_R2(p.MH_huimpacts25) # Bayesian R2
waic(p.MH_huimpacts25) # waic
saveRDS(p.MH_huimpacts25, "brms_models/p.MH_huimpacts25.rds")
p.MH_huimpacts25 <- readRDS("brms_models/p.MH_huimpacts25.rds")

#### Allelic Richness with priors ####

# Run Human Population Density Models
p.MAR_humanpopdens200 <-brm(m_alri_humanpopdens200, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens200km),
                            control = list(adapt_delta = 0.95, max_treedepth = 10), 
                            data =fishy.alri) 
print(summary(p.MAR_humanpopdens200)) 
plot(p.MAR_humanpopdens200)
p_direction(p.MAR_humanpopdens200) # Probability of Direction
pp_check(p.MAR_humanpopdens200, ndraws = 500)
bayes_R2(p.MAR_humanpopdens200) # Bayesian R2
waic(p.MAR_humanpopdens200) # waic
saveRDS(p.MAR_humanpopdens200, "brms_models/p.MAR_humanpopdens200.rds")
p.MAR_humanpopdens200<- readRDS("brms_models/p.MAR_humanpopdens200.rds")

p.MAR_humanpopdens100 <-brm(m_alri_humanpopdens100, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens100km),
                            control = list(adapt_delta = 0.95, max_treedepth = 10), 
                            data =fishy.alri) 
print(summary(p.MAR_humanpopdens100)) 
plot(p.MAR_humanpopdens100)
p_direction(p.MAR_humanpopdens100) # Probability of Direction
pp_check(p.MAR_humanpopdens100, ndraws = 500)
bayes_R2(p.MAR_humanpopdens100) # Bayesian R2
waic(p.MAR_humanpopdens100) # waic
saveRDS(p.MAR_humanpopdens100, "brms_models/p.MAR_humanpopdens100.rds")
p.MAR_humanpopdens100<- readRDS("brms_models/p.MAR_humanpopdens100.rds")

p.MAR_humanpopdens50 <-brm(m_alri_humanpopdens50, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens50km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data =fishy.alri) 
print(summary(p.MAR_humanpopdens50)) 
plot(p.MAR_humanpopdens50)
p_direction(p.MAR_humanpopdens50) # Probability of Direction
pp_check(p.MAR_humanpopdens50, ndraws = 500)
bayes_R2(p.MAR_humanpopdens50) # Bayesian R2
waic(p.MAR_humanpopdens50) # waic
saveRDS(p.MAR_humanpopdens50, "brms_models/p.MAR_humanpopdens50.rds")
p.MAR_humanpopdens50<- readRDS("brms_models/p.MAR_humanpopdens50.rds")

p.MAR_humanpopdens25 <-brm(m_alri_humanpopdens25, cores = 4, chains = 4, iter = 4000, warmup = 2000, 
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens25km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data =fishy.alri) 
print(summary(p.MAR_humanpopdens25)) 
plot(p.MAR_humanpopdens25)
p_direction(p.MAR_humanpopdens25) # Probability of Direction
pp_check(p.MAR_humanpopdens25, ndraws = 500)
bayes_R2(p.MAR_humanpopdens25) # Bayesian R2
waic(p.MAR_humanpopdens25) # waic
saveRDS(p.MAR_humanpopdens25, "brms_models/p.MAR_humanpopdens25.rds")
p.MAR_humanpopdens25<- readRDS("brms_models/p.MAR_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
p.MAR_huimpacts200 <-brm(m_alri_impacts200, cores = 4, chains = 4, iter = 4000, warmup = 2000, 
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts200km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data =fishy.alri) 
print(summary(p.MAR_huimpacts200))
plot(p.MAR_huimpacts200)
p_direction(p.MAR_huimpacts200) # Probability of Direction
pp_check(p.MAR_huimpacts200, ndraws = 500)
bayes_R2(p.MAR_huimpacts200) # Bayesian R2
waic(p.MAR_huimpacts200) # waic
saveRDS(p.MAR_huimpacts200, "brms_models/p.MAR_huimpacts200.rds")
p.MAR_huimpacts200<- readRDS("brms_models/p.MAR_huimpacts200.rds")

p.MAR_huimpacts100 <-brm(m_alri_impacts100, cores = 4, chains = 4, iter = 4000, warmup = 1000, 
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts100km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data =fishy.alri) 
print(summary(p.MAR_huimpacts100))
plot(p.MAR_huimpacts100)
p_direction(p.MAR_huimpacts100) # Probability of Direction
pp_check(p.MAR_huimpacts100, ndraws = 500)
bayes_R2(p.MAR_huimpacts100) # Bayesian R2
waic(p.MAR_huimpacts100) # waic
saveRDS(p.MAR_huimpacts100, "brms_models/p.MAR_huimpacts100.rds")
p.MAR_huimpacts100<- readRDS("brms_models/p.MAR_huimpacts100.rds")

p.MAR_huimpacts50 <- brm(m_alri_impacts50, cores = 4, chains = 4, iter = 4000, warmup = 2000, 
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts50km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data =fishy.alri) 
print(summary(p.MAR_huimpacts50))
plot(p.MAR_huimpacts50)
p_direction(p.MAR_huimpacts50) # Probability of Direction
pp_check(p.MAR_huimpacts50, ndraws = 500)
bayes_R2(p.MAR_huimpacts50) # Bayesian R2
waic(p.MAR_huimpacts50) # waic
saveRDS(p.MAR_huimpacts50, "brms_models/p.MAR_huimpacts50.rds")
p.MAR_huimpacts50<- readRDS("brms_models/p.MAR_huimpacts50.rds")

p.MAR_huimpacts25 <-brm(m_alri_impacts25, cores = 4, chains = 4, iter = 4000, warmup = 2000, 
                        prior = prior(normal(0,1), class = b, coef = scaleimpacts25km),
                        control = list(adapt_delta = 0.95, max_treedepth = 10), 
                        data =fishy.alri)
print(summary(p.MAR_huimpacts25))
plot(p.MAR_huimpacts25)
p_direction(p.MAR_huimpacts25) # Probability of Direction
pp_check(p.MAR_huimpacts25, ndraws = 500)
bayes_R2(p.MAR_huimpacts25) # Bayesian R2
waic(p.MAR_huimpacts25) # waic
saveRDS(p.MAR_huimpacts25, "brms_models/p.MAR_huimpacts25.rds")
p.MAR_huimpacts25<- readRDS("brms_models/p.MAR_huimpacts25.rds")

#### FST with priors ####

# Run Human Population Density Models
p.MFST_humanpopdens200 <-brm(m_fst_humanpopdens200, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                             prior = prior(normal(0,1), class = b, coef = scalehupopdens200km),
                             control = list(adapt_delta = 0.95, max_treedepth = 10), 
                             data =fishy.fst) 
print(summary(p.MFST_humanpopdens200)) 
plot(p.MFST_humanpopdens200)
p_direction(p.MFST_humanpopdens200) # Probability of Direction
pp_check(p.MFST_humanpopdens200, ndraws = 500)
bayes_R2(p.MFST_humanpopdens200) # Bayesian R2
waic(p.MFST_humanpopdens200) # waic
saveRDS(p.MFST_humanpopdens200, "brms_models/p.MFST_humanpopdens200.rds")
p.MFST_humanpopdens200<- readRDS("brms_models/p.MFST_humanpopdens200.rds")

p.MFST_humanpopdens100 <-brm(m_fst_humanpopdens100, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                             prior = prior(normal(0,1), class = b, coef = scalehupopdens100km),
                             control = list(adapt_delta = 0.95, max_treedepth = 10), 
                             data =fishy.fst) 
print(summary(p.MFST_humanpopdens100)) 
plot(p.MFST_humanpopdens100)
p_direction(p.MFST_humanpopdens100) # Probability of Direction
pp_check(p.MFST_humanpopdens100, ndraws = 500)
bayes_R2(p.MFST_humanpopdens100) # Bayesian R2
waic(p.MFST_humanpopdens100) # waic
saveRDS(p.MFST_humanpopdens100, "brms_models/p.MFST_humanpopdens100.rds")
p.MFST_humanpopdens100<- readRDS("brms_models/p.MFST_humanpopdens100.rds")

p.MFST_humanpopdens50 <-brm(m_fst_humanpopdens50, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens50km),
                            control = list(adapt_delta = 0.95, max_treedepth = 10), 
                            data =fishy.fst)  
print(summary(p.MFST_humanpopdens50)) 
plot(p.MFST_humanpopdens50)
p_direction(p.MFST_humanpopdens50) # Probability of Direction
pp_check(p.MFST_humanpopdens50, ndraws = 500)
bayes_R2(p.MFST_humanpopdens50) # Bayesian R2
waic(p.MFST_humanpopdens50) # waic
saveRDS(p.MFST_humanpopdens50, "brms_models/p.MFST_humanpopdens50.rds")
p.MFST_humanpopdens50<- readRDS("brms_models/p.MFST_humanpopdens50.rds")

p.MFST_humanpopdens25 <-brm(m_fst_humanpopdens25, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens25km),
                            control = list(adapt_delta = 0.95, max_treedepth = 10), 
                            data =fishy.fst) 
print(summary(p.MFST_humanpopdens25)) 
plot(p.MFST_humanpopdens25)
p_direction(p.MFST_humanpopdens25) # Probability of Direction
pp_check(p.MFST_humanpopdens25, ndraws = 500)
bayes_R2(p.MFST_humanpopdens25) # Bayesian R2
waic(p.MFST_humanpopdens25) # waic
saveRDS(p.MFST_humanpopdens25, "brms_models/p.MFST_humanpopdens25.rds")
p.MFST_humanpopdens25<- readRDS("brms_models/p.MFST_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
p.MFST_huimpacts200 <-brm(m_fst_impacts200, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                          prior = prior(normal(0,1), class = b, coef = scaleimpacts200km),
                          control = list(adapt_delta = 0.95, max_treedepth = 10), 
                          data =fishy.fst) 
print(summary(p.MFST_huimpacts200))
plot(p.MFST_huimpacts200)
p_direction(p.MFST_huimpacts200) # Probability of Direction
pp_check(p.MFST_huimpacts200, ndraws = 500)
bayes_R2(p.MFST_huimpacts200) # Bayesian R2
waic(p.MFST_huimpacts200) # waic
saveRDS(p.MFST_huimpacts200, "brms_models/p.MFST_huimpacts200.rds")
p.MFST_huimpacts200 <- readRDS("brms_models/p.MFST_huimpacts200.rds")

p.MFST_huimpacts100 <-brm(m_fst_impacts100, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                          prior = prior(normal(0,1), class = b, coef = scaleimpacts100km),
                          control = list(adapt_delta = 0.99, max_treedepth = 10), 
                          data =fishy.fst) 
print(summary(p.MFST_huimpacts100))
plot(p.MFST_huimpacts100)
p_direction(p.MFST_huimpacts100) # Probability of Direction
pp_check(p.MFST_huimpacts100, ndraws = 500)
bayes_R2(p.MFST_huimpacts100) # Bayesian R2
waic(p.MFST_huimpacts100) # waic
saveRDS(p.MFST_huimpacts100, "brms_models/p.MFST_huimpacts100.rds")
p.MFST_huimpacts100 <- readRDS("brms_models/p.MFST_huimpacts100.rds")

p.MFST_huimpacts50 <-brm(m_fst_impacts50, cores = 4, chains = 4, iter = 4000, warmup = 2000,
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts50km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data =fishy.fst) 
print(summary(p.MFST_huimpacts50))
plot(p.MFST_huimpacts50)
p_direction(p.MFST_huimpacts50) # Probability of Direction
pp_check(p.MFST_huimpacts50, ndraws = 500)
bayes_R2(p.MFST_huimpacts50) # Bayesian R2
waic(p.MFST_huimpacts50) # waic
saveRDS(p.MFST_huimpacts50, "brms_models/p.MFST_huimpacts50.rds")
p.MFST_huimpacts50 <- readRDS("brms_models/p.MFST_huimpacts50.rds")

p.MFST_huimpacts25 <-brm(m_fst_impacts25, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts25km),
                         control = list(adapt_delta = 0.99, max_treedepth = 10), 
                         data =fishy.fst) 
print(summary(p.MFST_huimpacts25))
plot(p.MFST_huimpacts25)
p_direction(p.MFST_huimpacts25) # Probability of Direction
pp_check(p.MFST_huimpacts25, ndraws = 500)
bayes_R2(p.MFST_huimpacts25) # Bayesian R2
waic(p.MFST_huimpacts25) # waic
saveRDS(p.MFST_huimpacts25, "brms_models/p.MFST_huimpacts25.rds")
p.MFST_huimpacts25 <- readRDS("brms_models/p.MFST_huimpacts25.rds")

### Ne with priors ###

# Run Human Population Density Models
p.MNE_humanpopdens200 <-brm(m_ne_humanpopdens200, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens200km),
                            control = list(adapt_delta = 0.99, max_treedepth = 10), 
                            data = fishy.ne)  
print(summary(p.MNE_humanpopdens200)) 
plot(p.MNE_humanpopdens200)
p_direction(p.MNE_humanpopdens200) # Probability of Direction
pp_check(p.MNE_humanpopdens200, ndraws = 500)
bayes_R2(p.MNE_humanpopdens200) # Bayesian R2
waic(p.MNE_humanpopdens200) # waic
saveRDS(p.MNE_humanpopdens200, "brms_models/p.MNE_humanpopdens200.rds")
p.MNE_humanpopdens200<- readRDS("brms_models/p.MNE_humanpopdens200.rds")

p.MNE_humanpopdens100 <-brm(m_ne_humanpopdens100, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                            prior = prior(normal(0,1), class = b, coef = scalehupopdens100km),
                            control = list(adapt_delta = 0.95, max_treedepth = 10), 
                            data = fishy.ne) 
print(summary(p.MNE_humanpopdens100)) 
plot(p.MNE_humanpopdens100)
p_direction(p.MNE_humanpopdens100) # Probability of Direction
pp_check(p.MNE_humanpopdens100, ndraws = 500)
bayes_R2(p.MNE_humanpopdens100) # Bayesian R2
waic(p.MNE_humanpopdens100) # waic
saveRDS(p.MNE_humanpopdens100, "brms_models/p.MNE_humanpopdens100.rds")
p.MNE_humanpopdens100<- readRDS("brms_models/p.MNE_humanpopdens100.rds")

p.MNE_humanpopdens50 <-brm(m_ne_humanpopdens50, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens50km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data = fishy.ne) 
print(summary(p.MNE_humanpopdens50)) 
plot(p.MNE_humanpopdens50)
p_direction(p.MNE_humanpopdens50) # Probability of Direction
pp_check(p.MNE_humanpopdens50, ndraws = 500)
bayes_R2(p.MNE_humanpopdens50) # Bayesian R2
waic(p.MNE_humanpopdens50) # waic
saveRDS(p.MNE_humanpopdens50, "brms_models/p.MNE_humanpopdens50.rds")
p.MNE_humanpopdens50<- readRDS("brms_models/p.MNE_humanpopdens50.rds")

p.MNE_humanpopdens25 <-brm(m_ne_humanpopdens25, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                           prior = prior(normal(0,1), class = b, coef = scalehupopdens25km),
                           control = list(adapt_delta = 0.95, max_treedepth = 10), 
                           data = fishy.ne) 
print(summary(p.MNE_humanpopdens25)) 
plot(p.MNE_humanpopdens25)
p_direction(p.MNE_humanpopdens25) # Probability of Direction
pp_check(p.MNE_humanpopdens25, ndraws = 500)
bayes_R2(p.MNE_humanpopdens25) # Bayesian R2
waic(p.MNE_humanpopdens25) # waic
saveRDS(p.MNE_humanpopdens25, "brms_models/p.MNE_humanpopdens25.rds")
p.MNE_humanpopdens25<- readRDS("brms_models/p.MNE_humanpopdens25.rds")

# Run Cumulative Human Impacts Models
p.MNE_huimpacts200 <-brm(m_ne_impacts200, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts200km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data = fishy.ne) 
print(summary(p.MNE_huimpacts200))
plot(p.MNE_huimpacts200)
p_direction(p.MNE_huimpacts200) # Probability of Direction
pp_check(p.MNE_huimpacts200, ndraws = 500)
bayes_R2(p.MNE_huimpacts200) # Bayesian R2
waic(p.MNE_huimpacts200) # waic
saveRDS(p.MNE_huimpacts200, "brms_models/p.MNE_huimpacts200.rds")
p.MNE_huimpacts200 <- readRDS("brms_models/p.MNE_huimpacts200.rds")

p.MNE_huimpacts100 <-brm(m_ne_impacts100, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                         prior = prior(normal(0,1), class = b, coef = scaleimpacts100km),
                         control = list(adapt_delta = 0.95, max_treedepth = 10), 
                         data = fishy.ne) 
print(summary(p.MNE_huimpacts100))
plot(p.MNE_huimpacts100)
p_direction(p.MNE_huimpacts100) # Probability of Direction
pp_check(p.MNE_huimpacts100, ndraws = 500)
bayes_R2(p.MNE_huimpacts100) # Bayesian R2
waic(p.MNE_huimpacts100) # waic
saveRDS(p.MNE_huimpacts100, "brms_models/p.MNE_huimpacts100.rds")
p.MNE_huimpacts100 <- readRDS("brms_models/p.MNE_huimpacts100.rds")

p.MNE_huimpacts50 <-brm(m_ne_impacts50, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                        prior = prior(normal(0,1), class = b, coef = scaleimpacts50km),
                        control = list(adapt_delta = 0.95, max_treedepth = 10), 
                        data = fishy.ne) # EFFECT
print(summary(p.MNE_huimpacts50))
plot(p.MNE_huimpacts50)
p_direction(p.MNE_huimpacts50) # Probability of Direction
pp_check(p.MNE_huimpacts50, ndraws = 500)
bayes_R2(p.MNE_huimpacts50) # Bayesian R2
waic(p.MNE_huimpacts50) # waic
saveRDS(p.MNE_huimpacts50, "brms_models/p.MNE_huimpacts50.rds")
p.MNE_huimpacts50 <- readRDS("brms_models/p.MNE_huimpacts50.rds")

p.MNE_huimpacts25 <-brm(m_ne_impacts25, cores = 4, chains = 4, iter = 2000, warmup = 1000,
                        prior = prior(normal(0,1), class = b, coef = scaleimpacts25km),
                        control = list(adapt_delta = 0.95, max_treedepth = 10), 
                        data = fishy.ne) 
print(summary(p.MNE_huimpacts25))
plot(p.MNE_huimpacts25)
p_direction(p.MNE_huimpacts25) # Probability of Direction
pp_check(p.MNE_huimpacts25, ndraws = 500)
bayes_R2(p.MNE_huimpacts25) # Bayesian R2
waic(p.MNE_huimpacts50) # waic
saveRDS(p.MNE_huimpacts25, "brms_models/p.MNE_huimpacts25.rds")
p.MNE_huimpacts25 <- readRDS("brms_models/p.MNE_huimpacts25.rds")