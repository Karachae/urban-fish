# R version 4.1.2
# Genetic data prep 

library(adegenet)
library(hierfstat)
library(raster)
library(rgdal)
library(readxl)

# List of raw STRUCTURE files
file_list <- list.files("/StructureData", full.names = FALSE)
species_list <- list()

setwd("/StructureData")

# Read in all files as genind objects
for (i in file_list){
  if(grepl(".stru", i)){
    junk <- read.table(i)
    s <- read.structure(i, n.ind = (nrow(junk)), n.loc = ((ncol(junk)-1)/2),
                        onerowperind = TRUE, col.lab = 0, col.pop = 1, row.marknames = 0, ask = FALSE)}
  else {
    k <- try(read.genepop(i, ncode = 2))
    s <- if(inherits(k, "try-error")){read.genepop(i, ncode = 3)}
    else {
      s <- (read.genepop(i, ncode = 2))}
  }
  
  species_list[[i]] <- s
  
}

# read in data

popmaster <- read_excel("GlobalMasterpopsheet.xlsx")

# List of names for all populations
allpoplist <- data.frame(unlist(lapply(species_list, popNames)))

## Gene diversity ##
# expected heterozygosity
gendiv <- data.frame(unlist(lapply(species_list, Hs)))
popgendiv <- cbind(allpoplist, gendiv)
names(popgendiv) <- c("pop","gene_diversity")

## Allelic richness ##
# Rarefied allelic richness averaged across loci per population
# minimum number of alleles = 10

arich <- data.frame(unlist(lapply(species_list, 
                                  function(x) colMeans(allelic.richness(x, min.n = 10, diploid = TRUE)$Ar, 
                                                       na.rm = TRUE))))
arich <- na.omit(arich)
arich <- cbind(allpoplist, arich)
names(arich) <- c("pop","allelic_richness")

## Population-specific FST ##
# requires >1 population for calculation

a <- unlist(lapply(species_list, function(x) length(levels(pop(x)))))
c <- which(a>1)
morethan1pop <- species_list[c]

# calculate FST
species_fst <- data.frame(unlist(lapply(morethan1pop, function (x) betas(x)$betaiovl)))
morethan1popids <- data.frame(unlist(lapply(morethan1pop, popNames)))
globfst <- cbind(morethan1popids, species_fst)
names(globfst) <- c("pop","global_fst")

## Number of individuals in population ##

t <- unname(unlist(lapply(species_list, function(x) summary(x)$n.by.pop)))
indivcount<- cbind(allpoplist,t)
names(indivcount) <- c("pop","num_individuals")

## Effective population size (Ne) ##

# estimates from NeEstimator v2
Ne <- read_excel("NeOnlyForMerge.xlsx")
Ne$Ne<-as.numeric(Ne$Ne) 

## Merge ##

popscoordsumm <- Reduce(function(x,y) merge(x,y, all=TRUE, incomparables = "NA"), 
                        list(popmaster, indivcount, popgendiv, arich, globfst, Ne))

# Remove populations with <5 individuals and populations with no coordinates
data.noindiv <- popscoordsumm[popscoordsumm$num_individuals>=5,]
data.noindiv.narm <- data.noindiv[!is.na(data.noindiv$lat),]

write.csv(data.noindiv.narm, "LoopedGlobalMasterPopSheet.csv", row.names = FALSE)
