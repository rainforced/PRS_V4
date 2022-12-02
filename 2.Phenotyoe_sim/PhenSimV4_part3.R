library(PhenotypeSimulator)
source("./utils/utilityFunctions.R")
source("./utils/PhenSim_Config.R")



# 3. simulate 4 different confounder effects
## * 1 binomial covariate effect 
## * 1 categorical (3 categories) covariate effects
## * 1 normally distributed independent and shared covariate effects
noiseFixed <- noiseFixedEffects(N = numOfPeople, P = numOfPheno, NrFixedEffects = 3,
                                NrConfounders = c(1, 1, 1),
                                distConfounders = c("bin", "cat_unif", "norm"),
                                probConfounders = 0.2,
                                catConfounders = c(3))

saveRDS(noiseFixed, file=paste(tempdir, "/noiseFixed_RAW.RData", sep=""))

#  4. simulate correlated effects with max correlation of 0.8
correlatedBg <- correlatedBgEffects(N = numOfPeople, P = numOfPheno, pcorr = 0.8)

saveRDS(correlatedBg, file=paste(tempdir, "/correlatedBg_RAW.RData", sep=""))

# 5. simulate observational noise effects
noiseBg <- noiseBgEffects(N = numOfPeople, P = numOfPheno)

saveRDS(noiseBg, file=paste(tempdir, "/noiseBg_RAW.RData", sep=""))