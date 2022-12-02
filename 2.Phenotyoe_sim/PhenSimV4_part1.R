library(PhenotypeSimulator)
source("./utils/utilityFunctions.R")
source("./utils/PhenSim_Config.R")


# 1. simulate genetic variant effects
causalSNPs <- getCausalSNPs(N=numOfPeople, NrCausalSNPs = numOfSNPs, chr = chrs,
                            genoFilePrefix = genoFilePrefix,
                            genoFileSuffix = genoFileSuffix,
                            delimiter = " ", format = "oxgen", 
                            verbose = TRUE)
genFixed <- geneticFixedEffects(N = numOfPeople, P = numOfPheno, 
                                X_causal = causalSNPs,
                                pIndependentGenetic = independent,
                                distBeta = "norm", mBeta = 0, sdBeta = 1,
                                verbose = TRUE)
## ADD non-linearity
genFixed$independent <- genFixed$independent^2

saveRDS(genBg, file=paste(tempdir, "/genFixed_RAW.RData", sep=""))
# readRDS(paste(tempdir, "/genFixed_RAW.RData", sep="")) - to load the file