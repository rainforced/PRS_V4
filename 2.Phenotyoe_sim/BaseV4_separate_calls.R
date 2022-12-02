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


# 2. kinship
## kinship estimate based on SNPs
kinship <- getKinship(N=numOfPeople, kinshipfile = kinshipfile,
                      sep = "\t", header = FALSE, verbose = TRUE)
## infinitesimal effect
genBg <- geneticBgEffects(N=numOfPeople, kinship = kinship, P = numOfPheno)

saveRDS(genBg, file=paste(tempdir, "/genBg_RAW.RData", sep=""))
# readRDS(paste(tempdir, "/genBg_RAW.RData", sep="")) - to load the file


# 3. simulate 4 different confounder effects
## * 1 binomial covariate effect 
## * 1 categorical (3 categories) covariate effects
## * 1 normally distributed independent and shared covariate effects
noiseFixed <- noiseFixedEffects(N = numOfPeople, P = numOfPheno, NrFixedEffects = 3,
                                NrConfounders = c(1, 1, 1),
                                pIndependentConfounders = c(1, 1, 1),
                                distConfounders = c("bin", "cat_unif", "norm"),
                                probConfounders = 0.2,
                                catConfounders = c(3))


#  4. simulate correlated effects with max correlation of 0.8
correlatedBg <- correlatedBgEffects(N = numOfPeople, P = numOfPheno, pcorr = 0.8)


# 5. simulate observational noise effects
noiseBg <- noiseBgEffects(N = numOfPeople, P = numOfPheno)


# rescale phenotype components
genFixed_shared_scaled <- rescaleVariance(genFixed$shared, shared * h2s *genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent,
                                               independent * h2s *genVar)
genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent,
                                            independent * (1-h2s) * genVar)
noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * phi* noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent,
                                              independent * phi* noiseVar)
correlatedBg_scaled <- rescaleVariance(correlatedBg$correlatedBg,
                                       shared * rho * noiseVar)
noiseFixed_shared_scaled <- rescaleVariance(noiseFixed$shared, shared * delta *
                                              noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent,
                                                 independent * delta * noiseVar)

# Total variance proportions have to add up yo 1
total <- shared * h2s *genVar + independent * h2s * genVar +
  shared * (1-h2s) * genVar + independent * (1-h2s) * genVar +
  shared * phi* noiseVar + independent * phi* noiseVar +
  rho * noiseVar + shared * delta * noiseVar + independent * delta * noiseVar

total == 1 


Y <- scale(genBg_independent_scaled$component + noiseBg_independent_scaled$component +
             genFixed_independent_scaled$component + noiseFixed_independent_scaled$component)
genFixed_independent_scaled

write.csv(Y, file = "Phenotype.csv") 
################################################################################