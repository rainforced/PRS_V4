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

saveRDS(genFixed, file=paste(tempdir, "/genFixed_200k_RAW.RData", sep=""))
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
saveRDS(noiseFixed, file=paste(tempdir, "/noiseFixed_RAW.RData", sep=""))


# 5. simulate observational noise effects
noiseBg <- noiseBgEffects(N = numOfPeople, P = numOfPheno)
saveRDS(noiseBg, file=paste(tempdir, "/noiseBg_RAW.RData", sep=""))

# rescale phenotype components
genFixed_independent_scaled <- rescaleVariance(genFixed$independent,
                                               independent * h2s *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent,
                                            independent * (1-h2s) * genVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent,
                                              independent * phi* noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent,
                                                 independent * delta * noiseVar)

# Total variance proportions have to add up yo 1
total <- independent * h2s * genVar + independent * (1-h2s) * genVar +
  independent * phi* noiseVar + independent * delta * noiseVar

total == 1 



Y <- scale(genBg_independent_scaled$component + 
             noiseBg_independent_scaled$component +
             genFixed_independent_scaled$component + 
             noiseFixed_independent_scaled$component)


write.csv(Y, file = paste(datadir, "/Phenotype_22k.csv", sep = ''))
write.csv(causalSNPs, file = paste(datadir, "/causalSNPs_22k.csv", sep = ''))
################################################################################
#writeStandardOutput(directory = datadir,
#                    phenotypes = Y,
#                    genotypes = causalSNPs, #??
#                    additionalPhenotypes = NULL,
#                    covariates = NULL,
#                    kinship = NULL,
#                    eval_kinship = NULL,
#                    evec_kinship = NULL,
                    #id_samples,
                    #id_snps,
                    # id_phenos,
#                    outstring = NULL,
#                    standardInput_samples = NULL,
#                    standardInput_genotypes = NULL,
#                    format = c("plink", "csv"),
#                    intercept_gemma = FALSE,
#                    nameAdditional = "_nonLinear",
#                    verbose = TRUE)
