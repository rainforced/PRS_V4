library(PhenotypeSimulator)
source("./utils/utilityFunctions.R")
source("./utils/PhenSim_Config.R")




genFixed <- readRDS(paste(tempdir, "/genFixed_RAW.RData", sep=""))
genBg <- readRDS(paste(tempdir, "/genBg_RAW.RData", sep=""))
noiseFixed <- readRDS(paste(tempdir, "/noiseFixed_RAW.RData", sep=""))
correlatedBg <- readRDS(paste(tempdir, "/correlatedBg_RAW.RData", sep=""))
noiseBg <- readRDS(paste(tempdir, "/noiseBg_RAW.RData", sep=""))



# Total variance proportions have to add up yo 1
total <- shared * h2s *genVar + independent * h2s * genVar +
  shared * (1-h2s) * genVar + independent * (1-h2s) * genVar +
  shared * phi* noiseVar + independent * phi* noiseVar +
  rho * noiseVar + shared * delta * noiseVar + independent * delta * noiseVar

total == 1 

# combine components into final phenotype
Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component +
             genBg_independent_scaled$component + noiseBg_independent_scaled$component +
             genFixed_shared_scaled$component + noiseFixed_shared_scaled$component +
             genFixed_independent_scaled$component + noiseFixed_independent_scaled$component +
             correlatedBg_scaled$component)

write.csv(Y, file = paste(datadir, "/Phenotype.csv", sep="")) 