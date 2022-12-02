library(PhenotypeSimulator)
source("./utils/utilityFunctions.R")
source("./utils/PhenSim_Config.R")


# 2. kinship
## kinship estimate based on SNPs
kinship <- getKinship(N=numOfPeople, kinshipfile = kinshipfile,
                      sep = "\t", header = FALSE, verbose = TRUE)
## infinitesimal effect
genBg <- geneticBgEffects(N=numOfPeople, kinship = kinship, P = numOfPheno)

saveRDS(genBg, file=paste(tempdir, "/genBg_RAW.RData", sep=""))
# readRDS(paste(tempdir, "/genBg_RAW.RData", sep="")) - to load the file