library(PhenotypeSimulator)
# directories
indir <- "/media/oohprivet/Storage Disk/Lab/LargeBatch" # genotypes dir
datadir <- "/media/oohprivet/Storage Disk/Lab/LargeBatch/PhenSimV4" # where phen data is saved
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)
tempdir = "/media/oohprivet/Storage Disk/Lab/LargeBatch/PhenSimV4/temp"
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)

intermediatedir = 
kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
                     sep="")
# settings
numOfPeople = 10000
numOfPheno = 1

genVar <- 0.6
noiseVar <- 1 - genVar
totalSNPeffect <- 0.01
h2s <- totalSNPeffect/genVar
phi <- 0.6
rho <- 0.1
delta <- 0.3
shared <- 0.8
independent <- 1 - shared

# kinship
## kinship estimate based on SNPs
kinship <- getKinship(N=numOfPeople, kinshipfile = kinshipfile,
                      sep = "\t", header = FALSE, verbose = TRUE)

## infinitesimal effect
genBg <- geneticBgEffects(N=numOfPeople, kinship = kinship, P = numOfPheno)


saveRDS(genBg, file=paste(tempdir, "/genBg_RAW.RData", sep=""))

genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent,
                                            independent * (1-h2s) * genVar)
