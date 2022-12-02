library(PhenotypeSimulator)

# directories
indir <- "/mnt/jack-5/vperelygin/data/GenSim/Hapgen_b36/60k" # genotypes dir
datadir <- "/mnt/jack-5/vperelygin/data/GenSim/Hapgen_b36/60k/PhenSimV4_200k" # where phen data is saved
#indir <- "/media/oohprivet/Storage Disk/Lab/LargeBatch" # genotypes dir
#datadir <- "/media/oohprivet/Storage Disk/Lab/LargeBatch/PhenSimV4_200k" # genotypes dir
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)
tempdir = "/mnt/jack-5/vperelygin/data/GenSim/Hapgen_b36/60k/PhenSimV4_200k/temp"
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)

kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
                     sep="")
genoFilePrefix <- paste(indir, "/genotypes_", sep="")
genoFileSuffix <- "_hapgen.controls.gen"

# Simulation settings
numOfPeople = 60000
numOfSNPs = 198000 # 22 chr*9000 snps/chr
numOfPheno = 1
chrs = 1:22 

# genetic effects
genVar <- 0.6
noiseVar <- 1 - genVar
totalSNPeffect <- 0.5 # infinitesimal effect = genVar - totalSNPeffect = 0.1
h2s <- totalSNPeffect/genVar 
theta = 0 # proportion of SNPs all traits will share
# non-genetic effect
phi = 0.25 # experimental noise contribution -> Y_noiseBG
delta = 0.75 # non-genetic covariates -> Y_noiseFixed
rho = 0 # correlated non-genetic -> Y_correlatedBg
shared <- 0.0
independent <- 1 - shared