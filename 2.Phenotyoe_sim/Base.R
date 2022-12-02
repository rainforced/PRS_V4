library(PhenotypeSimulator)
# directories
indir <- "/mnt/jack-5/vperelygin/data/GenSim/Hapgen_b36/100k" # genotypes dir
datadir <- "/mnt/jack-5/vperelygin/data/GenSim/Hapgen_b36/100k/PhenSim" # where phen data is saved
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)

kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
                     sep="")
genoFilePrefix <- paste(indir, "/genotypes_", sep="")
genoFileSuffix <- "_hapgen.controls.gen"

# Simulation settings
numOfPeople = 100000
numOfSNPs = 499400 # 22 chr*22700 snps/chr
numOfPheno = 1
chrs = 1:22 

# genetic effects
genVar <- 0.6
totalSNPeffect <- 0.5 # infinitesimal effect = genVar - totalSNPeffect = 0.1
h2s <- totalSNPeffect/genVar 
theta = 0 # proportion of SNPs all traits will share

# non-genetic effect
phi = 0.25 # experimental noise contribution -> Y_noiseBG
delta = 0.75 # non-genetic covariates -> Y_noiseFixed
rho = 0 # correlated non-genetic -> Y_correlatedBg

# simulate phenotype with one phenotype component
simulation <- runSimulation(N = numOfPeople, P = numOfPheno, cNrSNP=numOfSNPs,
                            format = "oxgen", 
                            genoFilePrefix = genoFilePrefix,
                            genoFileSuffix = genoFileSuffix,
                            genoDelimiter = " ",
                            chr = chrs, theta=theta,
                            # genetic effects:
                            genVar = genVar, h2s = h2s,
                            ## genetic effects distribution parameters
                            mBetaGenetic = 0, sdBetaGenetic = 0.2,
                            # noise effects:
                            phi = phi, delta = delta,
                            ## non-genetic covariate simulation
                            NrFixedEffects = 3, NrConfounders = c(1, 1, 1),
                            distConfounders = c("bin", "cat_unif", "norm"),
                            probConfounders = 0.2, catConfounders = c(3),
                            ## kinship info
                            kinshipfile = kinshipfile,
                            kinshipDelimiter="\t",
                            kinshipHeader=FALSE,
                            verbose = TRUE, seed=43)


outdirectory <- savePheno(simulation, directory = datadir,
                          format=c("csv", "plink"),
                          saveIntermediate=TRUE)