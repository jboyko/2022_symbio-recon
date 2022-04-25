setwd("2022_symbio-recon/")

require(castor)
require(corHMM)
require(hisse)

phy_EcMAM <- readRDS("data/phy_GBMB_mycEcMAM.RDS")
phy_ErM <- readRDS("data/phy_GBMB_mycErM.RDS")
phy_NMAM <- readRDS("data/phy_GBMB_mycNMAM.RDS")

dat_EcMAM <- read.csv("data/mycEcM-AM.csv")[,-1]
dat_ErM <- read.csv("data/mycErM.csv")[,-1]
dat_NMAM <- read.csv("data/mycNM-AM.csv")[,-1]

total_sample <- 353185 # based on Smith and Brown

# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# DATASET 1: EcMAM
# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# calculate the sampling fraction
sample_f <- length(phy_EcMAM$tip.label)/total_sample
# create the dataset
summary(as.factor(dat_EcMAM$myc))
tip_states <- factor(dat_EcMAM$myc, levels = c("AM", "NM", "EcM", "EcM-AM"))
names(tip_states) <- dat_EcMAM$taxa
tip_states <- tip_states[phy_EcMAM$tip.label]
tip_states <- as.numeric(tip_states)

# CHARACTER DEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
Nstates  = 4
NPstates = 4
proxy_map = c(1,2,3,4)
birth_rate_model = c(1,2,3,4)
death_rate_model = c(1,1,1,1)
transition_matrix = TransMatMakerMuHiSSE()
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_EcMAM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)

# CHARACTER INDEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
Nstates  = 8
NPstates = 4
proxy_map = c(1,2,3,4,1,2,3,4)
birth_rate_model = c(1,1,1,1,2,2,2,2)
death_rate_model = c(1,1,1,1,2,2,2,2)
transition_matrix = TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_EcMAM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)



# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# DATASET 2: ErM
# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# calculate the sampling fraction
sample_f <- length(phy_ErM$tip.label)/total_sample
# create the dataset
summary(as.factor(dat_ErM$myc))
tip_states <- factor(dat_ErM$myc, levels = c("AM", "NM", "EcM", "ErM"))
names(tip_states) <- dat_ErM$taxa
tip_states <- tip_states[phy_ErM$tip.label]
tip_states <- as.numeric(tip_states)

# CHARACTER DEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
Nstates  = 4
NPstates = 4
proxy_map = c(1,2,3,4)
birth_rate_model = c(1,2,3,4)
death_rate_model = c(1,1,1,1)
transition_matrix = TransMatMakerMuHiSSE()
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_ErM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)

# CHARACTER INDEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
# number of states and proxy states 
Nstates  = 8
NPstates = 4
proxy_map = c(1,2,3,4,1,2,3,4)
birth_rate_model = c(1,1,1,1,2,2,2,2)
death_rate_model = c(1,1,1,1,2,2,2,2)
transition_matrix = TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_ErM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)



# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# DATASET 3: NMAM
# - - - - - - - - - - - - - - # - - - - - - - - - - - - - -
# calculate the sampling fraction
sample_f <- length(phy_NMAM$tip.label)/total_sample
# create the dataset
summary(as.factor(dat_NMAM$myc))
tip_states <- factor(dat_NMAM$myc, levels = c("AM", "NM", "EcM", "NM-AM"))
names(tip_states) <- dat_NMAM$taxa
tip_states <- tip_states[phy_NMAM$tip.label]
tip_states <- as.numeric(tip_states)

# CHARACTER DEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
Nstates  = 4
NPstates = 4
proxy_map = c(1,2,3,4)
birth_rate_model = c(1,2,3,4)
death_rate_model = c(1,1,1,1)
transition_matrix = TransMatMakerMuHiSSE()
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_NMAM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)

# CHARACTER INDEPENDENT MODEL
# - - - - - - - - - - - - - -
# number of states and proxy states 
# number of states and proxy states 
Nstates  = 8
NPstates = 4
proxy_map = c(1,2,3,4,1,2,3,4)
birth_rate_model = c(1,1,1,1,2,2,2,2)
death_rate_model = c(1,1,1,1,2,2,2,2)
transition_matrix = TransMatMakerMuHiSSE(hidden.traits = 1, make.null = TRUE)
transition_matrix[is.na(transition_matrix)] <- 0
# run multiple trials to ensure global optimum
Ntrials <- Nthreads <- 1
# fit HiSSE model to tree & tip data
fit = fit_musse(phy_NMAM,
                Nstates                       = Nstates,
                NPstates                      = NPstates,
                proxy_map                     = proxy_map,
                tip_pstates                   = tip_states,
                birth_rate_model              = birth_rate_model,
                death_rate_model              = death_rate_model,
                transition_rate_model         = transition_matrix,
                sampling_fractions            = sample_f,
                Ntrials                       = Ntrials,
                include_ancestral_likelihoods = TRUE,
                Nthreads                      = Nthreads)





