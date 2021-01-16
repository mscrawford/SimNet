# Perfect Plasticity Approximation model (Strigul et al. 2008)
# Adapted from Rüger et al. 2020 (code written by Caroline Farrior, cfarrior@gmail.com; https://github.com/cfarrior/Ruger_etal_2020)

# Runs:
#   - Serially on a personal computer (PARALLEL = FALSE, CLUSTER = FALSE)
#   - In parallel on a personal computer with snow (PARALLEL = TRUE, CLUSTER = FALSE)
#   - In parallel on a personal computer with snow (PARALLEL = TRUE, CLUSTER = TRUE)

# This program uses:
# tidyverse (dplyr in parallelize function)
# data.table (fwrite)
# snow (parallelization)


# Global constants ----------------------------------------------------------------------------

# Note:
#   crown radius = 0.5*dbh^0.62; dbh in cm and crown radius in m
#   crown area = phi*dbh^theta

# Commonly used constants
phi <- round(pi * 0.5^2, 4)
theta <- 1.24

# Global constants
PA <- 10000 # plot size: 1 ha
dnot <- 1 # initial tree diameter: 1 cm
mincohortN <- 0.001 # minimum cohort size
maxT <- 2000 # simulation terminal year
deltaT <- 5 # model timestep. FIXED.
cutT <- deltaT


# Global mutable parameters -------------------------------------------------------------------

set.seed(1)

PARALLEL <- TRUE # Run parallel job using snow?
CLUSTER <- FALSE # Run on cluster computer?
printInternalSeedRain <- FALSE # Rather than normal results, print a data.frame of the internal seed rain

# Directories
base.directory <- getwd() # Ensure that your base directory is the PPA folder
species.data.file = "/input/PPA_vitalRates.csv"
out.directory = "/output"

parameterization <- NULL
communities <- NULL


# SimNet parameters ---------------------------------------------------------------------------

# SimNet constants
metacommunityStage <- maxT / 2
averageMonocultureInteralSeedRain <- 551 # This value is calculated through by running only
# monocultures with printInternalSeedRain = TRUE, dividing the N for each species by 5
# (to calculate yearly reproduction) and then averaging across all the monocultures.

# Parameterization variables
nRep = 64
Ninitial = c(1, 2, 4, 8, 16, 32, 64)
SeedRain = c(100)


# Run the PPA model -------------------------------------------------------
run <- function()
{
    parameterization <- parameterize()
    communities <- parameterization$communities
    parameterization <- parameterization$parameterization

    results <- NULL
    if (PARALLEL) {
        results <- run_parallel(parameterization, communities)
    }
    else {
        results <- run_serial(parameterization, communities)
    }

    saveRDS(object = results,
            file = paste0(base.directory, out.directory, "/PPA_output.rds"))
}

# Run PPA on a cluster ------------------------------------------------------------------------
run_parallel <- function(parameterization, communities)
{
    # Create MPI cluster
    cluster <- NULL
    if (CLUSTER) {
        cluster <- snow::getMPIcluster()
    } else {
        cluster <- snow::makeCluster(3)
    }

    ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))

    snow::clusterExport(cl = cluster, c(ex,
                                        "printInternalSeedRain",
                                        "parameterization", "communities",
                                        "phi", "theta",
                                        "PA", "dnot", "mincohortN", "maxT", "deltaT", "cutT",
                                        "metacommunityStage", "averageMonocultureInteralSeedRain"))

    results <- snow::parApply(cl = cluster,
                              X = as.array(1:dim(parameterization)[1]),
                              MARGIN = 1,
                              FUN = function(i)
                              {
                                  parameter_set <- parameterization[i, ]

                                  community <- communities[[parameter_set$Ninitial]][[parameter_set$Rep]]

                                  result = run_simulation(community, parameter_set)

                                  return(result)
                              }
    )

    snow::stopCluster(cluster)

    results <- do.call("rbind", results)

    return(results)
}


# Run PPA serially ----------------------------------------------------------------------------
run_serial <- function(parameterization, communities)
{
    results <- NULL
    for (i in seq(1, dim(parameterization)[1]))
    {
        parameter_set <- parameterization[i, ]

        community <- communities[[parameter_set$Ninitial]][[parameter_set$Rep]]

        results <- rbind(results, run_simulation(community, parameter_set))
    }

    return(results)
}


# Parameterization ----------------------------------------------------------------------------
parameterize = function()
{
    require(dplyr)

    speciesPool <- read.table(paste0(base.directory, species.data.file),
                              sep = ",",
                              header = TRUE)

    # Generate a list of list of data.frames, corresponding to nRep communities for each species pool size (Ninitial).

    communities <- list()

    for (i in Ninitial)
    {
        reps <- list()

        if (i == 1) # with 64 replicates, ensure that each monoculture is run once
        {
            monocultures <- speciesPool %>%
                sample_n(nRep, replace = FALSE)

            for (r in seq(1, nRep))
            {
                reps[[r]] <- monocultures[r, ]
            }

            communities[[1]] <- reps
        }

        else
        {
            for (r in seq(1, nRep))
            {
                community <- speciesPool %>%
                    sample_n(i, replace = FALSE)

                reps[[r]] <- community

                if (i == 64)
                {
                    break # only one 64 species community should be included
                }
            }
            communities[[i]] <- reps
        }
    }

    parameterization <- expand.grid(Ninitial = Ninitial,
                                    Rep = seq(1, nRep),
                                    SeedRain = SeedRain) %>%
        filter(!(Ninitial == 64 & Rep > 1)) # Ensure the presence of only one 64-species community.

    # Because PPA discards species with biomasses of 0, this creates a separate
    # document that contains the starting communities' mixtures.
    print.communities <- list()
    for (Ninitial in 1:length(communities))
    {
        Ninitial.list <- bind_rows(communities[[Ninitial]], .id = "index") %>%
            rename(Rep = index)

        print.communities[[Ninitial]] = Ninitial.list
    }

    print.communities <- bind_rows(print.communities, .id = "index") %>%
        rename(Ninitial = index) %>%
        select(Ninitial, Rep, SpeciesID)

    data.table::fwrite(x = print.communities,
                       file = paste0(base.directory, out.directory, "/PPA_initialCommunities.csv"),
                       sep = "\t")

    return(list(parameterization = parameterization,
                communities = communities))
}


# Run simulation ------------------------------------------------------------------------------
# runs the simulations, and calculates the statistics for each timestep
# community: data.frame of the vital rates of the functional groups (sp, G1, G2, G3, G4, mu1, mu2, mu3, mu4, F)
run_simulation <- function(community,
                           parameter_set)
{
    # Initialize the community
    data <- matrix(1, nrow = dim(community)[1], ncol = 4)

    data[,1] <- dnot;                           # All the trees are 1 cm (dnot) in diameter
    data[,2] <- PA / parameter_set$Ninitial;    # SimNet protocol is to initialize with 1 ind/m^2 & equal # of ind/species
    data[,3] <- 1;                              # By default, they are the in the crown
    data[,4] <- as.numeric(community$SpeciesID) # Identify them by their species

    if (parameter_set$SeedRain > 0)
    {
        externalSeedRain <- generateExternalSeedRain(community, parameter_set$SeedRain, parameter_set$Ninitial)
    }

    results <- NULL

    # Main loop. Remember data matrix has columns:
    # (1) diameter per individual
    # (2) # of individuals
    # (3) crown class (1: Canopy, 2: Understory, 3: 2nd understory, 4: 3rd understory)
    # (4) species
    for (t in seq(0, maxT, by = deltaT))
    {
        startingData <- data

        # Mortality matrix has columns:
        # (1) diameter per dead individual
        # (2) # of dead individuals
        # (3) crown class of dead individuals (not used)
        # (4) species
        mortality <- matrix(0, nrow = 0, ncol = 4)

        for (sp in unique(data[,4]))
        {
            # Step 1: Mortality
            muV <- c(community[community[,1] == sp, seq(6, 9)], 1)
            for (i in seq(1, dim(data)[1])[data[,4] == sp])
            {
                mortality <- rbind(mortality, c(data[i, 1],
                                                data[i, 2] * (1 - (1 - muV[[data[i, 3]]])^deltaT),
                                                data[i, 3],
                                                sp))

                data[i, 2] <- data[i, 2] * (1 - muV[[data[i, 3]]])^deltaT
            }
            data <- data[data[,2] > mincohortN, , drop = FALSE]

            # Step 2: Growth by layer
            GV <- community[community[,1] == sp, seq(2, 5)]
            for (layer in unique(data[data[,4] == sp, , drop = FALSE][,3]))
            {
                data[(data[,4] == sp) & (data[,3] == layer), 1] <-
                    data[(data[,4] == sp) & (data[,3] == layer), 1] + GV[[layer]] * deltaT
            }
        }

        data <- data[data[,2] > mincohortN, , drop = FALSE] # Remove any cohort with too few individuals
        data <- data[data[,1] > 0, , drop = FALSE]          # Remove any cohort with a negative avg. diameter

        # Step 3: Reproduce
        # Step 3a. Internal seed rain
        internalSeedRain <- generateInternalSeedRain(data, community)
        data <- rbind(data, internalSeedRain)

        # Step 3b. SimNet external seed rain
        if (parameter_set$SeedRain > 0 & t <= metacommunityStage)
        {
            data <- rbind(data, externalSeedRain)
        }

        # Step 4: Assign crown class
        CA <- sum(phi * data[,1]^theta * data[,2]) # calculate total crown area of individuals in the plot
        if (CA <= PA)
        {
            # if less than the ground area, no need for CCassign. Everyone is in the canopy.
            data[,3] = 1
        } else {
            # if greater than the ground area, go through CCassign
            data <- CCassign_manylayers_continuous(data, PA, deltaT)
        }

        data <- data[data[,3] < 5, , drop = FALSE] # kill plants that fall into a fifth layer immediately

        # Step 5: Record
        if (floor(t/cutT) == t/cutT & t != 0) # Records every cutT timesteps
        {
            if (printInternalSeedRain)
            {
                stopifnot(parameter_set$SeedRain == 0 & parameter_set$Ninitial == 1)

                # Print the internal seed dispersal that has occurred over the last time step
                # When calculating average yearly monoculture seed rain, remember to divide N by the time step!
                results <- rbind(results, calculate_output(internalSeedRain, startingData, mortality, t, parameter_set, community))
            }
            else
            {
                results <- rbind(results, calculate_output(data, startingData, mortality, t, parameter_set, community))
            }
        }
    }

    return(results)
}


generateInternalSeedRain <- function(data, community)
{
    template_internalSeedRain <- NULL
    for (sp in unique(data[,4]))
    {
        t_sp_cohorts <- data[data[,4] == sp, , drop = FALSE]

        nbaabg <- n_ba_agb(t_sp_cohorts, community)

        template_internalSeedRain <- rbind(template_internalSeedRain,
                                           c(dnot,
                                             nbaabg$ba * community[community[,1] == sp, , drop = FALSE]$F,
                                             4,
                                             sp))
    }

    template_internalSeedRain <- as.data.frame(template_internalSeedRain)
    colnames(template_internalSeedRain) <- c("diameter", "N", "layer", "SpeciesID")
    template_internalSeedRain <- merge(template_internalSeedRain, community[, c(1, 5)], by = "SpeciesID")

    # Generate deltaT cohorts, each at a different stage of growth
    internalSeedRain <- NULL
    for (i in seq(deltaT - 1, 0, -1))
    {
        baby <- template_internalSeedRain
        baby$diameter <- baby$diameter + baby$G4*i;
        # Because the census occurs on a five year time step, seedling mortality is implicitly considered
        baby <- as.matrix(baby[, c(2, 3, 4, 1)])

        internalSeedRain <- rbind(internalSeedRain, baby)
    }

    return(internalSeedRain)
}


generateExternalSeedRain <- function(community, SeedRain, Ninitial)
{
    stopifnot(SeedRain > 0)

    # Main matrix with columns:
    # [,1] cohort diameter. Initially, all the trees are 1 cm (dnot) in diameter.
    # [,2] number of individuals. External seed rain is a % of the average monoculture internal rain.
    # [,3] crown class. By default, they are in the lowest understory.
    # [,4] species
    template_externalSeedRain <- matrix(1, nrow = dim(community)[1], ncol = 4)

    template_externalSeedRain[,1] <- dnot;
    template_externalSeedRain[,2] <- (SeedRain/100) * (averageMonocultureInteralSeedRain / Ninitial);
    template_externalSeedRain[,3] <- 4;
    template_externalSeedRain[,4] <- as.numeric(community$SpeciesID)

    externalSeedRain <- NULL
    for (i in seq(deltaT - 1, 0, -1))
    {
        baby <- template_externalSeedRain;
        baby[,1] <- baby[,1] + community$G4*i;
        # Because the census occurs on a five year time step, seedling mortality is implicitly considered
        externalSeedRain <- rbind(externalSeedRain, baby)
    }

    return(externalSeedRain)
}


# Calculate simulation output -----------------------------------------------------------------
calculate_output <- function(data, startingData, mortality, year, parameter_set, community)
{
    out <- NULL

    for (species in unique(data[,4]))
    {
        ss_t1 <- data[data[,4] == species, , drop = FALSE]
        ss_t0 <- startingData[startingData[,4] == species, , drop = FALSE]
        ss_t0_mort <- mortality[mortality[,4] == species, , drop = FALSE]

        nbaabg_t1 <- n_ba_agb(ss_t1, community)
        nbaabg_t0 <- n_ba_agb(ss_t0, community)
        nbaabg_t0_mort <- n_ba_agb(ss_t0_mort, community)

        productivity <- nbaabg_t1$agb - nbaabg_t0$agb + nbaabg_t0_mort$agb
        productivity <- ifelse(productivity < 0, 0, productivity)

        sp.out <- data.frame(Model = "PPA",
                             Ninitial = parameter_set$Ninitial,
                             Rep = parameter_set$Rep,
                             SeedRain = parameter_set$SeedRain,
                             SpeciesID = species,
                             `F` = community[community[,1] == species, , drop = FALSE]$F,
                             Stage = ifelse(year <= metacommunityStage, "assembly", "disassembly"),
                             Year = year,
                             N = nbaabg_t1$n,
                             BasalArea = nbaabg_t1$ba,
                             Biomass = nbaabg_t1$agb,
                             Productivity = productivity)

        out <- rbind(out, sp.out)
    }

    return(out)
}


# Calculate abundance, basal area, and aboveground biomass per species -------------------------
n_ba_agb <- function(data, community)
{
    stopifnot(length(unique(data[,4])) %in% c(0, 1)) # This function should only be run on a per species basis

    n <- sum(data[,2])
    ba <- sum( (data[,1]/200)^2 * pi * data[,2] )

    wd <- community[community$SpeciesID == unique(data[,4]), ]$wd
    agb <- sum(wd *
                   exp(-1.499 + 2.148*log(data[,1]) +
                           0.207*log(data[,1])^2 -
                           0.0281*log(data[,1])^3) /
                   1000 * data[,2])

    return(list(n = n, ba = ba, agb = agb))
}


# Main cohorts function -----------------------------------------------------------------------
# assigns crown class based on the PPA assumption
# assumes all individuals have the same crown area and height allometries
# assumes that CAtot > PA
# works for any number of layers
CCassign_manylayers_continuous = function(data, PA, deltaT)
{
    CAv <- phi*data[,1]^theta
    data <- data[order(CAv,decreasing = TRUE), ]
    CAv <- CAv[order(CAv,decreasing = TRUE)]
    cohortCAv <- CAv*data[,2]
    cacaV <- cumsum(cohortCAv)

    data[,3] <- 1

    for (layers in seq(1, floor(sum(cohortCAv)/PA)))
    {
        # make a vector cacaV where the ith entry is the crown area of the i'th cohort plus all cohorts
        # with individuals of greater diameter.
        CAv <- phi*data[,1]^theta
        data <- data[order(CAv, decreasing = TRUE), ]
        CAv <- CAv[order(CAv, decreasing = TRUE)]
        cohortCAv <- CAv * data[,2]
        cacaV <- cumsum(cohortCAv)

        # pull out individuals that are definitely in the canopy and understory
        und <- data[cacaV > PA * layers, , drop = FALSE]
        can <- data[cacaV <= PA * layers, , drop = FALSE]

        # split the first cohort in the understory to fill the leftover open canopy space
        canCA <- max(0, sum(phi * can[,1]^theta * can[,2]))
        tosplit <- und[1, , drop = FALSE]
        opencan <- layers * PA - canCA
        splitind_incan <- opencan / (phi * tosplit[1, 1]^theta)
        und[1, 2] <- und[1, 2] - splitind_incan
        tosplit[,2] <- splitind_incan
        tosplit[,3] <- layers
        can <- rbind(can, tosplit)
        if (max(can[,3]) != layers)
        {
            print("error")
        }
        und[,3] <- layers + 1

        # piece the data back together
        data <- rbind(can, und)

        data <- data[data[,2] > mincohortN, , drop = FALSE]
    }

    return(data)
}


# Run -----------------------------------------------------------------------------------------
run()
