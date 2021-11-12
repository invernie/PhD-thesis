library(parallel)
library(momentuHMM)

filepath <- "path/to/folder/"

raw <- read.csv("rates.csv")
raw_ordered <- raw[with(raw, order(ID, timepoints)),]
raw_cols <- raw_ordered[, c("ID", "D", "stone_dens", "distance", "n_ants")]
#raw_cols$colony <- as.factor(raw_cols$colony)
View(raw_cols)

# check data distribution and estimate starting parameters
hist(raw_cols$D)
nrow(raw_cols[raw_cols$D == 0,]) # for zero-inflation
nrow(raw_cols)
plot(y = raw_ordered$D, x = raw_ordered$timepoints)
plot(y=raw_ordered[raw_ordered$ID == "R54","D"], x = raw_ordered[raw_ordered$ID == "R54", "timepoints"])

hmmdata <- prepData(raw_cols, coordNames = NULL, 
                 covNames = c("stone_dens", "distance", "n_ants"))

# Create cluster of size ncores
ncores <- detectCores() - 1
cl <- makeCluster(getOption("cl.cores", ncores))

niter <- 100 # n of iterations with different starting values

# by ROI analysis (local variables)
# formulaList <- c(~colony+n_ants*distance+stone_dens*distance, ~colony+distance+n_ants*stone_dens, ~colony+stone_dens+n_ants*distance, ~colony+n_ants+stone_dens*distance, ~colony+n_ants+stone_dens+distance, ~colony*distance+stone_dens+n_ants, ~colony*stone_dens+distance+n_ants, ~colony*I(distance^2)+stone_dens+n_ants+distance+I(distance^2), ~colony+stone_dens+distance+n_ants+I(stone_dens^2), ~colony+stone_dens+distance+I(stone_dens^2)+I(stone_dens^2)*distance, ~colony+stone_dens+distance+I(stone_dens^2)+I(stone_dens^2)*n_ants, ~colony+stone_dens, ~colony+stone_dens+I(stone_dens^2), ~colony+n_ants, ~colony+distance)
# subfixList <- c("nd-ds", "ns", "nd", "sd", "n-d-s", "cd", "cs", "cs2", "c-n-d-s-s2", "ds2", "ns2", "c-s", "c-s2", "c-n", "c-d")

# global variable analysis
formulaList <- c(~stone_dens, ~stone_dens+I(stone_dens^2), ~n_ants, ~distance, ~stone_dens+n_ants+distance, ~stone_dens+n_ants+distance+I(stone_dens^2), ~n_ants+stone_dens*distance, ~stone_dens+n_ants*distance)
subfixList <- c("s", "s-s2", "n", "d", "s-n-d", "s-n-d-s2", "n-ds", "s-nd")

for (f in (1:length(subfixList))) {
  
  # Export objects needed in parallelised function to cluster
  clusterExport(cl, list("hmmdata", "fitHMM", "getPar0", "formulaList", "subfixList", "f"))
  
  # Create list of starting values (state 1 = high activity, state 2 = low activity)
  
  allPar0 <- lapply(as.list(1:niter), function(x) {
    
    # D
    mean0 <- runif(2,
                   min = c(0.75, 0.05),
                   max = c(2, 1.5)
    )
    
    # P
    # mean0 <- runif(2,
    #                min = c(0.4, 0.05),
    #                max = c(1, 1.5)
    # )
    
    sd0 <- runif(2,
                 min = c(0.5, 0.5),
                 max = c(2.5, 2.5)
    )
    
    # when running the analysis on ROI subsets or on P, zeromass is needed
    zm0 <- runif(2,
                 min = c(0, 0.05),
                 max = c(0.1, 0.25)
    )
    
    
    return(list(D = c(mean0, sd0, zm0)))
    
  })

  
  #  Fit niter models in parallel
  
  allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
    
    dist <- list(D = "gamma")
    
    m0 <- fitHMM(data = hmmdata, 
                 nbStates = 2, 
                 dist = dist, 
                 Par0 = par0, 
                 nlmPar = list(print.level = 2))
    
    #DM <- list(P = c(mean = formulaList [f], sd = ~1, zeromass = ~1))
    DM <- list(D = c(mean = formulaList [f], sd = ~1))
    Par0 <- getPar0(m0, DM = DM)
    m <- fitHMM(data = hmmdata, 
                nbStates = 2, 
                dist = dist, 
                Par0 = Par0$Par, 
                DM = DM, 
                nlmPar = list(print.level = 2))
    
    return(m)
    
  })
  
  # save models as RDS file
  sn <- paste0(filepath, "allm_", subfixList [f], ".rds")
  saveRDS(c(allm_parallel, allPar0), sn) 
  
}


########################
## No covariate model ##
########################
 
clusterExport(cl, list("hmmdata", "fitHMM"))

# allPar0 <- lapply(as.list(1:niter), function(x) {
# 
#   mean0 <- runif(2,
#                  min = c(0.75, 0.05),
#                  max = c(2, 1.5)
#   )
# 
#   sd0 <- runif(2,
#                min = c(0.5, 0.5),
#                max = c(2, 2)
#   )
# 
#   zm0 <- runif(2,
#                min = c(0, 0.07),
#                max = c(0.1, 0.2)
#   )
# 
# 
#   return(list(C = c(mean0, sd0, zm = zm0)))
# 
# })


allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  
  dist <- list(D = "gamma")
  
  m0 <- fitHMM(data = hmmdata, 
               nbStates = 2, 
               dist = dist, 
               Par0 = par0, 
               nlmPar = list(print.level = 2))
  
  return(m0)
  
})

sn <- paste0(filepath, "allm_0.rds")
saveRDS(c(allm_parallel, allPar0), sn) 


##################################
## Transition probability model ##
##################################

# Export objects needed in parallelised function to cluster
clusterExport(cl, list("hmmdata", "fitHMM"))

allPar0.trans <- lapply(as.list(1:niter), function(x) {

  mean0 <- runif(2,
                 min = c(0.75, 0.05),
                 max = c(2, 1.5)
  )

  sd0 <- runif(2,
               min = c(0.5, 0.5),
               max = c(2, 2)
  )

  # zm0 <- runif(2,
  #              min = c(0, 0.07),
  #              max = c(0.1, 0.2)
  # )


  return(list(D = c(mean0, sd0, zm = zm0)))
  #return(list(P = c(mean0, sd0)))

})

#allPar0.trans <- allPar0

allm_parallel <- parLapply(cl = cl, X = allPar0.trans, fun = function(par0) {
  
  dist <- list(D = "gamma")
  
  m.trans <- fitHMM(data = hmmdata, 
               nbStates = 2, 
               dist = dist, 
               Par0 = par0,
               formula = ~distance, # this is where the covariate is included
               nlmPar = list(print.level = 2))
  
  return(m.trans)
  
})

sn <- paste0(filepath, "allm_trans_d.rds")
saveRDS(c(allm_parallel, allPar0.trans), sn)
