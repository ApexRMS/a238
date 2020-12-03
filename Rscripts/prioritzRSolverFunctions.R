## Source files for 6.2_RunPrioritizationScenarios.R 
## Custom defined prioritizR solvers for different scenario types
## C Tucker


rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}


## Set Solver for single layer objective problems
prob <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}

# Error checking solver for 1 feature layer prblems
genericSolver <- function(inputfeaturesAll, costLayerZ, NumSitesGoalZ){
  inputfeatures <- mask(inputfeaturesAll, costLayerZ) 

tryCatch({  
  tryCatch({
    tryCatch({prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=2) %>% 
        solve(.)}, 
        error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=5) %>% 
            solve(.)}
    )},
    error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=15) %>%
        solve(.)}) 
  },  
  error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=25) %>%
      solve(.)}) 
  
}

# Error checking solver for combine 1 feature layer prblems
genericSolver2 <- function(inputfeaturesAll, costLayerZ, NumSitesGoalZ){
  inputfeatures <- mask(inputfeaturesAll, costLayerZ)
  inputfeatures <- calc(inputfeatures, fun = rescaleR)
  
  tryCatch({  
    tryCatch({
      tryCatch({prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=2) %>% 
          solve(.)}, 
          error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=5) %>% 
              solve(.)}
      )},
      error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=15) %>%
          solve(.)}) 
  },  
  error=function(e){prob(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=25) %>%
      solve(.)}) 
  
}

# Set up minimum shortfall solver for single layer objective problems
probMS <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_min_shortfall_objective(z) %>% #minimize cost surface
    add_relative_targets(0.75) %>%
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}
MSSolver <- function(inputfeaturesAll, costLayerZ, NumSitesGoalZ){
  inputfeatures <- mask(inputfeaturesAll, costLayerZ) 
  
            tryCatch({  
              tryCatch({
                tryCatch({probMS(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=2) %>% 
                    solve(.)}, 
                    error=function(e){probMS(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=5) %>% 
                        solve(.)}
                )},
                error=function(e){probMS(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=15) %>%
                    solve(.)}) 
            },  
            error=function(e){probMS(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=25) %>%
                solve(.)}) 
            }

# Define multiobjective solver, max utility
probMU <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}
MUSolver <- function(inputfeaturesAll, costLayerZ, NumSitesGoalZ){
  inputfeatures <- mask(inputfeaturesAll, costLayerZ) 
  
  tryCatch({  
    tryCatch({
      tryCatch({probMU(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=2) %>% 
          solve(.)}, 
          error=function(e){probMU(costLayerZ, inputfeatures, NumSitesGoalZ, first=F, gap=5) %>% 
              solve(.)}
      )},
      error=function(e){probMU(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=15) %>%
          solve(.)}) 
  },  
  error=function(e){probMU(costLayerZ, inputfeatures, NumSitesGoalZ, first=T, gap=25) %>%
      solve(.)}) 
  }

