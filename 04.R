#  -----------------------------------------------------------------------------------

# file 04: simulate networks from the Bayesian ERGMs

# last updated: 23/04/2025

# ------------------------------------------------------------------------------------



# function to simulate networks from Bayesian models
simulator <- function(thetas, model, nsim){
  
  # model parameters
  theta <- colMeans(thetas)
  
  # wrapper function for igraph local average transitivity
  sims <- stats::simulate(
    model, # model specification
    coef = theta, # model parameters
    seed = 20240517, # Malone's birthday
    nsim = nsim, # number of simulations
    output = "network" # thing to simulate
  )
  
  # return 
  return(sims)
}
g_siren.sim <- simulator(
  thetas = bayes_01.siren$Theta,
  model = bayes_01.siren$formula,
  nsim = 10000
  )
g_togo.sim <- simulator(
  thetas = bayes_02.togo$Theta,
  model = bayes_02.togo$formula,
  nsim = 10000
  )
g_caviar.sim <- simulator(
  thetas = bayes_03.caviar$Theta,
  model = bayes_03.caviar$formula,
  nsim = 10000
  )
g_cielnet.sim <- simulator(
  thetas = bayes_04.cielnet$Theta,
  model = bayes_04.cielnet$formula,
  nsim = 10000
  )
g_cocaine.sim <- simulator(
  thetas = bayes_05.cocaine$Theta,
  model = bayes_05.cocaine$formula,
  nsim = 10000
  )
g_heroin.sim <- simulator(
  thetas = bayes_06.heroin$Theta,
  model = bayes_06.heroin$formula,
  nsim = 10000
  )
g_oversize.sim <- simulator(
  thetas = bayes_07.oversize$Theta,
  model = bayes_07.oversize$formula,
  nsim = 10000
  )
g_montagna.sim <- simulator(
  thetas = bayes_08.montagna$Theta,
  model = bayes_08.montagna$formula,
  nsim = 10000
  )



# close .R script


