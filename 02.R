#  -----------------------------------------------------------------------------------

# file 02: estimate the Bayesian exponential random graph models

# last updated: 15/04/2025

# ------------------------------------------------------------------------------------


# function to estimate the variance-covariance matrix for the priors
prior_var <- function(g){
  
  # require
  require(ergm); require(network)
  
  # set seed for replication
  set.seed(20110210) # Halle's birthday
  
  # formula
  formula <- g ~ edges + gwdegree(decay = 0.50, fixed = TRUE) + gwdsp(decay = 0.50, fixed = TRUE) + gwesp(decay = 0.50, fixed = TRUE) + degcor
  
  # prior means
  thetas <- c(-log(network::network.size(g)), -1.00, 0.50, 0.50, -1.00)
  
  # simulate networks
  sims <- simulate(
    formula, 
    coef = thetas, 
    nsim = 10,
    seed = 20110210, # Halle's birthday
    output = "network"
    )
  
  # number of simulations
  n <- length(sims)
  
  # empty vector to store the results 
  results <- list()
  
  # estimate the model for each of the simulations
  for(i in 1:n){
  results[[i]] <- Bergm::bergm(
    sims[[i]] ~ edges + gwdegree(decay = 1.00, fixed = TRUE) + gwdsp(decay = 1.00, fixed = TRUE) + gwesp(decay = 1.00, fixed = TRUE) + degcor,
    gamma = 0.1 # empirically, gamma ranges 0.1 - 1.5, where smaller gamma leads to better acceptance in big graphs
    )
  }
  
  # calculate the variance-covariance matrix
  vcov <- lapply(results, function(model) stats::cov(model$Theta))
  
  # 3D array of variance-covariance matrices
  vcov <- simplify2array(vcov) # number of parameters x number of parameters x number of models
  vcov <- apply(vcov, c(1, 2), mean)
  print(heatmap(vcov))
  return(vcov)
}
test <- prior_var(g_siren)





# posterior parameter estimation -----------------------------------------------
bayes <- function(y, x, mu, sigma){
  i <- nrow(y) * 2 # burn in iterations to begin the MCMC run 
  k <-    250 # sample iterations
  h <-    5 * 2 # chains in the posterior distribution ... approximately twice the number of model parameters 
  n <-    h * k # per Caimo & Friel (2011), auxiliary chain = # chains (h) * sample iterations (k)
  
  # load 'bergm' package
  require('Bergm')
  
  # set seed for replication
  set.seed(20110210) # Halle's birthday
  
  # prior means
  # c <- -log(network::network.size(y))
  d <- sna::gden(g_siren, mode = "graph")
  y0 <- log(d/(1-d))
  mu <- c(
         y0, # density
      -1.00, # centralization, negative terms 
       0.50, 
       0.50, 
    -100.00
    )
  
  # prior variance
  # sigma <- diag(10, nrow = 5, ncol = 5)
  # sigma[upper.tri(sigma)] <- 0
  # sigma[lower.tri(sigma)] <- 0
  
  # estimate the model
  bayes <- Bergm::bergm(
    x, # equation
    # prior.mean = mu,
    # prior.sigma = sigma,
    gamma = 0.1 # empirically, gamma ranges 0.1 - 1.5, where smaller gamma leads to better acceptance in big graphs
    )
  return(bayes)
}



# compute Bayes' factors -------------------------------------------------------
# https://cran.r-project.org/web/packages/BFpack/vignettes/vignette_BFpack.html
BF <- function(model, priors){ 
  
  # required packages
  require("BFpack")
  
  # for replication
  set.seed(20240517) # Bugsy's birthday
  
  # don't run
  # this function assumes equal prior probabilities of p = 1/3 = 0.33
  # bf <- BFpack::BF(model)
  
  # function to calculate Bayes' factor
  # https://github.com/jomulder/BFpack
  # https://bfpack.info
  # https://cran.r-project.org/web/packages/BFpack/vignettes/vignette_BFpack.html
  bf <- BFpack::BF(
    model,
    prior.hyp = priors
    )
  
  # hypothesis test
  bf <- bf$PHP_exploratory
  bf <- round(bf, digits = 3) # round posterior probabilities to three decimals
  
  # label coefficients
  row.names(bf) <- c("gwdegree", "gwdsp", "gwesp", "degcor")
  
  # interpretation of the hypothesis tests
  # see Mulder et al., p.48: https://www.sciencedirect.com/science/article/pii/S0378873323000801?via%3Dihub
  cat(message("Bayesian hypothesis test -- posterior probabilities of the hypotheses for each model parameter:"))
  cat("\n")
  print(bf)
  cat("\n")
  cat(message("The columns indicate the probability that the coefficient = 0, < 0, and > 0.")) # e.g., P(B > 0 | Y)
  cat("\n")
  cat(message("Note: Bayesian hypothesis tests not provided for the edges coefficient because an improper flat prior is used for this parameter."))
}



# goodness-of-fit --------------------------------------------------------------
GOF <- function(bayes){
  
  # required packages
  require("Bergm")
  
  # set seed for replication
  set.seed(20110211) # Halle's birthday
  
  # inputs
  n <- 10000 # graph simulations
  i <- 15 # degree distribution range
  j <- 15 # geodisstance range 
  k <- 15 # edgewise-shared partners range
  
  # goodness-of-fit tests
  gof <- Bergm::bgof(
    bayes,
    n.deg  = i,
    n.dist = j,
    n.esp  = k,
    directed = F,    # symmetric graph
    sample.size = n # random graph realizations
    )
  return(gof)
}



# estimate the model -----------------------------------------------------------
bayes_01.siren <- bayes(
  y = g_siren,
  x = g_siren ~ edges + 
    gwdegree(decay = 0.5, fixed = TRUE) + # gwdegree(decay = 0.5, fixed = TRUE) 
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_01.siren)

# hypothesis test
test <- BF(bayes_01.siren, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_01.siren <- GOF(bayes_01.siren)



# estimate the model -----------------------------------------------------------
bayes_02.togo <- bayes(
  y = g_togo,
  x = g_togo ~ edges + 
    gwdegree(decay = 0.70, fixed = TRUE) + # gwdegree(decay = 0.5, fixed = TRUE)
    gwdsp(decay = 1.0, fixed = TRUE) +
    gwesp(decay = 1.0, fixed = TRUE) +
    degcor
    )
summary(bayes_02.togo)

# hypothesis test
BF(bayes_02.togo, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_02.togo <- GOF(bayes_02.togo)



# estimate the model -----------------------------------------------------------
bayes_03.caviar <- bayes(
  y = g_caviar,
  x = g_caviar ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) + # gwdegree(decay = 2.0, fixed = TRUE)
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_03.caviar)

# hypothesis test
BF(bayes_03.caviar, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_03.caviar <- GOF(bayes_03.caviar)



# estimate the model -----------------------------------------------------------
bayes_04.cielnet <- bayes(
  y = g_cielnet,
  x = g_cielnet ~ edges + 
    gwdegree(decay = 1.0, fixed = TRUE) +
    gwdsp(decay = 0.2, fixed = TRUE) +
    gwesp(decay = 0.2, fixed = TRUE) +
    degcor
    )
summary(bayes_04.cielnet)

# hypothesis test
BF(bayes_04.cielnet, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_04.cielnet <- GOF(bayes_04.cielnet)



# estimate the model -----------------------------------------------------------
bayes_05.cocaine <- bayes(
  y = g_cocaine,
  x = g_cocaine ~ edges + 
    gwdegree(decay = 1.2, fixed = TRUE) + # gwdegree(decay = 1.5, fixed = TRUE)
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_05.cocaine)

# hypothesis test
BF(bayes_05.cocaine, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_05.cocaine <- GOF(bayes_05.cocaine)



# estimate the model -----------------------------------------------------------
bayes_06.heroin <- bayes(
  y = g_heroin,
  x = g_heroin ~ edges + 
    gwdegree(decay = 0.30, fixed = TRUE) +  # gwdegree(decay = 2.5, fixed = TRUE)
    gwdsp(decay = 0.3, fixed = TRUE) +
    gwesp(decay = 0.3, fixed = TRUE) +
    degcor
   )
summary(bayes_06.heroin)

# hypothesis test
BF(bayes_06.heroin, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_06.heroin <- GOF(bayes_06.heroin)



# estimate the model -----------------------------------------------------------
bayes_07.oversize <- bayes(
  y = g_oversize,
  x = g_oversize ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) +
    gwdsp(decay = 1.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_07.oversize)

# hypothesis test
BF(bayes_07.oversize, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_07.oversize <- GOF(bayes_07.oversize)



# estimate the model -----------------------------------------------------------
bayes_08.montagna <- bayes(
  y = g_montagna,
  x = g_montagna ~ edges + 
    gwdegree(decay = 2.0, fixed = TRUE) + # gwdegree(decay = 2.0, fixed = TRUE)
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_08.montagna)

# hypothesis test
BF(bayes_08.montagna, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_08.montagna <- GOF(bayes_08.montagna)





# function to calculate the goodness-of-fit statistics from the models 
gof_stats <- function(obs.data, sim.data){
  
  # required packages
  require("matrixStats"); require("dplyr"); require("tidyr"); require("magrittr")
  
  # number of simulations 
  n <- ncol(sim.data)
  
  # range for the distribution
  range <- 1:16
  
  # observed distribution, and restrict the range
  obs <- obs.data
  obs <- obs[range] # constrict range of the distribution
  obs <- as.data.frame(obs)
  colnames(obs) <- c("obs")
  obs <- dplyr::mutate(obs, stat = 0:15) # add column for statistic
  
  # simulated distribution, and restrict the range
  sim <- sim.data
  
  # compute quantile range for plot
  
  # 50% quantile i.e., median
  median <- matrixStats::rowQuantiles(sim, probs = 0.50, na.rm = T) # median <- matrixStats::rowMedians(sim)
  
  # 25% quantile
  Q1 <- matrixStats::rowQuantiles(sim, probs = 0.25, na.rm = T)
  
  # 75% quantile
  Q3 <- matrixStats::rowQuantiles(sim, probs = 0.75, na.rm = T)
  
  # compute lower quantile for the error bars
  ci.lo <- matrixStats::rowQuantiles(sim, probs = 0.025, na.rm = T)
  
  # compute higher quantile for the error bars
  ci.hi <- matrixStats::rowQuantiles(sim, probs = 0.975, na.rm = T)
  
  # join
  stats <- as.data.frame(cbind(median, Q1, Q3, ci.lo, ci.hi))
  stats <- stats[range, ] # constrict range of the distribution
  colnames(stats) <- c("median", "Q1", "Q3", "ci.lo", "ci.hi")
  stats <- dplyr::mutate(stats, stat = 0:15) # add column for statistic
  
  # simulated data points
  sim <- as.data.frame(sim)
  sim <- sim[range, ] # constrict range of the distribution
  
  # add statistic as first colum
  stat <- rownames(sim); stat <- as.numeric(stat); stat <- as.data.frame(stat)
  sim <- cbind(stat, sim)
  # colnames(sim)[2:n] <- paste("stat", colnames(sim)[2:n], sep = "")
  
  # reshape from wide to long
  n <- ncol(sim)
  `%>%` <- magrittr::`%>%`
  sim <- sim %>% 
    tidyr::pivot_longer(
      cols = 2:n,
      names_to = "simulation",
      values_to = "value"
      )
  
  # list
  data <- list(obs, sim, stats)
  
  # data
  return(data)
}

# goodness-of-fit for the degree distribution
gof_01.siren.stats.deg <- gof_stats(
  obs.data = gof_01.siren$obs.degree, 
  sim.data = gof_01.siren$sim.degree
  )
gof_02.togo.stats.deg <- gof_stats(
  obs.data = gof_02.togo$obs.degree, 
  sim.data = gof_02.togo$sim.degree
  )
gof_03.caviar.stats.deg <- gof_stats(
  obs.data = gof_03.caviar$obs.degree, 
  sim.data = gof_03.caviar$sim.degree
  )
gof_04.cielnet.stats.deg <- gof_stats(
  obs.data = gof_04.cielnet$obs.degree, 
  sim.data = gof_04.cielnet$sim.degree
  )
gof_05.cocaine.stats.deg <- gof_stats(
  obs.data = gof_05.cocaine$obs.degree, 
  sim.data = gof_05.cocaine$sim.degree
  )
gof_06.heroin.stats.deg <- gof_stats(
  obs.data = gof_06.heroin$obs.degree,
  sim.data = gof_06.heroin$sim.degree
  )
gof_07.oversize.stats.deg <- gof_stats(
  obs.data = gof_07.oversize$obs.degree, 
  sim.data = gof_07.oversize$sim.degree
  )
gof_08.montagna.stats.deg <- gof_stats(
  obs.data = gof_08.montagna$obs.degree, 
  sim.data = gof_08.montagna$sim.degree
  )

# goodness-of-fit for the edgewise shared partner distribution i.e., cliques of different sizes
gof_01.siren.stats.esp <- gof_stats(
  obs.data = gof_01.siren$obs.esp, 
  sim.data = gof_01.siren$sim.esp
  )
gof_02.togo.stats.esp <- gof_stats(
  obs.data = gof_02.togo$obs.esp, 
  sim.data = gof_02.togo$sim.esp
  )
gof_03.caviar.stats.esp <- gof_stats(
  obs.data = gof_03.caviar$obs.esp, 
  sim.data = gof_03.caviar$sim.esp
  )
gof_04.cielnet.stats.esp <- gof_stats(
  obs.data = gof_04.cielnet$obs.esp, 
  sim.data = gof_04.cielnet$sim.esp
  )
gof_05.cocaine.stats.esp <- gof_stats(
  obs.data = gof_05.cocaine$obs.esp, 
  sim.data = gof_05.cocaine$sim.esp
  )
gof_06.heroin.stats.esp <- gof_stats(
  obs.data = gof_06.heroin$obs.esp,
  sim.data = gof_06.heroin$sim.esp
  )
gof_07.oversize.stats.esp <- gof_stats(
  obs.data = gof_07.oversize$obs.esp, 
  sim.data = gof_07.oversize$sim.esp
  )
gof_08.montagna.stats.esp <- gof_stats(
  obs.data = gof_08.montagna$obs.esp, 
  sim.data = gof_08.montagna$sim.esp
  )



# construct goodness-of-fit plots ----------------------------------------------
gof_plot <- function(data, title){
  
  # set seed for replication
  set.seed(20190816) # Maeve's birthday
  
  # function to set the labels for the y-axis
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # observed data
  obs <- data[[1]]
  
  # simulated data
  sim <- data[[3]]
  
  # construct the boxplot
  fig <- ggplot2::ggplot(data = sim, ggplot2::aes(x = stat, y = median)) +
    # error bars to visualize the confidence intervals
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black", alpha = 1.0) +
    # boxplot
    ggplot2::geom_boxplot(
      data = sim,
      mapping = ggplot2::aes(lower = Q1, upper = Q3, middle = median, ymin = ci.lo, ymax = ci.hi, group = stat), 
      stat = "identity",
      fill = "white", color = "black", width = 0.5, alpha = 1.0
      ) +
    # line graph for the degree distribution
    ggplot2::geom_line(data = obs, ggplot2::aes(x = stat, y = obs), color = "red", linetype = "solid", size = 1, alpha = 0.8) +
    ggplot2::scale_y_continuous(limits = c(0.00, 1.00), breaks = seq(0, 1.00, 0.25), labels = scaleFUN) +
    # ggplot2::scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, 1)) +
    ggplot2::labs(
      title = title,
      y = NULL,
      x = NULL
      ) +
    ggthemes::theme_few() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 5, face = "plain"),
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      )
  return(fig)
}



# function to output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 2.5, 
    height = 2, 
    device = 'pdf', 
    dpi = 700
    )
}

# appendix fig 3 degree distribution
output(plot = gof_plot(
  data = gof_01.siren.stats.deg, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA3a.pdf"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.deg, 
  title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA3b.pdf"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.deg, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA3c.pdf"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.deg, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA3d.pdf"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.deg, 
  title = "(E) NEW YORK CITY COCAINE TRAFFICKERS"), 
  filename = "figA3e.pdf"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.deg, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"), 
  filename = "figA3f.pdf"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.deg, 
  title = "(G) 'NDRANGHETA - OPERATION OVERSIZE"), 
  filename = "figA3g.pdf"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.deg, 
  title = "(H) COSA NOSTRA - OPERATION MONTAGNA"), 
  filename = "figA3h.pdf"
  )

# appendix fig 4 edgewise shared partner distribution i.e., the distribution of cliques of different sizes
output(plot = gof_plot(
  data = gof_01.siren.stats.esp, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA4a.pdf"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.esp, title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA4b.pdf"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.esp, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA4c.pdf"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.esp, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA4d.pdf"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.esp, 
  title = "(E) NEW YORK CITY COCAINE TRAFFICKERS"), 
  filename = "figA4e.pdf"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.esp, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"), 
  filename = "figA4f.pdf"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.esp, 
  title = "(G) 'NDRANGHETA - OPERATION OVERSIZE"), 
  filename = "figA4g.pdf"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.esp, 
  title = "(H) COSA NOSTRA - OPERATION MONTAGNA"), 
  filename = "figA4h.pdf"
  )



# close .R script





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



# function to calculate the goodness-of-fit statistics from the models 
gof_stats <- function(obs.data, sim.data){
  
  # set seed for replication
  set.seed(20110210) # Halle's birthday
  
  # required packages
  require("matrixStats"); require("dplyr"); require("tidyr"); require("magrittr")
  
  # number of simulations 
  n <- length(sim.data)
  
  # range for the distribution
  range <- 1:16
  
  # number of nodes
  v <- network::network.size(obs.data)
  
  # observed distribution, and restrict range
  obs <- summary(obs.data ~ dsp(0:v))
  
  # frequency
  obs_total <- sum(obs)
  obs <- obs/obs_total
  
  # data frame
  obs <- as.data.frame(obs); colnames(obs) <- c("obs")
  obs <- obs[range, ] # restrict range of the distribution
  obs <- as.data.frame(obs); colnames(obs) <- c("obs")
  obs <- dplyr::mutate(obs, stat = 0:15) # add column for statistic
  
  # store goodness-of-fit simulations statistics in list
  # sim <- data.frame(matrix(ncol = 0, nrow = n))
  sim <- c()
  
  # loop for the simulations
  for(i in 1:n){
    
    # simulated distribution, and restrict range
    run <- summary(sim.data[[i]] ~ dsp(0:v))
    
    # frequency
    total <- sum(run)
    run <- run/total
    
    # store the clustering coefficients in the results vector
    sim <- cbind(sim, run)
  }
  sim <- as.matrix(sim)

  # compute quantile range for plot
  
  # 50% quantile i.e., median
  median <- matrixStats::rowQuantiles(sim, probs = 0.50, na.rm = T) # median <- matrixStats::rowMedians(sim)
  
  # 25% quantile
  Q1 <- matrixStats::rowQuantiles(sim, probs = 0.25, na.rm = T)
  
  # 75% quantile
  Q3 <- matrixStats::rowQuantiles(sim, probs = 0.75, na.rm = T)
  
  # compute lower quantile for the error bars
  ci.lo <- matrixStats::rowQuantiles(sim, probs = 0.025, na.rm = T)
  
  # compute higher quantile for the error bars
  ci.hi <- matrixStats::rowQuantiles(sim, probs = 0.975, na.rm = T)
  
  # join
  stats <- as.data.frame(cbind(median, Q1, Q3, ci.lo, ci.hi))
  stats <- stats[range, ] # constrict range of the distribution
  colnames(stats) <- c("median", "Q1", "Q3", "ci.lo", "ci.hi")
  stats <- dplyr::mutate(stats, stat = 0:15) # add column for statistic
  
  # list
  data <- list(obs, stats)
  
  # data
  return(data)
}
gof_01.siren.stats.dsp <- gof_stats(
  obs.data = g_siren, 
  sim.data = g_siren.sim
  )
gof_02.togo.stats.dsp <- gof_stats(
  obs.data = g_togo, 
  sim.data = g_togo.sim
  )
gof_03.caviar.stats.dsp <- gof_stats(
  obs.data = g_caviar, 
  sim.data = g_caviar.sim
  )
gof_04.cielnet.stats.dsp <- gof_stats(
  obs.data = g_cielnet, 
  sim.data = g_cielnet.sim
  )
gof_05.cocaine.stats.dsp <- gof_stats(
  obs.data = g_cocaine,
  sim.data = g_cocaine.sim
  )
gof_06.heroin.stats.dsp <- gof_stats(
  obs.data = g_heroin, 
  sim.data = g_heroin.sim
  )
gof_07.oversize.stats.dsp <- gof_stats(
  obs.data = g_oversize, 
  sim.data = g_oversize.sim
  )
gof_08.montagna.stats.dsp <- gof_stats(
  obs.data = g_montagna, 
  sim.data = g_montagna.sim
  )




# construct goodness-of-fit plots ----------------------------------------------
gof_plot <- function(data, title){
  
  # set seed for replication
  set.seed(20190816) # Maeve's birthday
  
  # function to set the labels for the y-axis
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # observed data
  obs <- data[[1]]
  
  # simulated data
  sim <- data[[2]]
  
  # construct the boxplot
  fig <- ggplot2::ggplot(data = sim, ggplot2::aes(x = stat, y = median)) +
    # error bars to visualize the confidence intervals
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black", alpha = 1.0) +
    # boxplot
    ggplot2::geom_boxplot(
      data = sim,
      mapping = ggplot2::aes(lower = Q1, upper = Q3, middle = median, ymin = ci.lo, ymax = ci.hi, group = stat), 
      stat = "identity",
      fill = "white", color = "black", width = 0.5, alpha = 1.0
    ) +
    # line graph for the degree distribution
    ggplot2::geom_line(data = obs, ggplot2::aes(x = stat, y = obs), color = "red", linetype = "solid", size = 1, alpha = 0.8) +
    ggplot2::scale_y_continuous(limits = c(0.00, 1.00), breaks = seq(0, 1.00, 0.25), labels = scaleFUN) +
    # ggplot2::scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, 1)) +
    ggplot2::labs(
      title = title,
      y = NULL,
      x = NULL
    ) +
    ggthemes::theme_few() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 5, face = "plain"),
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      )
  plot(fig)
  return(fig)
}
gof_plot(data = gof_01.siren.stats.dsp, title = "(A) SIREN AUTO THEFT RING")
gof_plot(data = gof_02.togo.stats.dsp, title = "(B) TOGO AUTO THEFT RING")
gof_plot(data = gof_03.caviar.stats.dsp, title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION")
gof_plot(data = gof_04.cielnet.stats.dsp, title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION")
gof_plot(data = gof_05.cocaine.stats.dsp, title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT") # good
gof_plot(data = gof_06.heroin.stats.dsp, title = "(F) NEW YORK CITY HEROIN TRAFFICKERS")
gof_plot(data = gof_07.oversize.stats.dsp, title = "(G) 'NDRANGHETA - OPERATION OVERSIZE")
gof_plot(data = gof_08.montagna.stats.dsp, title = "(H) COSA NOSTRA - OPERATION MONTAGNA")




