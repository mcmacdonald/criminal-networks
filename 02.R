#  -----------------------------------------------------------------------------------

# file 02: estimate the Bayesian exponential random graph models

# note: you must run this file before you run file named '03.R'

# ------------------------------------------------------------------------------------



# posterior parameter estimation -----------------------------------------------
bayes <- function(y, formula){
  
  # load 'bergm' package
  require('Bergm'); require('network')
  
  # inputs
  i <- network::network.size(y) * 2 # burn in iterations to begin the MCMC run 
  k <-  1000 # sample iterations
  h <- 5 * 2 # chains in the posterior distribution ... approximately twice the number of model parameters 
  n <- h * k # per Caimo & Friel (2011), auxiliary chain = # chains (h) * sample iterations (k)
  
  # set seed for replication
  set.seed(20110210) # Halle's birthday
  
  # estimate the model
  bayes <- Bergm::bergm(
    formula = formula, # equation
    burn.in = i, 
    main.iters = k, 
    nchains = h,
    aux.iters = n,
    gamma = 0.1 # empirically, gamma ranges 0.1 - 1.5, where smaller gamma leads to better acceptance in big graphs
    )
  return(bayes)
}



# goodness-of-fit --------------------------------------------------------------
GOF <- function(bayes){
  
  # required packages
  require("Bergm")
  
  # set seed for replication
  set.seed(20110211) # Halle's birthday
  
  # inputs
  n <- 100 # graph simulations
  i <- 15 # degree distribution range
  j <- 15 # geodisstance range 
  k <- 15 # edgewise-shared partners range
  
  # goodness-of-fit
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
  formula = g_siren ~ edges + 
    gwdegree(decay = 0.7, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_01.siren)

# goodness-of-fit
gof_01.siren <- GOF(bayes_01.siren)



# estimate the model -----------------------------------------------------------
bayes_02.togo <- bayes(
  y = g_togo,
  formula = g_togo ~ edges + 
    gwdegree(decay = 0.7, fixed = TRUE) +
    gwdsp(decay = 1.7, fixed = TRUE) +
    gwesp(decay = 1.0, fixed = TRUE) +
    degcor
    )
summary(bayes_02.togo)

# goodness-of-fit
gof_02.togo <- GOF(bayes_02.togo)



# estimate the model -----------------------------------------------------------
bayes_03.caviar <- bayes(
  y = g_caviar,
  formula = g_caviar ~ edges + 
    gwdegree(decay = 2.0, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_03.caviar)

# goodness-of-fit
gof_03.caviar <- GOF(bayes_03.caviar)



# estimate the model -----------------------------------------------------------
bayes_04.cielnet <- bayes(
  y = g_cielnet,
  formula = g_cielnet ~ edges + 
    gwdegree(decay = 1.7, fixed = TRUE) +
    gwdsp(decay = 0.2, fixed = TRUE) +
    gwesp(decay = 0.2, fixed = TRUE) +
    degcor
    )
summary(bayes_04.cielnet)

# goodness-of-fit
gof_04.cielnet <- GOF(bayes_04.cielnet)



# estimate the model -----------------------------------------------------------
bayes_05.cocaine <- bayes(
  y = g_cocaine,
  formula = g_cocaine ~ edges + 
    gwdegree(decay = 2.5, fixed = TRUE) +
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_05.cocaine)

# goodness-of-fit
gof_05.cocaine <- GOF(bayes_05.cocaine)



# estimate the model -----------------------------------------------------------
bayes_06.heroin <- bayes(
  y = g_heroin,
  formula = g_heroin ~ edges + 
    gwdegree(decay = 0.5, fixed = TRUE) +
    gwdsp(decay = 0.1, fixed = TRUE) +
    gwesp(decay = 1.0, fixed = TRUE) +
    degcor
   )
summary(bayes_06.heroin)

# goodness-of-fit
gof_06.heroin <- GOF(bayes_06.heroin)



# estimate the model -----------------------------------------------------------
bayes_07.oversize <- bayes(
  y = g_oversize,
  formula = g_oversize ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) +
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 1.7, fixed = TRUE) +
    degcor
    )
summary(bayes_07.oversize)

# goodness-of-fit
gof_07.oversize <- GOF(bayes_07.oversize)



# estimate the model -----------------------------------------------------------
bayes_08.montagna <- bayes(
  y = g_montagna,
  formula = g_montagna ~ edges + 
    gwdegree(decay = 1.2, fixed = TRUE) +
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_08.montagna)

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
    ggplot2::geom_line(data = obs, ggplot2::aes(x = stat, y = obs), color = "red", linetype = "solid", linewidth = 1, alpha = 0.8) +
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
      plot.title = ggplot2::element_text(size = 5, face = "plain", hjust = 0.5),
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
    device = 'png', 
    dpi = 700
    )
}

# appendix Fig 5 degree distribution
output(plot = gof_plot(
  data = gof_01.siren.stats.deg, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA3a.png"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.deg, 
  title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA3b.png"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.deg, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA3c.png"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.deg, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA3d.png"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.deg, 
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS"), 
  filename = "figA3e.png"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.deg, 
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT"), 
  filename = "figA3f.png"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.deg, 
  title = "(G) OVERSIZE - 'NDRANGHETA RACKETEERING"), 
  filename = "figA3g.png"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.deg, 
  title = "(H) MONTAGNA - COSA NOSTRA BID-RIGGING CONSPIRACY"), 
  filename = "figA3h.png"
  )

# appendix Fig 6 edgewise shared partner distribution
output(plot = gof_plot(
  data = gof_01.siren.stats.esp, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA4a.png"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.esp, title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA4b.png"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.esp, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA4c.png"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.esp, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA4d.png"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.esp, 
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS"), 
  filename = "figA4e.png"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.esp, 
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT"), 
  filename = "figA4f.png"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.esp, 
  title = "(G) OVERSIZE - 'NDRANGHETA RACKETEERING"), 
  filename = "figA4g.png"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.esp, 
  title = "(H) MONTAGNA - COSA NOSTRA BID-RIGGING CONSPIRACY"), 
  filename = "figA4h.png"
  )



# close .R script




