#  ---------------------------------------------------------------------------------------------------

# file 04: test for sample bias in the criminal networks (Appendix 1)

# ----------------------------------------------------------------------------------------------------



# function to estimate the latent space model
lsm <- function(g){
  
  # required packages
  require('statnet'); require('igraph'); require('intergraph'); require('latentnet'); require('statip')
  
  # set seed for replication purposes
  seed <- set.seed(20190812) # Hayes' birthday
  
  # calculate the modal shortest path length
  d <- igraph::distances(intergraph::asIgraph(g))
  d <- as.vector(d) 
  d <- statip::mfv(d) # mode
  
  # latent space model, controlling for the squared euclidean distance, where the dimension set to the modal shortest path
  # https://cran.r-project.org/web/packages/latentnet/latentnet.pdf
  model <- latentnet::ergmm(
    g ~ intercept + euclidean2(d = d) + sociality(),
    family = "Bernoulli",
    # tofit = "mcmc",
    control = latentnet::ergmm.control(
      burnin = 100000,
      interval = 100
      ),
    seed = seed
    )
  
  # print model parameters
  print(summary(model))
  
  # return
  return(model)
}
lsm_siren <- lsm(g_siren)
lsm_togo <- lsm(g_togo)
lsm_caviar <- lsm(g_caviar)
lsm_cielnet <- lsm(g_cielnet)
lsm_cocaine <- lsm(g_cocaine)
lsm_heroin <- lsm(g_heroin)
lsm_oversize <- lsm(g_oversize)
lsm_montagna <- lsm(g_montagna)



# goodness-of-fit of the graph parameters
lsm_gof <- function(model){
  set.seed(20110210) # Halle's birthday
  gof <- gof( # goodness-of-fit function
    model,
    nsim = 100, 
    GOF = ~distance + degree + dspart + espart, 
    verbose = FALSE
    )
  return(gof) # return
}
lsm_gof.siren <- lsm_gof(lsm_siren)
lsm_gof.togo <- lsm_gof(lsm_togo)
lsm_gof.caviar <- lsm_gof(lsm_caviar)
lsm_gof.cielnet <- lsm_gof(lsm_cielnet)
lsm_gof.cocaine <- lsm_gof(lsm_cocaine)
lsm_gof.heroin <- lsm_gof(lsm_heroin)
lsm_gof.oversize <- lsm_gof(lsm_oversize)
lsm_gof.montagna <- lsm_gof(lsm_montagna)



# function to calculate the goodness-of-fit statistics from the latent space models 
gof_lsm <- function(obs.data, sim.data){
  
  # required packages
  require('dplyr'); require('matrixStats')
  
  # range for the distribution
  range <- 1:6 # six degrees of separation
  
  # observed distribution
  obs <- obs.data
  
  # calculate relative frequency
  obs_total <- sum(obs)
  obs <- obs/obs_total
  
  # constrict range of the distribution
  obs <- obs[range]
  
  # construct data frame for observed distribution
  obs <- as.data.frame(obs)
  colnames(obs) <- c("obs")
  
  # add column for statistic
  n <- nrow(obs)
  obs <- dplyr::mutate(obs, stat = 1:n) 
  
  # simulated distribution
  sim <- sim.data
  
  # calculate relative frequency
  sim_total <- rowSums(sim)
  sim <- sim/sim_total
  
  # compute quantiles for plot
  
  # 50% quantile i.e., median
  median <- matrixStats::colQuantiles(sim, probs = 0.50, na.rm = T) # median <- matrixStats::colMedians(sim)
  
  # 25% quantile
  Q1 <- matrixStats::colQuantiles(sim, probs = 0.25, na.rm = T)
  
  # 75% quantile
  Q3 <- matrixStats::colQuantiles(sim, probs = 0.75, na.rm = T)
  
  # compute lower quantile for the error bars
  ci.lo <- matrixStats::colQuantiles(sim, probs = 0.025, na.rm = T)
  
  # compute higher quantile for the error bars
  ci.hi <- matrixStats::colQuantiles(sim, probs = 0.975, na.rm = T)
  
  # join
  stats <- as.data.frame(cbind(median, Q1, Q3, ci.lo, ci.hi))
  stats <- stats[range, ] # constrict range of the distribution
  colnames(stats) <- c("median", "Q1", "Q3", "ci.lo", "ci.hi")
  
  # add column for statistic
  n <- nrow(stats)
  stats <- dplyr::mutate(stats, stat = 1:n)
  
  # list
  data <- list(obs, stats)
  
  # data
  return(data)
}

# distribution of path lengths
gof_01.siren.stats.dist <- gof_lsm(
  obs.data = lsm_gof.siren$obs.dist, 
  sim.data = lsm_gof.siren$sim.dist
  )
gof_02.togo.stats.dist <- gof_lsm(
  obs.data = lsm_gof.togo$obs.dist, 
  sim.data = lsm_gof.togo$sim.dist
  )
gof_03.caviar.stats.dist <- gof_lsm(
  obs.data = lsm_gof.caviar$obs.dist, 
  sim.data = lsm_gof.caviar$sim.dist
  )
gof_04.cielnet.stats.dist <- gof_lsm(
  obs.data = lsm_gof.cielnet$obs.dist, 
  sim.data = lsm_gof.cielnet$sim.dist
  )
gof_05.cocaine.stats.dist <- gof_lsm(
  obs.data = lsm_gof.cocaine$obs.dist, 
  sim.data = lsm_gof.cocaine$sim.dist
  )
gof_06.heroin.stats.dist <- gof_lsm(
  obs.data = lsm_gof.heroin$obs.dist,
  sim.data = lsm_gof.heroin$sim.dist
  )
gof_07.oversize.stats.dist <- gof_lsm(
  obs.data = lsm_gof.oversize$obs.dist, 
  sim.data = lsm_gof.oversize$sim.dist
  )
gof_08.montagna.stats.dist <- gof_lsm(
  obs.data = lsm_gof.montagna$obs.dist, 
  sim.data = lsm_gof.montagna$sim.dist
  )



# construct goodness-of-fit plots for the path length distribution i.e., degrees of separation
gof_plot <- function(data, title){
  
  # set seed for replication
  set.seed(20190816) # Maeve's birthday
  
  # required packages
  require('ggplot2'); require('ggthemes')
  
  # function to set the labels for the y-axis
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # observed data
  obs <- data[[1]]
  
  # simulated data
  sim <- data[[2]]
  
  # construct the boxplot
  fig <- ggplot2::ggplot(sim, ggplot2::aes(x = stat, y = median)) +
    # error bars to visualize the confidence intervals
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black", alpha = 1.0) +
    # boxplot
    ggplot2::geom_boxplot(
      mapping = ggplot2::aes(lower = Q1, upper = Q3, middle = median, ymin = ci.lo, ymax = ci.hi, group = stat), 
      stat = "identity",
      fill = "white", color = "black", width = 0.5, alpha = 1.0
      ) +
    # line graph for the degree distribution
    ggplot2::geom_line(data = obs, ggplot2::aes(x = stat, y = obs), color = "red", linetype = "solid", size = 1, alpha = 0.8) +
    ggplot2::scale_y_continuous(limits = c(0.00, 1.00), breaks = seq(0, 1.00, 0.25), labels = scaleFUN) +
    # ggplot2::scale_x_continuous(limits = c(1, 6), breaks = seq(1, 6)) +
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

# appendix fig 1 for the goodness-of-fit of the path length distribution i.e., degrees of separation
output(plot = gof_plot(
  data = gof_01.siren.stats.dist, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA1a.png"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.dist, 
  title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA1b.png"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.dist, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA1c.png"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.dist, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA1d.png"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.dist, 
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS"), 
  filename = "figA1e.png"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.dist, 
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT"), 
  filename = "figA1f.png"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.dist, 
  title = "(G) OVERSIZE - 'NDRANGHETA RACKETEERING"), 
  filename = "figA1g.png"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.dist, 
  title = "(H) MONTAGNA - COSA NOSTRA BID-RIGGING CONSPIRACY"), 
  filename = "figA1h.png"
  )





# calculate measures to test for spotlight effects i.e.,, sampling bias
spotlight <- function(model, g){
  
  # required packages
  require('CINNA'); require('intergraph'); require('igraph')

  # sociality effects
  effects <- model
  
  # normalize
  effects <- (effects - min(effects))/(max(effects) - min(effects))

  # don't run
  # normalized closeness centrality
  # c <- igraph::closeness(intergraph::asIgraph(g), mode = "total", normalize = TRUE)
  
  # normalized harmonic closeness centrality
  c <- CINNA::harmonic_centrality(intergraph::asIgraph(g), mode = "all")
  c <- (c - min(c))/(max(c) - min(c))
  
  # join
  data <- as.data.frame(cbind(effects, c))
  colnames(data) <- c("effects", "c")
  
  # return
  return(data)
}
spotlight_siren <- spotlight(model = lsm_siren$mkl$beta, g = g_siren)
spotlight_togo <- spotlight(model = lsm_togo$mkl$beta, g = g_togo)
spotlight_caviar <- spotlight(model = lsm_caviar$mkl$beta, g = g_caviar)
spotlight_cielnet <- spotlight(model = lsm_cielnet$mkl$beta, g = g_cielnet)
spotlight_cocaine <- spotlight(model = lsm_cocaine$mkl$beta, g = g_cocaine)
spotlight_heroin <- spotlight(model = lsm_heroin$mkl$beta, g = g_heroin)
spotlight_oversize <- spotlight(model = lsm_oversize$mkl$beta, g = g_oversize)
spotlight_montagna <- spotlight(model = lsm_montagna$mkl$beta, g = g_montagna)



# r-squared for locally weighted regression
# https://fibosworld.wordpress.com/2012/11/04/loess-regression-with-r/
rsq <- function(data){
  
  # required packages
  require('stats')
    
  # model  
  model <-  stats::loess(effects ~ c, span = 1, data = data)
 
  # r-squared for loess
  y.obs <- data$effects
  y.hat <- stats::predict(model)
  rsq <- stats::cor(y.obs, y.hat)^2
  
  # return
  return(rsq)
}
rsq(spotlight_siren)
rsq(spotlight_togo)
rsq(spotlight_caviar)
rsq(spotlight_cielnet)
rsq(spotlight_cocaine)
rsq(spotlight_heroin)
rsq(spotlight_oversize)
rsq(spotlight_montagna)




# function to set the labels for the y-axis
scaleFUN <- function(x) sprintf("%.2f", x)

# construct the plots for fig a2
spotlight_plot <- function(data, title,rsqFUN){
  
  # required packages
  require('ggplot2'); require('ggthemes')
  
  # plot
  ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = c, y = effects)) +
    ggplot2::geom_point(shape = 21, size = 1, color = "black", fill = "white") +
    ggplot2::geom_smooth(method = "loess", span = 1, formula = y ~ x, se = FALSE, color = "black", alpha = 0.80) +
    ggthemes::theme_few() +
    ggplot2::scale_y_continuous(limits = c(0.00, 1.00), breaks = seq(0.00, 1.00, 0.25), labels = scaleFUN) +
    ggplot2::scale_x_continuous(limits = c(0.00, 1.00), breaks = seq(0.00, 1.00, 0.25), labels = scaleFUN) +
    ggplot2::labs(
      title = title,
      y = NULL,
      x = NULL,
      ) +
    ggplot2::annotate(
      "text", x = 0.25, y = 0.95, 
      label = bquote(italic(R)^2 == .(round(rsqFUN, digits = 2))), 
      hjust = 1, 
      size = 3
      ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 5, face = "plain", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 5),
      )
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

# appendix fig 2
output(plot = spotlight_plot(
  spotlight_siren, 
  title = "(A) SIREN AUTO THEFT RING", 
  rsqFUN = rsq(spotlight_siren)
  ), 
  filename = "figA2a.png"
  )
output(plot = spotlight_plot(
  spotlight_togo, 
  title = "(B) TOGO AUTO THEFT RING",
  rsqFUN = rsq(spotlight_togo)
  ), 
  filename = "figA2b.png"
  )
output(plot = spotlight_plot(
  spotlight_caviar, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION",
  rsqFUN = rsq(spotlight_caviar)
  ), 
  filename = "figA2c.png"
  )
output(plot = spotlight_plot(
  spotlight_cielnet, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION",
  rsqFUN = rsq(spotlight_cielnet)
  ), 
  filename = "figA2d.png"
  )
output(plot = spotlight_plot(
  spotlight_cocaine, 
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS",
  rsqFUN = rsq(spotlight_cocaine)
  ), 
  filename = "figA2e.png"
  )
output(plot = spotlight_plot(
  spotlight_heroin, 
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT",
  rsqFUN = rsq(spotlight_heroin)
  ), 
  filename = "figA2f.png"
  )
output(plot = spotlight_plot(
  spotlight_oversize, 
  title = "(G) OVERSIZE - 'NDRANGHETA RACKETEERING",
  rsqFUN = rsq(spotlight_oversize)
  ), 
  filename = "figA2g.png"
  )
output(plot = spotlight_plot(
  spotlight_montagna, 
  title = "(H) MONTAGNA - COSA NOSTRA BID-RIGGING CONSPIRACY",
  rsqFUN = rsq(spotlight_montagna)
  ), 
  filename = "figA2h.png"
  )



# close .R file


