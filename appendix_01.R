#  ---------------------------------------------------------------------------------------------------

# file 03: test for sample bias in the flows of communication

# note: you must run this file before you run file named '04.R'

# last updated: 27/02/2025

# ----------------------------------------------------------------------------------------------------




# function to generate the clustering coefficients of simulated graphs from the latent space model
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
    nsim = 10000, 
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



# function to calculate the goodness-of-fit statistics from the models 
gof_lsm <- function(obs.data, sim.data){
  
  # required packages
  
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

# appendix fig 1 for the goodness-of-fit of the path length distribution i.e., degrees of separation
output(plot = gof_plot(
  data = gof_01.siren.stats.dist, 
  title = "(A) SIREN AUTO THEFT RING"), 
  filename = "figA1a.pdf"
  )
output(plot = gof_plot(
  data = gof_02.togo.stats.dist, 
  title = "(B) TOGO AUTO THEFT RING"), 
  filename = "figA1b.pdf"
  )
output(plot = gof_plot(
  data = gof_03.caviar.stats.dist, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA1c.pdf"
  )
output(plot = gof_plot(
  data = gof_04.cielnet.stats.dist, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"), 
  filename = "figA1d.pdf"
  )
output(plot = gof_plot(
  data = gof_05.cocaine.stats.dist, 
  title = "(E) NEW YORK CITY COCAINE TRAFFICKERS"), 
  filename = "figA1e.pdf"
  )
output(plot = gof_plot(
  data = gof_06.heroin.stats.dist, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"), 
  filename = "figA1f.pdf"
  )
output(plot = gof_plot(
  data = gof_07.oversize.stats.dist, 
  title = "(G) 'NDRANGHETA - OPERATION OVERSIZE"), 
  filename = "figA1g.pdf"
  )
output(plot = gof_plot(
  data = gof_08.montagna.stats.dist, 
  title = "(H) COSA NOSTRA - OPERATION MONTAGNA"), 
  filename = "figA1h.pdf"
  )





# calcualte measures to test for spotlight effects
spotlight <- function(model, g){

  # sociality effects
  effects <- model
  
  # normalize
  effects <- (effects - min(effects))/(max(effects) - min(effects))

  # normalized closeness centrality
  # don't run
  # c <- igraph::closeness(intergraph::asIgraph(g), mode = "total", normalize = TRUE)
  
  # normalized harmonic closeness centrality
  c <- CINNA::harmonic_centrality(intergraph::asIgraph(g), mode = "all")
  c <- (c - min(c))/(max(c) - min(c))
  
  # join
  data <- as.data.frame(cbind(effects, c))
  colnames(data) <- c("effects", "c")
  
  # average the spotlight (sociality) effects by closeness centrality
  # data <- data %>% dplyr::group_by(c) %>% summarise(effects = mean(effects))
  
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
  require(stats)
    
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
      plot.title = ggplot2::element_text(size = 5, face = "plain"),
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
    device = 'pdf', 
    dpi = 700
    )
}

# appendix fig 2
output(plot = spotlight_plot(
  spotlight_siren, 
  title = "(A) SIREN AUTO THEFT RING", 
  rsqFUN = rsq(spotlight_siren)
  ), 
  filename = "figA2a.pdf"
  )
output(plot = spotlight_plot(
  spotlight_togo, 
  title = "(B) TOGO AUTO THEFT RING",
  rsqFUN = rsq(spotlight_togo)
  ), 
  filename = "figA2b.pdf"
  )
output(plot = spotlight_plot(
  spotlight_caviar, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION",
  rsqFUN = rsq(spotlight_caviar)
  ), 
  filename = "figA2c.pdf"
  )
output(plot = spotlight_plot(
  spotlight_cielnet, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION",
  rsqFUN = rsq(spotlight_cielnet)
  ), 
  filename = "figA2d.pdf"
  )
output(plot = spotlight_plot(
  spotlight_cocaine, 
  title = "(E) NEW YORK CITY COCAINE TRAFFICKERS",
  rsqFUN = rsq(spotlight_cocaine)
  ), 
  filename = "figA2e.pdf"
  )
output(plot = spotlight_plot(
  spotlight_heroin, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS",
  rsqFUN = rsq(spotlight_heroin)
  ), 
  filename = "figA2f.pdf"
  )
output(plot = spotlight_plot(
  spotlight_oversize, 
  title = "(G) 'NDRANGHETA - OPERATION OVERSIZE",
  rsqFUN = rsq(spotlight_oversize)
  ), 
  filename = "figA2g.pdf"
  )
output(plot = spotlight_plot(
  spotlight_montagna, 
  title = "(H) COSA NOSTRA - OPERATION MONTAGNA",
  rsqFUN = rsq(spotlight_montagna)
  ), 
  filename = "figA2h.pdf"
  )




####### extra



# function to calculate the cut-off for the number of co-conspirators 
deg <- function(g){
  
  # degree
  d <- sna::degree(g_caviar, gmode = "graph")
  
  # max value to cut-off the distribution
  cut <- max(d)
  return(cut)
}



# function to calculate the cut-off for the number of common co-conspirators
dsp <- function(g){
  
  # adjacency graph
  adj <- network::as.matrix.network(g, matrix.type = "adjacency")
  
  # shared partners for each dyad
  dsp <- adj %*% adj # matrix multiplication
  
  # lower triangle, excluding self-pairs
  dsp <- dsp[lower.tri(dsp)]
  
  # max value to cut-off the distribution
  cut <- max(dsp)
  return(cut)
}



# function to calculate the cut-off for the number of edgewise shared parners
esp <- function(g){
  
  # get all edges for all nodes 
  edges <- network::as.matrix.network.edgelist(g)
  
  # function to calcualte neighbors
  neighbors <- function(i){
    which(g[i, ] == 1)
  }
  
  # shared partners for each edge
  esp <- apply(edges, 1, function(dyad) {
    i <- dyad[1]
    j <- dyad[2]
    length(dplyr::intersect(neighbors(i), neighbors(j)))
    })
  
  # max value to cut-off the distribution
  cut <- max(esp)
  return(cut)
}
esp(g_caviar)



# function to calculate the cut-off for the path length or geodesic distance
dis <- function(g){
  # calculate the modal shortest path length
  d <- igraph::distances(intergraph::asIgraph(g))
  d <- as.vector(d)
  
  # ignore infinite values
  d[is.infinite(d)] <- 0 
  cut <- max(d)
  return(cut)
}



# model parameters
results <- function(g, model){
  
  # model inside function
  model <- model
  
  # posterior sample matrix
  b <- model$sample$beta
  b <- apply(b, 2, function(x) { # calculate coefficients and 95% credible intervals
    c(
      mu = mean(x),
      lo = stats::quantile(x, 0.025),
      hi = stats::quantile(x, 0.975)
    )
  }
    )
  b <- as.data.frame(t(b)) # transpose and transform into data frame
  
  # coefficient names
  terms <- lsm_siren$model$coef.names

  # centrality measures
  degree <- sna::degree(g, gmode = "graph"); eigen <- sna::evcent(g_siren, gmode = "graph")
  
  # drop reference 
  degree <-degree[-1]; eigen <- eigen[-1]
  
  # include 0 for intercept
  degree <- c(0, degree); eigen <- c(0, eigen)
  
  # return
  b <- cbind(terms, b, degree, eigen)
  colnames(b) <- c("term", "mu", "lo", "hi", "freeman", "eigen")
  return(b)
}
results(g = g_siren, model = lsm_siren)
results(g = g_togo, model = lsm_togo)
results(g = g_caviar, model = lsm_caviar)
results(g = g_cielnet, model = lsm_cielnet)
results(g = g_cocaine, model = lsm_cocaine)
results(g = g_heroin, model = lsm_heroin)
results(g = g_oversize, model = lsm_oversize)
results(g = g_montagna, model = lsm_montagna)



# function to transform the goodness-of-fit statistics to data frame
gof2df <- function(gof, cutoff){
  
  # transform into data frame
  gof <- as.data.frame(gof)
  
  # statistic as column
  gof$bin <- 0 : (nrow(gof) - 1) # nrow() counts the column names as a row, so subtract it from the calculation
  
  # rename columns
  colnames(gof) <- c("n", "lo", "mu", "hi", "p", "bin")
  
  # drop all test statistics outside the maximum value
  gof <- dplyr::filter(gof, bin <= cutoff)
  
  # normalize x-axis for plotting
  gof$bin <- gof$bin/cutoff
  
  # recale y-axis for plotting
  d <- max(gof[, 1:4]) # largest statistic in the observed and simulated values
  
  # normalize y-axis
  gof$n  <- gof$n/d
  gof$mu <- gof$mu/d
  gof$lo <- gof$lo/d
  gof$hi <- gof$hi/d
  
  # return data frame
  return(gof)
}



# function to plot goodness-of-fit statistics --------------------------------------------------------
gof_plot <- function(stat, xlab, group, title){

# plot data  
ggplot2::ggplot(stat, ggplot2::aes(x = bin, y = mu, group = 1)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey", alpha = 0.2) +
  ggplot2::geom_line(color = "grey", alpha = 0.8, size = 1.5) +
  # ggplot2::geom_point(ggplot2::aes(y = n), color = "red", size = 2) +
  ggplot2::geom_segment(ggplot2::aes(x = bin, xend = bin, y = 0, yend = n), color = "black", size = 1) +
  # ggplot2::geom_col(ggplot2::aes(x = bin, y = n), color = "black", fill = "white", size = 1) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = scales::label_number(accuracy = 0.01)  # optional: format to 2 decimal places
    ) +
  ggthemes::theme_clean() + 
  ggplot2::labs(
    title = paste(group),
    subtitle = paste(title),
    x = xlab,
    y = "NORMALIZED FREQUENCY"
    ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.50),
    plot.subtitle = ggplot2::element_text(size = 10, face = "plain")
    )
}
<- gof_plot(
  stat = gof2df(gof = lsm_gof.siren$summary.dist, cutoff = dis(g_siren)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "SIREN AUTO THEFT RING",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.togo$summary.dist, cutoff = dis(g_togo)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "TOGO AUTO THEFT RING",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.caviar$summary.dist, cutoff = dis(g_caviar)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "CAVIAR DRUG TRAFFICKING ORGANIZATION",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.cielnet$summary.dist, cutoff = dis(g_cielnet)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "CIELNET DRUG TRAFFICKING ORGANIZATION",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.cocaine$summary.dist, cutoff = dis(g_cocaine)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "NATARAJAN COCAINE TRAFFICKING ORGANIZATION",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.heroin$summary.dist, cutoff = dis(g_heroin)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "NATARAJAN HEROIN TRAFFICKING ORGANIZATION",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.oversize$summary.dist, cutoff = dis(g_oversize)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "N'DRANGHETA OVERSIZE DRUG TRAFFICKING ORGANIZATION",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )
gof_plot(
  stat = gof2df(gof = lsm_gof.montagna$summary.dist, cutoff = dis(g_montagna)), 
  xlab = "NORMALIZED DEGREES OF SEPARTION", 
  group = "MONTAGNA BID-RIGGING CONSPIRACY, COSA NOSTRA",  
  title = "GOODNESS-OF-FIT FOR DEGREES OF SEPARATION"
  )





# function to plot goodness-of-fit statistics --------------------------------------------------------
gof_plot <- function(stat, xlab, group, title){
  
  # plot data  
  ggplot2::ggplot(stat, ggplot2::aes(x = bin, y = mu, group = 1)) +
    # observed distribution
    ggplot2::geom_segment(ggplot2::aes(x = bin, xend = bin, y = 0, yend = n), color = "black", linewidth = 1.5, size = 1, alpha = 0.5) +
    # simulated distributiona
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey", alpha = 0.2) +
    ggplot2::geom_line(color = "black", alpha = 0.8, size = 1) +
     # scale x-axis
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      labels = scales::label_number(accuracy = 0.01)  # optional: format to 2 decimal places
      ) +
    # scale y-axis
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0.25, 0.50, 0.75, 1.00),
      labels = scales::label_number(accuracy = 0.01)  # optional: format to 2 decimal places
      ) +
    # plot theme
    ggthemes::theme_clean() + 
    # axis labels
    ggplot2::labs(
      title = paste(group),
      subtitle = paste(title),
      x = xlab,
      y = "NORMALIZED FREQUENCY"
      ) +
    # edit the text size
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 5, face = "bold", hjust = 0.50),
      plot.subtitle = ggplot2::element_text(size = 5, face = "plain"),
      axis.text.x = ggplot2::element_text(size = 5),
      axis.text.y = ggplot2::element_text(size = 5),
      axis.title.x = ggplot2::element_text(size = 5),
      axis.title.y = ggplot2::element_text(size = 5)
      )
}

# Siren auto theft ring
figA2a <- gof_plot(
  stat = gof2df(gof = lsm_gof.siren$summary.deg, cutoff = deg(g_siren)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2b <- gof_plot(
  stat = gof2df(gof = lsm_gof.siren$summary.dspart, cutoff = dsp(g_siren)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "SIREN AUTO THEFT RING GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2c <- gof_plot(
  stat = gof2df(gof = lsm_gof.siren$summary.espart, cutoff = esp(g_siren)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Togo auto theft ring
figA2d <- gof_plot(
  stat = gof2df(gof = lsm_gof.togo$summary.deg, cutoff = deg(g_togo)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2e <- gof_plot(
  stat = gof2df(gof = lsm_gof.togo$summary.dspart, cutoff = dsp(g_togo)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "TOGO AUTO THEFT RING GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2f <- gof_plot(
  stat = gof2df(gof = lsm_gof.togo$summary.espart, cutoff = esp(g_togo)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Caviar drug trafficking organization
figA2g <- gof_plot(
  stat = gof2df(gof = lsm_gof.caviar$summary.deg, cutoff = deg(g_caviar)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2h <- gof_plot(
  stat = gof2df(gof = lsm_gof.caviar$summary.dspart, cutoff = dsp(g_caviar)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "CAVIAR DRUG TRAFFICKING ORGANIZATION GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2i <- gof_plot(
  stat = gof2df(gof = lsm_gof.caviar$summary.espart, cutoff = esp(g_caviar)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Cielnet drug trafficking organization
figA2j <- gof_plot(
  stat = gof2df(gof = lsm_gof.cielnet$summary.deg, cutoff = deg(g_cielnet)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2k <- gof_plot(
  stat = gof2df(gof = lsm_gof.cielnet$summary.dspart, cutoff = dsp(g_cielnet)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "CIELNET DRUG TRAFFICKING ORGANIZATION GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2l <- gof_plot(
  stat = gof2df(gof = lsm_gof.cielnet$summary.espart, cutoff = esp(g_cielnet)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Natarajan cocaine drug trafficking organization
figA2m <- gof_plot(
  stat = gof2df(gof = lsm_gof.cocaine$summary.deg, cutoff = deg(g_cocaine)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2n <- gof_plot(
  stat = gof2df(gof = lsm_gof.cocaine$summary.dspart, cutoff = dsp(g_cocaine)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "NATARAJAN COCAINE TRAFFICKING ORGANIZATION GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2o <- gof_plot(
  stat = gof2df(gof = lsm_gof.cocaine$summary.espart, cutoff = esp(g_cocaine)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Natarajan heroin drug trafficking organization
figA2p <- gof_plot(
  stat = gof2df(gof = lsm_gof.heroin$summary.deg, cutoff = deg(g_heroin)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2q <- gof_plot(
  stat = gof2df(gof = lsm_gof.heroin$summary.dspart, cutoff = dsp(g_heroin)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "NATARAJAN HEROIN TRAFFICKING ORGANIZATION GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2r <- gof_plot(
  stat = gof2df(gof = lsm_gof.heroin$summary.espart, cutoff = esp(g_heroin)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# oversize network,'Ndrangheta
figA2s <- gof_plot(
  stat = gof2df(gof = lsm_gof.oversize$summary.deg, cutoff = deg(g_oversize)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2t <- gof_plot(
  stat = gof2df(gof = lsm_gof.oversize$summary.dspart, cutoff = dsp(g_oversize)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "'NDRANGHETA OVERSIZE NETWORK GOF STATISITICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2u <- gof_plot(
  stat = gof2df(gof = lsm_gof.oversize$summary.espart, cutoff = esp(g_oversize)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )

# Montagna bid-rigging conspiracy, Sicily
figA2v <- gof_plot(
  stat = gof2df(gof = lsm_gof.montagna$summary.deg, cutoff = deg(g_montagna)), 
  xlab = "NUMBER OF CO-CONSPIRATORS", 
  group = "",  
  title = "(A) DEGREE DISTRIBUTION"
  )
figA2w <- gof_plot(
  stat = gof2df(gof = lsm_gof.montagna$summary.dspart, cutoff = dsp(g_montagna)), 
  xlab = "NUMBER OF COMMON CO-CONSPIRATORS", 
  group = "MONTAGNA BID-RIGGING CONSPIRACY GOF STATISTICS",  
  title = "(B) DYADWISE SHARED PARTNERS"
  )
figA2x <- gof_plot(
  stat = gof2df(gof = lsm_gof.montagna$summary.espart, cutoff = esp(g_montagna)), 
  xlab = "NUMBER OF CLIQUES CONSPIRATORS BELONG TO", 
  group = "",  
  title = "(C) EDGEWISE SHARED PARTNERS"
  )



# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 2.5, 
    height = 2.5, 
    device = 'png', 
    dpi = 700
    )
}
output(plot = figA2a, filename = "figA2a.png")
output(plot = figA2b, filename = "figA2b.png")
output(plot = figA2c, filename = "figA2c.png")
output(plot = figA2d, filename = "figA2d.png")
output(plot = figA2e, filename = "figA2e.png")
output(plot = figA2f, filename = "figA2f.png")
output(plot = figA2g, filename = "figA2g.png")
output(plot = figA2h, filename = "figA2h.png")
output(plot = figA2i, filename = "figA2i.png")
output(plot = figA2j, filename = "figA2j.png")
output(plot = figA2k, filename = "figA2k.png")
output(plot = figA2l, filename = "figA2l.png")
output(plot = figA2m, filename = "figA2m.png")
output(plot = figA2n, filename = "figA2n.png")
output(plot = figA2o, filename = "figA2o.png")
output(plot = figA2p, filename = "figA2p.png")
output(plot = figA2q, filename = "figA2q.png")
output(plot = figA2r, filename = "figA2r.png")
output(plot = figA2s, filename = "figA2s.png")
output(plot = figA2t, filename = "figA2t.png")
output(plot = figA2u, filename = "figA2u.png")
output(plot = figA2v, filename = "figA2v.png")
output(plot = figA2w, filename = "figA2w.png")
output(plot = figA2x, filename = "figA2x.png")


