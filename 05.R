#  -----------------------------------------------------------------------------------

# file 04: simulate the clustering coefficients from the Bayesian ERGMs

# last updated: 15/04/2025

# ------------------------------------------------------------------------------------



# calculate the clustering coefficient of each of the criminal networks
clustering_obs <- function(g){
  
  # required packages
  require("statnet"); require("igraph"); require("intergraph")
  
  # clustering coefficient
  cc <- igraph::transitivity(
    intergraph::asIgraph(g),
    type = "localaverageundirected",
    isolates = "zero"
    )
  
  # return
  return(cc)
}
cc_siren <- clustering_obs(g_siren)
cc_togo <- clustering_obs(g_togo)
cc_caviar <- clustering_obs(g_caviar)
cc_cielnet <- clustering_obs(g_cielnet)
cc_cocaine <- clustering_obs(g_cocaine)
cc_heroin <- clustering_obs(g_heroin)
cc_oversize <- clustering_obs(g_oversize)
cc_montagna <- clustering_obs(g_montagna)



# distribution of clustering coefficients from the simulations
clustering_sim <- function(simulations){
  
  # required packages
  require("statnet"); require("igraph"); require("intergraph")

  # number of simulations
  samples <- length(simulations)
  
  # store results in list
  results <- c()
  
  # loop for the simulations
  for(i in 1:samples){
    
    # calculate the clustering coefficient for each of the simulations
    cc <- igraph::transitivity(
      graph = intergraph::asIgraph(simulations[[i]]),
      type = "localaverageundirected",
      isolates = "zero"
      )
    
    # store the clustering coefficients in the results vector
    results <- c(results, cc)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("cc")
  
  # return
  return(results)
}
cc_siren.sim <- clustering_sim(g_siren.sim)
cc_togo.sim <- clustering_sim(g_togo.sim)
cc_caviar.sim <- clustering_sim(g_caviar.sim)
cc_cielnet.sim <- clustering_sim(g_cielnet.sim)
cc_cocaine.sim <- clustering_sim(g_cocaine.sim)
cc_heroin.sim <- clustering_sim(g_heroin.sim)
cc_oversize.sim <- clustering_sim(g_oversize.sim)
cc_montagna.sim <- clustering_sim(g_montagna.sim)



# figure 3. plot histogram of the clustering coefficients for the real and simulated graphs
plot_clustering <- function(obs, sim, title){
  
  # required packages
  require('ggplot2'); require('scales'); library('ggplot2')
  
  # mean and standard deviation of the latent space models
  mu <- mean(sim$cc); sd <- sd(sim$cc)
  
  # compute z-score
  z <- (obs - mu)/sd
  
  # two-tailed significance test
  p <- 2 * (1 - stats::pnorm(abs(z)))
  
  # flag statistical significance
  sig <- if (p < 0.001) {"***"} 
  else{
    if (p < 0.1) {"**"}
    else{
      if (p < 0.5){"*"}
      else{
        # unicode
        # "\u2099\u00B7\u209B\u00B7"
        " n.s."
      }
    }
  }
  
  # label to annotate the mean coefficient of the simulated networks
  label1 <- mean(sim$cc); label1 <- round(label1, digits = 2)
  
  # label to annotate the coefficient for the criminal networks
  label2 <- round(obs, digits = 2)
  
  # distribution of correlation coefficients
  histogram <- ggplot2::ggplot(sim, ggplot2::aes( x = cc ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 20, color = "black", fill = "white") +
    # line markers for the clustering coefficients
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(cc)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = obs), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
    # annotation for the mean clustering coefficient for the random networks
    ggplot2::annotate(geom = "label", x = label1 - 0.02, y = 0.35, label = sprintf('%0.2f', label1), size = 2) +
    # annotation for the clustering coefficient for the actual networks
    ggplot2::annotate(geom = "label", x = label2 + 0.02, y = 0.30, label = paste0(sprintf('%0.2f', label2), sig), size = 2) +
    # transform y-axis to percentage scale
    ggplot2::aes(y = after_stat(count)/sum(after_stat(count))) + 
    ggplot2::scale_y_continuous(
      name = "PROBABILITY DENSITY FUNCTION (PDF)", 
      labels = scales::percent_format(accuracy = 1L), # 2L to round to one decimal place, 3L to round to two decimal places, etc.
      limits = c(0.00, 0.40), 
      breaks = c(0.00, 0.10, 0.20, 0.30, 0.40)
      ) +
    # ggplot2::scale_x_continuous(
      # name = "CLUSTERING COEFFICIENT",
      # labels = scales::label_number(accuracy = 0.01),
      # limits = c(-0.001, 1.00),
      # breaks = c(-0.001, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
      # ) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "CLUSTERING COEFFICIENT") +
    ggplot2::ggtitle(title) +
    ggthemes::theme_few() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 5, face = "plain"),
      axis.text.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain"),
      axis.text.y = ggplot2::element_text(color = "black", size = 5, hjust = 1.0, vjust = 0.0, face = "plain"),  
      axis.title.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.0, face = "plain"),
      axis.title.y = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain")
      )
  plot(histogram)
  return(histogram)
}
clustering_siren <- plot_clustering(
  obs = cc_siren,
  sim = cc_siren.sim, 
  title = "(A) SIREN AUTO THEFT RING"
  )
clustering_togo <- plot_clustering(
  obs = cc_togo, 
  sim = cc_togo.sim, 
  title = "(B) TOGO AUTO THEFT RING"
  )
clustering_caviar <- plot_clustering(
  obs = cc_caviar, 
  sim = cc_caviar.sim, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
  )
clustering_cielnet <- plot_clustering(
  obs = cc_cielnet, 
  sim = cc_cielnet.sim, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
  )
clustering_cocaine <- plot_clustering(
  obs = cc_cocaine, 
  sim = cc_cocaine.sim, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
  )
clustering_heroin <- plot_clustering(
  obs = cc_heroin, 
  sim = cc_heroin.sim, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
  )
clustering_oversize <- plot_clustering(
  obs = cc_oversize, 
  sim = cc_oversize.sim, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
  )
clustering_montagna <- plot_clustering(
  obs = cc_montagna, 
  sim = cc_montagna.sim, 
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
  )



# output high resolution images
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
output(plot = clustering_siren, filename = "fig3a.pdf")
output(plot = clustering_togo, filename = "fig3b.pdf")
output(plot = clustering_caviar, filename = "fig3c.pdf")
output(plot = clustering_cielnet, filename = "fig3d.pdf")
output(plot = clustering_cocaine, filename = "fig3e.pdf")
output(plot = clustering_heroin, filename = "fig3f.pdf")
output(plot = clustering_oversize, filename = "fig3g.pdf")
output(plot = clustering_montagna, filename = "fig3h.pdf")



# close .r script



#  -----------------------------------------------------------------------------------------------------------------------------

# file 05: estimate the clustering coefficient of each graph and compare them to random graphs

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------------------------------------------------



# calculate the clustering coefficient of each of the criminal networks ------------------------------
clustering <- function(g){
  # required packages
  require("statnet"); require("igraph"); require("intergraph")
  cc <- igraph::transitivity(
    intergraph::asIgraph(g),
    type = "localaverageundirected",
    isolates = "zero"
  )
  print(cc)
}
cc_siren <- clustering(g_siren)
cc_togo <- clustering(g_togo)
cc_caviar <- clustering(g_caviar)
cc_cielnet <- clustering(g_cielnet)
cc_cocaine <- clustering(g_cocaine)
cc_heroin <- clustering(g_heroin)
cc_oversize <- clustering(g_oversize)
cc_montagna <- clustering(g_montagna)
cc_tfc <- clustering(g_tfc)
cc_super <- clustering(g_super)





# function to generate the clustering coefficients of simulated graphs from the latent space model
simulator <- function(g){
  
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
    g ~ intercept + euclidean2(d = d),
    family = "Bernoulli",
    control = ergmm.control(
      burnin = 100000,
      interval = 100
    ),
    seed = seed
  )
  
  # goodness-of-fit of the geodesic distances or path lengths
  gof <- gof(
    model,
    nsim = 10000, 
    GOF = ~distance, 
    verbose = FALSE
  )
  
  # print goodness-of-fit
  plot(gof); print(gof)
  
  # simulate large number of networks from the latent space model
  simulations <- simulate(model, nsim = 10000, seed = seed)
  simulations <- simulations$networks # retain the simulated networks
  
  # number of simulations
  samples <- length(simulations)
  
  # store results in list
  results <- c()
  
  # loop for the simulations
  for(i in 1:samples){
    
    # calculate the clustering coefficient for each of the simulations
    cc <- igraph::transitivity(
      graph = intergraph::asIgraph(simulations[[i]]),
      type = "localaverageundirected",
      isolates = "zero"
    )
    
    # store the clustering coefficients in the results vector
    results <- c(results, cc)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("cc")
  
  # return
  return(results)
}
cc_siren.lsm <- simulator(g_siren)
cc_togo.lsm <- simulator(g_togo)
cc_caviar.lsm <- simulator(g_caviar)
cc_cielnet.lsm <- simulator(g_cielnet)
cc_cocaine.lsm <- simulator(g_cocaine)
cc_heroin.lsm <- simulator(g_heroin)
cc_oversize.lsm <- simulator(g_oversize)
cc_montagna.lsm <- simulator(g_montagna)
cc_tfc.lsm <- simulator(g_tfc)





# plot histogram of the clustering coefficients for the real and simulated graphs
plot_clustering <- function(coeff, lsm, title){
  
  # required packages
  require('ggplot2'); require('scales'); library('ggplot2')
  
  # mean and standard deviation of the latent space models
  mu <- mean(lsm$cc); sd <- sd(lsm$cc)
  
  # compute z-score
  z <- (coeff - mu)/sd
  
  # two-tailed significance test
  p <- 2 * (1 - stats::pnorm(abs(z)))
  
  # flag statistical significance
  sig <- if (p < 0.001) {"***"} 
  else{
    if (p < 0.1) {"**"}
    else{
      if (p < 0.5){"*"}
      else{
        NULL
      }
    }
  }
  
  # label to annotate the mean coefficient of the simulated networks
  label1 <- mean(lsm$cc); label1 <- round(label1, digits = 2)
  
  # label to annotate the coefficient for the criminal networks
  label2 <- round(coeff, digits = 2)
  
  # distribution of correlation coefficients
  histogram <- ggplot2::ggplot(lsm, ggplot2::aes( x = cc ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
    # line markers for the clustering coefficients
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(cc)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = coeff), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
    # annotation for the mean clustering coefficient for the random networks
    ggplot2::annotate(geom = "label", x = label1 - 0.02, y = 0.20, label = sprintf('%0.2f', label1), size = 3) +
    # annotation for the clustering coefficient for the actual networks
    ggplot2::annotate(geom = "label", x = label2 + 0.02, y = 0.18, label = paste0(sprintf('%0.2f', label2), sig), size = 3) +
    # transform y-axis to percentage scale
    ggplot2::aes(y = after_stat(count)/sum(after_stat(count))) + 
    ggplot2::scale_y_continuous(
      name = "PROBABILITY DENSITY FUNCTION (PDF)", 
      labels = scales::percent_format(accuracy = 1L), # 2L to round to one decimal place, 3L to round to two decimal places, etc.
      limits = c(0.00, 0.20), 
      breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)
    ) +
    ggplot2::scale_x_continuous(
      name = "CLUSTERING COEFFICIENT",
      labels = scales::label_number(accuracy = 0.01),
      limits = c(0.00, 1.00),
      breaks = c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
    ) +
    ggplot2::ggtitle(title) +
    ggthemes::theme_clean() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "plain"),
      axis.text.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain"),
      axis.text.y = ggplot2::element_text(color = "black", size = 8, hjust = 1.0, vjust = 0.0, face = "plain"),  
      axis.title.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.0, face = "plain"),
      axis.title.y = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain")
    )
  plot(histogram)
  return(histogram)
}
clustering_siren <- plot_clustering(
  coeff = cc_siren,
  lsm = cc_siren.lsm, 
  title = "(A) SIREN AUTO THEFT RING"
)
clustering_togo <- plot_clustering(
  coeff = cc_togo, 
  lsm = cc_togo.lsm, 
  title = "(B) TOGO AUTO THEFT RING"
)
clustering_caviar <- plot_clustering(
  coeff = cc_caviar, 
  lsm = cc_caviar.lsm, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
)
clustering_cielnet <- plot_clustering(
  coeff = cc_cielnet, 
  lsm = cc_cielnet.lsm, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
)
clustering_cocaine <- plot_clustering(
  coeff = cc_cocaine, 
  lsm = cc_cocaine.lsm, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
)
clustering_heroin <- plot_clustering(
  coeff = cc_heroin, 
  lsm = cc_heroin.lsm, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
)
clustering_oversize <- plot_clustering(
  coeff = cc_oversize, 
  lsm = cc_oversize.lsm, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
)
clustering_montagna <- plot_clustering(
  coeff = cc_montagna, 
  lsm = cc_montagna.lsm, 
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
)
clustering_tfc <- plot_clustering(
  coeff = cc_tfc, 
  lsm = cc_tfc.lsm, 
  title = "(I) THE FRENCH CONNECTION - FEDERAL BUREAU OF NARCOTICS"
)



# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 5, 
    height = 5, 
    device = 'pdf', 
    dpi = 700
  )
}
output(plot = clustering_siren, filename = "fig4a.pdf")
output(plot = clustering_togo, filename = "fig4b.pdf")
output(plot = clustering_caviar, filename = "fig4c.pdf")
output(plot = clustering_cielnet, filename = "fig4d.pdf")
output(plot = clustering_cocaine, filename = "fig4e.pdf")
output(plot = clustering_heroin, filename = "fig4f.pdf")
output(plot = clustering_oversize, filename = "fig4g.pdf")
output(plot = clustering_montagna, filename = "fig4h.pdf")
output(plot = clustering_tfc, file = "fig4i.pdf")





# close .r script





# function to generate the clustering coefficients of of random graphs the have similar degree distribution 
sampler <- function(g){
  
  # required packages
  require('statnet'); require('igraph'); require('intergraph')
  
  # set seed for replication purposes
  set.seed(20190812) # Hayes' birthday
  
  # dimensions of the graph
  n <- length(igraph::V(intergraph::asIgraph(g))) # 'n' nodes
  m <- length(igraph::E(intergraph::asIgraph(g))) # 'm' edges
  
  # don't run
  # calculate the modal number of neighbors from the degree distribution
  # d <- igraph::degree(intergraph::asIgraph(g), mode = "total", loops = FALSE)
  # d <- mode(d)
  
  # don't run
  # calculate the number of neighbors from the degree distribution
  # d <- igraph::degree(intergraph::asIgraph(g), mode = "total", loops = FALSE)
  
  # function to generate random graphs from the Erdos-Renyi model
  samples <- 10000 # number of random graphs to generate
  results <- c() # results vector
  for(i in 1:samples){
    
    # generate random graphs
    g_random = igraph::sample_gnm(
      n = n,
      m = m, # m edges taken from the uniform random distribution from the set of all possible edges 
      directed = FALSE, 
      loops = FALSE
    )
    
    # randomly select the number of neighbors from the degree distribution
    # nei <- sample(d, 1)
    
    # randomly select the probability to rewire edges
    # p <- runif(1, 0, 1)
    
    # don't run
    # small world sampler 
    # g_random <- igraph::sample_smallworld(
    # dim = 1, # = 1 for a lattice
    # size = n,
    # nei = nei,
    # p = p
    # )
    
    # calculate clustering coefficient for the random graphs
    cc <- igraph::transitivity(
      g_random,
      type = "localaverageundirected",
      isolates = "zero"
    )
    
    # store the clustering coefficients in the results vector
    results <- c(results, cc)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("cc")
  return(results)
}
cc_siren.random <- sampler(g_siren)
cc_togo.random <- sampler(g_togo)
cc_caviar.random <- sampler(g_caviar)
cc_cielnet.random <- sampler(g_cielnet)
cc_cocaine.random <- sampler(g_cocaine)
cc_heroin.random <- sampler(g_heroin)
cc_oversize.random <- sampler(g_oversize)
cc_montagna.random <- sampler(g_montagna)
cc_tfc.random <- sampler(g_tfc)
cc_super.random <- sampler(g_super)



# compare the clustering coefficients from the laternt space model and erdos-renyi random graph models
compare_distributions <- function(x, y){
  
  # t-test
  result <- stats::t.test(x, y)
  
  # test statistic
  t <- result$statistic
  message("t statistic:")
  cat(t); cat("\n"); cat("\n")
  
  # p-value
  p <- result$p.value
  p <- round(p, digits = 4)
  message("p-value:")
  cat(p); cat("\n"); cat("\n")
  
  # 
  x <- result$estimate[1]
  x <- round(x, digits = 2)
  message("mean clustering coefficient of Erdos-Renyi random graph model:")
  cat(x); cat("\n"); cat("\n")
  
  # 
  y <- result$estimate[2]
  y <- round(y, digits = 2)
  message("mean clustering coefficient of the latent space model:")
  cat(y); cat("\n"); cat("\n")
}
compare_distributions(x = cc_siren.random, cc_siren.lsm)
compare_distributions(x = cc_togo.random, cc_togo.lsm)
compare_distributions(x = cc_cavair.random, cc_caviar.lsm)
compare_distributions(x = cc_cielnet.random, cc_cielnet.lsm)
compare_distributions(x = cc_cocaine.random, cc_cocaine.lsm)
compare_distributions(x = cc_heroin.random, cc_heroin.lsm)
compare_distributions(x = cc_oversize.random, cc_oversize.lsm)
compare_distributions(x = cc_montagna.random, cc_montagna.lsm)
compare_distributions(x = cc_tfc.random, cc_tfc.lsm)



