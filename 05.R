#  -----------------------------------------------------------------------------------------------------------------------------

# file 05: estimate the clustering coefficient of each graph and compare them to random graphs

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------------------------------------------------



# calculate the degree assortativity coefficient of each of the criminal networks ------------------------------
clustering <- function(g){
  # required packages
  require("statnet"); require("igraph"); require("intergraph")
   cc <- igraph::transitivity(
     intergraph::asIgraph(g),
     type = "globalundirected"
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



# function to generate the clustering coefficients of of random graphs the have similar degree distribution 
sampler <- function(g){
  
  # required packages
  require('statnet'); require('igraph'); require('intergraph')
  
  # set seed for replication purposes
  set.seed(20190812) # Hayes' birthday
  
  # dimensions of the graph
  n <- length(igraph::V(intergraph::asIgraph(g))) # 'n' nodes
  m <- length(igraph::E(intergraph::asIgraph(g))) # 'm' edges
  
  # function to generate random graphs from the Erdos-Renyi model
  samples <- 10000 # number of random graphs to generate
  results <- c() # results vector
  for(i in 1:samples){
    # random graphs
    g_random = igraph::sample_gnm(
      n = n,
      m = m, # m edges taken from the uniform random distribution from the set of all possible edges 
      directed = FALSE, 
      loops = FALSE
      )
    
    # calculate clustering coefficient for the random graphs
    cc <- igraph::transitivity(
      g_random,
      type = "globalundirected"
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



# plot histogram of the clustering coefficients for the real and random graphs
plot_clustering <- function(coeff, random, title){
  require('ggplot2'); require('scales'); library('ggplot2')
  label1 <- mean(random$cc); label1 <- round(label1, digits = 2) # label to annotate the mean coefficient of the random networks
  label2 <- round(coeff, digits = 2) # label to annotate the coefficient for the criminal networks
  histogram <- ggplot2::ggplot(random, ggplot2::aes( x = cc ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
    # line marker for the mean clustering coefficient for the random networks
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(cc)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
    ggplot2::annotate(geom = "label", x = label1, y = 0.25, label = as.character(label1), size = 5) +
    # line marker for the clustering coefficient for the actual networks
    ggplot2::geom_vline(ggplot2::aes(xintercept = coeff), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
    ggplot2::annotate(geom = "label", x = label2, y = 0.25, label = as.character(label2), size = 5) +
    # transform y-axis to percentage scale
    ggplot2::aes(y = after_stat(count)/sum(after_stat(count))) + 
    ggplot2::scale_y_continuous(
      name = "PROBABILITY DENSITY FUNCTION (PDF)", 
      labels = scales::percent_format(accuracy = 1L), # 2L to round to one decimal place, 3L to round to two decimal places, etc.
      limits = c(0.00, 0.29), 
      breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25)
    ) +
    ggplot2::scale_x_continuous(
      name = "CLUSTERING COEFFICIENT",
      labels = scales::label_number(accuracy = 0.01),
      limits = c(0.00, 0.50),
      breaks = c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50)
    ) +
    ggplot2::ggtitle(title) +
    ggthemes::theme_clean()
  return(histogram)
}
clustering_siren <- plot_clustering(
  coeff = cc_siren,
  random = cc_siren.random, 
  title = "(A) SIREN AUTO THEFT RING"
  )
clustering_togo <- plot_clustering(
  coeff = cc_togo, 
  random = cc_togo.random, 
  title = "(B) TOGO AUTO THEFT RING"
  )
clustering_caviar <- plot_clustering(
  coeff = cc_caviar, 
  random = cc_caviar.random, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
  )
clustering_cielnet <- plot_clustering(
  coeff = cc_cielnet, 
  random = cc_cielnet.random, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
  )
clustering_cocaine <- plot_clustering(
  coeff = cc_cocaine, 
  random = cc_cocaine.random, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
  )
clustering_heroin <- plot_clustering(
  coeff = cc_heroin, 
  random = cc_heroin.random, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
  )
clustering_oversize <- plot_clustering(
  coeff = cc_oversize, 
  random = cc_oversize.random, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
  )
clustering_montagna <- plot_clustering(
  coeff = cc_montagna, 
  random = cc_montagna.random, 
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
  )
clustering_tfc <- plot_clustering(
  coeff = cc_tfc, 
  random = cc_tfc.random, 
  title = "(I) THE FRENCH CONNECTION - FEDERAL BUREAU OF NARCOTICS"
  )
clustering_super <- plot_clustering(
  coeff = cc_super, 
  random = cc_super.random, 
  title = "(J) SUPER POPULATION OF CRIMINAL NETWORKS"
  )



# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop", 
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
output(plot = clustering_super, filename = "fig4j.pdf")



# close .r script



