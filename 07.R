#  -----------------------------------------------------------------------------------

# file 04: simulate the degree assortativity coefficients from the Bayesian ERGMs

# last updated: 15/04/2025

# ------------------------------------------------------------------------------------



# calculate the degree assortativity coefficient of each of the criminal networks
assortativity_obs <- function(g){
  
  # required packages
  require("statnet"); require("igraph"); require("intergraph")
  
  # degree assortativity
  coef <- igraph::assortativity_degree(
    intergraph::asIgraph(g),
    directed = FALSE
  )
  
  # return
  return(coef)
}
degcor_siren <- assortativity_obs(g_siren)
degcor_togo <- assortativity_obs(g_togo)
degcor_caviar <- assortativity_obs(g_caviar)
degcor_cielnet <- assortativity_obs(g_cielnet)
degcor_cocaine <- assortativity_obs(g_cocaine)
degcor_heroin <- assortativity_obs(g_heroin)
degcor_oversize <- assortativity_obs(g_oversize)
degcor_montagna <- assortativity_obs(g_montagna)



# distribution of degree assortativity coefficients from the simulations
assortativity_sim <- function(simulations){
  
  # required packages
  require("statnet"); require("igraph"); require("intergraph")
  
  # number of simulations
  samples <- length(simulations)
  
  # store results in list
  results <- c()
  
  # loop for the simulations
  for(i in 1:samples){
    
    # calculate the degree assortativity coefficient for each of the simulations
    coef <- igraph::assortativity_degree(
      intergraph::asIgraph(simulations[[i]]),
      directed = FALSE
    )
    
    # store the degree assortativity coefficients in the results vector
    results <- c(results, coef)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("coef")
  
  # return
  return(results)
}
degcor_siren.sim <- assortativity_sim(g_siren.sim)
degcor_togo.sim <- assortativity_sim(g_togo.sim)
degcor_caviar.sim <- assortativity_sim(g_caviar.sim)
degcor_cielnet.sim <- assortativity_sim(g_cielnet.sim)
degcor_cocaine.sim <- assortativity_sim(g_cocaine.sim)
degcor_heroin.sim <- assortativity_sim(g_heroin.sim)
degcor_oversize.sim <- assortativity_sim(g_oversize.sim)
degcor_montagna.sim <- assortativity_sim(g_montagna.sim)



# figure 4. plot histogram of the degree assortativity coefficients for the real and simulated graphs
plot_assortativity <- function(obs, sim, title){
  
  # required packages
  require('ggplot2'); require('scales'); library('ggplot2')
  
  # mean and standard deviation of the latent space models
  mu <- mean(sim$coef); sd <- sd(sim$coef)
  
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
  label1 <- mean(sim$coef); label1 <- round(label1, digits = 2)
  
  # label to annotate the coefficient for the criminal networks
  label2 <- round(obs, digits = 2)
  
  # distribution of correlation coefficients
  histogram <- ggplot2::ggplot(sim, ggplot2::aes( x = coef ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 20, color = "black", fill = "white") +
    # line markers for the clustering coefficients
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(coef)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
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
    ggplot2::coord_cartesian(xlim = c(-1, 0)) +
    ggplot2::labs(x = "DEGREE ASSORTATIVITY COEFFICIENT") +
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
assortativity_siren <- plot_assortativity(
  obs = degcor_siren,
  sim = degcor_siren.sim, 
  title = "(A) SIREN AUTO THEFT RING"
)
assortativity_togo <- plot_assortativity(
  obs = degcor_togo, 
  sim = degcor_togo.sim, 
  title = "(B) TOGO AUTO THEFT RING"
)
assortativity_caviar <- plot_assortativity(
  obs = degcor_caviar, 
  sim = degcor_caviar.sim, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
)
assortativity_cielnet <- plot_assortativity(
  obs = degcor_cielnet, 
  sim = degcor_cielnet.sim, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
)
assortativity_cocaine <- plot_assortativity(
  obs = degcor_cocaine, 
  sim = degcor_cocaine.sim, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
)
assortativity_heroin <- plot_assortativity(
  obs = degcor_heroin, 
  sim = degcor_heroin.sim, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
)
assortativity_oversize <- plot_assortativity(
  obs = degcor_oversize, 
  sim = degcor_oversize.sim, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
)
assortativity_montagna <- plot_assortativity(
  obs = degcor_montagna, 
  sim = degcor_montagna.sim, 
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
output(plot = assortativity_siren, filename = "fig4a.pdf")
output(plot = assortativity_togo, filename = "fig4b.pdf")
output(plot = assortativity_caviar, filename = "fig4c.pdf")
output(plot = assortativity_cielnet, filename = "fig4d.pdf")
output(plot = assortativity_cocaine, filename = "fig4e.pdf")
output(plot = assortativity_heroin, filename = "fig4f.pdf")
output(plot = assortativity_oversize, filename = "fig4g.pdf")
output(plot = assortativity_montagna, filename = "fig4h.pdf")





# close .r script





#  ----------------------------------------------------------------------------------------------------------------------------

# file 04: estimate the degree assortativity of each graph and compare them to random graphs with the same degree distribution

# last updated: 27/02/2025

#  ----------------------------------------------------------------------------------------------------------------------------



# calculate the degree assortativity coefficient of each of the criminal networks ------------------------------
assortativity <- function(g){
  # required packages
  require('statnet'); require('igraph'); require('intergraph')
  degcor <- igraph::assortativity_degree(
    intergraph::asIgraph(g), 
    directed = FALSE)
  print(degcor)
}
b_siren <- assortativity(g_siren)
b_togo <- assortativity(g_togo)
b_caviar <- assortativity(g_caviar)
b_cielnet <- assortativity(g_cielnet)
b_cocaine <- assortativity(g_cocaine)
b_heroin <- assortativity(g_heroin)
b_oversize <- assortativity(g_oversize)
b_montagna <- assortativity(g_montagna)



# function to generate the degree assortativity coefficients of of random graphs the have similar degree distribution 
sampler <- function(g, method){
  
  # required packages
  require('statnet'); require('igraph'); require('intergraph')
  
  # set seed for replication purposes
  set.seed(20190812) # Hayes' birthday
  
  # degree distribution
  degrees <- igraph::degree(
    graph = intergraph::asIgraph(g),
    mode = "total",
    loops = FALSE,
    normalized = FALSE
    )
  degrees <- degrees[order(degrees, decreasing = FALSE)] # order from smallest to largest 
  
  # function to generate random graphs that have the same degree distribution
  samples <- 10000 # number of random graphs to generate
  results <- c() # results vector
  for(i in 1:samples){
    # sampling method of specific degree distribution
    g_random <- igraph::sample_degseq( # https://igraph.org/r/html/1.2.5/sample_degseq.html
      degrees, # degree distribution for the actual criminal networks
      method = method # the method by which the function generates the random graphs
      )
    
    # don't run
    # confirm that the object is a simple, random graph
    # print(igraph::is_simple(g_random)) # 'vl' method always generates simple, undirected random graphs
    
    # calculate assortativity coefficient for the random graphs
    b <- igraph::assortativity_degree(g_random, directed = FALSE)
    
    # store the assortativity coefficients in the results vector
    results <- c(results, b)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("b")
  return(results)
}
b_siren.random <- sampler(g_siren, method = "vl")
b_togo.random <- sampler(g_togo, method = "vl")
b_caviar.random <- sampler(g_caviar, method = "vl")
b_cielnet.random <- sampler(g_cielnet, method = "vl")
b_cocaine.random <- sampler(g_cocaine, method = "vl")
b_heroin.random <- sampler(g_heroin, method = "vl")
b_oversize.random <- sampler(g_oversize, method = "vl")
b_montagna.random <- sampler(g_montagna, method = "vl")
b_tfc.random <- sampler(g_tfc, method = "vl")
b_super.random <- sampler(g_super, method = "vl")



# plot histogram of the assortativity coefficients for the real and random graphs
plot_assortativity <- function(coeff, random, title){
  require('ggplot2'); require('scales'); require("ggthemes"); library('ggplot2')
  label1 <- mean(random$b); label1 <- round(label1, digits = 2) # label to annotate the mean coefficient of the random networks
  label2 <- round(coeff, digits = 2) # label t0o annotate the coefficient for the criminal networks
  histogram <- ggplot2::ggplot(random, ggplot2::aes( x = b ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
    # line marker for the mean assortativity coefficient for the random networks
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(b)), color = "skyblue2", linewidth = 1, linetype = "solid") +
    ggplot2::annotate(geom = "label", x = label1, y = 0.50, label = as.character(label1), size = 3) +
    # line marker for the assortativty coefficient for the actual networks
    ggplot2::geom_vline(ggplot2::aes(xintercept = coeff), colour = "firebrick1", linewidth = 1, linetype = "dashed") +
    ggplot2::annotate(geom = "label", x = label2, y = 0.50, label = as.character(label2), size = 3) +
    # transform y-axis to percentage scale
    ggplot2::aes(y = after_stat(count)/sum(after_stat(count))) + 
    ggplot2::scale_y_continuous(
      name = "PROBABILITY DENSITY FUNCTION (PDF)", 
      labels = scales::percent_format(accuracy = 1L), # 2L to round to one decimal place, 3L to round to two decimal places, etc.
      limits = c(0.00, 0.50), 
      breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
      ) +
    ggplot2::scale_x_continuous(
      name = "DEGREE ASSORTATIVITY COEFFICIENT",
      labels = scales::label_number(accuracy = 0.01),
      limits = c(-0.61, 0.10),
      breaks = c(-0.50, -0.40, -0.30, -0.20, -0.10, 0.00)
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
  print(histogram)
  return(histogram)
}
assortativity_siren <- plot_assortativity(
  coeff = b_siren, 
  random = b_siren.random, 
  title = "(A) SIREN AUTO THEFT RING"
  )
assortativity_togo <- plot_assortativity(
  coeff = b_togo, 
  random = b_togo.random, 
  title = "(B) TOGO AUTO THEFT RING"
  )
assortativity_caviar <- plot_assortativity(
  coeff = b_caviar, 
  random = b_caviar.random, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
  )
assortativity_cielnet <- plot_assortativity(
  coeff = b_cielnet, 
  random = b_cielnet.random, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
  )
assortativity_cocaine <- plot_assortativity(
  coeff = b_cocaine, 
  random = b_cocaine.random, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
  )
assortativity_heroin <- plot_assortativity(
  coeff = b_heroin, 
  random = b_heroin.random, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
  )
assortativity_oversize <- plot_assortativity(
  coeff = b_oversize, 
  random = b_oversize.random, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
  )
assortativity_montagna <- plot_assortativity(
  coeff = b_montagna, 
  random = b_montagna.random, 
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
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
output(plot = assortativity_siren, filename = "fig3a.pdf")
output(plot = assortativity_togo, filename = "fig3b.pdf")
output(plot = assortativity_caviar, filename = "fig3c.pdf")
output(plot = assortativity_cielnet, filename = "fig3d.pdf")
output(plot = assortativity_cocaine, filename = "fig3e.pdf")
output(plot = assortativity_heroin, filename = "fig3f.pdf")
output(plot = assortativity_oversize, filename = "fig3g.pdf")
output(plot = assortativity_montagna, filename = "fig3h.pdf")
output(plot = assortativity_tfc, filename = "fig3i.pdf")
output(plot = assortativity_super, filename = "fig3j.pdf")



# close .r script





# don't run
# this function rewires the edges of each graph instead of randomly generating one through sampler() function 

# function to rewire edges of graph, but keeping the original degree distribution intact
rewire <- function(g){
  
  # required packages
  require('statnet'); require('igraph'); require('intergraph'); require("magrittr")
  
  # set seed for replication purposes
  set.seed(20190812) # Hayes' birthday
  
  # call pipe to string functions together
  `%>%` <- magrittr::`%>%`
  # function to rewire the graph while keeping the original degree distribution intact
  samples <- 10000 # number of edges to rewire
  results <- c() # results vector
  for(i in 1:samples){
  g_rewire <- intergraph::asIgraph(g) %>% 
    igraph::rewire(igraph::keeping_degseq(niter = samples, loops = FALSE))
  # print_all(igraph::rewire(g, with = igraph::keeping_degseq(niter = igraph::vcount(g) * 10)))
  
  # calculate assortativity coefficient for the random graphs
  b <- igraph::assortativity_degree(g_rewire, directed = FALSE)
  
  # store the assortativity coefficients in the results vector
  results <- c(results, b)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("b")
  return(results)
}
b_siren.rewire <- rewire(g_siren)
b_togo.rewire <- rewire(g_togo)
b_caviar.rewire <- rewire(g_caviar)
b_cielnet.rewire <- rewire(g_cielnet)
b_cocaine.rewire <- rewire(g_cocaine)
b_heroin.rewire <- rewire(g_heroin)
b_oversize.rewire <- rewire(g_oversize)
b_montagna.rewire <- rewire(g_montagna)
b_super.rewire <- rewire(g_super)



# test if the assortativity coefficients calculated from the different methods considerably differ from one another
test <- function(x, y){
  x <- x[, 1]
  y <- y[, 1]
  # 'independence test' of the data generating process
  message("Pearson correlation (r) of the assortativity coefficients calculated from the different methods:")
  r <- cor(x, y, method = "pearson"); cat(r); cat("\n"); cat("\n")
  message("... r should be approximately 0.00 if each method to generate is independent of the other.")
  # calculate the sample means
  message("Average assortativity coefficient for the sampler function:")
  cat(mean(x)); cat("\n"); cat("\n")
  message("Average degree assortativity coefficient for the rewire function:")
  cat(mean(y)); cat("\n"); cat("\n")
}
test(x = b_siren.random, y = b_siren.rewire)
test(x = b_togo.random, y = b_togo.rewire)
test(x = b_cavair.random, y = b_caviar.rewire)
test(x = b_cielnet.random, y = b_cielnet.rewire)
test(x = b_cocaine.random, y = b_cocaine.rewire)
test(x = b_heroin.random, y = b_heroin.rewire)
test(x = b_oversize.random, y = b_oversize.rewire)
test(x = b_montagna.random, y = b_montagna.rewire)
test(x = b_super.random, y = b_super.rewire)










# figure 4. plot histogram of the degree assortativity coefficients for the real and simulated graphs
plot_assortativity <- function(obs, sim, title){
  
  # required packages
  require('ggplot2'); require('scales'); library('ggplot2')
  
  # mean and standard deviation of the latent space models
  mu <- mean(sim$coef); sd <- sd(sim$coef)
  
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
  label1 <- mean(sim$coef); label1 <- round(label1, digits = 2)
  
  # label to annotate the coefficient for the criminal networks
  label2 <- round(obs, digits = 2)
  
  # distribution of correlation coefficients
  histogram <- ggplot2::ggplot(sim, ggplot2::aes( x = coef ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 20, color = "black", fill = "white") +
    # line markers for the clustering coefficients
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(coef)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
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
    ggplot2::coord_cartesian(xlim = c(-1, 0)) +
    ggplot2::labs(x = "DEGREE ASSORTATIVITY COEFFICIENT") +
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
assortativity_siren <- plot_assortativity(
  obs = degcor_siren,
  sim = degcor_siren.sim, 
  title = "(A) SIREN AUTO THEFT RING"
)
assortativity_togo <- plot_assortativity(
  obs = degcor_togo, 
  sim = degcor_togo.sim, 
  title = "(B) TOGO AUTO THEFT RING"
)
assortativity_caviar <- plot_assortativity(
  obs = degcor_caviar, 
  sim = degcor_caviar.sim, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
)
assortativity_cielnet <- plot_assortativity(
  obs = degcor_cielnet, 
  sim = degcor_cielnet.sim, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
)
assortativity_cocaine <- plot_assortativity(
  obs = degcor_cocaine, 
  sim = degcor_cocaine.sim, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
)
assortativity_heroin <- plot_assortativity(
  obs = degcor_heroin, 
  sim = degcor_heroin.sim, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
)
assortativity_oversize <- plot_assortativity(
  obs = degcor_oversize, 
  sim = degcor_oversize.sim, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
)
assortativity_montagna <- plot_assortativity(
  obs = degcor_montagna, 
  sim = degcor_montagna.sim, 
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
output(plot = assortativity_siren, filename = "fig4a.pdf")
output(plot = assortativity_togo, filename = "fig4b.pdf")
output(plot = assortativity_caviar, filename = "fig4c.pdf")
output(plot = assortativity_cielnet, filename = "fig4d.pdf")
output(plot = assortativity_cocaine, filename = "fig4e.pdf")
output(plot = assortativity_heroin, filename = "fig4f.pdf")
output(plot = assortativity_oversize, filename = "fig4g.pdf")
output(plot = assortativity_montagna, filename = "fig4h.pdf")



# close .r script

