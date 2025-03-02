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
b_super <- assortativity(g_super)



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
b_cavair.random <- sampler(g_caviar, method = "vl")
b_cielnet.random <- sampler(g_cielnet, method = "vl")
b_cocaine.random <- sampler(g_cocaine, method = "vl")
b_heroin.random <- sampler(g_heroin, method = "vl")
b_oversize.random <- sampler(g_oversize, method = "vl")
b_montagna.random <- sampler(g_montagna, method = "vl")
b_super.random <- sampler(g_super, method = "vl")



# plot histogram of the assortativity coefficients for the real and random graphs
plot_assortativity <- function(coeff, random, title){
  require('ggplot2'); require('scales'); library('ggplot2')
  label1 <- mean(random$b); label1 <- round(label1, digits = 2) # label to annotate the mean coefficient of the random networks
  label2 <- round(coeff, digits = 2) # label t0o annotate the coefficient for the criminal networks
  histogram <- ggplot2::ggplot(random, ggplot2::aes( x = b ) ) + 
    # ggplot2::geom_histogram(ggplot2::aes(y = stat(density)), bins = 25, color = "black", fill = "white") +
    ggplot2::geom_histogram(bins = 100, color = "black", fill = "white") +
    # line marker for the mean assortativity coefficent for the random networks
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(b)), color = "skyblue2", linewidth = 1, linetype = "dashed") +
    ggplot2::annotate(geom = "label", x = label1, y = 0.25, label = as.character(label1), size = 5) +
    # line marker for the assortativty coefficient for the actual networks
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
      name = "DEGREE ASSORTATIVITY COEFFICIENT",
      labels = scales::label_number(accuracy = 0.01),
      limits = c(-0.52, 0.00),
      breaks = c(-0.50, -0.40, -0.30, -0.20, -0.10, 0.00)
      ) +
    ggplot2::ggtitle(title) +
    ggthemes::theme_clean()
  return(histogram)
}
assortativity_siren <- plot_assortativity(coeff = b_siren, random = b_siren.random, title = "(A) SIREN AUTO THEFT RING")
assortativity_togo <- plot_assortativity(coeff = b_togo, random = b_togo.random, title = "(B) TOGO AUTO THEFT RING")
assortativity_caviar <- plot_assortativity(coeff = b_caviar, random = b_cavair.random, title = "(C) CAVAIR DRUG TRAFFICKING ORGANIZATION")
assortativity_cielnet <- plot_assortativity(coeff = b_cielnet, random = b_cielnet.random, title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION")
assortativity_cocaine <- plot_assortativity(coeff = b_cocaine, random = b_cocaine.random, title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT")
assortativity_heroin <- plot_assortativity(coeff = b_heroin, random = b_heroin.random, title = "(F) NEW YORK CITY HEROIN TRAFFICKERS")
assortativity_oversize <- plot_assortativity(coeff = b_oversize, random = b_oversize.random, title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE")
assortativity_montagna <- plot_assortativity(coeff = b_montagna, random = b_montagna.random, title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA")
assortativity_super <- plot_assortativity(coeff = b_super, random = b_super.random, title = "(I) SUPER POPULATION OF CRIMINAL NETWORKS")



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
output(plot = assortativity_super, filename = "fig3i.pdf")



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




