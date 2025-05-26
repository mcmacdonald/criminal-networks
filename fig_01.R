#  ---------------------------------------------------------------------------------------------------

# file 01: construct archetypes of criminal networks

# ----------------------------------------------------------------------------------------------------

# require packages
require('igraph'); require('ggraph')



# set parameters to construct graphs
set.seed(20240517) # for replication
n <- 36  # number of nodes
m <- 2  # each node connects to m number of nodes, with probability proportional to their degree
kn <- 3  # each node connects to kn nearest neighbors
p <- 0.10  # probability to rewire edges



# generate the centralized graph
g1 <- igraph::sample_pa( # generative model for scale-free networks
  n = n, 
  power = 2, 
  m = m, 
  directed = FALSE
  )

# generate the decentralized graph
g2 <- igraph::sample_smallworld( # generative model for cliquish networks
  dim = 1, 
  size = n, 
  nei = kn,
  p = p
  )


# see 'ggraph' tutorial for plotting
# https://r-graph-gallery.com/package/ggraph.html



# plot the centralized organizational structure
fig1a <- ggraph::ggraph(tidygraph::as_tbl_graph(g1), layout = "fr") +  # Kamada-Kawai layout for hierarchy
  ggraph::geom_edge_link(color = "grey", alpha = 0.50, ggplot2::aes(width = 2)) +
  ggraph::geom_node_point(size = igraph::degree(g1), color = "black", fill = "white", shape = 21, stroke = 2) +  # Size by degree
  ggraph::theme_graph(background = "white") +
  ggplot2::labs(title = "(A) CENTRALIZED ORGANIZATIONAL STRUCTURE") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = 'plain', hjust = 0.50))



# plot the decentralized organizational structure
fig1b <- ggraph::ggraph(tidygraph::as_tbl_graph(g2), layout = "fr") +
    ggraph::geom_edge_link(color = "grey", alpha = 0.50, ggplot2::aes(width = 2)) +
    ggraph::geom_node_point(size = igraph::degree(g2), color = "black", fill = "white", shape = 21, stroke = 2) +
    ggraph::theme_graph(background = "white") +
    ggplot2::labs(title = "(B) DECENTRALIZED ORGANIZATIONAL STRUCTURE") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = 'plain', hjust = 0.50))

  

# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 10, 
    height = 10, 
    device = 'png', 
    dpi = 700
  )
}
output(fig1a, "fig1a.png")
output(fig1b, "fig1b.png")


