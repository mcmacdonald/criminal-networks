#  -----------------------------------------------------------------------------

# file 01: import and manipulate data used to estimate the hierarchical model

# note: you must run this file before you run file named '02.R'

# ------------------------------------------------------------------------------



# first load data into R


# label rows and columns with unique names
nodes <- function(x, prefix){
  rownames(x) <- c(1:nrow(x)); rownames(x) <- paste0(prefix, rownames(x)) # prefix = "N"
  colnames(x) <- c(1:nrow(x)); colnames(x) <- paste0(prefix, colnames(x)) # prefix = "N"
  return(x)
}
caviar   <- nodes(caviar, prefix = "caviar")
cielnet  <- nodes(cielnet, prefix = "cielnet")
cocaine  <- nodes(cocaine, prefix = "cocaine")
heroin   <- nodes(heroin, prefix = "heroin")
siren    <- nodes(siren, prefix = "siren")
togo     <- nodes(togo, prefix = "togo")
oversize <- nodes(oversize, prefix = "oversize")
montagna <- nodes(montagna, prefix = "montagna")



# symmetrize the relational data -----------------------------------------------
symmetrize = function(adj){
  
  # required packages
  require("igraph")
  
  # matrix manipulation
  adj[is.na(adj)] <- 0 # set missing cells = 0 i.e., no tie (in case of missings)
  adj = as.matrix(adj) # coerce to matrix
  adj = igraph::graph_from_adjacency_matrix( # adjacency matrix
    adjmatrix = adj,
    mode = "max", 
    weighted = NULL, 
    diag = FALSE
    )
  adj = igraph::simplify( # drop edge weights, transform to binary adjacency matrix
    graph = adj, 
    remove.loops = TRUE, 
    remove.multiple = TRUE
    )
  adj = igraph::as_adjacency_matrix( # symmetric adjacency matrix to symmetrize
    adj, 
    type = "both", 
    names = TRUE
    )
  return(adj)
}
siren    = symmetrize(siren)
togo     = symmetrize(togo)
caviar   = symmetrize(caviar)
cielnet  = symmetrize(cielnet)
cocaine  = symmetrize(cocaine)
heroin   = symmetrize(heroin)
oversize = symmetrize(oversize)
montagna = symmetrize(montagna)



# function to delete isolates --------------------------------------------------
delete_isolates <- function(adj){
  require('igraph')
  g <- igraph::graph_from_adjacency_matrix(
    adj, 
    mode = "undirected", 
    weighted = NULL, 
    diag = FALSE
    )
  g <- igraph::delete_vertices( # delete isolates
    g, 
    which(igraph::degree(g, v = igraph::V(g), mode = "total", loops = F)==0) # degree = 0
    )
  el = igraph::as_edgelist( # into edgelist
    graph = g, 
    names = TRUE
    )
  el = as.data.frame( # back to dataframe
    x = el, 
    stringsAsFactors = FALSE
    ) 
  el = dplyr::rename( # rename columns
    el, 
    i = V1, 
    j = V2
    )
  return(el)
}
siren    = delete_isolates(siren)
togo     = delete_isolates(togo)
caviar   = delete_isolates(caviar)
cielnet  = delete_isolates(cielnet)
cocaine  = delete_isolates(cocaine)
heroin   = delete_isolates(heroin)
oversize = delete_isolates(oversize)
montagna = delete_isolates(montagna)



# function to generate nodelists -----------------------------------------------
nl <- function(el, net){
  el$weight <- NULL
  nl <- stack(el)
  nl$ind <- NULL
  nl <- dplyr::rename(nl, name = values)
  nl <- unique(nl)
  nl <- dplyr::mutate(nl, group = net) # attach name of networks to nodelist -- need them to assign 'group' membership in the graph
  return(nl)
}
v_siren    <- nl(siren, net = "siren")
v_togo     <- nl(togo, net = "togo")
v_caviar   <- nl(caviar, net = "caviar")
v_cielnet  <- nl(cielnet, net = "cielnet")
v_cocaine  <- nl(cocaine, net = "cocaine")
v_heroin   <- nl(heroin, net = "heroin")
v_oversize <- nl(oversize, net = "oversize")
v_montagna <- nl(montagna, net = "montagna")



# transform igraph objects to network objects to estimate models
graph <- function(e, v){
  g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = v)
  g <- igraph::simplify(g) # force the graph to not be of type 'multiple'
  g <- intergraph::asNetwork(g)
  return(g)
}
g_siren  <- graph(siren, v = v_siren)
g_togo   <- graph(togo, v = v_togo)
g_caviar <- graph(caviar, v = v_caviar)
g_cielnet  <- graph(cielnet, v = v_cielnet)
g_cocaine  <- graph(cocaine, v = v_cocaine)
g_heroin   <- graph(heroin, v = v_heroin)
g_oversize <- graph(oversize, v = v_oversize)
g_montagna <- graph(montagna, v = v_montagna)



# close .R script

