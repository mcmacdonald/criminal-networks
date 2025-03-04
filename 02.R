#  ---------------------------------------------------------------------------------------------------

# file 02: import and manipulate data used to estimate the hierarchical exponential random graph model

# note: you must run this file before you run file named '03.R'

# last updated: 27/02/2025

# ----------------------------------------------------------------------------------------------------



# function to import data ------------------------------------------------------
import = function(file){
  
  #  required packages
  require("utils"); require("readr")
  
  # read.csv
  data = utils::read.csv(
    file   = file, 
    header = FALSE, 
    stringsAsFactors = FALSE, 
    fileEncoding = "UTF-8-BOM"
    ); 
  data = data[-1, -1] # drop first row, column (i.e, vertex labels)
  return(data)
}
caviar   = import(file = "~/criminal-networks/data/morselli_caviar.csv") # caviar drug trafficking network - Morselli
cielnet  = import(file = "~/criminal-networks/data/morselli_cielnet.csv") # CielNet (Montreal drug runner) network - Morselli
cocaine  = import(file = "~/criminal-networks/data/natajaran_cocaine.csv") # NY cocaine trafficking network - Natarajan
heroin   = import(file = "~/criminal-networks/data/natajaran_heroin.csv") # NY heroin trafficking network - Natarajan
siren    = import(file = "~/criminal-networks/data/morselli_siren.csv") # siren auto theft - Morselli
togo     = import(file = "~/criminal-networks/data/morselli_togo.csv") # togo auto theft - Morselli
oversize = import(file = "~/criminal-networks/data/berlusconietal_oversize.csv") # 'Ndrangheta -- wiretap records for those named in court records
montagna = import(file = "~/criminal-networks/data/luciaetal_montagna.csv") # 'Cosa Nostra -- wiretap records for those arrested and indicted on criminal charges
tfc      = import(file = "~/criminal-networks/data/macdonald_tfc.csv") # The French Connection international heroin trafficking organization


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
# tfc     <- nodes(tfc, prefix = "tfc) # keep names 


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
tfc      = symmetrize(tfc)



# function to transform adjacency matrix to edgelist ---------------------------
to_edgelist <- function(adj){
  g <- igraph::graph_from_adjacency_matrix(adj)
  e <- igraph::as_data_frame(g)
  e = dplyr::rename( # rename columns
    e, 
    i = from, 
    j = to
    )
  return(e)
}
siren <- to_edgelist(siren)
togo <- to_edgelist(togo)
caviar <- to_edgelist(caviar)
cielnet <- to_edgelist(cielnet)
cocaine <- to_edgelist(cocaine)
heroin <- to_edgelist(heroin)
oversize <- to_edgelist(oversize)
montagna <- to_edgelist(montagna)
tfc <- to_edgelist(tfc)



# function to delete isolates --------------------------------------------------
delete_isolates <- function(el){
  require('igraph')
  g <- igraph::graph_from_data_frame(el, directed = FALSE)
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
tfc = delete_isolates(tfc)



# function to generate nodelists -----------------------------------------------
nl <- function(el, net){
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
v_tfc      <- nl(tfc, net = "tfc")



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
g_tfc <- graph(tfc, v = v_tfc)



# append edgelists to create 'super network' -----------------------------------
g_super <- rbind(
  siren,
  togo,
  caviar,
  cielnet,
  cocaine,
  heroin,
  oversize,
  montagna,
  tfc
  )

# join into attribute data for super network
v_super <- rbind(
  v_siren,
  v_togo,
  v_caviar,
  v_cielnet,
  v_cocaine,
  v_heroin,
  v_oversize,
  v_montagna,
  v_tfc
  )

# construct super network with names and 'group membership'
g_super <- igraph::graph_from_data_frame( # 'group membership' automatically assigned as a node attribute
  g_super, 
  directed = FALSE, 
  vertices = v_super
  )
g_super <- intergraph::asNetwork(g_super) # network object
network::is.network(g_super) # check that graph is of type 'network'

rm( # drop attribute data for individual networks
  v_siren,
  v_togo,
  v_caviar,
  v_cielnet,
  v_cocaine,
  v_heroin,
  v_oversize,
  v_montagna,
  v_tfc
)
rm( # drop edgelists for individual networks
  siren,
  togo,
  caviar,
  cielnet,
  cocaine,
  heroin,
  oversize,
  montagna,
  tfc
  )





# close .r script


