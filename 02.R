#  ---------------------------------------------------------------------------------------------------

# file 02: import and manipulate data used to estimate the hierarchical exponential random graph model

# note: you must run this file before you run file named '03.R'

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
oversize = import(file = "~/criminal-networks/data/berlusconietal_oversize.csv") # oversize -- wiretap records



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



# symmetrize the network data --------------------------------------------------

# function has three steps:
# 1. transform matrices into adjacency matrix
# 2. transform into binary adjacency matrix, or adj matrix with no edge weights
# 3. transform into symmetric edgelist
symmetrize = function(adj){
  
  # required packages
  require("igraph")
  
  # matrix manipulation
  adj[is.na(adj)] <- 0 # set missing cells = 0 i.e., no tie (in case of missings)
  adj = as.matrix(adj) # coerce to matrix
  adj = igraph::graph_from_adjacency_matrix( # adjacency matrix
    adjmatrix = adj,
    mode = "undirected", 
    weighted = NULL, 
    diag = FALSE
    )
  adj = igraph::simplify( # into binary adjacency mat
    graph = adj, 
    remove.loops = TRUE, 
    remove.multiple = TRUE
    )
  adj = igraph::delete_vertices( # delete isolates
    adj, 
    which(igraph::degree(adj, v = igraph::V(adj), mode = "total", loops = F)==0)
    )
  el = igraph::as_edgelist( # into edgelist
    graph = adj, 
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

# symmetrize
siren    = symmetrize(siren)
togo     = symmetrize(togo)
caviar   = symmetrize(caviar)
cielnet  = symmetrize(cielnet)
cocaine  = symmetrize(cocaine)
heroin   = symmetrize(heroin)
oversize = symmetrize(oversize)



# join into super network
g_super <- rbind(
  siren,
  togo,
  caviar,
  cielnet,
  cocaine,
  heroin,
  oversize
  )


# function to generate nodelists -----------------------------------------------
nl <- function(el, net){
  nl <- stack(el)
  nl$ind <- NULL
  nl <- dplyr::rename(nl, name = values)
  nl <- unique(nl)
  nl <- dplyr::mutate(nl, group = net) # attach name of networks to nodelist -- need them to assign 'group' membership in the graph
  return(nl)
}
d_siren   <- nl(siren, net = "siren")
d_togo    <- nl(togo, net = "togo")
d_caviar  <- nl(caviar, net = "caviar")
d_cielnet <- nl(cielnet, net = "cielnet")
d_cocaine <- nl(cocaine, net = "cocaine")
d_heroin  <- nl(heroin, net = "heroin")
d_oversize<- nl(oversize, net = "oversize")

# join into attribute data for super network
d_super <- rbind(
  d_siren,
  d_togo,
  d_caviar,
  d_cielnet,
  d_cocaine,
  d_heroin,
  d_oversize
  )
rm( # drop attribute data for individual networks
  d_siren,
  d_togo,
  d_caviar,
  d_cielnet,
  d_cocaine,
  d_heroin,
  d_oversize
  )
rm( # drop edgelists for individual networks
  siren,
  togo,
  caviar,
  cielnet,
  cocaine,
  heroin,
  oversize
  )



# construct super network with names and 'group membership'
g_super <- igraph::graph_from_data_frame( # 'group membership' automatically assigned as a node attribute
  g_super, 
  directed = FALSE, 
  vertices = d_super
  )
g_super <- intergraph::asNetwork(g_super) # network object



# close .r script


