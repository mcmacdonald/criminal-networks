# ---------------------------------------------------------------------------------------------------------------------

# File 01: preparing the network data

# call packages
library('utils')
library('readr')
library('igraph')
# ---------------------------------------------------------------------------------------------------------------------

# UCInet has 'covert networks' for download, including many classic criminal networks
# source: https://www.sites.google.com/site/ucinetsoftware/datasets/covert-networks?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1

# cross-referenced from UC Boulder's index of complex networks:
# https://icon.colorado.edu/#!/networks

# A series of networks of Cosa Nostra were published either shared or published ...
# ... from previous confidential records
# source: https://github.com/lcucav/criminal-nets/tree/master/dataset





# Function to import data --------------------------------------------------
import = function(file){
data = utils::read.csv(
  file   = file, 
  header = FALSE, 
  stringsAsFactors = FALSE, 
  fileEncoding = "UTF-8-BOM"
); 
data = data[-1, -1] # drop first row, colum (i.e, vertex labels)
return(data)
}
# Below, the prefix indicates network type:
# r = ringing (auto theft) network
# d = drug trafficking network
# g = gang network
# m = mafia i.e.,  "Cosa Nostra" or 'Ndrangheta
r_siren    = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/CAVIAR_FULL.csv") # siren auto theft
r_togo     = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/TOGO.csv" ) # togo auto theft
d_caviar   = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/CAVIAR_FULL.csv") # caviar drug trafficking network - Morselli
d_cocaine  = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/COCAINE_DEALING.csv") # NY cocaine trafficking network - Natarajan
d_heroin   = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/HEROIN_DEALING.csv") # NY heroin trafficking network - Natarajan
d_cielnet  = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/CIELNET.csv") # CielNet (Montreal drug runner) network - Morselli
g_ity      = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/ITALIAN_GANGS.csv") # Italian gangs
g_ldn      = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/LONDON_GANG.csv") # London gangs
g_mtl      = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/MONTREALGANG.csv") # Montreal gangs
m_infinito = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/NDRANGHETAMAFIA_2M.csv") # Ndrangheta - meetings --- from Operation Infinito

    # transpose the infinto meeting network from two-mode to one mode projection
    m_infinito = igraph::graph.incidence(
      incidence = m_infinito, # graph by incidence/event i.e., meeting attendance
      directed = FALSE,
      mode = "total",
      weighted = NULL
      )

    # bipariate projection of the incidence/event two-mode network
    m_infinito = igraph::bipartite.projection(graph = m_infinito)
    m_infinito = m_infinito$proj1 # keep the one-mode projection
    m_infinito = igraph::get.adjacency( # convert graph to adj matrix
      graph = m_infinito,
      type = "both"
      ) 
    m_infinito = as.matrix(m_infinito)     # into matrix format
    m_infinito = as.data.frame(m_infinito) # into data frame format

    
# You can work with netork data in matrix format in various packages igraph(), sna() ...
# ... but I convert matrices into edgelist format


  
# symmetrize the network data --------------------------------------------------

# Function has three steps:
# 1. transform matrices into graphs, 
# 2. get binary graphs, or graphs with no edge weights
# 3. transform to edgelist to get symmetry
symmetrize = function(mat){
  
  # step-by-step transform into symmetrized edgelist
  mat[is.na(mat)] <- 0 # set missing cells = 0 i.e., no tie
  mat = as.matrix(mat) # coerce to matrix
  mat = igraph::graph_from_adjacency_matrix( # Adjacency matrix
    adjmatrix = mat, 
    mode = "undirected", 
    weighted = NULL, 
    diag = FALSE
    )
  mat = igraph::simplify( # into binary network
    graph = mat, 
    remove.loops = TRUE, 
    remove.multiple = TRUE
    )
  mat = igraph::as_edgelist( # into edgelist
    graph = mat, 
    names = TRUE
    )
  mat = as.data.frame( # back to dataframe
    x = mat, 
    stringsAsFactors = FALSE
    ) 
  mat = dplyr::rename( # rename columns
    mat, 
    i = V1, 
    j = V2
    )
  mat = igraph::graph_from_data_frame( # return graph object
    d = mat, 
    directed = FALSE
    )
  
  # End function
}
# symmetrize graphs 
r_siren    = symmetrize(mat = r_siren)
r_togo     = symmetrize(mat = r_togo)
d_caviar   = symmetrize(mat = d_caviar)
d_cocaine  = symmetrize(mat = d_cocaine)
d_heroin   = symmetrize(mat = d_heroin)
d_cielnet  = symmetrize(mat = d_cielnet)
g_ity      = symmetrize(mat = g_ity)
g_ldn      = symmetrize(mat = g_ldn)
g_mtl      = symmetrize(mat = g_mtl)
m_infinito = symmetrize(mat = m_infinito)
    


# ... close #.R script