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
    header = TRUE, 
    stringsAsFactors = FALSE, 
    fileEncoding = "UTF-8-BOM"
  ); 
  # data = data[-1, -1] # drop first row, colum (i.e, vertex labels)
  return(data)
}
# import matrix
# m = mafia i.e., 'Ndrangheta
m_infinito = import(file = "https://raw.githubusercontent.com/mcmacdonald/criminal-networks/master/data/NDRANGHETAMAFIA_2M.csv") # Ndrangheta meeting participation i.e., Operation Infinito

# transpose the infinto meeting network from two-mode to one-mode projection
m_1 = igraph::graph.incidence(
  incidence = m_infinito, # graph by incidence/event i.e., meeting attendance
  directed = FALSE,
  mode = "total",
  weighted = NULL
)

# two-mode graph
m_2 = igraph::graph_from_incidence_matrix(
  incidence = m_infinito, # graph by incidence/event i.e., meeting attendance
  directed = FALSE,
  mode = c("total"),
  multiple = FALSE,
  weighted = NULL,
  add.names = NULL
  )

# check whether the graph is two-mode
igraph::is_bipartite(m_2) # = TRUE when two-mode


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
# 1. transform matrices into graphs
# 2. get binary graphs, or graphs with no edge weights
# 3. transform into symmetric edgelist
symmetrize = function(mat){
  mat[is.na(mat)] <- 0 # set missing cells = 0 i.e., no tie (in case of missings)
  mat = as.matrix(mat) # coerce to matrix
  mat = igraph::graph_from_adjacency_matrix( # adj matrix
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
  # close function
}
# symmetrize graphs 
m_infinito = symmetrize(mat = m_infinito)



# ... close #.R script