#  -----------------------------------------------------------------------------------

# file 01: import and manipulate relational data for the 'oversize' network of mafias

# note: you must run this file before you run file named '02.R'

# link to download the original data:
# Berlusconi et al.: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0154244
# most suspects under investigation were members or associates of the ‘Ndrangheta of the region of Calabria

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------



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
oversize = import(file = "~/Desktop/oversize.csv") # oversize -- wiretap records


# label rows and columns with unique names
nodes <- function(x, prefix){
  rownames(x) <- c(1:nrow(x)); rownames(x) <- paste0(prefix, rownames(x)) # prefix = "N"
  colnames(x) <- c(1:nrow(x)); colnames(x) <- paste0(prefix, colnames(x)) # prefix = "N"
  return(x)
}
oversizeWR <- nodes(oversizeWR, prefix = "oversize")
oversizeJU <- nodes(oversizeJU, prefix = "oversize")



# symmetrize the relational data ------------------------------------------------

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
  # adj = igraph::delete_vertices( # delete isolates
    # adj, 
    # which(igraph::degree(adj, v = igraph::V(adj), mode = "total", loops = F)==0)
    # )
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
oversizeWR = symmetrize(oversizeWR)
oversizeJU = symmetrize(oversizeJU)



# delete nodes from oversize wiretap records not named in court records

  # nodelist
  v <- stack(oversizeWR)
  v$ind <- NULL
  v <- unique(v$values)
  v <- as.data.frame(v)
  v <- dplyr::rename(v, name = v)

  # graph
  g <- igraph::graph_from_data_frame(
    oversizeJU, 
    directed = FALSE, 
    vertices = v
    )
  
  # degree centrality 
  degree <- igraph::degree(
    g, 
    mode = "total", 
    loops = FALSE
    )
  degree0 <- degree[degree == 0] # isolates

  # drop anyone not named in court records
  oversize <- igraph::graph_from_data_frame(oversizeWR, directed = FALSE, vertices = v)
  oversize <- igraph::delete_vertices( # delete isolates
    oversize, 
    names(degree0)
    )
  rm(degree, degree0, oversizeWR, oversizeJU, v)

# into adjacency matrix format
oversize <- igraph::as_adjacency_matrix(oversize, names = TRUE)
oversize <- as.matrix(oversize)
oversize <- as.data.frame(oversize)

# don't run
# setwd("~/Desktop")
# write.csv(oversize, file = "oversize.csv")
<<<<<<< HEAD



# close .r script
=======
>>>>>>> aa72bc6595c755ddd1eb690dcef6eea99ac85398



# close .r script


