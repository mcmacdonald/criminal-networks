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
oversizeWR = import(file = "~/Desktop/oversizeWR.csv") # wiretap records
overwizeJU = import(file = "~/Desktop/oversizeJU.csv") # court records

# label rows and columns with unique names
nodes <- function(x, prefix){
  rownames(x) <- c(1:nrow(x)); rownames(x) <- paste0(prefix, rownames(x)) # prefix = "N"
  colnames(x) <- c(1:nrow(x)); colnames(x) <- paste0(prefix, colnames(x)) # prefix = "N"
  return(x)
}
oversizeWR <- nodes(oversizeWR, prefix = "oversize")
oversizeJU <- nodes(oversizeJU, prefix = "oversize")



# symmetrize the relational data ------------------------------------------------
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

# symmetrize
oversizeWR = symmetrize(oversizeWR)
oversizeJU = symmetrize(oversizeJU)



# function to transform adjacency matrix to edgelist ---------------------------
to_edgelist <- function(adj){
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
oversizeWR <- to_edgelist(oversizeWR)
oversizeJU <- to_edgelist(oversizeJU)



# delete nodes from oversize wiretap records not named in court records

  # nodelist
  v <- stack(oversizeWR) # names of all nodes in the wire tap records
  v$ind <- NULL
  v <- unique(v$values)
  v <- as.data.frame(v)
  v <- dplyr::rename(v, name = v)
  
  # undirected, binary graph object for the wire tap records 
  oversizeWR <- igraph::graph_from_data_frame(
    oversizeWR,
    directed = FALSE, 
    vertices = v # nodes
    )

  # undirected, binary graph object for the court records
  oversizeJU <- igraph::graph_from_data_frame(
    oversizeJU,
    directed = FALSE, 
    vertices = v # nodes
    )
  
  # degree centrality for the court records
  degreeJU <- igraph::degree(
    oversizeJU, 
    mode = "total", 
    loops = FALSE
    )
  degreeJU0 <- degreeJU[degreeJU == 0] # isolates

  # drop anyone from the wiretap record not named in court records
  oversizeWR <- igraph::delete_vertices( # delete isolates
    oversizeWR, 
    names(degreeJU0)
    )
  rm(degreeJU, 
     degreeJU0,
     oversizeJU, 
     v
     )

# transform into adjacency matrix format
oversizeWR <- igraph::as_adjacency_matrix(oversizeWR, names = TRUE)
oversizeWR <- as.matrix(oversizeWR)
oversizeWR <- as.data.frame(oversizeWR)

# don't run
# setwd("~/Desktop")
# write.csv(oversizeWR, file = "oversize.csv")



# close .r script



