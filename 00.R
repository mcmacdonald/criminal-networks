


# import relational data
montagna <- readr::read_csv("~/Desktop/montagna_phone.csv") # Montagna wiretap phone records

# import node attribute data
montagna_attributes <- read_csv("~/Desktop/montagna_attributes.csv")
montagna_attributes <- dplyr::select(montagna_attributes, node, arrest)

# nodes to delete
suspects <- dplyr::filter(montagna_attributes, arrest == "N")
suspects <- dplyr::select(suspects, node)

# graph
montagna <- igraph::graph_from_data_frame(
  montagna,
  directed = TRUE,
  vertices = montagna_attributes$node
  )

# drop anyone from the wiretap record not named in court records
montagna <- igraph::delete_vertices( # delete isolates
  montagna, 
  suspects$node
  )
rm(montagna_attributes, suspects) # drop



# symmetric adjacency matrix to symmetrize
montagna <- igraph::as_adjacency_matrix(montagna, attr = "weight")
montagna <- as.matrix(montagna)
montagna <- as.data.frame(montagna)

# don't run
# setwd("~/Desktop")
# write.csv(montagna, file = "luciaetal_montagna.csv")

