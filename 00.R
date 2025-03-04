


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





# import the relational data for the french connection 
fbn <- readxl::read_excel("~/Desktop/frenchconnection.xlsx", sheet = "Sheet1")
fbn$i <- stringr::str_trim(fbn$i, side = "both")
fbn$j <- stringr::str_trim(fbn$j, side = "both")

# import nodelist and attribute
fbn_nodelist <- readxl::read_excel("~/Desktop/frenchconnection.xlsx", sheet = "Sheet2")
fbn_nodelist$name <- stringr::str_trim(fbn_nodelist$name)

# retain members and associates of the french connection

# nodelist does not contain everyone in the edgelist, so construct full nodelist
v_fbn <- stack(fbn)
v_fbn <- v_fbn[, 1]
v_fbn <- unique(v_fbn)
v_fbn <- as.data.frame(v_fbn)
colnames(v_fbn) <- "name"

# drop anyone not in the crime family

# graph
fbn <- igraph::graph_from_data_frame(
  fbn,
  directed = FALSE,
  vertices = v_fbn
  )

# nodes part of or in direct contact with members of french connection
v_fc <- fbn_nodelist[, 1]
v_fc <- as.data.frame(v_fc)
colnames(v_fc) <- "name"

# nodes to delete
`%>%` <- magrittr::`%>%`
delete <- v_fbn %>% dplyr::anti_join(v_fc)

# delete nodes
fbn <- igraph::delete_vertices( # delete isolates
  fbn, 
  delete$name
  )
rm(fbn_nodelist, # drop
   v_fbn,
   v_fc,
   delete
   )

# symmetric adjacency matrix to symmetrize
fbn <- igraph::as_adjacency_matrix(fbn)
fbn <- as.matrix(fbn)
fbn <- as.data.frame(fbn)

# don't run
# setwd("~/Desktop")
# write.csv(fbn, file = "macdonald_tfc.csv")




