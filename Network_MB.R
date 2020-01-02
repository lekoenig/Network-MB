# The objective of this script is to initialize a river-network mass balance model
# Last updated 5 November 2019 by LE Koenig

library(igraph)        # network analysis; create river network graphs
library(dplyr)         # general data manipulation and cleaning
library(nhdplusTools)  # R package for interafacing with NHDPlus
library(dataRetrieval) # R package for interfacing with NWIS data
library(sf)            # spatial analysis

source("./R/Network_MB_Analysis_Functions.R")


## ---------- Sample an NHD network to use as a heuristic river network structure ---------- ##

  # Bunnell Brook, CT (as an example):
  site.number <- "01186500"
  site.info <- readNWISsite(site.number)

  # identify comid that overlaps usgs gaging station at Bunnell Brook:
  start_point <- st_sfc(st_point(c(site.info$dec_long_va, site.info$dec_lat_va)), crs = 4269)
  start_comid <- discover_nhdplus_id(start_point)

  # sample upstream tributaries:
  flowline <- navigate_nldi(list(featureSource = "comid", 
                                 featureID = start_comid), 
                            mode = "upstreamTributaries", 
                            data_source = "")

  # save subset NHD network as a gpkg:
  subset_gpkg <-subset_nhdplus(comids = flowline$nhdplus_comid,
                               output_file = tempfile(fileext = ".gpkg"),
                               nhdplus_data = "download")
  flowline <- sf::read_sf(subset_gpkg, "NHDFlowline_Network")
  catchment <- sf::read_sf(subset_gpkg, "CatchmentSP")
  #waterbody <- sf::read_sf(subset_gpkg, "NHDWaterbody")  # not worrying about waterbodies for now

  # Plot the NHD network being modeled:
  plot(sf::st_geometry(flowline), col = "blue")
  plot(sf::st_geometry(catchment), add = TRUE)

  # Select variables of interest from flowline table, and cast to points/vertices along linestrings:
  #flowline.test <- flowline %>%
  #  select(.,comid,lengthkm,ftype,reachcode,streamleve,streamorde,fromnode,tonode,areasqkm,totdasqkm) %>%
  #  mutate(edgeID = c(1:n())) %>%
  #  st_cast(.,"POINT") %>%
  #  mutate(X = st_coordinates(.)[,2],
  #         Y = st_coordinates(.)[,1]) 

  # split flowlines at vertices:
  #flowline.test2 <- st_collection_extract(lwgeom::st_split(flowline, flowline.test),"LINESTRING") %>%
  #                  select(.,comid,streamorde,streamcalc,fromnode,tonode)


## ---------- Prep river network as an igraph object ---------- ##

  # Note that for the network graphs, NHD reaches are represented as nodes, and edges describe the flow relationships between NHD reaches.

  # Use hydraulic geometry scaling based on Raymond et al. 2012: https://aslopubs.onlinelibrary.wiley.com/doi/full/10.1215/21573689-1597669
  a=12.88
  b=0.42
  c=0.408
  d=0.294
  
  # Prep flowline data frame: 
  flowline.sub <- flowline %>%
                  select(.,comid,lengthkm,ftype,reachcode,streamcalc,streamorde,fromnode,tonode,areasqkm,totdasqkm,q0001e) %>%
                  filter(streamcalc == streamorde) %>%
                  mutate(edgeID = c(1:n()),
                         tocomid = NA,
                         X_centroid = NA,
                         Y_centroid = NA,
                         local_C_in = 5,
                         # Uptake/decay rate (m day-1):
                         vf = 0.29,
                         runoff_mday = (q0001e*0.0283168*86400)/(totdasqkm*10^6),
                         width = (a*(q0001e*0.0283168)^b),
                         depth = (c*(q0001e*0.0283168)^d))

  # For each reach, find the downstream COMID and lat/lon of the flowline centroid (i.e., the middle of each respective NHD reach):
  for(i in 1:length(flowline.sub$comid)){
    tocomid <- flowline.sub$comid[which(flowline.sub$fromnode==flowline.sub$tonode[i])]
    if(length(tocomid)==0){tocomid <- NA}
    
    flowline.sub$tocomid[i] <- tocomid
    flowline.sub$X_centroid[i] <- get_coords_reach_mid(flowline.sub[i,])[,1]
    flowline.sub$Y_centroid[i] <- get_coords_reach_mid(flowline.sub[i,])[,2]
  }

  # Create nodes and links objects and define igraph object:
  nodes <- data.frame(flowline.sub$comid,flowline.sub$X_centroid,flowline.sub$Y_centroid,flowline.sub$lengthkm,flowline.sub$areasqkm,flowline.sub$totdasqkm,
                      flowline.sub$streamorde,flowline.sub$local_C_in,flowline.sub$vf,flowline.sub$runoff_mday,flowline.sub$width,flowline.sub$depth)
  names(nodes) <- c("comid","x","y","lengthkm","areasqkm","totdasqkm","streamorder","local_C_in","vf","runoff_mday","width","depth")
  
  links <- flowline.sub[-which(is.na(flowline.sub$tocomid)),c("comid","tocomid")]  # remove network outlet from links (tocomid is empty)
  
  net <- graph.data.frame(d=links,vertices=nodes,directed = T)
  
  # Plot igraph object:
  plot(net,vertex.label=NA,vertex.shape='circle',vertex.size=5,edge.width=1,edge.arrow.size=0.4,edge.color="darkgray")
  
  # Inspect igraph object:
  get.vertex.attribute(net)

  
## ---------- Apply mass balance model to all reaches within the network ---------- ##
  
  # How many upstream reaches flow into each reach?
  for(i in 1:length(V(net))){
    up.all <- ego(net,order=length(V(net)),nodes=V(net)[i],mode=c("in"),mindist=0)
    V(net)$up.all[i] <- length(unlist(up.all))
  }
  
  # Re-arrange network graph so that model starts with headwaters, and moves down the network to larger river reaches:
  net2 <- igraph::graph_from_data_frame(d=igraph::as_data_frame(net,what="edges"),
                                   vertices=igraph::as_data_frame(net,what="vertices") %>%
                                   arrange(up.all))
  # Visually check new network graph: 
  plot(net2,vertex.label=NA,vertex.shape='circle',vertex.size=5,edge.width=1,edge.arrow.size=0.4,edge.color="darkgray")
  
  # Run the mass-balance model:
  net2.mod <- solveMB(network=net2)
  print(net2.mod$Prop_C_lost)
  

        
    
    
    