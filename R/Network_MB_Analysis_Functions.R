## Functions used in the Network-MB analysis script
## Last updated 7 November 2019 by LE Koenig


## ---------- Functions for interfacing with NHD network structures --------- ##

  # Function to find the X,Y coordinates of the center point of each NHD flowline:
    get_coords_reach_mid <- function(x){
      # Project to albers equal area and get coordinates for the center point along the flowline:
      fline <- x %>% st_transform(5070)
      mid <- suppressWarnings(st_centroid(fline) %>% st_transform(4269))   
      X <- st_coordinates(mid)[,1]
      Y <- st_coordinates(mid)[,2]
      out <- data.frame(X=X,Y=Y)
      return(out)
    }


## ---------- Functions used in the mass balance model --------- ##

  # Function to implement the C mass-balance model:
    solveMB <- function(network){
      
      # Define input parameters:
      
      # Inflow discharge from local catchment (m3 d-1):
      V(network)$Qlocal <- 0
      # Inflow discharge from upstream reaches (m3 d-1):
      V(network)$Qnet <- 0
      # Outflow discharge from each reach (m3 d-1):
      V(network)$Qout <- NA
      # C load from upstream catchment (g d-1):
      V(network)$Clocal <- 0
      # C load from upstream reaches (g d-1):
      V(network)$Cnet <- 0
      # Exported C load from each reach (g d-1):
      V(network)$Cout <- NA
      # Uptake/decay rate (m day-1) - change vf within analysis file:
      #V(network)$vf <- 0
      # Fraction of C load lost during transport within each reach (unitless):
      V(network)$Clost <- NA
      
      # Calculate mass-balance for each reach moving down the network from headwaters to mouth:
      for(i in 1:length(V(network))){
        
        # Find neighboring reaches upstream that flow in to the reach:
        up <- igraph::neighbors(network,i,mode=c("in"))
        up.all.nodes <- unlist(ego(network,order=length(V(network)),nodes=V(network)[i],mode=c("in"),mindist=0))[-1]
        
        # Define hydrologic inflows/outflows for each reach (m3 d-1)
        # Discharge inflow from local catchment (m3 d-1):
        V(network)$Qlocal[i] <- V(network)$runoff_mday[i] * (V(network)$areasqkm[i]*10^6)
        
        # Discharge inflow from upstream network (m3 d-1):
        if(length(up)>0){
          V(network)$Qnet[i] <- sum(V(network)$Qout[up])
        }  
        
        # Discharge outflow to downstream reach (m3 d-1):
        V(network)$Qout[i] <- V(network)$Qlocal[i] + V(network)$Qnet[i]
        
        # Calculate reach hydraulic load (m d-1):
        HL <- V(network)$Qout[i]/(V(network)$width[i]*(V(network)$lengthkm[i]*1000))
        
        # Define carbon inflows/outflows for each reach (g d-1)
        # Carbon inflow from local catchment (g d-1):
        V(network)$Clocal[i] <- V(network)$Qlocal[i]*V(network)$local_C_in[i]
        
        # Carbon inflow from upstream network (g d-1):
        if(length(up)>0){
          V(network)$Cnet[i] <- sum(V(network)$Cout[up])
        } 
        
        # Exported carbon load to downstream reach (g d-1):
        V(network)$Cout[i] <- V(network)$Clocal[i] + V(network)$Cnet[i] - ((V(network)$Clocal[i] + V(network)$Cnet[i])*(1-exp(-V(network)$vf[i]/HL)))
        
        # Calculate fraction C lost
        V(network)$Clost[i] <- 1-(V(network)$Cout[i]/(V(network)$Clocal[i] + V(network)$Cnet[i]))
      }
      
      # Get list with attributes
      out <- get.vertex.attribute(network)
      
      # Calculate proportion C removed in the network:
      out$Prop_C_lost <- 1-(V(network)$Cout[length(V(network))]/sum(V(network)$Clocal))
      
      # Export network:
      return(out)
      
    }
