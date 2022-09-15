# Load packages
library(igraph)
library(doParallel)
library(foreach)
library(pROC)
library(splitstackshape)

# Set seed
set.seed(20190909)

###############
## Functions ##
###############

#' Generate graph
#' 
#' This function generates a graph.
#' The vertex names are the numbers 1 through N,
#' but as characters, not integers or floats.
#' 
#' @param type This is the type of graph.
#' Options are "PA" for preferential attachment,
#' "ER" for G(n,p),
#' "complete", "normal" for vertices drawn
#' from a bivariate standard normal distribution
#' (and edge weights corresponding to distance),
#' and "SBM" for starting with a "normal" graph
#' and then adding to each coordinate
#' 1 times the sign of that coordinate
#' (thus separating the four quadrants from
#' each other).
#' @param N Number of vertices in the graph.
#' @param m Number of new edges emanating
#' from each new vertex in a preferential
#' attachment graph.
#' @param d Distribution to draw from
#' for the edge weights.
#' Options are "uniform", "gamma",
#' "folded normal", "exponential",
#' and "1".
#' The "1" means that every edge has weight 1.
#' No matter which distribution is chosen,
#' the edge weights are independent and identically
#' distributed.
#' @param r Number of digits
#' to round the edge weights to.
#' Any resulting zeroes are changed to half
#' the minimum edge weight.
#' 
#' @return A graph.
genGraph <- function(type,
                     N = 1e2,
                     m = 3,
                     d = "uniform",
                     r = NULL){
  # Generate random graph
  if(type == "PA"){
    g <- sample_pa(n = N,
                   m = m,
                   directed = FALSE)
  } else if(type == "ER"){
    g <- sample_gnp(n = N,
                    p = 0.5)
  } else if(type == "complete"){
    g <- make_full_graph(n = N)
  } else if(type == "normal" | type == "SBM"){
    coords <- matrix(data = rnorm(n = N * 2),
                     nrow = N,
                     ncol = 2)
    if(type == "SBM"){
      coords <- coords + sign(coords)
    }
    distance <- as.matrix(dist(coords))
    g <- graph_from_adjacency_matrix(adjmatrix = distance,
                                     mode = "undirected",
                                     weighted = TRUE,
                                     diag = FALSE)
    V(g)$coord1 <- coords[ ,1]
    V(g)$coord2 <- coords[ ,2]
  }
  
  # The vertex names are characters
  V(g)$name <- as.character(1:N)
  
  if(d == "uniform"){
    E(g)$weight <- runif(n = length(E(g)))
  } else if(d == "gamma"){
    E(g)$weight <- rgamma(n = length(E(g)),
                          shape = 5,
                          rate = 1)
  } else if(d == "folded normal"){
    E(g)$weight <- abs(rnorm(n = length(E(g))))
  } else if(d == "exponential"){
    E(g)$weight <- rexp(n = length(E(g)))
  } else if(d == "1"){
    E(g)$weight <- 1
  }
  if(!is.null(r)){
    E(g)$weight <- round(x = E(g)$weight,
                         digits = r)
    zeros <- E(g)$weight == 0
    E(g)$weight[zeros] <- min(E(g)$weight[!zeros]) / 2
  }
  
  return(g)
}

#' Get subgraph
#' 
#' Sample the vertices of a graph and return
#' either the vertices or the induced subgraph. 
#' 
#' @param g A graph.
#' @param p Proportion of vertices to sample.
#' @param V1 Vector of vertex names.
#' Should be in order such that there is
#' an edge between the ith element of V1
#' and the ith element of V2.
#' In other words,
#' the edges in g are given by rows of cbind(V1,V2).
#' E(g) must be in the same order as cbind(V1,V2).
#' @param V2 
#' @param type Type of sampling.
#' Options are "uniform"
#' (sample vertices uniformly at random);
#' "far" (if the graph is complete,
#' sample vertices in proportion to their strength,
#' or total edge weight;
#' if the graph is not complete,
#' sample vertices in approximately inverse
#' proportion to their degree);
#' "near" (if the graph is complete,
#' sample vertices in roughly inverse proportion to their strength,
#' or total edge weight;
#' if the graph is not complete,
#' sample vertices in approximate
#' proportion to their degree);
#' "walk" (random walk,
#' as described in manuscript);
#' and "quadrant" (uses points
#' in specified quadrants of the Cartesian plane).
#' @param returnVertices 
#' 
#' @return If returnVertices is TRUE,
#' a character vector of vertex names.
#' Otherwise, the subgraph induced by sampling.
getSubGraph <- function(g,
                        p,
                        V1,
                        V2,
                        type,
                        returnVertices = FALSE){
  # Get number of vertices
  N <- gorder(g)
  
  # Number of vertices in sampled subgraph
  n <- round(p * N)
  
  # Sample from graph
  if(type == "uniform"){
    v <- sample(x = V(g)$name,
                size = n)
  } else if(type == "far"){
    if(gsize(g) == choose(N,2)){
      v <- sample(x = V(g)$name,
                  size = n,
                  prob = strength(g))
    } else {
      v <- sample(x = V(g)$name,
                  size = n,
                  prob = max(degree(g)) - degree(g) + max(1,
                                                          min(degree(g))))
    }
  } else if(type == "near"){
    if(gsize(g) == choose(N,2)){
      v <- sample(x = V(g)$name,
                  size = n,
                  prob = max(strength(g)) - strength(g) + min(strength(g)))
    } else {
      v <- sample(x = V(g)$name,
                  size = n,
                  prob = degree(g) + (min(degree(g)) == 0))
    }
  } else if(type == "walk"){
    v <- rep(x = NA,
             times = n)
    v[1] <- sample(x = V(g)$name,
                   size = 1)
    for(i in 2:n){
      indices <- which(V1 == v[i - 1] | V2 == v[i - 1])
      if(length(indices) < 1){
        v[i] <- sample(x = V(g)$name[!(V(g)$name %in% v)],
                       size = 1)
      } else {
        iNeighbors <- as.vector(t(as.matrix(cbind(V1,V2)[indices, ])))
        iNeighbors <- iNeighbors[iNeighbors != v[i - 1]]
        # This next line relies on E(g) being in the same order
        # as cbind(V1,V2)
        iEdges <- E(g)$weight[indices][!(iNeighbors %in% v)]
        iNeighbors <- iNeighbors[!(iNeighbors %in% v)]
        if(length(iNeighbors) < 1){
          v[i] <- sample(x = V(g)$name[!(V(g)$name %in% v)],
                         size = 1)
        } else if(length(iNeighbors) == 1){
          # This next line requires that V1 and V2 be of the same class
          # as V(g)
          v[i] <- iNeighbors
        } else {
          v[i] <- sample(x = iNeighbors,
                         size = 1,
                         prob = max(iEdges) - iEdges + min(iEdges))
        }
      }
    }
  } else if(type == "quadrant"){
    if(p >= 0 & p <= 0.25){
      v <- V(g)$name[V(g)$coord1 >= 0 & V(g)$coord2 >= 0]
    } else if(p <= 0.5){
      v <- V(g)$name[V(g)$coord1 >= 0]
    } else if(p <= 0.75){
      v <- V(g)$name[V(g)$coord1 >= 0 | V(g)$coord2 >= 0]
    } else if(p <= 1){
      v <- V(g)$name
    }
  }
  
  if(returnVertices){
    return(v)
  } else {
    return(induced_subgraph(graph = g,
                            vids = v))
  }
}

#' Sample graph and check edges of MST
#' 
#' This function takes in a graph g,
#' samples the vertices to get a subgraph,
#' finds the MST of the subgraph,
#' and returns information about which edges
#' from the original graph were either
#' (a) in the subgraph's MST or
#' (b) in the subgraph but not in
#' the subgraph's MST.
#' 
#' @param g A graph.
#' @param p Proportion of vertices to sample.
#' @param V1 Vector of vertex names.
#' Should be in order such that there is
#' an edge between the ith element of V1
#' and the ith element of V2.
#' In other words,
#' the edges in g are given by rows of cbind(V1,V2).
#' E(g) must be in the same order as cbind(V1,V2).
#' @param V2 
#' @param type Type of sampling.
#' Options are described in getSubGraph.
#' 
#' @return A numeric vector twice the length of V1
#' (i.e., with twice as many elements as there are edges
#' in g).
#' For the first half of the returned vector,
#' a 1 indicates that the corresponding edge
#' in g was in the subgraph's MST;
#' a 0 indicates that the corresponding edge
#' in g was not in the subgraph's MST.
#' For the second half of the returned vector,
#' a 1 indicates that the corresponding edge
#' in g was in the subgraph but not in
#' the subgraph's MST;
#' a 0 indicates that the corresponding edge
#' in g was either not in the subgraph at all
#' or was in the subgraph's MST.
randMST <- function(g,p,V1,V2,type){
  
  # Sample from graph
  h <- getSubGraph(g = g,
                   p = p,
                   V1 = V1,
                   V2 = V2,
                   type = type)
  
  # Get MST
  mstH <- as_edgelist(mst(h))
  h <- as_edgelist(h)
  
  # Get MST info
  result1 <- rep(x = 0,
                 times = length(V1))
  if(nrow(mstH) > 0){
    for(i in 1:nrow(mstH)){
      result1[V1 %in% mstH[i, ] & V2 %in% mstH[i, ]] <- 1
    }
  }
  
  # Get non-MST info
  result2 <- rep(x = 0,
                 times = length(V1))
  if(nrow(h) > 0){
    for(i in 1:nrow(h)){
      result2[V1 %in% h[i, ] & V2 %in% h[i, ]] <- 1
    }
  }
  result2 <- result2 * (1 - result1)
  
  return(c(result1,
           result2))
}

#' Single graph simulation
#' 
#' This function:
#' (1) Simulates a graph g or accepts it as an argument;
#' (2) Finds the MST of g;
#' (3) For each p in (0.25,0.50,0.75) and
#' each of a variety of sampling types,
#' samples a subgraph h;
#' (4) Finds the MST of h; and
#' (5) Bootstraps from h.
#' The process is described in the "Simulation Study"
#' subsection of the "Methods" section of the manuscript.
#' 
#' @param typeOfGraph Either a graph
#' or a type of graph to pass to genGraph.
#' @param N Number of vertices
#' to pass to genGraph.
#' @param m Preferential attachment parameter
#' to pass to genGraph.
#' @param d Distribution to pass to genGraph.
#' @param r Number of digits for rounding
#' to pass to genGraph.
#' 
#' @return A numeric vector with 7 sections:
#' (1) Length: For normal or SBM, 15; for anything else, 12.
#' Contents: The PPV_i from step 5 of the simulation study.
#' (2) Length: Same as for (1).
#' Contents: The average BPPV_i from step 7 of the simulation study.
#' (3) Length: Same as for (1).
#' Contents: The AUC_i from step 8 of the simulation study.
#' (4) Length: N.
#' Contents: The degree distribution of g.
#' (5) Length: Same as for (1).
#' Contents: The NPV_i from step 5 of the simulation study
#' (cut from the most recent version of the manuscript).
#' (6) Length: Same as for (1).
#' Contents: The average BNPV_i from step 7 of the simulation study
#' (cut from the most recent version of the manuscript).
#' (7) Length: Same as for (1).
#' Contents: The NPV AUC_i from step 8 of the simulation study
#' (cut from the most recent version of the manuscript).
sim <- function(typeOfGraph,
                N = 1e2,
                m = 3,
                d = "uniform",
                r = NULL){
  
  # Generate graph
  if(class(typeOfGraph) == "igraph"){
    g <- typeOfGraph
    N <- vcount(g)
  } else if (class(typeOfGraph) == "character"){
    g <- genGraph(type = typeOfGraph,
                  N = N,
                  m = m,
                  d = d,
                  r = r)
  }
  
  # Get degree distribution
  dd <- degree_distribution(graph = g)
  if(length(dd) < N){
    dd <- c(dd,
            rep(x = 0,
                times = N - length(dd)))
  }
  
  # Make data frame
  data <- as.data.frame(x = as_edgelist(g),
                        stringsAsFactors = FALSE)
  
  # Get MST
  mstG <- as_edgelist(mst(g))
  
  # Add MST info to data frame
  data$mstG <- 0
  for(i in 1:nrow(mstG)){
    data$mstG[data$V1 %in% mstG[i, ] & data$V2 %in% mstG[i, ]] <- 1
  }
  
  # Sample, calculate MST, and record edges
  test1 <- class(typeOfGraph) == "character" && typeOfGraph == "normal"
  test2 <- class(typeOfGraph) == "igraph" && !is.null(typeOfGraph$coord1)
  if(test1 | test2){
    result <- matrix(data = 0,
                     nrow = length(data$V1),
                     ncol = 15)
    colnames(result) <- c("uniform25",
                          "far25",
                          "near25",
                          "walk25",
                          "quadrant25",
                          "uniform50",
                          "far50",
                          "near50",
                          "walk50",
                          "quadrant50",
                          "uniform75",
                          "far75",
                          "near75",
                          "walk75",
                          "quadrant75")
    typeOfSampling <- c("uniform",
                        "far",
                        "near",
                        "walk",
                        "quadrant")
  } else {
    result <- matrix(data = 0,
                     nrow = length(data$V1),
                     ncol = 12)
    colnames(result) <- c("uniform25",
                          "far25",
                          "near25",
                          "walk25",
                          "uniform50",
                          "far50",
                          "near50",
                          "walk50",
                          "uniform75",
                          "far75",
                          "near75",
                          "walk75")
    typeOfSampling <- c("uniform",
                        "far",
                        "near",
                        "walk")
  }
  result2 <- result
  b <- rep(x = NA,
           times = ncol(result))
  names(b) <- colnames(result)
  AUCs <- b
  b2 <- b
  AUCs2 <- AUCs
  p <- c(0.25,0.5,0.75)
  for(j in 1:length(p)){
    for(k in 1:length(typeOfSampling)){
      index <- (j - 1) * length(typeOfSampling) + k
      
      # Sample from graph
      h <- getSubGraph(g = g,
                       p = p[j],
                       V1 = data$V1,
                       V2 = data$V2,
                       type = typeOfSampling[k])
      
      # Make data frame
      dataH <- as.data.frame(x = as_edgelist(h),
                             stringsAsFactors = FALSE)
      
      # Which edges from G are in H?
      for(i in 1:nrow(dataH)){
        result2[data$V1 %in% dataH[i, ] & data$V2 %in% dataH[i, ],index] <- 1
      }
      
      # Get MST
      mstH <- as_edgelist(mst(h))
      
      # Get MST info
      dataH$mstH <- 0
      if(nrow(mstH) > 0){
        for(i in 1:nrow(mstH)){
          # Which edges in G are in the MST of H?
          result[data$V1 %in% mstH[i, ] & data$V2 %in% mstH[i, ],index] <- 1
          result2[data$V1 %in% mstH[i, ] & data$V2 %in% mstH[i, ],index] <- 0
          
          # Which edges in H are in the MST of H?
          dataH$mstH[dataH$V1 %in% mstH[i, ] & dataH$V2 %in% mstH[i, ]] <- 1
        }
      }
      
      # Merge
      dataH <- merge(x = dataH,
                     y = data,
                     all.x = TRUE,
                     sort = FALSE)
      
      # Perform bootstrap
      temp <- replicate(n = 1e2,
                        expr = randMST(g = h,
                                       p = p[j],
                                       V1 = dataH$V1,
                                       V2 = dataH$V2,
                                       type = typeOfSampling[k]))
      # temp2 is "in the sample but not in the sample MST"
      temp2 <- temp[(nrow(temp) / 2 + 1):nrow(temp), ]
      # temp is "in the sample MST"
      temp <- temp[1:(nrow(temp) / 2), ]
      if(dim(temp)[1] > 1){
        b[index] <- mean(x = colSums(temp * dataH$mstH) / colSums(temp),
                         na.rm = TRUE)
      } else {
        b[index] <- mean(temp * dataH$mstH / temp)
      }
      if(mean(dataH$mstG) == 0 | mean(dataH$mstG) == 1){
        AUCs[index] <- NA
      } else {
        AUCs[index] <- auc(response = dataH$mstG,
                           predictor = rowSums(temp),
                           direction = "<")
      }
      
      if(dim(temp2)[1] > 1){
        b2[index] <- mean(x = colSums(temp2 * (!dataH$mstH)) / colSums(temp2),
                         na.rm = TRUE)
      } else {
        b2[index] <- mean(temp2 * (!dataH$mstH) / temp2)
      }
      if(mean(dataH$mstG) == 0 | mean(dataH$mstG) == 1){
        AUCs2[index] <- NA
      } else {
        AUCs2[index] <- auc(response = 1 - dataH$mstG,
                            predictor = rowSums(temp2),
                            direction = "<")
      }
    }
  }
  
  return(c(colSums(data$mstG * result) / colSums(result),
           b,
           AUCs,
           dd,
           colSums((!data$mstG) * result2) / colSums(result2),
           b2,
           AUCs2))
}

#' Make table
#' 
#' This function takes the results of running many iterations
#' of the sim function and produces a table in LaTeX code.
#' 
#' @param result A matrix where each row is the result of
#' a single iteration of the function sim.
#' @param n A vector of three elements denoting
#' the number of vertices sampled.
#' For the simulated graphs, the number of vertices in g is always 100,
#' and the proportions of sampled nodes are always (0.25,0.50,0.75),
#' so this is always (25,50,75).
#' @param symbol This is always "theta".
#' In a previous version of the manuscript, another Greek letter
#' was used to denote the NPV as opposed to the PPV,
#' but that was cut.
makeTable <- function(result,
                      n,
                      symbol = "theta"){
  mVec <- colMeans(x = result,
                   na.rm = TRUE)
  sVec <- apply(X = result,
                MARGIN = 2,
                FUN = sd,
                na.rm = TRUE)
  lVec <- mVec - qnorm(0.975) * sVec / sqrt(NROW(result))
  uVec <- mVec + qnorm(0.975) * sVec / sqrt(NROW(result))
  mVec <- format(x = mVec,
                 digits = 3)
  lVec <- format(x = lVec,
                 digits = 3)
  uVec <- format(x = uVec,
                 digits = 3)
  if(ncol(result) == 36){
    cat(paste0("& $",n[1],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[1]," (",lVec[1],"-",uVec[1],") & ",
               mVec[3]," (",lVec[3],"-",uVec[3],") & ",
               mVec[2]," (",lVec[2],"-",uVec[2],") & ",
               mVec[4]," (",lVec[4],"-",uVec[4],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[13]," (",lVec[13],"-",uVec[13],") & ",
               mVec[15]," (",lVec[15],"-",uVec[15],") & ",
               mVec[14]," (",lVec[14],"-",uVec[14],") & ",
               mVec[16]," (",lVec[16],"-",uVec[16],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[25]," (",lVec[25],"-",uVec[25],") & ",
               mVec[27]," (",lVec[27],"-",uVec[27],") & ",
               mVec[26]," (",lVec[26],"-",uVec[26],") & ",
               mVec[28]," (",lVec[28],"-",uVec[28],") \\\\\n"))
    cat(paste0("& $",n[2],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[5]," (",lVec[5],"-",uVec[5],") & ",
               mVec[7]," (",lVec[7],"-",uVec[7],") & ",
               mVec[6]," (",lVec[6],"-",uVec[6],") & ",
               mVec[8]," (",lVec[8],"-",uVec[8],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[17]," (",lVec[17],"-",uVec[17],") & ",
               mVec[19]," (",lVec[19],"-",uVec[19],") & ",
               mVec[18]," (",lVec[18],"-",uVec[18],") & ",
               mVec[20]," (",lVec[20],"-",uVec[20],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[29]," (",lVec[29],"-",uVec[29],") & ",
               mVec[31]," (",lVec[31],"-",uVec[31],") & ",
               mVec[30]," (",lVec[30],"-",uVec[30],") & ",
               mVec[32]," (",lVec[32],"-",uVec[32],") \\\\\n"))
    cat(paste0("& $",n[3],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[9]," (",lVec[9],"-",uVec[9],") & ",
               mVec[11]," (",lVec[11],"-",uVec[11],") & ",
               mVec[10]," (",lVec[10],"-",uVec[10],") & ",
               mVec[12]," (",lVec[12],"-",uVec[12],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[21]," (",lVec[21],"-",uVec[21],") & ",
               mVec[23]," (",lVec[23],"-",uVec[23],") & ",
               mVec[22]," (",lVec[22],"-",uVec[22],") & ",
               mVec[24]," (",lVec[24],"-",uVec[24],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[33]," (",lVec[33],"-",uVec[33],") & ",
               mVec[35]," (",lVec[35],"-",uVec[35],") & ",
               mVec[34]," (",lVec[34],"-",uVec[34],") & ",
               mVec[36]," (",lVec[36],"-",uVec[36],") \\\\\n"))
  } else if(ncol(result) == 45){
    cat(paste0("& $",n[1],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[1]," (",lVec[1],"-",uVec[1],") & ",
               mVec[3]," (",lVec[3],"-",uVec[3],") & ",
               mVec[2]," (",lVec[2],"-",uVec[2],") & ",
               mVec[4]," (",lVec[4],"-",uVec[4],") & ",
               mVec[5]," (",lVec[5],"-",uVec[5],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[16]," (",lVec[16],"-",uVec[16],") & ",
               mVec[18]," (",lVec[18],"-",uVec[18],") & ",
               mVec[17]," (",lVec[17],"-",uVec[17],") & ",
               mVec[19]," (",lVec[19],"-",uVec[19],") & ",
               mVec[20]," (",lVec[20],"-",uVec[20],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[31]," (",lVec[31],"-",uVec[31],") & ",
               mVec[33]," (",lVec[33],"-",uVec[33],") & ",
               mVec[32]," (",lVec[32],"-",uVec[32],") & ",
               mVec[34]," (",lVec[34],"-",uVec[34],") & ",
               mVec[35]," (",lVec[35],"-",uVec[35],") \\\\\n"))
    cat(paste0("& $",n[2],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[6]," (",lVec[6],"-",uVec[6],") & ",
               mVec[8]," (",lVec[8],"-",uVec[8],") & ",
               mVec[7]," (",lVec[7],"-",uVec[7],") & ",
               mVec[9]," (",lVec[9],"-",uVec[9],") & ",
               mVec[10]," (",lVec[10],"-",uVec[10],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[21]," (",lVec[21],"-",uVec[21],") & ",
               mVec[23]," (",lVec[23],"-",uVec[23],") & ",
               mVec[22]," (",lVec[22],"-",uVec[22],") & ",
               mVec[24]," (",lVec[24],"-",uVec[24],") & ",
               mVec[25]," (",lVec[25],"-",uVec[25],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[36]," (",lVec[36],"-",uVec[36],") & ",
               mVec[38]," (",lVec[38],"-",uVec[38],") & ",
               mVec[37]," (",lVec[37],"-",uVec[37],") & ",
               mVec[39]," (",lVec[39],"-",uVec[39],") & ",
               mVec[40]," (",lVec[40],"-",uVec[40],") \\\\\n"))
    cat(paste0("& $",n[3],"$ & $\\bar{\\hat{\\",symbol,"}}$ & ",
               mVec[11]," (",lVec[11],"-",uVec[11],") & ",
               mVec[13]," (",lVec[13],"-",uVec[13],") & ",
               mVec[12]," (",lVec[12],"-",uVec[12],") & ",
               mVec[14]," (",lVec[14],"-",uVec[14],") & ",
               mVec[15]," (",lVec[15],"-",uVec[15],") \\\\\n"))
    cat(paste0("& & $\\bar{\\hat{\\",symbol,"}}^{(B)}$ & ",
               mVec[26]," (",lVec[26],"-",uVec[26],") & ",
               mVec[28]," (",lVec[28],"-",uVec[28],") & ",
               mVec[27]," (",lVec[27],"-",uVec[27],") & ",
               mVec[29]," (",lVec[29],"-",uVec[29],") & ",
               mVec[30]," (",lVec[30],"-",uVec[30],") \\\\\n"))
    cat(paste0("& & $\\overline{AUC}$ & ",
               mVec[41]," (",lVec[41],"-",uVec[41],") & ",
               mVec[43]," (",lVec[43],"-",uVec[43],") & ",
               mVec[42]," (",lVec[42],"-",uVec[42],") & ",
               mVec[44]," (",lVec[44],"-",uVec[44],") & ",
               mVec[45]," (",lVec[45],"-",uVec[45],") \\\\\n"))
  }
}

#' Get PPV
#' 
#' This function takes a graph g and a vector v of vertices in the graph.
#' It finds the subgraph h induced by selecting only the vertices in v.
#' It returns the proportion of edges in the MST of h that are also in
#' the MST of g.
#' This function is used for the HIV data set.
#' 
#' @param g A graph.
#' @param v A vector of vertices in g to sample.
#' 
#' @return The PPV.
getProp <- function(g,v){
  
  # Make data frame
  data <- as.data.frame(x = as_edgelist(g),
                        stringsAsFactors = FALSE)
  
  # Get MST
  mstG <- as_edgelist(mst(g))
  
  # Add MST info to data frame
  data$mstG <- 0
  for(i in 1:nrow(mstG)){
    data$mstG[data$V1 %in% mstG[i, ] & data$V2 %in% mstG[i, ]] <- 1
  }
  
  # Get subgraph
  h <- induced_subgraph(graph = g,
                        vids = v)
  
  # Get MST
  mstH <- as_edgelist(mst(h))
  
  # Get MST info
  data$mstH <- 0
  if(nrow(mstH) > 0){
    for(i in 1:nrow(mstH)){
      data$mstH[data$V1 %in% mstH[i, ] & data$V2 %in% mstH[i, ]] <- 1
    }
  }
  return(sum(data$mstG * data$mstH) / sum(data$mstH))
}

#' Get mean
#' 
#' Computes the expected number of unique values when sampling M times
#' with replacement from a set of $m$ values.
#' This is from a previous draft of the manuscript.
#' 
#' @param m The number of values you're drawing from.
#' @param M The number of draws (with replacement).
#' 
#' @return Mean.
getMean <- function(m,M){
  result <- M * ((m - 1) / m)^(M - 1)
  return(result)
}

#' Get SD
#' 
#' Computes the SD of the number of unique values when sampling M times
#' with replacement from a set of $m$ values.
#' This is from a previous draft of the manuscript.
#' 
#' @param m The number of values you're drawing from.
#' @param M The number of draws (with replacement).
#' 
#' @return SD.
getSD <- function(m,M){
  term1 <- getMean(m,M)
  term2 <- M * (M - 1) * (m - 1) / m * ((m - 2) / m)^(M - 2)
  term3 <- M^2 * ((m - 1) / m)^(2 * (M - 1))
  result <- sqrt(term1 + term2 - term3)
  return(result)
}

#' Randomize edge order, by weight
#' 
#' This function takes the edge weights of a graph and assigns each
#' edge a rank.
#' The smallest edge is assigned rank 1, the second-smallest edge
#' is assigned rank 2, and so forth.
#' If any edges have the same weight, the ranking is randomized
#' among those edges.
#' If, for example, the edge weights are (2,1,2),
#' the function is equally likely to return (3,1,2) and (2,1,3).
#' 
#' @param edgeWeights The edge weights of a graph.
#' 
#' @return Ranks of the edges, with 1 indicated the smallest edge weight.
#' Ranks for edges with the same weight are randomized.
randEdgeOrder <- function(edgeWeights){
  d <- as.data.frame(x = table(edgeWeights),
                     stringsAsFactors = FALSE)
  class(d$edgeWeights) <- class(edgeWeights)
  d$cFreq <- cumsum(d$Freq)
  newVec <- rep(x = NA,
                times = length(edgeWeights))
  lower <- 1
  for(i in 1:length(d$cFreq)){
    if(d$Freq[i] == 1){
      r <- d$cFreq[i]
    } else {
      r <- sample(x = lower:d$cFreq[i],
                  size = d$Freq[i])
    }
    newVec[edgeWeights == d$edgeWeights[i]] <- r
    lower <- d$cFreq[i] + 1
  }
  return(newVec)
}

#' Randomize edge order, get MST
#' 
#' This function takes a graph g and a vector of vertices v.
#' Then it:
#' (1) Uses randEdgeOrder to randomize the ordering of any edges
#' with identical weight.
#' (2) Gets the subgraph h of g induced from the vertices in v.
#' (3) Finds the MST of g and the MST of h.
#' (4) Returns information abut which edges of g are in the MSTs
#' of g and h.
#' 
#' @param g A graph.
#' @param v Vertices of g to use to create a subgraph h.
#' @param V1 Vector of vertex names.
#' Should be in order such that there is
#' an edge between the ith element of V1
#' and the ith element of V2.
#' In other words,
#' the edges in g are given by rows of cbind(V1,V2).
#' E(g) must be in the same order as cbind(V1,V2).
#' @param V2 
#' 
#' @return A vector twice the length of V1.
#' In the first half, a 1 indicates that the corresponding edge of g
#' is in the MST of h.
#' In the second half, a 1 indicates that the corresponding edge of g
#' is in the MST of g.
#' All other elements are 0.
randMSTOrd <- function(g,v,V1,V2){
  E(g)$weight <- randEdgeOrder(E(g)$weight)
  h <- induced_subgraph(graph = g,
                        vids = v)
  mstG <- as_edgelist(mst(g))
  mstH <- as_edgelist(mst(h))
  x <- rep(x = 0,
           times = length(V1))
  y <- x
  for(i in 1:nrow(mstH)){
    x[V1 %in% mstH[i, ] & V2 %in% mstH[i, ]] <- 1
  }
  for(i in 1:nrow(mstG)){
    y[V1 %in% mstG[i, ] & V2 %in% mstG[i, ]] <- 1
  }
  return(c(x,y))
}

#' Get regression slope
#' 
#' This function is used with the HIV data
#' when sampling by 3-digit zip code prefix.
#' It takes a graph g and a vector of vertices v.
#' It runs randMSTOrd 100 times on g and returns the slope of a
#' simple linear regression of the number of population MSTs
#' an edge is in on the number of sample MSTs an edge is in.
#' I don't think this was used in the manuscript.
#' 
#' @param g A graph.
#' @param v Vertices of g to use to create a subgraph h.
#' @param V1 Vector of vertex names.
#' Should be in order such that there is
#' an edge between the ith element of V1
#' and the ith element of V2.
#' In other words,
#' the edges in g are given by rows of cbind(V1,V2).
#' E(g) must be in the same order as cbind(V1,V2).
#' @param V2 
#' 
#' @return Slope of simple linear regression.
getB <- function(g,
                 v,
                 V1,
                 V2){
  x <- rowSums(replicate(n = 1e2,
                         expr = randMSTOrd(g = g,
                                           v = v,
                                           V1 = V1,
                                           V2 = V2)))
  y <- x[(length(V1) + 1):length(x)]
  x <- x[1:length(V1)]
  b <- sum(x * y) / sum(x^2)
  
  return(b)
}

################
## Simulation ##
################

# Set size of original graph
N <- 1e2

# Set number of replications
nrep <- 1e3

# Register number of cores
registerDoParallel(cores = 10)

# Unique edge weights
resultPA <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = "PA")
resultER <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = "ER")
resultCo <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = "complete")
resultNo <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = "normal")
resultC1 <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = "complete",
                                                      d = "1")

###################
## Distance Data ##
###################

# Load data
distances <- read.csv(file = "all.csv",
                      stringsAsFactors = FALSE)

# Standardize spelling
distances$ID1 <- gsub(pattern = "plamsa",
                      replacement = "plasma",
                      x = distances$ID1)
distances$ID1 <- gsub(pattern = "Virologic",
                      replacement = "ViroLogic",
                      x = distances$ID1)
distances$ID2 <- gsub(pattern = "plamsa",
                      replacement = "plasma",
                      x = distances$ID2)
distances$ID2 <- gsub(pattern = "Virologic",
                      replacement = "ViroLogic",
                      x = distances$ID2)

# Check for rows where ID1 == ID2
same <- distances$ID1 == distances$ID2
sum(same)
rm(same)

# Check for duplicate rows
check <- do.call(what = paste0,
                 args = distances[ ,1:2])
length(check) == length(unique(check))
rm(check)

# Count unique elements
length(unique(distances$ID1))
length(unique(distances$ID2))
nUniqueID <- length(unique(c(distances$ID1,
                             distances$ID2)))
nUniqueID
choose(nUniqueID,2)
nrow(distances)
rm(nUniqueID)

# Create key
idStrings <- unique(c(distances$ID1,
                      distances$ID2))
idNumber <- 1:length(idStrings)
key <- data.frame(idStrings,
                  idNumber)
rm(idStrings,
   idNumber)

# Change strings to numbers
distances <- merge(x = distances,
                   y = key,
                   by.x = "ID1",
                   by.y = "idStrings",
                   all.x = TRUE)
distances$ID1 <- distances$idNumber
distances$idNumber <- NULL
distances <- merge(x = distances,
                   y = key,
                   by.x = "ID2",
                   by.y = "idStrings",
                   all.x = TRUE)
distances$ID2 <- distances$idNumber
distances$idNumber <- NULL

# How many edges have zero distance?
zeros <- distances$Distance == 0
sum(zeros)

# Assign very small positive distance to edges with zero distance
distances$Distance[zeros] <- min(distances$Distance[distances$Distance > 0]) / 2

# How many edges have zero distance?
zeros <- distances$Distance == 0
sum(zeros)

# Clean up
rm(zeros)

# Sort
distances <- distances[order(distances$ID1,
                             distances$ID2), ]

# Change "Distance" to "weight"
names(distances)[names(distances) == "Distance"] <- "weight"

# Remove edges with weight > 1.5%
distances <- distances[distances$weight <= 0.015, ]

# Create graph
gRD <- graph_from_data_frame(d = distances,
                             directed = FALSE)

# Analyze data
resultRD <- foreach(i = 1:nrep,
                    .combine = rbind,
                    .packages = "igraph") %dopar% sim(typeOfGraph = gRD)
gRDFix <- gRD
E(gRDFix)$weight <- randEdgeOrder(E(gRDFix)$weight)
resultRDFix <- foreach(i = 1:nrep,
                       .combine = rbind,
                       .packages = "igraph") %dopar% sim(typeOfGraph = gRDFix)


# How many unique participants are in the network?
ids <- cSplit(indt = key,
              splitCols = "idStrings",
              sep = "|")
ids$idStrings_1 <- as.character(ids$idStrings_1)
ids <- ids[order(ids$idStrings_1), ]
length(unique(ids$idStrings_1))
length(unique(ids$idNumber))
length(unique(c(distances$ID1,
                distances$ID2)))
sum(ids$idNumber %in% c(distances$ID1,
                        distances$ID2))

###################
## Zip Code Data ##
###################

# Load data
zips <- read.csv(file = "zips_20200117.csv",
                 stringsAsFactors = FALSE)

# Summary
table(zips$zip3)
sum(is.na(zips$zip3))

# Merge zip code data
ids <- merge(x = ids,
             y = zips,
             by.x = "idStrings_1",
             by.y = "pid",
             all.x = TRUE)
ids <- ids[order(ids$idNumber), ]

# Remove ids that aren't in graph
ids <- ids[as.character(ids$idNumber) %in% V(gRD)$name, ]

# Remove ids that have missing zip code
ids <- ids[!is.na(ids$zip3), ]

# Make data frame
data <- as.data.frame(x = as_edgelist(gRD),
                      stringsAsFactors = FALSE)

# Get slopes
prod(factorial(table(E(gRD)$weight)))
sum(ids$zip3 == 919)
sum(ids$zip3 == 919) / gorder(gRD)
temp <- as.character(ids$idNumber[ids$zip3 == 919])
prod(factorial(table(E(induced_subgraph(graph = gRD,
                                        vids = temp))$weight)))
getB(g = gRD,
     v = temp,
     V1 = data$V1,
     V2 = data$V2)
getProp(g = gRD,
        v = temp)
sum(ids$zip3 == 920)
sum(ids$zip3 == 920) / gorder(gRD)
temp <- as.character(ids$idNumber[ids$zip3 == 920])
prod(factorial(table(E(induced_subgraph(graph = gRD,
                                        vids = temp))$weight)))
getB(g = gRD,
     v = temp,
     V1 = data$V1,
     V2 = data$V2)
getProp(g = gRD,
        v = temp)
sum(ids$zip3 == 921)
sum(ids$zip3 == 921) / gorder(gRD)
temp <- as.character(ids$idNumber[ids$zip3 == 921])
prod(factorial(table(E(induced_subgraph(graph = gRD,
                                        vids = temp))$weight)))
getB(g = gRD,
     v = temp,
     V1 = data$V1,
     V2 = data$V2)
getProp(g = gRD,
        v = temp)
sum(c(ids$zip3 == 919,
      ids$zip3 == 920,
      ids$zip3 == 921))
sum(c(ids$zip3 == 919,
      ids$zip3 == 920,
      ids$zip3 == 921)) / gorder(gRD)

#################
## Make tables ##
#################

makeTable(result = resultCo[ ,1:36],
          n = c(25,50,75))
makeTable(result = resultER[ ,1:36],
          n = c(25,50,75))
makeTable(result = resultNo[ ,1:45],
          n = c(25,50,75))
makeTable(result = resultPA[ ,1:36],
          n = c(25,50,75))

###############
## Make plot ##
###############

# Put data in matrix
dd <- matrix(data = NA,
             nrow = 4,
             ncol = N - 1)
temp <- colMeans(x = resultCo[ ,38:136])
dd[1, ] <- 1 - cumsum(temp)
temp <- colMeans(x = resultER[ ,38:136])
dd[2, ] <- 1 - cumsum(temp)
temp <- colMeans(x = resultPA[ ,38:136])
dd[3, ] <- 1 - cumsum(temp)
temp <- degree_distribution(graph = gRD)
if(length(temp) < N){
  temp <- c(temp,
            rep(x = 0,
                times = N - length(temp)))
}
dd[4, ] <- 1 - cumsum(temp[1:(N - 1)])

# Clean up
rm(temp)

# Only go until there are zeroes
dd[dd <= 0] <- NA

# Make plot
pdf(file = "draft06a.pdf",
    width = 7.5,
    height = 4.635)
plot(x = 1:ncol(dd),
     y = dd[1, ],
     type = "l",
     log = "xy",
     xlab = "k",
     ylab = "Proportion of Nodes with Degree > k",
     ylim = c(min(dd,
                  na.rm = TRUE),
              max(dd,
                  na.rm = TRUE)))
lines(x = 1:ncol(dd),
      y = dd[2, ],
      pch = 0,
      type = "b")
lines(x = 1:ncol(dd),
      y = dd[3, ],
      pch = 1,
      type = "b")
lines(x = 1:ncol(dd),
      y = dd[4, ],
      pch = 2,
      type = "b")
legend(x = "bottomleft",
       legend = c("Complete and Normal",
                  "G(n,p)",
                  "B-A",
                  "PIRC"),
       lty = 1,
       pch = c(NA,0:2))
dev.off()

# Plot of number of unique edge weights
x <- 1:1e4
y <- getMean(m = x,
             M = choose(100,2))
sds <- getSD(m = x,
             M = choose(100,2))
pdf(file = "draft06b.pdf",
    width = 7.5,
    height = 4.635)
plot(x = x,
     y = y,
     type = "l",
     xlab = "m",
     ylab = "E(X)")
points(x = x,
       y = y - 2 * sds,
       type = "l",
       lty = 2)
points(x = x,
       y = y + 2 * sds,
       type = "l",
       lty = 2)
dev.off()

# Duplicate edge weights
table(table(distances$weight))

# Plot displaying duplicate edge weights
tab <- table(distances$weight)
pdf(file = "draft06c.pdf",
    width = 7.5,
    height = 4.635)
plot(x = as.numeric(names(tab)),
     y = as.numeric(unname(tab)),
     log = "y",
     xlab = "Weight",
     ylab = "Number of Edges with This Weight")
dev.off()

# Resolution of edge weights
temp <- sort(distances$weight)
temp <- temp[-1] - temp[-length(temp)]
min(temp[temp > 0])
rm(temp)

# How many components?
cc <- count_components(gRD)
cd <- component_distribution(gRD)
cd * cc
hist(x = rep(x = 0:(length(cd) - 1),
             times = cd * cc))
pdf(file = "draft06d.pdf",
    width = 7.5,
    height = 4.635)
plot(x = 1:(length(cd) - 1),
     y = (cd * cc)[-1],
     log = "y",
     xlab = "Number of Nodes",
     ylab = "Number of Components")
dev.off()

##################
## Save results ##
##################

save.image("draft06.RData")










