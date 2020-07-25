topo.error <- function(somobj, type = c("nodedist", "bmu"), data) {
 # if (somobj$method != "kohonen")
#    stop("Topographic error measures as yet only defined for class 'som'")
  
  type = match.arg(type)
  dmat <- switch(type,
                 nodedist = {
                   as.matrix(dist(somobj$codes[[1]]))
                 },
                 bmu = {
                   if (missing(data) & !is.null(somobj$data[[1]]))
                     data <- somobj$data[[1]]
                   
                   ## following code adapted from map.kohonen, may be put in a
                   ## separate auxiliary function lateron
                   nd <- nrow(data)
                   ncodes <- nrow(somobj$codes[[1]])
                   np <- ncol(somobj$codes[[1]])
                   distances <- .C("mapKohonen",
                                   as.double(data),
                                   as.double(somobj$codes[[1]]), 
                                   as.integer(ncodes),
                                   as.integer(nd),
                                   as.integer(np),
                                   dists = double(nd * ncodes),
                                   NAOK = TRUE, PACKAGE = "kohonen")$dists
                   
                   matrix(distances, nd, ncodes, byrow = TRUE)
                 })
  
  ## here we assume that there are no cases where more than one zero
  ## exists, something that would be _very_ strange
  dmat.ordered <- t(apply(dmat, 1, order))
  
  dists.2D <- unit.distances(somobj$grid, somobj$grid$toroidal)
  mean(mapply(function(i, j) dists.2D[i,j],
              dmat.ordered[,1], dmat.ordered[,2]))
}

