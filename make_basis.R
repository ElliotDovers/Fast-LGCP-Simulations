make_basis <- function (nodes.on.long.edge, data, radius.type = c("diag", "limiting"), 
          coord.names = c("x", "y")) 
{
  if (!all(coord.names %in% colnames(data))) {
    stop("at least one of 'coord.names' not found in the data provided")
  }
  radius.type <- match.arg(radius.type)
  xrange <- range(data[, coord.names[1]], na.rm = T)
  yrange <- range(data[, coord.names[2]], na.rm = T)
  dx <- diff(xrange)
  dy <- diff(yrange)
  big.axis.id <- which.max(c(dx, dy))
  small.axis.id <- c(1, 2)[!c(1, 2) %in% big.axis.id]
  big.axis.jitter <- c(dx, dy)[big.axis.id]/(nodes.on.long.edge * 
                                               4)
  big.axis.nodes <- seq(c(xrange[1], yrange[1])[big.axis.id] - 
                          big.axis.jitter, c(xrange[2], yrange[2])[big.axis.id] + 
                          big.axis.jitter, length.out = nodes.on.long.edge)
  small.axis.nodes <- seq(c(xrange[1], yrange[1])[small.axis.id] - 
                            big.axis.jitter, c(xrange[2], yrange[2])[small.axis.id] + 
                            big.axis.jitter, by = diff(big.axis.nodes)[1])
  node.list <- list()
  node.list[[big.axis.id]] <- big.axis.nodes
  node.list[[small.axis.id]] <- small.axis.nodes
  bf.info <- as.data.frame(expand.grid(node.list[[1]], node.list[[2]]))
  colnames(bf.info) <- coord.names
  bf.info$scale <- switch(radius.type, diag = diff(big.axis.nodes)[1] * 
                            sqrt(2), limiting = sqrt(dx * dy)/log(nrow(bf.info)))
  bf.info$res <- 1
  return(bf.info)
}