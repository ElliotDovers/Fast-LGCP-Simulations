# Function to convert a vector to a spatstat image object via locations ######
vec2im <- function(vec, x.loc, y.loc){
  ux <- sort(unique(x.loc))
  uy <- sort(unique(y.loc))
  nx <- length(ux)
  ny <- length(uy)
  col.ref <- match(x.loc, ux)
  row.ref <- match(y.loc, uy)
  Vec <- rep(NA, max(row.ref)*max(col.ref))
  vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
  Vec[vec.ref] <- vec
  grid.mask <- matrix(Vec, max(row.ref), max(col.ref),dimnames = list(uy, ux))
  vec.im <- spatstat.geom::im(grid.mask, xcol = ux, yrow = uy)

  return(vec.im)
}