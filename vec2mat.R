vec2mat <- function (vec, x.loc, y.loc)
{
  ux <- sort(unique(x.loc))
  uy <- sort(unique(y.loc))
  nx <- length(ux)
  ny <- length(uy)
  row.ref <- match(x.loc, ux)
  col.ref <- match(y.loc, uy)
  Vec <- rep(NA, max(row.ref) * max(col.ref))
  vec.ref <- (col.ref - 1) * max(row.ref) + row.ref
  Vec[vec.ref] <- vec
  grid.mask <- matrix(Vec, max(row.ref), max(col.ref), dimnames = list(ux,
                                                                       uy))
  return(grid.mask)
}