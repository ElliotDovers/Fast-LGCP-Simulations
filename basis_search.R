basis_search <- function (sim.list, approx_type = c("VA", "Lp")) {
  # function for obtaining the (approx.) marginal log-likelihood
  source("get_mll.R")
  
  approx_type <- match.arg(approx_type)
  
  # initialise storage vectors
  sqrt_node_number <- NULL
  lls <- NULL
  # hard code the first two searches so we can perform a while loop
  search_result <- get_mll(sim.list, 7, approx_type)
  if (is.na(search_result["mll"])) {
    return(7)
    stop("basis search failed at start")
  }
  search_result <- rbind(search_result, get_mll(sim.list, 14, approx_type))
  if (is.na(search_result[nrow(search_result), "mll"])) {
    return(7)
    stop("basis search failed beyond start")
  }
  # create the lcv
  lcv <- 14 + 7
  # perform the search loop
  while ((search_result[nrow(search_result), "mll"] >= search_result[nrow(search_result) - 1, "mll"]) &
         (search_result[nrow(search_result), "sqrt_nodes"]^2 < (0.25*sum(sim.list$dframe$pres))) &
         (!is.na(search_result[nrow(search_result), "mll"]))) {
    search_result <- rbind(search_result, get_mll(sim.list, lcv, approx_type))
    lcv <- lcv + 7
  }
  return(search_result)
}