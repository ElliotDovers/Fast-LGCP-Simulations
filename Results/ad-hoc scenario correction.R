# ad-hoc scenario correction
setwd("C:/Users/iamel/Desktop/Work/Post Doc/LGCP/HPC Simulations/Results/SS_N_1000_exp")

res.list <- list.files()
for (i in res.list) {
  load(i)
  sim.res$scenario <- "SS_N_1000_exp"
  save(list = "sim.res", file = i)
}