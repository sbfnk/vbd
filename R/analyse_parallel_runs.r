library('RBi')
library('RBi.helpers')
library('coda')
library('dplyr')

traces <- list()
Z_traj <- list()
for (i in 1:20)
{
    filename <-
        paste("~/Data/Zika/vbd_beta_sero", i, sep = "_")
    res_filename <- paste(filename, "rds", sep = ".")
    model_filename <- paste(filename, "bi", sep = ".")
    if (file.exists(res_filename) && file.exists(model_filename)) 
    {
      res <- readRDS(res_filename)
      model <- bi_model(model_filename)
      traces[[i]] <- mcmc(get_traces(res, model = model))
      Z_traj[[i]] <- res[["Z_h"]]
    }
}

res <- list(Z_h = bind_rows(Z_traj))

plot_libbi(res, model, extra.aes = c(color = "disease", linetype = "setting"))

mcmc <- mcmc.list(traces[!sapply(traces, is.null)])
