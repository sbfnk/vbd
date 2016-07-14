library('RBi')
library('RBi.helpers')
library('coda')
library('dplyr')

output_dir <- path.expand("~/Data/Zika")

mcmc_files <- list.files(output_dir, pattern = "vbd_poisson_sero_patch_fnh_[0-9]+\\.rds", full.names = TRUE)

traces <- list()

for (file in mcmc_files)
{
    traces[[file]] <- readRDS(file)
}

np_translate <- NULL
combined <- lapply(names(traces[[1]]), function(x) {
  z <- rbindlist(lapply(seq_along(traces), function(y) {
    data.table(traces[[y]][[x]])[, unique_np := paste(np, y, sep = "_")]
  }))
  if (is.null(np_translate)) {
    np_translate <-
      data.table(unique_np = unique(z$unique_np),
                 new_np = seq_along(unique(z$unique_np)) - 1)
  }
  z <- merge(z, np_translate, by = "unique_np")
  z[, np := NULL]
  z[, unique_np := NULL]
  setnames(z, "new_np", "np")
  setkey(z, np)
  z
})


model <- bi_model("~/code/vbd/bi/vbd_patch.bi")
names(combined) <- names(traces[[1]])

res <- combined

p <- plot_libbi(combined, model, extra.aes = c(color = "disease", linetype = "setting"))
