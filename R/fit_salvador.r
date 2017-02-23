library('rbi.helpers')
library('magrittr')
library('cowplot')
library('scales')

###################
## model fitting ##
###################

## number of samples in serology
n_serology=633

## read observation data
obs <- readRDS("fit_data.rds")

## set model directory
model_dir <- path.expand("~/code/vbd/bi")

## load model and fix population size
vbd_model <- bi_model(paste(model_dir, "vbd.bi", sep="/")) %>%
  fix(N=2.675e+6, p_p_immune=0.06, p_p_risk=1)

## fit
bi <- libbi(vbd_model,
            input=list(serology_sample=n_serology),
            obs=obs,
            end_time=max(obs$Sero$time))

bi_prior <- sample(bi, target="prior", nsamples=100)

bi %<>%
  optimise() %>%
  sample(proposal="prior", nsamples=1000) %>%
  adapt_proposal(min=0.1, max=0.4) %>%
  sample(sample_obs=TRUE, nsamples=5000, thin=5)

save_libbi(bi, "salvador.rds")

## bi <- read_libbi("~/Research/Analysis/Zika/salvador_poisson_over_fullN.rds")
## bi$model[48] <- "Z <- 0"
pred <- predict(bi, end_time=104, noutputs=104, sample_obs=TRUE, verbose=TRUE)

res <- bi_read(pred)

## calc R_eff
## turn Serology into %

common_plot_options <-
  list(x=pred,
       all.times=TRUE,
       select=list(time=1:104),
       hline=c(`R[eff]`=1),
       labels=c(beta_track="beta", Reff="R[eff]"),
       date.origin=as.Date("2015-01-05") - 7, date.unit="week")

p <- list()
p[["A"]] <- do.call(plot_libbi, c(common_plot_options, list(type="obs", obs="Incidence")))$trajectories + ggtitle("Incidence")
p[["B"]] <- do.call(plot_libbi, c(common_plot_options, list(type="obs", obs="Serology")))$trajectories + ggtitle("Serology")
p[["C"]] <- do.call(plot_libbi, c(common_plot_options, list(type="state", state="beta_track")))$trajectories + ggtitle("Seasonal variation in transmission")
p[["D"]] <- do.call(plot_libbi, c(common_plot_options, list(type="state", state="Reff")))$trajectories + ggtitle("Reproduction number")

plot <- do.call(plot_grid, c(p, list(labels=names(p), ncol=2)))

p <- plot(pred, prior=bi_prior, type=c("param", "logeval"))

save_libbi(pred, "salvador_prediction.rds")
