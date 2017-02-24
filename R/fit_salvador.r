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

pred <- predict(bi, end_time=80, noutputs=80, sample_obs=TRUE)

common_plot_options <-
  list(x=pred,
       model=vbd_model,
       all.times=TRUE,
       select=list(time=1:104),
       hline=c(`Reff`=1),
       labels=c(beta_track="beta", Reff="R[eff]"),
       date.origin=as.Date("2015-01-05") - 7, date.unit="week")

## calculate R0 at the beginning of the outbreak
R0_calc <- bi_read(pred)[c("beta_track", "p_d_inf_h")]
R0_calc <- lapply(names(R0_calc), function(x) {setnames(R0_calc[[x]], "value", x)})
R0_df <- data.table(merge(R0_calc[[1]], R0_calc[[2]]))
R0_df <- R0_df[time == 0]
R0_df <- R0_df[, value := beta_track * p_d_inf_h]
res$

## plot figure
p <- list()
p[["A"]] <- do.call(plot_libbi, c(common_plot_options, list(type="obs", obs="Incidence")))$trajectories + scale_y_continuous("Weekly incidence") + scale_x_date("", date_labels="%b %Y")
p[["B"]] <- do.call(plot_libbi, c(common_plot_options, list(type="obs", obs="Serology")))$trajectories + scale_y_continuous("Percent immune", labels=percent) + scale_x_date("", date_labels="%b %Y")
p[["C"]] <- do.call(plot_libbi, c(common_plot_options, list(type="state", state="beta_track")))$trajectories + scale_y_continuous("Transmission rate") + scale_x_date("", date_labels="%b %Y")
p[["D"]] <- do.call(plot_libbi, c(common_plot_options, list(type="state", state="Reff")))$trajectories + scale_y_continuous("Reproduction number") + scale_x_date("", date_labels="%b %Y")

plot <- do.call(plot_grid, c(p, list(labels=names(p), ncol=2)))

p <- plot(pred, prior=bi_prior, type=c("param", "logeval"))

save_libbi(pred, "salvador_prediction.rds")
