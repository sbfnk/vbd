library('rbi.helpers')
library('magrittr')
library('cowplot')
library('scales')
library('dplyr')
library('data.table')

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

bi_prior <- sample(bi, target="prior", nsamples=5000)
save_libbi(bi_prior, "salvador_prior.rds")

bi %<>%
  optimise() %>%
  sample(proposal="prior", nsamples=5000) %>%
  adapt_proposal(min=0.1, max=0.4) %>%
  sample(sample_obs=TRUE, nsamples=500000, thin=50)
save_libbi(bi, "salvador_posterior.rds")
date()

pred <- predict(bi, end_time=80, noutputs=80, sample_obs=TRUE)
save_libbi(pred, "salvador_prediction.rds")

res <- bi_read(pred)

## calculate R0 at the beginning of the outbreak
R0_calc <- copy(res[c("beta_track", "p_d_inf_h")])
R0_calc <- lapply(names(R0_calc), function(x) {setnames(R0_calc[[x]], "value", x)})
R0_df <- data.table(merge(R0_calc[[1]], R0_calc[[2]]))
R0_df <- R0_df[time == 0]
R0_df <- R0_df[, value := beta_track * p_d_inf_h]

mean(R0_df$value)
quantile(R0_df$value, c(0.025, 0.975))
mean(res$p_p_rep$value)
quantile(res$p_p_rep$value, c(0.025, 0.975))

res$Serology %<>%
  mutate(value=value/n_serology)
pred %<>% attach_file(file="output", data=res, force=TRUE)

pred_obs <- bi_read(pred, file="obs")
pred_obs$Serology %<>% mutate(value=value/n_serology)
pred %<>% attach_file(file="obs", data=pred_obs, force=TRUE)

common_plot_options <-
  list(x=pred,
       plot=FALSE,
       all.times=TRUE,
       select=list(time=1:104),
       hline=c(`Reff`=1),
       labels=c(beta_track="beta", Reff="R[eff]"),
       date.origin=as.Date("2015-01-05") - 7, date.unit="week")

## plot figure
p <- list()
p[["A"]] <- do.call(plot, c(common_plot_options, list(type="obs", obs="Incidence")))$trajectories + scale_y_continuous("Weekly incidence", labels=comma) + scale_x_date("", date_labels="%b %Y")
p[["B"]] <- do.call(plot, c(common_plot_options, list(type="obs", obs="Serology")))$trajectories + scale_y_continuous("Percent immune", labels=percent) + scale_x_date("", date_labels="%b %Y")
p[["C"]] <- do.call(plot, c(common_plot_options, list(type="state", state="beta_track")))$trajectories + scale_y_continuous("Transmission rate") + scale_x_date("", date_labels="%b %Y")
p[["D"]] <- do.call(plot, c(common_plot_options, list(type="state", state="Reff")))$trajectories + scale_y_continuous("Reproduction number") + scale_x_date("", date_labels="%b %Y")

pg <- do.call(plot_grid, c(p, list(labels=names(p), ncol=2)))
ggsave("salvador_trajectories.pdf", pg)

p <- plot(pred,
          prior=bi_prior,
          type=c("param", "logeval"),
          labels=c(p_d_inc_h="D[inc,h]",
                   p_d_inf_h="D[inf,h]",
                   p_R0="R[0]",
                   p_p_rep="r",
                   p_p_over="phi",
                   p_s_peak="Peak~week",
                   p_s_amp="Peak~amplitude",
                   initI="I[0]"
                   ),
          pairs=FALSE, 
          plot=FALSE)
ggsave("salvador_densities.pdf", p$densities)
ggsave("salvador_traces.pdf", p$traces)
ggsave("salvador_correlations.pdf", p$correlations)
