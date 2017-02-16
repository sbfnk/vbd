## install_github("sbfnk/rbi.helpers")
library('rbi.helpers')
library('magrittr')
library('cowplot')

###################
## model fitting ##
###################

## read observation data
obs <- readRDS("fit_data.rds")

## set model directory
model_dir <- path.expand("~/code/vbd/bi")

## sample size for serology
serology_sample <- data.frame(time=0, value=633)

## load model and fix population size
vbd_model <- bi_model(paste(model_dir, "vbd.bi", sep="/")) %>%
  fix(N=2.675e+6, p_p_immune=0.06, p_p_risk=1)

## fit
bi <- libbi(vbd_model, input=list(serology_sample=serology_sample),
            obs=obs,
            end_time=max(obs$Sero$time))

bi_prior <- sample(bi, target="prior", nsamples=10000)

bi %<>%
  optimise() %>%
  sample(proposal="prior", nsamples=1000) %>%
  adapt_proposal(min=0.1, max=0.4) %>%
  sample(sample_obs=TRUE, nsamples=500000, thin=50)

## plot
p <- plot(bi, prior=bi_prior, date.origin=as.Date("2015-01-05") - 7, date.unit="week", obs=c("Serology", "Incidence"), state=c("beta_track"), verbose=TRUE, type=c("obs", "param", "logeval", "state"))
model_name <- "poisson_over_fullN"

## save
save_libbi(bi, paste0("salvador_", model_name, ".rds"))

ggsave(paste0("salvador_", model_name, "_states.pdf"), p$trajectories)
ggsave(paste0("salvador_", model_name, "_traces.pdf"), p$traces)
ggsave(paste0("salvador_", model_name, "_densities.pdf"), p$densities)
ggsave(paste0("salvador_", model_name, "_pairs.pdf"), p$pairs)
ggsave(paste0("salvador_", model_name, "_correlations.pdf"), p$correlations)
ggsave(paste0("salvador_", model_name, "_logevals.pdf"), p$logevals)
