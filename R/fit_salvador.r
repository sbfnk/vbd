## install_github("libbi/rbi")
library('rbi')
## install_github("sbfnk/rbi.helpers")
library('rbi.helpers')

library('magrittr')

###################
## model fitting ##
###################


obs <- readRDS("fit_data.rds")
model_dir <- path.expand("~/code/vbd/bi")

vbd_model <- bi_model(paste(model_dir, "vbd.bi", sep="/")) %>%
  fix(N=2.675e+6)

serology_sample <- data.frame(time=max(obs$Serology$time), value=633)

## fit

bi <- libbi(vbd_model, input=list(serology_sample=serology_sample),
            obs=obs,
            end_time=max(obs$Sero$time)) %>%
  sample(proposal="prior", nsamples=1000) %>%
  adapt_proposal(min=0.1, max=0.4) %>%
  sample(sample_obs=TRUE, nsamples=10000)

p <- plot(bi, date.origin=as.Date("2015-01-05") - 7)

save_results(bi, "salvador.rds")
