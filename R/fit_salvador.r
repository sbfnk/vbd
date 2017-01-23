## install_github("libbi/rbi")
library('rbi')
## install_github("sbfnk/rbi.helpers")
library('rbi.helpers')

library('magrittr')

library('cowplot')

###################
## model fitting ##
###################


obs <- readRDS("fit_data.rds")
model_dir <- path.expand("~/code/vbd/bi")

vbd_model <- bi_model(paste(model_dir, "vbd.bi", sep="/")) %>%
  fix(N=2.675e+6,
      p_p_risk=1)

serology_sample <- data.frame(time=0, value=633)

## fit with serology

bi <- libbi(vbd_model, input=list(serology_sample=serology_sample),
            obs=obs,
            end_time=max(obs$Sero$time)) %>%
  optimise() %>%
  sample(proposal="prior", nsamples=1000) %>%
  adapt_proposal(min=0.1, max=0.4) %>%
  sample(sample_obs=TRUE, nsamples=100000, thin=10)

p <- plot(bi, date.origin=as.Date("2015-01-05") - 7, date.unit="week", states=c("Serology", "Incidence"), verbose=TRUE)
ggsave("sero.pdf", p$states)

save_results(bi, "salvador.rds")

