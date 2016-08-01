library('RBi')
library('RBi.helpers')
library('coda')
library('dplyr')
library('tidyr')
library('msm')
library('cowplot')

p_tau <- 7

output_dir <- path.expand("~/Data/Zika")
code_dir <- path.expand("~/code/vbd/")

models <- c("vbd_sero_fnh_yap_zika",
            "vbd_sero_fnh_earlier_yap_zika",
            "vbd_sero_patch_fnh_yap_zika",
            "vbd_sero_patch_fnh_earlier_yap_zika",
            "vbd_fnh_fais_dengue",
            "vbd_fnh_earlier_fais_dengue",
            "vbd_sero_fnh_yap_dengue",
            "vbd_sero_fnh_earlier_yap_dengue",
            "vbd_sero_patch_fnh_yap_dengue",
            "vbd_sero_patch_fnh_earlier_yap_dengue",
            "vbd_sero_fnh",
            "vbd_sero_fnh_earlier",
            "vbd_sero_patch_fnh",
            "vbd_sero_patch_fnh_earlier")

obs_id_levels <- c("yap_dengue", "fais_dengue", "yap_zika")

traces <- list()
bim <- list()
params <- list()

for (model in models)
{

    cat(date(), model, "\n")
    mcmc_pattern <- paste0(model, "_[0-9]+")
    m <- merge_parallel_runs(output_dir, mcmc_pattern, concatenate = TRUE)
    cat(date(), "merged\n")
    traces[[model]] <- m$traces
    bim[[model]] <- m$model
  ## l <- merge_parallel_runs(output_dir, mcmc_pattern, concatenate = FALSE)
  ## mc <- mcmc.list(lapply(l$traces, function(x) { mcmc(get_traces(x, l$model))}))

    ## p <- plot_libbi(m$traces, m$model, extra.aes = c(color = "setting", linetype = "disease"))

    if (length(m$traces) > 0)
    {

        cat(date(), "get parameters and states\n")
        l <- lapply(names(m$traces), function(x) {
            m$traces[[x]] %>% mutate(state = x)
        })

        found_disease <- sub("^.*(zika|dengue).*$", "\\1", model)
        if (found_disease != model) {
            l <- lapply(l, function(x) {
                if ("disease" %in% names(x))
                {
                    x %>% filter(disease == found_disease)
                } else
                {
                    x
                }
            })
        }

        found_setting <- sub("^.*(fais|yap).*$", "\\1", model)
        if (found_setting != model) {
            l <- lapply(l, function(x) {
                if ("setting" %in% names(x))
                {
                    x %>% filter(setting == found_setting)
                } else
                {
                    x
                }
            })
        }

        names(l) <- names(traces[[model]])
        params[[model]] <- bind_rows(l)
        if ("time" %in% names(params[[model]]))
        {
            params[[model]] <- params[[model]] %>%
                filter(is.na(time)) %>%
                select(-time)
        }

        params_disease <-
            lapply(levels(factor(params[[model]]$disease)),
                   function(x) { params[[model]] %>% filter(is.na(disease)) %>% mutate(disease = x) })

        params_all <-
            rbind(params[[model]] %>%
                  filter(!is.na(disease)),
                  bind_rows(params_disease))

        params_setting <-
            lapply(levels(factor(params[[model]]$setting)),
                   function(x) { params_all %>% filter(is.na(setting)) %>% mutate(setting = x) })

        params_all <-
            rbind(params_all %>%
                  filter(!is.na(setting)),
                  bind_rows(params_setting))

        if (grepl("earlier", model))
        {
            p_d_life_m <- 1
        } else
        {
            p_d_life_m <- 2
        }

        r0 <- params_all %>%
            spread(state, value) %>%
            mutate(R0 = (p_tau * p_d_life_m)**2 * 
                       p_b_h * p_b_m * 10**(p_lm) * p_d_inf_h /
                       (p_d_life_m + p_d_inc_m),
                   GI = (p_d_inc_h + p_d_inf_h + p_d_inc_m + p_d_life_m)) %>%
            gather(state, value, loglikelihood:GI) %>%
          filter(state %in% c("R0", "GI"))

        params[[model]] <- rbind(params[[model]], r0) %>%
            mutate(disease = ifelse(is.na(disease), "n/a", as.character(disease)),
                   setting = ifelse(is.na(setting), "n/a", as.character(setting)))

        states <- bind_rows(l) %>%
            filter(state == "Z_h") %>%
            group_by(time, disease, np, state, setting) %>%
            summarise(value = sum(value)) %>%
            ungroup() %>%
            spread(state, value)

        rep_params <- params_all %>%
            filter(state %in% c("p_rep")) %>%
            spread(state, value)

        for (rep_param in setdiff(c("p_rep"), colnames(rep_params)))
        {
            rep_params[[rep_param]] <- 1
        }

        if ("patch" %in% colnames(rep_params))
        {
            rep_params <- rep_params %>%
                select(-patch)
        }

        states <- states %>%
            left_join(rep_params, by = c("disease", "np", "setting")) %>%
            mutate(obs_id =
                       factor(paste(setting, disease, sep = "_"),
                              levels = obs_id_levels))  %>%
            filter(!is.na(obs_id))

        if (grepl("(yap|dengue)", model))
        {
            model_obs_id <- sub("^.*((yap|fais)_[^_]*)($|_.*$)", "\\1", model)
            states <- states %>%
                filter(obs_id == model_obs_id)
        }

        cat(date(), "sample\n")
        l$Cases <- states %>%
            mutate(Z_h = ifelse(Z_h < 0, 0, Z_h)) %>%
            mutate(mean = p_rep * Z_h,
                   sd = sqrt(p_rep * Z_h)) %>%
            mutate(sd = ifelse(sd < 1, 1, sd)) %>%
            mutate(value = rtnorm(n = nrow(states), lower = 0,
                                  mean = mean,
                                  sd = sd))

        traces[[model]] <- l

    }
}

dic <- c()
obs <- list()
joint_params <- list()

m1 <- merge(data.table(traces[["vbd_fnh_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_fnh_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")
m2 <- merge(m1, data.table(traces[["vbd_sero_fnh_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]
dic["separate"] <- compute_DIC(list(loglikelihood = m2))
obs[["separate"]] <- rbind(traces[["vbd_fnh_fais_dengue"]]$Cases,
                           traces[["vbd_sero_fnh_yap_dengue"]]$Cases,
                           traces[["vbd_sero_fnh_yap_zika"]]$Cases)
joint_params[["separate"]] <- rbind(params[["vbd_fnh_fais_dengue"]],
                                    params[["vbd_sero_fnh_yap_dengue"]],
                                    params[["vbd_sero_fnh_yap_zika"]])


m1 <- merge(data.table(traces[["vbd_fnh_earlier_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_fnh_earlier_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")
m2 <- merge(m1, data.table(traces[["vbd_sero_fnh_earlier_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]
dic["separate_earlier"] <- compute_DIC(list(loglikelihood = m2))
obs[["separate_earlier"]] <-
  rbind(traces[["vbd_fnh_earlier_fais_dengue"]]$Cases,
        traces[["vbd_sero_fnh_earlier_yap_dengue"]]$Cases,
        traces[["vbd_sero_fnh_earlier_yap_zika"]]$Cases)
joint_params[["separate_earlier"]] <-
  rbind(params[["vbd_fnh_earlier_fais_dengue"]],
        params[["vbd_sero_fnh_earlier_yap_dengue"]],
        params[["vbd_sero_fnh_earlier_yap_zika"]])



m1 <- merge(data.table(traces[["vbd_fnh_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_patch_fnh_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")
m2 <- merge(m1, data.table(traces[["vbd_sero_patch_fnh_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]
dic["separate_patch"] <- compute_DIC(list(loglikelihood = m2))
obs[["separate_patch"]] <-
  rbind(traces[["vbd_fnh_fais_dengue"]]$Cases,
        traces[["vbd_sero_patch_fnh_yap_dengue"]]$Cases,
        traces[["vbd_sero_patch_fnh_yap_zika"]]$Cases)
joint_params[["separate_patch"]] <-
  rbind(params[["vbd_fnh_fais_dengue"]],
        params[["vbd_sero_patch_fnh_yap_dengue"]],
        params[["vbd_sero_patch_fnh_yap_zika"]])

m1 <- merge(data.table(traces[["vbd_fnh_earlier_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_patch_fnh_earlier_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")
m2 <- merge(m1, data.table(traces[["vbd_sero_patch_fnh_earlier_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]
dic["separate_patch_earlier"] <- compute_DIC(list(loglikelihood = m2))
obs[["separate_patch_earlier"]] <-
  rbind(traces[["vbd_fnh_earlier_fais_dengue"]]$Cases,
        traces[["vbd_sero_patch_fnh_earlier_yap_dengue"]]$Cases,
        traces[["vbd_sero_patch_fnh_earlier_yap_zika"]]$Cases)
joint_params[["separate_patch_earlier"]] <-
  rbind(params[["vbd_fnh_earlier_fais_dengue"]],
        params[["vbd_sero_patch_fnh_earlier_yap_dengue"]],
        params[["vbd_sero_patch_fnh_earlier_yap_zika"]])


for (model in grep(paste0("(", paste(obs_id_levels, collapse = "|"), ")"), models, value = TRUE, invert = TRUE))
{
    dic[model] <- compute_DIC(traces[[model]])
    obs[[model]] <- traces[[model]]$Cases
    joint_params[[model]] <- params[[model]]
}

## read data
ts <- list()
analyses <- data.frame(setting = c("yap", "yap", "fais"), disease = c("dengue", "zika", "dengue"))

for (i in 1:nrow(analyses))
{
    this_setting <- analyses[i, "setting"]
    this_disease <- analyses[i, "disease"]
    this_filename <-
      paste(code_dir, "data",
            paste(this_setting, this_disease, "data.rds", sep = "_"),
            sep = "/")
    this_ts <- readRDS(this_filename) %>%
      mutate(setting = this_setting, disease = this_disease,
             week = ceiling(nr / 7))
    ts <- c(ts, list(this_ts))
}

dt_ts <- bind_rows(ts) %>%
    group_by(week, setting, disease) %>%
    summarize(value = sum(value)) %>%
    ungroup() %>%
    mutate(obs_id = factor(paste(setting, disease, sep = "_"),
                           levels = c("yap_dengue", "fais_dengue", "yap_zika"))) %>%
    arrange(week, obs_id) %>%
    select(week, obs_id, value) ## %>%
    ## complete(week, obs_id, fill = list(value = 0))

data <- dt_ts %>%
  mutate(time = week, state = "Cases")
first_obs <- data %>%
  filter(value > 0) %>%
  slice(which.min(time)) %>%
  .[["time"]]
last_obs <- data %>%
  filter(value > 0) %>%
  slice(which.max(time)) %>%
  .[["time"]]
data <- data %>%
  filter(time >= first_obs & time <= last_obs)

p_obs <- list()
p_r0 <- list()
p_ll <- list()
for (model in names(obs))
{
  model_name <- "vbd_sero"
  if (grepl("patch", model)) {
    model_name <- paste(model_name, "patch", sep = "_")
  }
  model_name <- paste(model_name, "fnh", sep = "_")
  if (grepl("earlier", model)) {
    model_name <- paste(model_name, "earlier", sep = "_")
  }

  temp_plot <- plot_libbi(read = list(Cases = obs[[model]]),
                          model = bim[[model_name]], 
                          data = data %>% filter(value > 0),
                          density_args = list(adjust = 2),
                          extra.aes = list(color = "obs_id"),
                          states = "Cases", trend = "mean", plot = FALSE,
                      limit.to.data = TRUE)
  p_obs[[model]] <- temp_plot$states + facet_wrap(~ obs_id, scales = "free")
  p_ll[[model]] <- temp_plot$likelihoods 
  p_r0[[model]] <-
    ggplot(joint_params[[model]] %>% filter(state == "R0"),
           aes(x = value, color = disease, linetype = setting)) +
    geom_line(stat = "density", adjust = 2, lwd = 2) +
    scale_color_brewer(palette = "Set1")
}

## zika patch estimated parameters
params[["vbd_sero_patch_fnh_yap_zika"]] %>%
  group_by(state) %>%
  summarise(value = median(value))

