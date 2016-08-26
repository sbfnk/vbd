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

models <-
    sub("_0\\.bi$", "", list.files(path = output_dir, pattern = ".*_0\\.bi"))
models <-
    intersect(models,
              c("vbd_fnh", "vbd_fnh_earlier",
                "vbd_sero_fnh", "vbd_sero_fnh_earlier",
                "vbd", "vbd_sero", "vbd_sero_earlier"))

obs_id_levels <- c("yap_dengue", "fais_dengue", "yap_zika")

traces <- list()
params <- list()
priors <- list()
bim <- list()

for (model in models)
{

    cat(date(), model, "\n")
    traces[[model]] <- list(prior = list(), posterior = list())
    params[[model]] <- list(prior = list(), posterior = list())
    posterior_pattern <- paste0(model, "_[0-9]+")
    posterior <-
        merge_parallel_runs(output_dir, posterior_pattern, concatenate = TRUE)
    cat(date(), "merged\n")
    traces[[model]][["posterior"]] <- posterior$traces
    bim[[model]] <- posterior$model
    prior_pattern <- paste0(model, "_[0-9]+_prior")
    prior <- merge_parallel_runs(output_dir, prior_pattern, concatenate = TRUE)
    traces[[model]][["prior"]] <- prior$traces

    if (length(traces[[model]][["posterior"]]) > 0)
    {
        for (dist in c("prior", "posterior"))
        {
            cat(date(), "get", dist, "parameters and states\n")
            l <- lapply(names(traces[[model]][[dist]]), function(x) {
              traces[[model]][[dist]][[x]] %>%
                mutate(state = x)
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

            names(l) <- names(traces[[model]][[dist]])
            params[[model]][[dist]] <- bind_rows(l)
            if ("time" %in% names(params[[model]][[dist]]))
            {
                params[[model]][[dist]] <- params[[model]][[dist]] %>%
                    filter(is.na(time)) %>%
                    select(-time)
            }

            params_disease <-
                lapply(levels(factor(params[[model]][[dist]]$disease)),
                       function(x) { params[[model]][[dist]] %>%
                                         filter(is.na(disease)) %>%
                                         mutate(disease = x) })

            params_all <-
                rbind(params[[model]][[dist]] %>%
                      filter(!is.na(disease)),
                      bind_rows(params_disease))

            params_setting <-
                lapply(levels(factor(params[[model]][[dist]]$setting)),
                       function(x) {
                           params_all %>%
                               filter(is.na(setting)) %>%
                               mutate(setting = x)
                       })

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
                       GI = (p_d_inc_h + p_d_inf_h +
                             p_d_inc_m + p_d_life_m))
            first_col <-
                min(which(!(colnames(r0) %in% c("disease", "setting",
                                                "patch", "np", "run"))))
            r0 <- r0 %>% gather(state, value, first_col:ncol(.)) %>%
                filter(state %in% c("R0", "GI"))

            l[["R0"]] <- r0 %>% filter(state == "R0")
            l[["GI"]] <- r0 %>% filter(state == "GI")

            params[[model]][[dist]] <- rbind(params[[model]][[dist]], r0) %>%
                mutate(disease = ifelse(is.na(disease), "n/a",
                                        as.character(disease)),
                       setting = ifelse(is.na(setting), "n/a",
                                        as.character(setting)))

            l$final_size <- l$Z_h %>%
              group_by(disease, setting, np, time) %>%
              summarise(value = sum(value))
            if (dist == "prior")
            {
              l$final_size <- l$final_size %>%
                  group_by(disease, setting, np) %>%
                  summarise(value = last(value))
            } else
            {
              l$final_size <- l$final_size %>%
                  group_by(disease, setting, np) %>%
                  summarise(value = sum(value))
            }
           l$proportion_infected <- l$final_size %>%
               left_join(l$p_N_h %>%
                         rename(N = value) %>%
                         select(-state)) %>%
               left_join(l$p_initial_susceptible_yap %>%
                         rename(init = value) %>%
                         mutate(setting = factor("yap", levels = levels(l$final_size$setting))) %>%
                         select(-state)) %>%
               mutate(init = ifelse(is.na(init), 1, init),
                      value = value / (N * init),
                      state = "proportion_infected") %>%
              select(-N, -init)
            l$final_size <- l$final_size %>%
              left_join(l$p_N_h %>%
                          rename(N = value) %>%
                          select(-state)) %>%
              mutate(value = value / N,
                     state = "final_size") %>%
              select(-N)

            if (dist == "posterior")
            {
                cat(date(), "observations\n")

                o <- lapply(names(posterior$obs), function(x) {
                    posterior$obs[[x]] %>%
                        mutate(obs_id =
                                   factor(obs_id, levels = obs_id_levels))  %>%
                        mutate(obs = x)
                })
                obs <- bind_rows(o) %>%
                    mutate(state = ifelse(obs == "Cases", "Z_h",
                                   ifelse(obs == "Sero", "final_size",
                                          NA_character_))) %>%
                    rename(time = week, data = value)

                for (state in c("Z_h", "final_size"))
                {
                    l[[state]] <- l[[state]] %>%
                        mutate(obs_id =
                                   factor(paste(setting, disease, sep = "_"),
                                          levels = obs_id_levels))
                }

                if ("Sero" %in% obs$obs)
                {
                    sero_time <- obs %>%
                        filter(obs == "Sero") %>%
                        select(time, obs_id)
                    l[["final_size"]] <- l[["final_size"]] %>%
                        left_join(sero_time)
                }

                states <- bind_rows(l) %>%
                    filter(state %in% c("Z_h", "final_size")) %>%
                    filter(!is.na(time)) %>% 
                    group_by(time, disease, np, state, setting, obs_id) %>%
                    summarise(value = sum(value)) %>%
                    ungroup() %>%
                    spread(state, value)

                rep_params <- params_all %>%
                    filter(state %in% c("p_rep", "p_N_h")) %>%
                    spread(state, value)

                for (rep_param in setdiff(c("p_rep", "p_N_h"),
                                          colnames(rep_params)))
                {
                    rep_params[[rep_param]] <- 1
                }

                if ("patch" %in% colnames(rep_params))
                {
                    rep_params <- rep_params %>%
                        select(-patch)
                }

                states <- states %>%
                    left_join(rep_params,
                              by = c("disease", "np", "setting")) %>%
                    filter(!is.na(obs_id))

                if (grepl("(yap|dengue)", model))
                {
                    model_obs_id <-
                        sub("^.*((yap|fais)_[^_]*)($|_.*$)", "\\1", model)
                    states <- states %>%
                        filter(obs_id == model_obs_id)
                }

                first_column <-
                    min(which(colnames(states) %in% c("final_size", "Z_h")))
                last_column <-
                    max(which(colnames(states) %in% c("final_size", "Z_h")))
                 obsdens <- states %>%
                    gather(state, value, first_column:last_column) %>%
                    left_join(obs, by = c("time", "obs_id", "state")) %>%
                    filter(!is.na(data)) %>%
                    mutate(value = ifelse(value < 0, 0, value)) %>%
                    mutate(mean = ifelse(state == "Z_h", p_rep * value,
                                  ifelse(state == "final_size", value,
                                         NA_real_)),
                           sd = ifelse(state == "Z_h",
                                ifelse(p_rep * value < 1, 1, sqrt(p_rep * value)),
                                ifelse(state == "final_size", 0.09 / 3.98, NA_real_))) %>%
                    mutate(sample = rtnorm(n = n(), lower = 0, mean = mean, sd = sd),
                           density = dtnorm(x = data, mean = mean, sd = sd,
                                            lower = 0, log = TRUE))

                l$Cases <- obsdens %>%
                    select(-value) %>%
                    rename(value = sample)

                l$pointll <- obsdens %>%
                    select(-value) %>%
                    rename(value = density)
            }
            traces[[model]][[dist]] <- l
        }
    }
}

## saveRDS(traces, "vbd_traces.rds")

dic <- c()
waic <- c()
obs <- list()
joint_params <- list()

## gather model strings
model_strings <- rep("", length(obs_id_levels))
names(model_strings) <- obs_id_levels
separate_models <- list()
for (earlier_string in c("", "_earlier"))
{
    model_strings["fais_dengue"] <-
        paste0("vbd_fnh", earlier_string, "_fais_dengue")
    for (patch_string in c("", "_patch"))
    {
        model_strings["yap_zika"] <-
            paste0("vbd_sero", patch_string, "_fnh", earlier_string,
                   "_yap_zika")
        move_strings <- ""
        if (nchar(patch_string) > 0) move_strings <- c(move_strings, "_move")
        for (move_string in move_strings)
        {
            reverse_strings <- ""
            if (nchar(move_string) > 0)
            {
              reverse_strings <- c(reverse_strings, "_reverse")
            }
            for (reverse_string in reverse_strings)
            {
                model_strings["yap_dengue"] <-
                    paste0("vbd", move_string, patch_string, reverse_string,
                           "_fnh", earlier_string, "_yap_dengue")
                full_model_string <-
                    paste0("vbd_sero", move_string, patch_string,
                           reverse_string, "_fnh", earlier_string)
                separate_models[[full_model_string]] <- model_strings
            }
        }
    }
}

available_models <-
    separate_models[sapply(separate_models, function(x) {all(x %in% models)})]

for (model_name in names(available_models))
{
    temp_ll <- list()
    model <- available_models[[model_name]]
    for (ll in c("loglikelihood", "pointll"))
    {
        by_string <-
            intersect(c("np", "time"),
                      colnames(traces[[model["fais_dengue"]]][["posterior"]][[ll]]))
        temp_traces <- list()
        temp_traces[["fd"]] <-
            data.table(traces[[model["fais_dengue"]]][[ll]])
        temp_traces[["yd"]] <-
            data.table(traces[[model["yap_dengue"]]][[ll]])
        temp_traces[["yz"]] <-
            data.table(traces[[model["yap_zika"]]][[ll]])

        for (name in names(temp_traces))
        {
            setnames(temp_traces[[name]], "value", name)
        }

        temp_ll[[ll]] <-
            merge(temp_traces[["fd"]], temp_traces[["yd"]], by = by_string)
        temp_ll[[ll]] <-
            merge(temp_ll[[ll]], temp_traces[["yz"]], by = by_string)
        temp_ll[[ll]][, value := fd+yd+yz]
    }
    sep_name <- paste(model_name, "separate", sep = "_")
    dic[[sep_name]] <- compute_DIC(temp_ll["loglikelihood"])
    mean_lik <-
        temp_ll[["pointll"]][, list(mean = mean(exp(value))), by = time]
    lppd <- mean_lik[, sum(log(mean))]
    var_ll <-
        temp_ll[["pointll"]][, list(var = var(value)), by = time]
    pwaic2 <- var_ll[, sum(var)]
    waic[sep_name] <- -2 * (lppd - pwaic2)
    obs[[sep_name]] <- do.call(rbind, lapply(traces[model], function(x) {
        x[["Cases"]]
    }))
    joint_params[[sep_name]] <- do.call(rbind, params[model])
}

for (model in grep(paste0("(", paste(obs_id_levels, collapse = "|"), ")"), models, value = TRUE, invert = TRUE))
{
    dic[model] <- compute_DIC(traces[[model]])
    mean_lik <- data.table(traces[[model]][["posterior"]][["pointll"]]) [, list(mean = mean(exp(value))), by = time]
    lppd <- mean_lik[, sum(log(mean))]
    var_ll <-
        data.table(traces[[model]][["posterior"]][["pointll"]])[, list(var = var(value)), by = time]
    pwaic2 <- var_ll[, sum(var)]
    waic[model] <- -2 * (lppd - pwaic2)
    obs[[model]] <- traces[[model]][["posterior"]]$Cases
    joint_params[[model]] <- params[[model]]
}

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
p_r0_sqrt <- list()
p_libbi <- list()
p_prior <- list()
for (model in names(obs))
{
  model_name <- "vbd_sero"
  if (grepl("patch", model)) {
    model_name <- paste(model_name, "patch", sep = "_")
  }
  if (grepl("stoch", model)) {
    model_name <- paste(model_name, "stoch", sep = "_")
  }
  if (grepl("fnh", model)) {
    model_name <- paste(model_name, "fnh", sep = "_")
  }
  if (grepl("earlier", model)) {
    model_name <- paste(model_name, "earlier", sep = "_")
  }

  temp_plot <- plot_libbi(read = list(Cases = obs[[model]]),
                          model = bim[[model_name]],
                          data = data %>% filter(value > 0),
                          density_args = list(adjust = 2),
                          ## densities = "histogram", 
                          extra.aes = list(color = "obs_id"),
                          states = "Cases", trend = "mean", plot = FALSE,
                          limit.to.data = TRUE,
                          quantiles = c(0.5, 0.72, 0.95))
  p_obs[[model]] <- temp_plot$states + facet_wrap(~ obs_id, scales = "free")
  if (model %in% names(traces))
  {
      p_libbi[[model]] <- list()
      for (type in names(traces[[model]]))
      {
          p_libbi[[model]][[type]] <-
              plot_libbi(read = traces[[model]][[type]],
                         prior = traces[[model]][["prior"]], 
                         model = bim[[model_name]],
                         ## density_args = list(bins = 20), 
                         ## densities = "histogram", 
                         density_args = list(adjust = 2, alpha = 0.5),
                         extra.aes = list(color = "disease",
                                          linetype = "setting"),
                         trend = "mean", plot = FALSE,
                         quantiles = c(0.5, 0.72, 0.95))
      }
  }
  p_r0[[model]] <-
      ggplot(traces[[model]][["posterior"]][["R0"]], 
             aes(x = value, color = disease, linetype = setting)) +
      geom_line(stat = "density", adjust = 2, lwd = 2) +
      scale_color_brewer(palette = "Set1")
  p_r0_sqrt[[model]] <-
      ggplot(traces[[model]][["posterior"]][["R0"]], 
             aes(x = sqrt(value), color = disease, linetype = setting)) +
      geom_line(stat = "density", adjust = 2, lwd = 2) +
      scale_color_brewer(palette = "Set1")
}

## plot R0 vs final size
prior_R0 <- lapply(traces, function(x) {x[["prior"]][["R0"]]})
prior_R0 <- lapply(names(prior_R0), function(x) {
    prior_R0[[x]] %>%
        mutate(model = x) %>%
        select(-patch) %>%
        select(-state) %>%
        rename(R0 = value)
})
prior_R0 <- bind_rows(prior_R0)

prior_GI <- lapply(traces, function(x) {x[["prior"]][["GI"]]})
prior_GI <- lapply(names(prior_GI), function(x) {
    prior_GI[[x]] %>%
        mutate(model = x) %>%
        select(-patch) %>%
        select(-state) %>%
        rename(GI = value)
})
prior_GI <- bind_rows(prior_GI)

prior_proportion_infected <- lapply(traces, function(x) {x[["prior"]][["proportion_infected"]]})
prior_proportion_infected <- lapply(names(prior_proportion_infected), function(x) {
    prior_proportion_infected[[x]] %>%
        mutate(model = x) %>%
        select(-state) %>%
        rename(proportion_infected = value)
})
prior_proportion_infected <- bind_rows(prior_proportion_infected)

final_sizes <- prior_final_size %>%
    left_join(prior_R0) %>%
    left_join(prior_GI) %>%
    left_join(prior_proportion_infected)

p <- ggplot(final_sizes %>% filter(R0 < 20 & proportion_infected <= 1),
            aes(x = R0, y = proportion_infected)) +
    geom_jitter()

## zika patch estimated parameters
params[["vbd_sero_patch_fnh_yap_zika"]] %>%
  group_by(state) %>%
  summarise(value = median(value))


all_params <- rbind(params[["vbd_fnh"]][["posterior"]] %>% spread(state, value) %>% mutate(model = "vbd_fnh"), params[["vbd_fnh_earlier"]][["posterior"]] %>% spread(state, value) %>% mutate(model = "vbd_fnh_earlier")) %>%
    filter(!(disease == "n/a" | setting == "n/a")) %>%
    mutate(data = paste(disease, setting, sep = "_"))

p <- ggplot(all_params %>% filter(setting == "fais" & disease == "dengue" & model == "vbd_sero_fnh"), aes(x = GI, y = R0)) +
    geom_jitter() +
    scale_color_brewer(palette = "Dark2")

all_params <- params[["vbd_sero"]][["posterior"]] %>% spread(state, value) %>%
    filter(!(disease == "n/a" | setting == "n/a")) %>%
    mutate(data = paste(disease, setting, sep = "_"))

p <- ggplot(all_params, aes(x = GI, y = R0, color = data)) +
    geom_jitter() +
    scale_color_brewer(palette = "Dark2")
