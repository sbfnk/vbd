library('RBi')
library('RBi.helpers')
library('dplyr')
library('tidyr')
library('msm')
library('cowplot')
library('ggExtra')
library('stringi')

p_tau <- 7 ## fixed biting rate, 1 per day (7 per week)

## set directory names for libbi output files and code
output_dir <- path.expand("~/Data/Zika")
code_dir <- path.expand("~/code/vbd/")

models <- sub("\\.cmd$", "", list.files(output_dir, "^.*\\.cmd$"))

obs_id_levels <- c("yap_dengue", "fais_dengue", "yap_zika")
## ordered in time
ordered_obs_id_levels <- c("yap_zika", "fais_dengue", "yap_dengue")

data_labels <- ordered_obs_id_levels
data_labels <- sub("^(.*)_(.*)$", "\\2 \\1", data_labels)
data_labels <- sub(" ", " in ", stri_trans_totitle(data_labels))
names(data_labels) <- ordered_obs_id_levels

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
             week = floor(nr / 7))
    ts <- c(ts, list(this_ts))
}

dt_ts <- bind_rows(ts) %>%
    group_by(week, setting, disease) %>%
    summarize(value = sum(value), onset_date = min(onset_date)) %>%
    ungroup() %>%
    mutate(obs_id = factor(paste(setting, disease, sep = "_"),
                           levels = c("yap_dengue", "fais_dengue", "yap_zika"))) %>%
    arrange(week, obs_id) %>%
    select(week, obs_id, value, onset_date) ## %>%
    ## complete(week, obs_id, fill = list(value = 0))

peak <- dt_ts %>%
    group_by(obs_id) %>%
    slice(which.max(value)) %>%
    ungroup() %>%
    select(peak_week = week, obs_id) 

growth <- dt_ts %>%
    left_join(peak) %>%
    filter(week <= peak_week & value > 0)

## estimate initial growth rate
r <- c()

for (id in levels(growth$obs_id))
{
    inc <- growth %>%
        filter(obs_id == id) %>%
        .$value
    time <- seq_along(inc)
    fit <- nls(inc ~ a * exp(time * b), start = c(a = 5, b = 0.1))
    new_el <- coef(fit)[["b"]]
    names(new_el) <- id
    r <- c(r, new_el)
}

r_tb <- tibble(obs_id = names(r), r = r) %>%
    mutate(setting = sub("_.*$", "", obs_id),
           disease = sub("^.*_", "", obs_id)) %>%
    select(-obs_id)

traces <- list()
params <- list()
priors <- list()
bi_models <- list()

for (model in models)
{

    cat(date(), model, "\n")
    traces[[model]] <- list(prior = list(), posterior = list())
    params[[model]] <- list(prior = list(), posterior = list())
    posterior <- readRDS(paste0(output_dir, "/", model, ".rds"))
    obs <- readRDS(paste0(output_dir, "/", model, "_obs.rds"))
    traces[[model]][["posterior"]] <- posterior
    bi_models[[model]] <- bi_model(paste0(output_dir, "/", model, ".bi"))
    prior <- readRDS(paste0(output_dir, "/", model, "_prior.rds"))
    traces[[model]][["prior"]] <- prior

    if (length(traces[[model]][["posterior"]]) > 0)
    {
        for (dist in c("prior", "posterior"))
        {
            cat(date(), "get", dist, "parameters and states\n")
            state_param_names <- setdiff(names(traces[[model]][[dist]]), "Cases")
            l <- lapply(state_param_names,
                        function(x)
                        {
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

            names(l) <- state_param_names
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

            littler <- r_tb %>%
                mutate(setting = factor(setting,
                                        levels = levels(params_all$setting)),
                       disease = factor(disease,
                                        levels = levels(params_all$disease)))

            r0 <- params_all %>%
                left_join(littler) %>%
                spread(state, value) %>%
                mutate(R0 = p_tau**2 *
                           p_b_h * p_b_m * 10**(p_lm) * p_d_inf_h  * p_d_life_m * p_d_life_m /
                           (p_d_life_m + p_d_inc_m),
                       GI = p_d_inc_h + p_d_inc_m + p_d_inf_h + p_d_life_m)

            first_col <-
                min(which(!(colnames(r0) %in% c("disease", "setting",
                                                "patch", "np"))))
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

                o <- lapply(names(obs), function(x) {
                  obs[[x]] %>%
                        mutate(obs_id =
                                 factor(obs_id, levels = obs_id_levels))  %>%
                        mutate(obs = x)
                })
                obs_states <- bind_rows(o) %>%
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

                if ("Sero" %in% obs_states$obs)
                {
                    sero_time <- obs_states %>%
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
                   left_join(obs_states,
                             by = c("time", "obs_id", "state")) %>%
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
        ## not fitted
        traces[[model]][["posterior"]][["p_t_start"]] <-
          traces[[model]][["posterior"]][["p_t_start"]] %>%
          filter(!(disease == "zika" & setting == "fais"))
    }
}

data <- dt_ts %>%
    rename(time = week) %>%
    mutate(state = "Cases")
first_obs <- data %>%
    group_by(obs_id) %>%
    filter(value > 0) %>%
    slice(which.min(time)) %>%
    select(time, obs_id) %>%
    rename(first_obs = time)
last_obs <- data %>%
    group_by(obs_id) %>%
    filter(value > 0) %>%
    slice(which.max(time)) %>%
    select(time, obs_id) %>%
    rename(last_obs = time)
data <- data %>%
    left_join(first_obs, by = "obs_id") %>%
    left_join(last_obs, by = "obs_id") %>%
    filter(time >= first_obs & time <= last_obs) %>%
    mutate(obs_id = factor(obs_id, levels = ordered_obs_id_levels,
                           labels = data_labels))

dic <- c()
waic <- c()
for (model in models)
{
    dic[model] <- compute_DIC(traces[[model]])
    mean_lik <- data.table(traces[[model]][["posterior"]][["pointll"]]) [, list(mean = mean(exp(value))), by = time]
    lppd <- mean_lik[, sum(log(mean))]
    var_ll <-
        data.table(traces[[model]][["posterior"]][["pointll"]])[, list(var = var(value)), by = time]
    pwaic2 <- var_ll[, sum(var)]
    waic[model] <- -2 * (lppd - pwaic2)
}

## plots
p_obs <- list()
p_r0gi <- list()
p_r0 <- list()
p_r0_sqrt <- list()
p_libbi <- list()
p_prior <- list()

labels <- c(p_d_inc_h = "italic(D)[plain(inc,H)]",
            p_d_inc_m = "italic(D)[plain(inc,M)]",
            p_d_inf_h = "italic(D)[plain(inf,H)]",
            p_lm = "log[10](italic(m))",
            p_initial_susceptible_yap = "italic(q)",
            p_rep = "italic(r)",
            p_b_h = "italic(b)[H]",
            p_b_m = "italic(b)[M]",
            p_t_start = "italic(t[0])",
            R0 = "italic(R)[H %->% H]",
            GI = "italic(G)",
            zika = "Zika",
            yap = "Yap",
            fais = "Fais")

for (model in models)
{
  max_r0 <- max(traces[[model]][["posterior"]][["R0"]]$value)
  traces[[model]][["prior"]][["R0"]] <- traces[[model]][["prior"]][["R0"]] %>%
    filter(value < max_r0)
  obs <- traces[[model]][["posterior"]][["Cases"]] %>%
      mutate(obs_id = factor(obs_id, levels = ordered_obs_id_levels,
                             labels = data_labels))
  temp_plot <-
    plot_libbi(read = list(Cases = obs),
               model = bi_models[[model]],
               data = data %>% filter(value > 0),
               density_args = list(adjust = 2),
               ## densities = "histogram", 
               extra.aes = list(group = "obs_id"),
               data.colour = "black", 
               states = "Cases", trend = "mean", plot = FALSE,
               limit.to.data = TRUE,
               quantiles = c(0.5, 0.72, 0.95))
  obs_states <- temp_plot$data$states %>%
      inner_join(data %>% select(time, obs_id, onset_date), by = c("time", "obs_id"))
  p_obs[[model]] <- ggplot(obs_states, aes(x = onset_date)) +
      geom_point(data = data, mapping = aes(y = value)) +
      facet_wrap(~ obs_id, scales = "free") +
      scale_x_date("Week", labels = scales::date_format("%e %b %Y")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~ obs_id, scales = "free") +
      scale_y_continuous("Disease incidence") +
      geom_line(aes(y = value)) +
      geom_ribbon(aes(ymin = min.1, ymax = max.1), alpha = 0.5) +
      geom_ribbon(aes(ymin = min.2, ymax = max.2), alpha = 0.25) +
      geom_ribbon(aes(ymin = min.3, ymax = max.3), alpha = 0.125)

  if (model %in% names(traces))
  {
      p_libbi[[model]] <- list()
      for (type in names(traces[[model]]))
      {
        ## rename
        param_names <- names(traces[[model]][[type]])

        p_libbi[[model]][[type]] <-
              plot_libbi(read = traces[[model]][[type]],
                         prior = traces[[model]][["prior"]], 
                         model = bi_models[[model]],
                         ## density_args = list(bins = 20), 
                         ## densities = "histogram", 
                         density_args = list(adjust = 2, alpha = 0.5),
                         extra.aes = list(color = "disease",
                                          linetype = "setting"),
                         trend = "mean", plot = FALSE,
                         quantiles = c(0.5, 0.95),
                         labels = labels, brewer.palette = "Set1")
      }
  }
  p_r0gi[[model]] <-
    plot_libbi(read = traces[[model]][[type]],
               prior = traces[[model]][["prior"]], 
               model = bi_models[[model]],
               ## density_args = list(bins = 20), 
               ## densities = "histogram", 
               density_args = list(adjust = 2, alpha = 0.5),
               extra.aes = list(color = "disease",
                                linetype = "setting"),
               trend = "mean", plot = FALSE,
               quantiles = c(0.5, 0.95),
               labels = labels,
               states = c(), 
               params = c("R0", "GI"),
               noises = c(),
               brewer.palette = "Set1")
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

ggsave("posterior_densities.pdf",
       p_libbi[["vbd_fnh"]][["posterior"]]$densities + scale_x_continuous(""),
       width = 7.5, height = 6.5)
ggsave("posterior_r0gi_densities.pdf",
       p_r0gi[["vbd_fnh"]]$densities + scale_x_continuous(""),
       width = 7, height = 4)

ggsave("GI_R0.pdf", p, width = 7, height = 2.3)

two_panels <- plot_grid(p_r0gi[["vbd_fnh"]]$densities + xlab(""), p,
                        ncol = 2, labels = c("A", "B"))

ggsave("dengue_zika_obs.pdf", p_obs[["vbd_fnh"]], width = 7, height = 3)

param_estimates <- list()
r0_estimates <- list()
r0_gi <- list()
r0_gi_summary <- list()
all_params <- list()
p_r0vgi <- list()

for (model in models)
{
  param_estimates[[model]] <-
    p_libbi[[model]][["posterior"]]$data$params %>%
    group_by(distribution, disease, setting, parameter) %>%
    summarise(mean = mean(value),
              median = median(value),
              min.1 = quantile(value, 0.25),
              max.1 = quantile(value, 0.75),
              min.2 = quantile(value, 0.025),
              max.2 = quantile(value, 0.975))

    r0_estimates[[model]] <- p_r0gi[[model]]$data$params %>%
      filter(distribution == "posterior" &
             parameter == "italic(R)[H %->% H]") %>% 
      group_by(disease, setting) %>%
      summarise(mean = mean(value),
                median = median(value), 
                min.1 = quantile(value, 0.25),
                max.1 = quantile(value, 0.75),
                min.2 = quantile(value, 0.025),
                max.2 = quantile(value, 0.975))
}


for (model in grep("_earlier$", models, invert = TRUE, value = TRUE))
{
    r0_gi[[model]] <- p_r0gi[[model]]$data$params %>%
        mutate(model = model,
               mosquito.lifespan = "2 weeks")
    all_params[[model]] <- params[[model]][["posterior"]] %>%
        spread(state, value) %>%
        mutate(model = model,
               mosquito.lifespan = "2 weeks")

    if (paste0(model, "_earlier") %in% names(p_r0gi))
    {
        r0_gi[[model]] <-
            rbind(r0_gi[[model]],
                  p_r0gi[[paste0(model, "_earlier")]]$data$params %>%
                  mutate(model = paste0(model, "_earlier"),
                         mosquito.lifespan = "1 week"))
        all_params[[model]] <-
            rbind(all_params[[model]],
                  params[[paste0(model, "_earlier")]][["posterior"]] %>%
                  spread(state, value) %>%
                  mutate(model = model,
                         mosquito.lifespan = "1 week"))
    }
    r0_gi[[model]] <- r0_gi[[model]] %>%
      filter(distribution == "posterior") %>%
      spread(parameter, value) %>%
      mutate(`italic(G)` = cut(`italic(G)`,
                               breaks = c(2, seq(2 + 3/7, 5, 1)),
                               labels = seq(2, 4, 1))) %>%
      filter(!is.na(`italic(G)`))

    all_params[[model]] <- all_params[[model]] %>%
        filter(!(disease == "n/a" | setting == "n/a")) %>%
        mutate(data = paste(setting, disease, sep = "_")) %>%
        filter(data != "fais_zika") %>%
        mutate(data = factor(data, levels = ordered_obs_id_levels,
                             labels = data_labels))

    r0_gi_summary[[model]] <- r0_gi[[model]] %>%
        group_by(`italic(G)`, disease, setting) %>%
        summarise(mean = mean(`italic(R)[H %->% H]`, na.rm = TRUE),
                  median = median(`italic(R)[H %->% H]`, na.rm = TRUE), 
                  min.1 = quantile(`italic(R)[H %->% H]`, 0.25, na.rm = TRUE),
                  max.1 = quantile(`italic(R)[H %->% H]`, 0.75, na.rm = TRUE),
                  min.2 = quantile(`italic(R)[H %->% H]`, 0.025, na.rm = TRUE),
                  max.2 = quantile(`italic(R)[H %->% H]`, 0.975, na.rm = TRUE))

    p_r0vgi[[model]] <- ggplot(all_params[[model]] %>% filter(data != "Zika in Fais"),
                               aes(x = GI, y = R0, color = mosquito.lifespan)) +
        geom_jitter() +
        facet_grid(~ data) +
        scale_x_continuous("Equilibrium generation interval (weeks)") +
        scale_y_continuous(expression(R[0])) +
        scale_color_brewer("Mosquito life span", palette = "Dark2") +
        theme(legend.position = "top")
}

## saveRDS(list(r0_gi = r0_gi,
##              all_params = all_params,
##              r0_gi_summary = r0_gi_summary,
##              p_r0vgi = p_r0vgi,
##              p_libbi = p_libbi,
##              p_r0 = p_r0,
##              p_r0gi = p_r0gi,
##              p_r0_sqrt = p_r0_sqrt),
##         "results.rds")

save_plot("r0gi_two.pdf", two_panels, ncol = 2)
