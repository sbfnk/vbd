analyse_traces <- function(models, dir)
{

    require('dplyr')
    require('tidyr')
    require('RBi')
    require('msm')

    res <- list()

    p_tau <- 7 ## mosquito biting rate fixed to 1/day
    obs_id_levels <- c("yap_dengue", "fais_dengue", "yap_zika") ## observations

    for (model in models)
    {
        message(model)
        ## read traces
        res[[model]] <- list(trace = list())
        res[[model]][["trace"]] <- list(prior = list(), posterior = list())
        res[[model]][["trace"]][["posterior"]] <- readRDS(paste0(dir, "/", model, ".rds"))
        res[[model]][["trace"]][["prior"]] <- readRDS(paste0(dir, "/", model, "_prior.rds"))
        res[[model]][["model"]] <- bi_model(paste0(dir, "/", model, ".bi"))
        obs <- readRDS(paste0(dir, "/", model, "_obs.rds"))
        res[[model]][["obs"]] <- bind_rows(lapply(names(obs), function(x) {obs[[x]] %>% mutate(state = x)}))

        for (dist in c("prior", "posterior"))
        {
            ## prepare variables for joining into one table
            var_names <- setdiff(names(res[[model]][["trace"]][[dist]]), "Cases")
            l <- lapply(var_names,
                        function(x)
                        {
                            res[[model]][["trace"]][[dist]][[x]] %>%
                                mutate(state = x)
                        })

            ## set disease and setting if they are contained in the model name (model only for one outbreak)
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

            names(l) <- var_names

            ## combine all
            params <- bind_rows(l)
            ## filter out time-independent parameters
            if ("time" %in% names(params))
            {
                params <- params %>%
                    filter(is.na(time)) %>%
                    select(-time)
            }

            ## disease-specific parameters
            params_disease <-
                lapply(levels(factor(params$disease)), function(x)
                {
                    params %>%
                        filter(is.na(disease)) %>%
                        mutate(disease = x)
                })

            ## join location- and disease-independent parameters with disease-specific parameters
            params_all <-
                rbind(params %>%
                      filter(!is.na(disease)),
                      bind_rows(params_disease))

            ## setting-specific parameters
            params_setting <-
                lapply(levels(factor(params$setting)),
                       function(x) {
                           params_all %>%
                               filter(is.na(setting)) %>%
                               mutate(setting = x)
                       })

            ## join setting-specific with other parameters
            params_all <-
                rbind(params_all %>%
                      filter(!is.na(setting)),
                      bind_rows(params_setting))

            ## check mosquito lifetime
            if (grepl("shorter", model))
            {
                p_d_life_m <- 1
            } else
            {
                p_d_life_m <- 2
            }

            ## calculate r0 and generation interval
            r0 <- params_all %>%
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

            ## add back to traces
            l[["R0"]] <- r0 %>% filter(state == "R0")
            l[["GI"]] <- r0 %>% filter(state == "GI")

            ## calculate final size
            l$final_size <- l$Z_h %>%
                group_by(disease, setting, np, time) %>%
                summarise(value = sum(value))

            ## in the prior runs, incidence is already cumulative (it's reset at every data point)
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

            ## if we're looking at the posterior, sample observations
            if (dist == "posterior")
            {
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
                    mutate(sample = rtnorm(n = n(), lower = 0, mean = mean, sd = sd))

                ## merge back into traces
                l$Cases <- obsdens %>%
                  select(-value) %>%
                  rename(value = sample) %>%
                  mutate(obs_id = factor(obs_id))
            }
            ## assign
            res[[model]][["trace"]][[dist]] <- l
        }
        ## not fitted
        res[[model]][["trace"]][["posterior"]][["p_t_start"]] <-
            res[[model]][["trace"]][["posterior"]][["p_t_start"]] %>%
            filter(!(disease == "zika" & setting == "fais"))
    }

    ## join models with different mosquito lifetimes

    shorter_models <-
        intersect(sub("_shorter", "", grep("_shorter", models, value = TRUE)),
                  models)

    for (model in shorter_models)
    {
        all_model <- paste(model, "all", sep = "_")
        message(all_model)
        res[[all_model]][["trace"]] <- list()
        short_model <- paste(model, "shorter", sep = "_")
        res[[all_model]][["model"]] <- res[[model]][["model"]]
        res[[all_model]][["obs"]] <- res[[model]][["obs"]]
        for (dist in names(res[[model]][["trace"]]))
        {
            res[[all_model]][["trace"]][[dist]] <- list()
            for (var in names(res[[model]][["trace"]][[dist]]))
            {
                max_np <- max(res[[model]][["trace"]][[dist]][[var]][["np"]])
                res[[all_model]][["trace"]][[dist]][[var]] <-
                    rbind(res[[model]][["trace"]][[dist]][[var]] %>%
                          mutate(p_d_life_m = 2),
                          res[[short_model]][["trace"]][[dist]][[var]] %>%
                          mutate(p_d_life_m = 1,
                                 np = np + 1 + max_np))
            }
        }
    }
    return(res)
}
