library('RBi')
library('RBi.helpers')
library('coda')
library('dplyr')
library('tidyr')
library('msm')

p_tau <- 7
p_d_life_m <- 2

output_dir <- path.expand("~/Data/Zika")

models <- c("vbd_fnh_earlier_fais_dengue",
            "vbd_sero_fnh_earlier_yap_zika",
            "vbd_sero_patch_fnh",
            "vbd_sero_patch_fnh_yap_dengue",
            "vbd_fnh_fais_dengue",
            "vbd_sero_fnh_yap_dengue",
            "vbd_sero_patch_fnh_earlier_yap_dengue",
            "vbd_sero_patch_fnh_yap_zika",
            "vbd_sero_fnh_earlier_yap_dengue",
            "vbd_sero_fnh_yap_zika",
            "vbd_sero_patch_fnh_earlier_yap_zika")

traces <- list()
bim <- list()
params <- list()

for (model in models)
{

    mcmc_pattern <- paste0(model, "_[0-9]+")
    m <- merge_parallel_runs(output_dir, mcmc_pattern, concatenate = TRUE)
    traces[[model]] <- m$traces
    bim[[model]] <- m$model
  ## l <- merge_parallel_runs(output_dir, mcmc_pattern, concatenate = FALSE)
  ## mc <- mcmc.list(lapply(l$traces, function(x) { mcmc(get_traces(x, l$model))}))

    ## p <- plot_libbi(m$traces, m$model, extra.aes = c(color = "setting", linetype = "disease"))

    if (length(m$traces) > 0)
    {

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

        r0 <- params_all %>%
            spread(state, value) %>%
            mutate(R0 = (p_tau * p_d_life_m)**2 * 
                       p_b_h * p_b_m * 10**(p_lm) * p_d_inf_h /
                       (p_d_life_m + p_d_inc_m)) %>%
            gather(state, value, loglikelihood:R0) %>%
            filter(state == "R0")

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
            filter(state %in% c("p_rep", "p_phi_mult", "p_phi_add")) %>%
            spread(state, value)

        for (rep_param in setdiff(c("p_rep", "p_phi_mult"), colnames(rep_params)))
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
                              levels = c("yap_dengue", "fais_dengue",
                                         "yap_zika"))) %>%
            filter(obs_id %in% dt_ts$obs_id)

        l$Cases <- states %>%
            mutate(value = rtnorm(n = nrow(states), lower = 0,
                                  mean = p_rep * Z_h,
                                  sd = sqrt((p_rep * Z_h + 1) / p_phi_mult)))

        traces[[model]] <- l

    }
}

m1 <- merge(data.table(traces[["vbd_fnh_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_fnh_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")

m2 <- merge(m1, data.table(traces[["vbd_sero_fnh_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]

compute_DIC(list(loglikelihood = m2))

m1 <- merge(data.table(traces[["vbd_fnh_earlier_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_fnh_earlier_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")

m2 <- merge(m1, data.table(traces[["vbd_sero_fnh_earlier_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]

compute_DIC(list(loglikelihood = m2))

m1 <- merge(data.table(traces[["vbd_patch_fnh_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_patch_fnh_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")

m2 <- merge(m1, data.table(traces[["vbd_sero_patch_fnh_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]

compute_DIC(list(loglikelihood = m2))

m1 <- merge(data.table(traces[["vbd_patch_fnh_earlier_fais_dengue"]]$loglikelihood)[, list(fd = value, np)],
            data.table(traces[["vbd_sero_patch_fnh_earlier_yap_dengue"]]$loglikelihood)[, list(yd = value, np)],
            by = "np")

m2 <- merge(m1, data.table(traces[["vbd_sero_patch_fnh_earlier_yap_zika"]]$loglikelihood[, list(yz = value, np)]),
            by = "np")
m2[, value := fd+yd+yz]

compute_DIC(list(loglikelihood = m2))

