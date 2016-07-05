############################################################################
## Script for running libbi analysis                                      ##
############################################################################

library('docopt')

"Script for fitting the dengue/zika model to the yap/fais data

Usage: dengue_zika_mcmc.r [options]

Options:
  -o --output=<output.file>         output file name
  -n --nsamples=<samples>           number of samples to obtain
  -c --nparticles=<num.particles>   number of particles
  -p --pre_samples=<pre.samples>    number of preparatory samples to obtain
  -e --seed=<seed>                  random seed
  -t --threads=<num.threads>        number of threads
  -w --erlang-v=<erlang.v>          number of Erlang compartments in the vector
  -u --erlang-h=<erlang.h>          number of Erlang compartments in the human
  -d --beta                         include stochastic beta
  -i --thin=<thin>                  thin
  -r --sample-prior                 sample prior
  -l --sample-observations          sample observations
  -s --sero                         include sero data
  -g --human-only                   only model humans
  -x --fix-natural-history          fix the natural history of the mosquito
  -m --model-file=<model.file>      given model file (means there will be no adaptation step)
  -k --keep                         keep working directory
  -f --force                        force overwrite
  -q --patch                        patch model
  -a --poisson                      poisson noise
  -y --setting=<setting>            only fit this setting
  -z --disease=<disease>            only fit this disease
  -v --verbose                      be verbose
  -b --parallel-number=<number>     parallel number
  -h --help                         show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

## read command line arguments
num_samples <- as.integer(opts[["nsamples"]])
num_particles <- as.integer(opts[["nparticles"]])
pre_samples <- as.integer(opts[["pre_samples"]])
num_threads <- as.integer(opts[["threads"]])
erlang_vector <- as.integer(opts[["erlang-v"]])
erlang_human <- as.integer(opts[["erlang-h"]])
seed <- as.integer(opts[["seed"]])
thin <- as.integer(opts[["thin"]])
output_file_name <- opts[["output"]]
model_file <- opts[["model-file"]]
beta <- opts[["beta"]]
sample_obs <- opts[["sample-observations"]]
sample_prior <- opts[["sample-prior"]]
sero <- opts[["sero"]]
poisson <- opts[["poisson"]]
human_only <- opts[["human-only"]]
fix_natural_history <- opts[["fix-natural-history"]]
patch <- opts[["patch"]]
force <- opts[["force"]]
keep <- opts[["keep"]]
verbose <- opts[["verbose"]]
par_nb <- as.integer(opts[["parallel-number"]])
analysis_setting <- opts[["setting"]]
analysis_disease <- opts[["disease"]]

stoch <- beta

library('dplyr')
library('tidyr')
library('RBi')
library('RBi.helpers')
library('cowplot')
library('stringi')
library('truncnorm')

code_dir <- path.expand("~/code/vbd/")
data_dir <-  path.expand("~/Data/Zika/")

## read data
analyses <- data.frame(setting = c("yap", "yap", "fais"), disease = c("dengue", "zika", "dengue"))

if (length(analysis_setting) > 0)
{
    analyses <- analyses %>% filter(setting == analysis_setting)
}

if (length(analysis_disease) > 0)
{
    analyses <- analyses %>% filter(disease == analysis_disease)
}

tend <- c()
ts <- list()

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
    ## set end time
    tend <- c(tend, this_ts %>% select(week) %>% max)
}

tend <- max(tend)
dt_ts <- bind_rows(ts) %>%
    group_by(week, setting, disease) %>%
    summarize(value = sum(value)) %>%
    ungroup() %>%
    mutate(obs_id = factor(paste(setting, disease, sep = "_"),
                           levels = c("yap_dengue", "fais_dengue", "yap_zika"))) %>%
    arrange(week, obs_id) %>%
    select(week, obs_id, value) ## %>%
    ## complete(week, obs_id, fill = list(value = 0))

init <- list(p_N_h = data.frame(setting = c("yap", "fais"), value = c(7391, 294)))

if (length(thin) == 0) thin <- 1

## get model
if (length(model_file) == 0)
{
    model_file_name <- paste(code_dir, "bi", "vbd", sep = "/")
    if (human_only)
    {
        model_file_name <- paste(model_file_name, "human", sep = "_")
    }
    if (patch)
    {
        model_file_name <- paste(model_file_name, "patch", sep = "_")
    }
    model_file_name <- paste(model_file_name, "bi", sep = ".")
} else
{
    model_file_name <- model_file
}

model <- bi_model(model_file_name)
if (verbose) ## all states
{
  no_output_pattern <- "has_output[[:space:]]*=[[:space:]]*0"
  no_output <- grep(no_output_pattern, model$get_lines())
  updated_lines <- sub(no_output_pattern, "",
                       model$get_lines()[no_output])
  model$update_lines(no_output, updated_lines)
}

## update model lines for requested number of Erlang compartments
erlang <- c(h = "erlang_human", m = "erlang_vector")

for (comp in names(erlang))
{
    if (length(opts[[erlang[comp]]]) > 0)
    {
        erlang_line_no <- grep(paste0("const e_delta_", comp, "[[:space:]]*="), model$get_lines())
        if (length(erlang_line_no) > 0)
        {
            erlang_line <- model$get_lines()[erlang_line_no]
            model$update_lines(erlang_line_no, sub("=.*$", paste("=", opts[[erlang[comp]]]), erlang_line))
        }
    }
}

if (!patch)
{
  model$fix(p_p_patch_yap = 1,
            p_lr_patch_yap = 0,
            e_patch = 1)
}

if (!beta)
{
  model$fix(n_transmission = 1)
}

## if (!human_only)
## {
##     model$fix(p_tau = 7)
## }

if (fix_natural_history)
{
    model$fix(p_d_inc_h = 5.9 / 7,
              p_d_inc_m = 9.8 / 7,
              p_d_life_m = 2)
}

if (poisson)
{
    model$fix(p_phi_mult = 1)
}

## set output file name
if (length(output_file_name) == 0)
{
    filebase <- "vbd"
    for (comp in names(erlang))
    {
        if (length(opts[[erlang[comp]]]) > 0)
        {
            filebase <- paste(filebase, paste0(comp, opts[[erlang[comp]]]), sep = "_")
        }
    }
    output_file_name <- paste0(data_dir, "/", filebase, ifelse(human_only, "_human", ""), ifelse(beta, "_beta", ""), ifelse(poisson, "_poisson", ""), ifelse(sero, "_sero", ""), ifelse(patch, "_patch", ""), ifelse(nrow(analyses) == 1, paste("", as.character(analyses[1, "setting"]), as.character(analyses[1, "disease"]), sep = "_"), ""),  ifelse(length(par_nb) == 0, "", paste0("_", par_nb)))
}
cat("Output: ",  output_file_name, "\n")

if (!force && file.exists(paste0(output_file_name, "_params.rds")))
{
    stop("exists")
}
## working folder
working_folder <- path.expand(output_file_name)
unlink(working_folder, recursive = TRUE)
suppressWarnings(dir.create(working_folder))

## set global options
global_options <-  list(nsamples = pre_samples,
                        "end-time" = tend,
                        noutputs = tend)

## set seed
if (length(seed) == 0) {
    seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
}

cat("Seed: ", seed, "\n")

set.seed(seed)

if (sample_prior)
{
    cat(date(), "Sampling from the prior distribution.\n")
    libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
    global_options[["seed"]] <- libbi_seed
    ## sample prior
    prior <- libbi(model = model, run = TRUE,
                   global_options = global_options, client = "sample",
                   working_folder = working_folder, target = "prior",
                   dims = list(disease = c("dengue", "zika")),
                   init = init, verbose = verbose)
    ## reading
    res_prior <- bi_read(prior, vars = model$get_vars("param"),
                         verbose = verbose)
    saveRDS(res_prior, paste(output_file_name, "prior.rds", sep = "_"))
    prior_model_file <- paste(output_file_name, "prior.bi", sep = "_")
    prior$model$write_model_file(prior_model_file)
}

## sample prior with likelihoods
cat(date(), "Sampling from the posterior distribution with prior = proposal.\n")
libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
global_options[["seed"]] <- libbi_seed
global_options[["nsamples"]] <- pre_samples
if (length(num_particles) > 0)
{
  global_options[["nparticles"]] <- num_particles
} else
{
  ## number of data points as number of particles
  if (stoch)
  {
    ## global_options[["nparticles"]] <- 2**floor(log(nrow(dt_ts), 2))
    global_options[["nparticles"]] <- 4
  }
}
obs <- list(Cases = dt_ts)
if (sero)
{
    obs[["Sero"]] <- data.frame(week = dt_ts %>% filter(obs_id == "yap_zika") %>% .$week %>% max, obs_id = "yap_zika", value = 0.73)
}

model_prior <- model$propose_prior()
bi_wrapper_prior <- libbi(model = model_prior, run = TRUE,
                          obs = obs, global_options = global_options, client = "sample",
                          dims = list(disease = c("dengue", "zika")),
                          working_folder = working_folder,
                          init = init, verbose = verbose)

cat(date(), "Running the model.\n")
bi_wrapper_adapted <- adapt_mcmc(bi_wrapper_prior, min = 0, max = 1)

if (stoch)
{
    bi_wrapper_stoch <- bi_wrapper_adapted
    if (length(num_particles) > 0)
    {
        bi_wrapper_adapted <- bi_wrapper_stoch
    } else
    {
        cat(date(), "Starting adaptation of the number of particles.\n")
        bi_wrapper_particle_adapted <-
            adapt_particles(bi_wrapper_stoch,
                            ## min = bi_wrapper_stoch$global_options[["nparticles"]]
                            min = 1024, max = 33000
                            )
        bi_wrapper_adapted <- bi_wrapper_particle_adapted
    }
}

if ("nparticles" %in% names(bi_wrapper_adapted$global_options))
{
    nparticles <- bi_wrapper_adapted$global_options[["nparticles"]]
} else
{
    nparticles <- 1
}

if (length(model_file) == 0)
{
    cat(date(), "Adapting the proposal distribution.\n")
    bi_wrapper_adapted <-
        adapt_mcmc(bi_wrapper_adapted, min = 0.1, max = 0.5, max_iter = 10, scale = 2)
} else
{
    bi_wrapper_adapted$model <- model
}

cat(date(), "Sampling from the posterior distribution of the full model.\n")

## if (sample_obs)
## {
##   for (state in c("Z_h", "C_h"))
##   {
##     state_line_no <- grep(paste0("state[[:space:]]*", state),
##                           bi_wrapper$model$get_lines())
##     state_line <- bi_wrapper$model$get_lines()[state_line_no]
##     bi_wrapper$model$update_lines(state_line_no, sub("has_output[[:space:]]*=[[:space:]]*0", "has_output = 1", state_line))
##   }
## }

libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
bi_wrapper <- bi_wrapper_adapted
bi_wrapper$run(add_options = list("init-np" = pre_samples - 1,
                                  nsamples = num_samples,
                                  seed = libbi_seed),
               init = bi_wrapper_adapted, verbose = verbose)

if (length(model_file) == 0)
{
    model_file <- paste(output_file_name, "bi", sep = ".")
    bi_wrapper$model$write_model_file(model_file)
}

command_file <- paste(output_file_name, "cmd", sep = ".")
cat(bi_wrapper$result$command, file = command_file)

final_model <- bi_wrapper$model

res <- bi_read(read = bi_wrapper, thin = thin, verbose = verbose)

## 25% burn-in
burn <- num_samples * 0.25
dic <- compute_DIC(read = res, burn = burn)

cat("DIC: ", dic, "\n")

saveRDS(list(dic = dic, nparticles = nparticles),
        file = paste0(output_file_name, "_dic.rds"))

if (length(par_nb) == 0)
{
    cat(date(), "Plotting.\n")

    plot_args <- list(read = res, model = final_model,
                      density_args = list(adjust = 2), burn = burn,
                      extra.aes = c(color = "disease", linetype = "setting"),
                      plot = FALSE)
    if (nrow(analyses) == 1)
    {
        plot_args[["select"]] <-
            c(setting = as.character(analyses[1, "setting"]),
              disease = as.character(analyses[1, "disease"]))
    }
    if (sample_prior)
    {
        plot_args[["prior"]] <- res_prior
    }
    p_param <- do.call(plot_libbi, plot_args)
    saveRDS(p_param$data, paste0(output_file_name, "_param_fits.rds"))

    if (!is.null(p_param[["states"]]))
    {
        ggsave(paste(output_file_name, "states.pdf", sep = "_"), p_param$states)
    }
    if (!is.null(p_param[["densities"]]))
    {
        ggsave(paste(output_file_name, "densities.pdf", sep = "_"), p_param$densities)
    }
    if (!is.null(p_param[["traces"]]))
    {
        ggsave(paste(output_file_name, "traces.pdf", sep = "_"), p_param$traces)
    }
    if (!is.null(p_param[["correlations"]]))
    {
        ggsave(paste(output_file_name, "correlations.pdf", sep = "_"),
               p_param$correlations)
    }
    if (!is.null(p_param[["noises"]]))
    {
        ggsave(paste(output_file_name, "noises.pdf", sep = "_"), p_param$noises)
    }
    if (!is.null(p_param[["likelihoods"]]))
    {
        ggsave(paste(output_file_name, "likelihoods.pdf", sep = "_"), p_param$likelihoods)
    }

    cat(date(), "..parameters.\n")
    l <- lapply(names(res), function(x) {
        res[[x]] %>% mutate(state = x)
    })

    params <- bind_rows(l)
    if ("time" %in% names(params))
    {
        params <- params %>%
            filter(is.na(time)) %>%
            select(-time)
    }

    params_disease <-
        lapply(levels(factor(params$disease)),
               function(x) { params %>% filter(is.na(disease)) %>% mutate(disease = x) })

    params_all <-
        rbind(params %>%
              filter(!is.na(disease)),
              bind_rows(params_disease))

    params_setting <-
        lapply(levels(factor(params$setting)),
               function(x) { params_all %>% filter(is.na(setting)) %>% mutate(setting = x) })

    params_all <-
        rbind(params_all %>%
              filter(!is.na(setting)),
              bind_rows(params_setting))

    r0 <- params_all %>%
        spread(state, value) %>%
        mutate(R0 = p_d_life_m * 
                   sqrt(p_b_h * p_b_m * 10**(p_lm) * p_d_inf_h /
                        (p_d_life_m + p_d_inc_m))) %>%
        gather(state, value, loglikelihood:R0) %>%
        filter(state == "R0")

    params <- rbind(params, r0) %>%
        mutate(disease = ifelse(is.na(disease), "n/a", as.character(disease)),
               setting = ifelse(is.na(setting), "n/a", as.character(setting)))

    saveRDS(params, paste0(output_file_name, "_params.rds"))

    cat(date(), "..R0.\n")
    p_R0 <- ggplot(params %>% filter(state == "R0" & disease %in% analyses$disease & setting %in% analyses$setting) %>%
                   mutate(disease = stri_trans_totitle(disease),
                          setting = stri_trans_totitle(setting)),
                   aes(x = value, color = disease, linetype = setting)) +
        geom_line(stat = "density", lwd = 2, adjust = 2) +
        scale_x_continuous(expression(R[0]), limits = c(0, 25)) +
        scale_color_brewer(palette = "Set1") +
        theme(legend.position = "top", legend.key.size=unit(1.5, "cm"))
    save_plot(paste(output_file_name, "r0.pdf", sep = "_"), p_R0)

    if (sample_obs)
    {
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

        res$Cases <- states %>%
            mutate(value = rtruncnorm(n = nrow(states), a = 0,
                                      mean = p_rep * Z_h,
                                      sd = sqrt((p_rep * Z_h + 1) / p_phi_mult)))

        ## manipulate data to match sampled observations
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

        ## p_obs <- plot_libbi(read = res, model = final_model, data = data,
        ##                     density_args = list(adjust = 2), burn = burn,
        ##                     extra.aes = list(color = "obs_id"),
        ##                     states = "Cases", params = NULL, noises = NULL,
        ##                     trend = "mean", plot = FALSE)

        p_obs <- plot_libbi(read = res, model = final_model, data = data,
                            density_args = list(adjust = 2), burn = burn,
                            extra.aes = list(color = "obs_id"),
                            states = "Cases", params = NULL, noises = NULL,
                            trend = "mean", plot = FALSE)

        p_obs_grid <- p_obs$states + facet_wrap(~ obs_id, scales = "free_y")

        ggsave(paste(output_file_name, "obs.pdf", sep = "_"), p_obs_grid)

        saveRDS(p_obs$data, paste0(output_file_name, "_obs_fits.rds"))
    }
}

saveRDS(res, paste0(output_file_name, ".rds"))

if (!keep) unlink(working_folder, recursive = TRUE)

quit()
