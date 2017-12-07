############################################################################
## Script for running libbi analysis                                      ##
############################################################################

library('docopt')

"Script for fitting the dengue/zika model to the yap/fais data

Usage: dengue_zika_mcmc.r [options]

Options:
  -o --output=<output.file>         output file name
  -n --nsamples=<samples>           number of samples to obtain
  -p --pre_samples=<pre.samples>    number of preparatory samples to obtain
  -c --nparticles=<particles>       number of particles
  -e --seed=<seed>                  random seed
  -t --threads=<num.threads>        number of threads
  -i --thin=<thin>                  thin
  -r --sample-prior                 sample prior
  -l --sample-observations          sample observations
  -d --shorter-mosquito-lifespan    earlier death of mosquitoes
  -s --sero                         include sero data
  -j --pop                          reduce population size
  -a --stoch                        stochastic model
  -q --patch                        patch model
  -u --reverse                      reverse the patches between zika and dengue
  -g --fix-move                     fix movement
  -x --fix-natural-history          fix the natural history of the mosquito
  -m --model-file=<model.file>      given model file (means there will be no adaptation step)
  -N --noise                        add multiplicative noise
  -k --keep                         keep working directory
  -f --force                        force overwrite
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

pop_size <- c(yap = 7370, fais = 294)

## read command line arguments
num_samples <- as.integer(opts[["nsamples"]])
num_particles <- as.integer(opts[["nparticles"]])
pre_samples <- as.integer(opts[["pre_samples"]])
num_threads <- as.integer(opts[["threads"]])
seed <- as.integer(opts[["seed"]])
thin <- as.integer(opts[["thin"]])
output_file_name <- opts[["output"]]
model_file <- opts[["model-file"]]
earlier_death <- opts[["shorter-mosquito-lifespan"]]
sample_obs <- opts[["sample-observations"]]
sample_prior <- opts[["sample-prior"]]
sero <- opts[["sero"]]
pop <- opts[["pop"]]
stoch <- opts[["stoch"]]
fix_natural_history <- opts[["fix-natural-history"]]
fix_move <- opts[["fix-move"]]
patch <- opts[["patch"]]
reverse <- opts[["reverse"]]
force <- opts[["force"]]
keep <- opts[["keep"]]
verbose <- opts[["verbose"]]
par_nb <- as.integer(opts[["parallel-number"]])
analysis_setting <- opts[["setting"]]
analysis_disease <- opts[["disease"]]
noise <- opts[["noise"]]

library('dplyr')
library('tidyr')
library('rbi')
library('rbi.helpers')
library('cowplot')
library('stringi')
library('truncnorm')
library('magrittr')

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
             week = floor(nr / 7))
    ts <- c(ts, list(this_ts))
    ## set end time
    tend <- c(tend, this_ts %>% select(week) %>% max)
}

tend <- max(tend)
dt_ts <- bind_rows(ts) %>%
    group_by(week, setting, disease) %>%
    summarize(value = sum(value)) %>%
    ungroup() ## %>%
    ## mutate(obs_id = factor(paste(setting, disease, sep = "_"),
    ##                        levels = c("yap_dengue", "fais_dengue", "yap_zika"))) %>%
    ## arrange(week, obs_id) %>%
    ## select(week, obs_id, value) ## %>%
    ## complete(week, obs_id, fill = list(value = 0))

if (length(thin) == 0) thin <- 1

## get model
if (length(model_file) == 0)
{
    model_file_name <- paste(code_dir, "bi", "vbd_yap_fais", sep = "/")
    model_file_name <- paste(model_file_name, "bi", sep = ".")
} else
{
    model_file_name <- model_file
}

model <- bi_model(model_file_name)
if (verbose) ## all states
{
  no_output_pattern <- "has_output[[:space:]]*=[[:space:]]*0"
  no_output <- grep(no_output_pattern, model)
  updated_lines <- sub(no_output_pattern, "", model[no_output])
  model[no_output] <-  updated_lines
}

if (!patch)
{
  model %<>% fix(p_p_patch_yap = 1,
                 p_red_foi_yap = 1)
}

if (!stoch)
{
  model %<>% fix(p_sd_lm = 0)
  model %<>% fix(n_lm = 0)
}

if (fix_natural_history)
{
    p_tau <- 7
    if (earlier_death)
    {
        p_d_life_m <- 1
    } else
    {
        p_d_life_m <- 2
    }
    model %<>% fix(p_d_life_m = p_d_life_m,
              p_tau = p_tau)
}

if (fix_move)
{
    p_p_patch_yap <- 5.674918e-01
    p_red_foi_yap <- 3.548147e-03
    model %<>% fix(p_p_patch_yap = p_p_patch_yap,
              p_red_foi_yap = p_red_foi_yap)
}

if (reverse)
{
    initial_line <- grep("^[[:space:]]*S_h\\[patch", model)
    new_line <- sub("patch == 0", "(patch == 0 && disease == 0)", model[initial_line])
    model[new_line] <- initial_line
}

if (!sero || pop)
{
  initial_susceptible_parameter_line <- grep("p_initial_susceptible_yap.*~", model)
  model %<>% insert_lines("p_initial_susceptible_yap[1] <- 1",
                          after = initial_susceptible_parameter_line)
}

if (!pop)
{
    model %<>% fix(p_pop_yap = 1)
}

if (!noise)
{
    model %<>% fix(p_phi = 0)
}

## set output file name
if (length(output_file_name) == 0)
{
    filebase <- "vbd"
    output_file_name <- paste0(data_dir, "/", filebase, ifelse(stoch, "_stoch", ""), ifelse(noise, "_noise", ""), ifelse(fix_move, "_move", ""), ifelse(sero, "_sero", ""), ifelse(pop, "_pop", ""), ifelse(patch, "_patch", ""), ifelse(reverse, "_reverse", ""), ifelse(fix_natural_history, "_fnh", ""), ifelse(earlier_death, "_shorter", ""), ifelse(nrow(analyses) == 1, paste("", as.character(analyses[1, "setting"]), as.character(analyses[1, "disease"]), sep = "_"), ""),  ifelse(length(par_nb) == 0, "", paste0("_", par_nb)))
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
                        end_time = tend,
                        noutputs = tend)

## set seed
if (length(seed) == 0) {
    seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
}

init <- list(p_N_h = data.frame(setting = factor(c("yap", "fais"),
                                                 levels = c("yap", "fais")),
                                value = pop_size[c("yap", "fais")]))
set.seed(seed)

obs <- list(Cases = dt_ts)
if (sero)
{
    if ("yap_zika" %in% dt_ts$obs_id)
    {
        obs[["Sero"]] <- data.frame(week = dt_ts %>%
                                        filter(obs_id == "yap_zika") %>%
                                        .$week %>% max,
                                    obs_id = "yap_zika", value = 0.73)
    }
}

if (sample_prior)
{
    cat(date(), "Sampling from the prior distribution.\n")
    libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
    global_options[["seed"]] <- libbi_seed
    ## sample prior
    prior <- sample(model = model, options = global_options,
                    working_folder = working_folder, target = "prior",
                    init = init, verbose = verbose, obs=obs)
    ## reading
    save_libbi(prior, paste(output_file_name, "prior.rds", sep = "_"))
}

## sample prior with likelihoods
cat(date(), "Sampling from the posterior distribution with prior = proposal.\n")
libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
global_options[["seed"]] <- libbi_seed
global_options[["nsamples"]] <- pre_samples
if (stoch)
{
  if (length(num_particles) == 0) {
    global_options[["nparticles"]] <- 2**floor(log2(nrow(dt_ts)))
  } else
  {
    global_options[["nparticles"]] <- num_particles
  }
} else
{
  global_options[["nparticles"]] <- 1
}

bi <- sample(model = model, proposal="prior", obs = obs, options = global_options,
             working_folder = working_folder,
             init = init, verbose = verbose)

if (stoch && length(num_particles) == 0)
{
  bi %<>% adapt_particles
}

if (length(model_file) == 0)
{
    cat(date(), "Adapting the proposal distribution.\n")
    bi %<>% adapt_proposal(min = 0.1, max = 0.4)
}

cat(date(), "Sampling from the posterior distribution of the full model.\n")

libbi_seed <- ceiling(runif(1, -1, .Machine$integer.max - 1))
bi %<>% sample(nsamples = num_samples, seed = libbi_seed, thin = thin, sample_obs=TRUE)

cat(date(), "Done.\n")

save_libbi(bi, paste0(output_file_name, ".rds"))

if (!keep) unlink(working_folder, recursive = TRUE)

quit()

## 25% burn-in
burn <- floor(num_samples / thin * 0.25)
dic <- compute_DIC(read = res, burn = burn)
cat("DIC: ", dic, "\n")

saveRDS(dic, file = paste0(output_file_name, "_dic.rds"))

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
        mutate(R0 = p_tau * p_d_life_m * 
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
                              levels = c("yap_dengue", "fais_dengue",
                                         "yap_zika"))) %>%
            filter(obs_id %in% dt_ts$obs_id)

        res$Cases <- states %>%
            mutate(value = rtruncnorm(n = nrow(states), a = 0,
                                      mean = p_rep * Z_h,
                                      sd = sqrt((p_rep * Z_h + 1))))

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

        p_obs <- plot_libbi(read = res, model = final_model, data = data,
                            density_args = list(adjust = 2), burn = burn,
                            extra.aes = list(color = "obs_id"),
                            states = "Cases", params = NULL, noises = NULL,
                            trend = "mean", plot = FALSE)

        p_obs_grid <- p_obs$states + facet_wrap(~ obs_id, scales = "free_y")

        ggsave(paste(output_file_name, "obs.pdf", sep = "_"), p_obs_grid, width = 10, height = 5)

        saveRDS(p_obs$data, paste0(output_file_name, "_obs_fits.rds"))
    }
}

saveRDS(res, paste0(output_file_name, ".rds"))

if (!keep) unlink(working_folder, recursive = TRUE)

quit()
