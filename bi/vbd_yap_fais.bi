model vbd_yap_fais {

  const e_setting = 2 // 0 = yap, 1 = fais
  const e_disease = 2 // 0 = dengue, 1 = zika, 2 = chik
  const e_patch = 2 // 2 patches

  dim setting(e_setting)
  dim disease(e_disease)
  dim patch(e_patch)

  param p_d_inc_h[disease]
  param p_d_inc_m[disease]
  param p_d_inf_h[disease]
  param p_d_life_m

  param p_tau[setting]

  param p_lm[setting] // number of female vectors per human (log base 10)
  param p_sd_lm[setting] // sd in number of vectors

  param p_N_h[setting]
  param p_initial_susceptible_yap[disease] // proportion initially susceptible for dengue in Yap
  param p_pop_yap

  param p_phi[disease,setting]

  param p_rep[disease] // reporting rate

  param p_b_h[disease] // probability that a bite on a human leads to infection
  param p_b_m[disease] // probability that a bite on a vector leads to infection

  param p_red_foi_yap // importation rate between two patches
  param p_p_patch_yap // proportion of individuals in each patch

  param p_t_start[setting,disease]

  // humans
  state S_h[patch,setting,disease](has_output = 0) // susceptible
  state E_h[patch,setting,disease](has_output = 0) // incubating
  state I_h[patch,setting,disease](has_output = 0) // infectious
  state R_h[patch,setting,disease](has_output = 0) // recovered
  state Z_h[patch,setting,disease] // incidence

  // vectors
  state S_m[patch,setting,disease](has_output = 0) // susceptible
  state E_m[patch,setting,disease](has_output = 0) // incubating
  state I_m[patch,setting,disease](has_output = 0) // infectious

  state next_obs[setting,disease](has_output = 0) // time of next observation
  state started[setting,disease](has_output = 0) // outbreak start switch
  state lm[patch,setting](has_output = 0)

  obs Cases[setting,disease]
  obs Sero[setting,disease]

  noise n_lm[patch,setting](has_output = 0)

  sub parameter {
    // 95% approximately 2 * std away from the mean
    // 1.5 see Ferguson et al., science
    p_d_inc_h[disease] ~ truncated_gaussian(mean = (5.9 - 1.5)/7, std = 0.25/7, lower = 0)
    p_d_inc_m[disease] ~ truncated_gaussian(mean = 6.5/7, std = 1.15/7, lower = 0)

    p_d_life_m ~ uniform(lower = 1, upper = 4)
    p_d_inf_h[disease] ~ truncated_gaussian(mean = 4.5/7, std = 1.75/7, lower = 0)

    p_rep[disease] ~ uniform(lower = 0, upper = 1)

    p_b_h[disease] ~ uniform(lower = 0, upper = 1)
    p_b_m[disease] ~ uniform(lower = 0, upper = 1)

    p_lm[setting] ~ uniform(lower = -1, upper = 2)
    p_sd_lm[setting] ~ uniform(lower = 0, upper = 1)

    p_tau[setting] ~ uniform(lower = 0.3 * 7, upper = 1 * 7)
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 9)

    p_initial_susceptible_yap[disease] ~ uniform(lower = 0, upper = 1)
    p_pop_yap ~ uniform(lower = 0, upper = 1)

    p_red_foi_yap ~ uniform(lower = 0, upper = 1)
    p_p_patch_yap ~ uniform(lower = 0.5, upper = 1)

    p_phi[disease,setting] ~ uniform(lower = 0, upper = 5)
  }

  sub initial {
    S_h[patch,setting,disease] <- p_N_h[setting] * (setting == 0 ? p_initial_susceptible_yap[disease] * p_pop_yap * (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 1 - patch)
    E_h[patch,setting,disease] <- 0
    I_h[patch,setting,disease] <- 0
    R_h[patch,setting,disease] <- 0
    Z_h[patch,setting,disease] <- 0
    E_m[patch,setting,disease] <- 0
    S_m[patch,setting,disease] <- 1
    I_m[patch,setting,disease] <- 0
    next_obs[setting,disease] <- 0
    started[setting,disease] <- 0
    lm[patch,setting] <- p_lm[setting]
  }

  sub transition {

    inline r_death_m = 1 / p_d_life_m
    inline r_births_m = 1 / p_d_life_m

    n_lm[patch,setting] ~ gaussian()
    lm[patch,setting] <- n_lm[patch,setting] * p_sd_lm[setting] + p_lm[setting]

    Z_h[patch,setting,disease] <- (t_next_obs > next_obs[setting,disease] ? 0 : Z_h[patch,setting,disease])
    next_obs[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? t_next_obs : next_obs[setting,disease])

    ode {
      dS_h[patch,setting,disease]/dt =
      - p_tau[setting] * p_b_h[disease] * pow(10, lm[patch,setting]) * I_m[patch,setting,disease] * S_h[patch,setting,disease]

      dE_h[patch,setting,disease]/dt =
      + p_tau[setting] * p_b_h[disease] * pow(10, lm[patch,setting]) * I_m[patch,setting,disease] * S_h[patch,setting,disease]
      - (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease]

      dI_h[patch,setting,disease]/dt =
      + (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease]
      - (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]

      dR_h[patch,setting,disease]/dt =
      + (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]

      dZ_h[patch,setting,disease]/dt =
      + (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease]

      dS_m[patch,setting,disease]/dt =
      + r_births_m
      - p_tau[setting] * p_b_m[disease] * (I_h[patch,setting,disease] + (setting == 0 ? p_red_foi_yap : 1) * I_h[1 - patch,setting,disease]) / (p_N_h[setting] * (setting == 0 ? p_pop_yap : 1)) * S_m[patch,setting,disease]
      - r_death_m * S_m[patch,setting,disease]

      dE_m[patch,setting,disease]/dt =
      + p_tau[setting] * p_b_m[disease] * (I_h[patch,setting,disease] + (setting == 0 ? p_red_foi_yap : 1) * I_h[1 - patch,setting,disease]) / (p_N_h[setting] * (setting == 0 ? p_pop_yap : 1)) * S_m[patch,setting,disease]
      - (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease]
      - r_death_m * E_m[patch,setting,disease]

      dI_m[patch,setting,disease]/dt =
      + (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease]
      - r_death_m * I_m[patch,setting,disease]
    }

    I_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : I_h[0,setting,disease])
    S_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? S_h[0,setting,disease] - 1 : S_h[0,setting,disease])
    Z_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : Z_h[0,setting,disease])

    started[setting,disease] <- ((t_now + 1) >= p_t_start[setting,disease] ? 1 : 0 * Z_h[0,setting,disease]) // put Z_h in so libbi doesn't rearrange

  }

  sub observation {
    Cases[setting,disease] ~ truncated_gaussian(mean = p_rep[disease] * (Z_h[0,setting,disease] + Z_h[1,setting,disease]), std = max(sqrt(p_rep[disease] * (Z_h[0,setting,disease] + Z_h[1,setting,disease]) + (p_rep[disease] ** 2) * ((Z_h[0,setting,disease] + Z_h[1,setting,disease]) ** 2) * (p_phi[setting, disease] ** 2)), 1), lower = 0)
    Sero[setting,disease] ~ gaussian(mean = (R_h[0,setting,disease] + R_h[1,setting,disease]) / p_N_h[setting], std = 0.09 / 3.98)
  }

}
