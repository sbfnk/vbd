model vbd_patch {

  const e_delta_h = 1
  const e_delta_m = 1

  const e_setting = 2 // 0 = yap, 1 = fais
  const e_disease = 2 // 0 = dengue, 1 = zika
  const e_obs_id = 3 // 0 = yap/dengue, 1 = fais/dengue, 2 = yap/zika
                     // setting = obs_id % 2, disease = obs_id / 2
  const e_patch = 2 // 2 patches

  dim delta_erlang_h(e_delta_h)
  dim delta_erlang_m(e_delta_m)
  dim setting(e_setting)
  dim disease(e_disease)
  dim obs_id(e_obs_id)
  dim patch(e_patch)

  param p_d_inc_h[disease]
  param p_d_inc_m[disease]
  param p_d_inf_h[disease]
  param p_d_life_m

  param p_tau[setting]

  param p_vol_transmission[setting]

  param p_p_asymptomatic[disease] // proportion of infections that are asymptomatic
  param p_lm[setting] // number of female vectors per human (log base 10)
  param p_N_h[setting]
  param p_initial_susceptible // proportion initially susceptible for dengue in Yap

  param p_rep[disease] // reporting rate

  param p_b_h[disease] // probability that a bite on a human leads to infection
  param p_b_m[disease] // probability that a bite on a vector leads to infection

  param p_lr_patch_yap // importation rate between two patches
  param p_p_patch_yap // proportion of individuals in each patch

  param p_t_start[setting,disease]

  param p_phi_mult[disease]
  // param p_phi_add[disease]

  // humans
  state S_h[patch,setting,disease](has_output = 0) // susceptible
  state E_h[patch,setting,disease,delta_erlang_h](has_output = 0) // incubating
  state I_h[patch,setting,disease](has_output = 0) // infectious
  state R_h[patch,setting,disease] // recovered
  state Z_h[patch,setting,disease] // incidence

  // vectors
  state S_m[patch,setting,disease](has_output = 0) // susceptible
  state E_m[patch,setting,disease,delta_erlang_m](has_output = 0) // incubating
  state I_m[patch,setting,disease](has_output = 0) // infectious

  state next_obs[setting,disease](has_output = 0) // time of next observation
  state started[setting,disease](has_output = 0) // outbreak start switch

  noise n_transmission[setting,disease](has_output = 0)

  obs Cases[obs_id]
  obs Sero[obs_id](has_output = 0)

  sub parameter {
    p_d_inc_h[disease] ~ log_gaussian(mean = log(5.9/7), std = 0.07/7)
    p_d_inc_m[disease] ~ log_gaussian(mean = log(9.8/7), std = 0.36/7)

    p_d_life_m ~ uniform(lower = 2, upper = 4)
    p_d_inf_h[disease] ~ truncated_gaussian(mean = 4.5/7, std = 1.78/7, lower = 0)

    p_p_asymptomatic[disease] ~ uniform(lower = 0, upper = 1)

    p_rep[disease] ~ uniform(lower = 0, upper = 1)

    p_b_h[disease] ~ uniform(lower = 0, upper = 1)
    p_b_m[disease] ~ uniform(lower = 0, upper = 1)

    p_lm[setting] ~ uniform(lower = -1, upper = 2)
    p_tau[setting] ~ uniform(lower = 0.3 * 7, upper = 1 * 7)
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 9)

    p_phi_mult[disease] ~ uniform(lower = 0, upper = 1)
    // p_phi_add[disease] ~ uniform(lower = 0, upper = 5)

    p_initial_susceptible ~ uniform(lower = 0, upper = 1)

    p_lr_patch_yap ~ uniform(lower = -5, upper = -2)
    p_p_patch_yap ~ uniform(lower = 0, upper = 1)

    p_vol_transmission[setting] ~ uniform(lower = 0, upper = 3)
  }

  sub initial {
    S_h[patch,setting,disease] <- (setting == 0 && disease == 0 ? p_initial_susceptible : 1) * p_N_h[setting] * (setting == 0 ? (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 1 - patch)
    E_h[patch,setting,disease,delta_erlang_h] <- 0
    I_h[patch,setting,disease] <- 0
    R_h[patch,setting,disease] <- 0
    Z_h[patch,setting,disease] <- 0
    E_m[patch,setting,disease,delta_erlang_m] <- 0
    S_m[patch,setting,disease] <- 1
    I_m[patch,setting,disease] <- 0
    next_obs[setting,disease] <- 0
    started[setting,disease] <- 0
  }

  sub transition {

    inline r_death_m = 1 / p_d_life_m
    inline r_births_m = 1 / p_d_life_m

    Z_h[patch,setting,disease] <- (t_next_obs > next_obs[setting,disease] ? 0 : Z_h[patch,setting,disease])
    next_obs[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? t_next_obs : next_obs[setting,disease])

    n_transmission[setting,disease] ~ gamma(shape = pow(10, p_vol_transmission[setting]), scale = 1 / pow(10, p_vol_transmission[setting]))

    ode {
      dS_h[patch,setting,disease]/dt =
      - (p_tau[setting] * p_b_h[disease] * pow(10, p_lm[setting])) * I_m[patch,setting,disease] * S_h[patch,setting,disease] * n_transmission[setting,disease]
      - (setting == 0 ? pow(10, p_lr_patch_yap) * S_h[patch,setting,disease] * (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)
      + (setting == 0 ? pow(10, p_lr_patch_yap) * S_h[patch - 1,setting,disease] * (patch == 1 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)


      dE_h[patch,setting,disease,delta_erlang_h]/dt =
      + (delta_erlang_h == 0 ? (p_tau[setting] * p_b_h[disease] * pow(10, p_lm[setting])) * I_m[patch,setting,disease] * S_h[patch,setting,disease] * n_transmission[setting,disease] : e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h - 1])
      - e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h]
      - (setting == 0 ? pow(10, p_lr_patch_yap) * E_h[patch,setting,disease,delta_erlang_h] * (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)
      + (setting == 0 ? pow(10, p_lr_patch_yap) * E_h[patch - 1,setting,disease,delta_erlang_h] * (patch == 1 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)

      dI_h[patch,setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      - (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]
      - (setting == 0 ? pow(10, p_lr_patch_yap) * I_h[patch,setting,disease] * (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)
      + (setting == 0 ? pow(10, p_lr_patch_yap) * I_h[patch - 1,setting,disease] * (patch == 1 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)


      dR_h[patch,setting,disease]/dt =
      + p_p_asymptomatic[disease] * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      + (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]
      - (setting == 0 ? pow(10, p_lr_patch_yap) * R_h[patch,setting,disease] * (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)
      + (setting == 0 ? pow(10, p_lr_patch_yap) * R_h[patch - 1,setting,disease] * (patch == 1 ? p_p_patch_yap : 1 - p_p_patch_yap) : 0)

      dZ_h[patch,setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]

      dS_m[patch,setting,disease]/dt =
      + r_births_m
      - p_tau[setting] * p_b_m[disease] * I_h[patch,setting,disease] / p_N_h[setting] * S_m[patch,setting,disease]
      - r_death_m * S_m[patch,setting,disease]

      dE_m[patch,setting,disease,delta_erlang_m]/dt =
      + (delta_erlang_m == 0 ? p_tau[setting] * p_b_m[disease] * I_h[patch,setting,disease] / p_N_h[setting] * S_m[patch,setting,disease] : e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,delta_erlang_m - 1])
      - e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,delta_erlang_m]
      - r_death_m * E_m[patch,setting,disease,delta_erlang_m]

      dI_m[patch,setting,disease]/dt =
      + e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,e_delta_m - 1]
      - r_death_m * I_m[patch,setting,disease]
    }

    I_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : I_h[0,setting,disease])
    S_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? S_h[0,setting,disease] - 1 : S_h[0,setting,disease])
    Z_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : Z_h[0,setting,disease])

    started[setting,disease] <- ((t_now + 1) >= p_t_start[setting,disease] ? 1 : 0 * Z_h[0,setting,disease]) // put Z_h in so libbi doesn't rearrange

  }

  sub observation {
    Cases[obs_id] ~ truncated_gaussian(mean = p_rep[obs_id / 2] * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]), std = max(sqrt(p_rep[obs_id / 2] * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]) / p_phi_mult[obs_id / 2]), 1), lower = 0)
    Sero[obs_id] ~ gaussian(mean = (R_h[0,obs_id % 2,obs_id / 2] + R_h[1,obs_id % 2,obs_id / 2]) / p_N_h[obs_id % 2], std = 0.09 / 3.98)
  }

}
