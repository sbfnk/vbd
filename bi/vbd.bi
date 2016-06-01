model vbd {

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

  param p_p_asymptomatic[disease] // proportion of infections that are asymptomatic
  param p_lm[setting] // number of female vectors per human (log base 10)
  param p_N_h[setting]
  param p_initial_susceptible // proportion initially susceptible for dengue in Yap

  param p_rep[disease] // reporting rate

  param p_b_h[disease] // probability that a bite on a human leads to infection
  param p_b_m[disease] // probability that a bite on a vector leads to infection

  param p_r_patch_yap // importation rate between two patches
  param p_p_patch_yap // proportion of individuals in each patch

  param p_tau[setting]

  param p_t_start[setting,disease]

  param p_phi_mult[disease]
  param p_phi_add[disease]

  // humans
  state S_h[patch,setting,disease](has_output = 0) // susceptible
  state E_h[patch,setting,disease,delta_erlang_h](has_output = 0) // incubating
  state I_h[patch,setting,disease](has_output = 0) // infectious
  state R_h[patch,setting,disease](has_output = 0) // recovered
  state Z_h[patch,setting,disease] // incidence
  state C_h[patch,setting,disease] // cumulative incidence


  // vectors
  state S_m[patch,setting,disease](has_output = 0) // susceptible
  state E_m[patch,setting,disease,delta_erlang_m](has_output = 0) // incubating
  state I_m[patch,setting,disease](has_output = 0) // infectious

  state next_obs[setting,disease](has_output = 0) // time of next observation
  state started[setting,disease](has_output = 0) // outbreak start switch

  state S_h_move[patch,disease]
  state E_h_move[patch,disease,delta_erlang_h]
  state I_h_move[patch,disease]
  state R_h_move[patch,disease]

  noise n_S_move[patch,disease]
  noise n_E_move[patch,disease,delta_erlang_h]
  noise n_I_move[patch,disease]
  noise n_R_move[patch,disease]

  obs Cases[obs_id]
  obs Sero[obs_id]

  sub parameter {
    p_d_inc_h[disease] ~ log_gaussian(mean = log(5.9), std = 0.07)
    p_d_inc_m[disease] ~ log_gaussian(mean = log(9.8), std = 0.36)

    p_d_life_m ~ uniform(lower = 4, upper = 30)
    p_d_inf_h[disease] ~ truncated_gaussian(mean = 4.5, std = 1.78, lower = 0)

    p_p_asymptomatic[disease] ~ uniform(lower = 0, upper = 1)

    p_rep[disease] ~ uniform(lower = 0, upper = 1)

    p_b_h[disease] ~ uniform(lower = 0, upper = 1)
    p_b_m[disease] ~ uniform(lower = 0, upper = 1)

    p_lm[setting] ~ uniform(lower = -1, upper = 2)
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 64)

    p_phi_mult[disease] ~ uniform(lower = 0, upper = 1)
    p_phi_add[disease] ~ uniform(lower = 1, upper = 5)

    p_initial_susceptible ~ uniform(lower = 0, upper = 1)

    p_tau[setting] ~ uniform(lower = 0.3, upper = 1)

    p_r_patch_yap ~ uniform(lower = 0, upper = 0.1)
    p_p_patch_yap ~ uniform(lower = 0, upper = 1)
  }

  sub initial {
    S_h[patch,setting,disease] <- (setting == 0 && disease == 0 ? p_initial_susceptible : 1) * p_N_h[setting] * (setting == 0 ? (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap) : 1 - patch)
    E_h[patch,setting,disease,delta_erlang_h] <- 0
    I_h[patch,setting,disease] <- 0
    R_h[patch,setting,disease] <- 0
    C_h[patch,setting,disease] <- 0
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

    // movement on yap
    n_S_move[patch,disease] ~ gaussian(p_r_patch_yap * S_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap), sqrt(p_r_patch_yap * (1 - p_r_patch_yap) * S_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap)))
    n_E_move[patch,disease,delta_erlang_h] ~ gaussian(p_r_patch_yap * E_h[patch,0,disease,delta_erlang_h] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap), sqrt(p_r_patch_yap * (1 - p_r_patch_yap) * E_h[patch,0,disease,delta_erlang_h] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap)))
    n_I_move[patch,disease] ~ gaussian(p_r_patch_yap * I_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap), sqrt(p_r_patch_yap * (1 - p_r_patch_yap) * I_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap)))
    n_R_move[patch,disease] ~ gaussian(p_r_patch_yap * R_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap), sqrt(p_r_patch_yap * (1 - p_r_patch_yap) * R_h[patch,0,disease] / (patch == 0 ? p_p_patch_yap : 1 - p_p_patch_yap)))

    S_h_move[patch,disease] <- - max(0, min(floor(n_S_move[patch,disease] + 0.5), S_h[patch,0,disease])) + max(0, min(floor(n_S_move[1 - patch,disease] + 0.5), S_h[1 - patch,0,disease]))
    E_h_move[patch,disease,delta_erlang_h] <- - max(0, min(floor(n_E_move[patch,disease,delta_erlang_h] + 0.5), E_h[patch,0,disease,delta_erlang_h])) + max(0, min(floor(n_E_move[1 - patch,disease,delta_erlang_h] + 0.5), E_h[1 - patch,0,disease,delta_erlang_h]))
    I_h_move[patch,disease] <- - max(0, min(floor(n_I_move[patch,disease] + 0.5), I_h[patch,0,disease])) + max(0, min(floor(n_I_move[1 - patch,disease] + 0.5), I_h[1 - patch,0,disease]))
    R_h_move[patch,disease] <- - max(0, min(floor(n_R_move[patch,disease] + 0.5), R_h[patch,0,disease])) + max(0, min(floor(n_R_move[1 - patch,disease] + 0.5), R_h[1 - patch,0,disease]))

    S_h[patch,0,disease] <- S_h[patch,0,disease] + S_h_move[patch,disease]
    E_h[patch,0,disease,delta_erlang_h] <- E_h[patch,0,disease,delta_erlang_h] + E_h_move[patch,disease,delta_erlang_h]
    I_h[patch,0,disease] <- I_h[patch,0,disease] + I_h_move[patch,disease]
    R_h[patch,0,disease] <- R_h[patch,0,disease] + R_h_move[patch,disease]

    ode {
      dS_h[patch,setting,disease]/dt =
      - p_tau[setting] * p_b_h[disease] * pow(10, p_lm[setting])* I_m[patch,setting,disease] * S_h[patch,setting,disease]

      dE_h[patch,setting,disease,delta_erlang_h]/dt =
      + (delta_erlang_h == 0 ? p_tau[setting] * p_b_h[disease] * pow(10, p_lm[setting])* I_m[patch,setting,disease] * S_h[patch,setting,disease] : e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h - 1])
      - e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h]

      dI_h[patch,setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      - (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]

      dR_h[patch,setting,disease]/dt =
      + (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]

      dZ_h[patch,setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]

      dC_h[patch,setting,disease]/dt =
      + e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]

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
    Cases[obs_id] ~ truncated_gaussian(mean = p_rep[obs_id / 2] * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]), std = sqrt((p_rep[obs_id / 2] * (1 - p_rep[obs_id / 2]) * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]) + p_phi_add[obs_id / 2]) / p_phi_mult[obs_id / 2]), lower = 0)
    Sero[obs_id] ~ gaussian(mean = C_h[0,obs_id % 2,obs_id / 2] / p_N_h[obs_id % 2], std = 0.09 / 3.98)
  }

}
