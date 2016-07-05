model vbd {

  const e_delta_h = 1
  const e_delta_m = 1

  dim delta_erlang_h(e_delta_h)
  dim delta_erlang_m(e_delta_m)
  dim setting(e_setting)
  dim disease(e_disease)
  dim obs_id(e_obs_id)

  param p_d_inc_h[disease]
  param p_d_inc_m[disease]
  param p_d_inf_h[disease]
  param p_d_life_m

  param p_tau[setting]

  param p_vol_transmission[setting]

  param p_p_asymptomatic[disease] // proportion of infections that are asymptomatic
  param p_lm[setting] // number of female vectors per human (log base 10)
  param p_initial_susceptible // proportion initially susceptible for dengue in Yap

  param p_rep[disease] // reporting rate

  param p_b_m[disease] // probability that a bite on a vector leads to infection
  param p_b_h[disease] // probability that a bite by a vector leads to infection

  param p_t_start[setting,disease]

  param p_phi_mult[disease]
  // param p_phi_add[disease]

  // humans
  state S_h[setting,disease](has_output = 0) // susceptible
  state E_h[setting,disease,delta_erlang_h](has_output = 0) // incubating
  state I_h[setting,disease](has_output = 0) // infectious
  state R_h[setting,disease] // recovered
  state Z_h[setting,disease] // incidence

  // vectors
  state S_m[setting,disease](has_output = 0) // susceptible
  state E_m[setting,disease,delta_erlang_m](has_output = 0) // incubating
  state I_m[setting,disease](has_output = 0) // infectious

  state next_obs[setting,disease](has_output = 0) // time of next observation
  state started[setting,disease](has_output = 0) // outbreak start switch

  state R0_m_h[setting,disease] // human-to-mosquito
  state R0_h_m[setting,disease] // mosquito-to-human

  noise n_transmission[setting,disease](has_output = 0)

  input p_N_h[setting]

  obs Cases[obs_id]
  obs Sero[obs_id]

  sub parameter {
    p_d_inc_h[disease] ~ log_gaussian(mean = log(5.9/7), std = 0.07/7)
    p_d_inc_m[disease] ~ log_gaussian(mean = log(9.8/7), std = 0.36/7)

    p_d_life_m ~ gaussian(mean = 2, std = 1)
    p_d_inf_h[disease] ~ truncated_gaussian(mean = 4.5/7, std = 1.78/7, lower = 0)

    p_p_asymptomatic[disease] ~ uniform(lower = 0, upper = 1)

    p_rep[disease] ~ uniform(lower = 0, upper = 1)

    p_b_h[disease] ~ uniform(lower = 0, upper = 1)
    p_b_m[disease] ~ uniform(lower = 0, upper = 1)

    p_lm[setting] ~ uniform(lower = -1, upper = 2)
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 9)

    p_phi_mult[disease] ~ uniform(lower = 0, upper = 1)
    // p_phi_add[disease] ~ uniform(lower = 0, upper = 5)

    p_initial_susceptible ~ uniform(lower = 0, upper = 1)

    p_vol_transmission[setting] ~ uniform(lower = 0, upper = 3)
  }

  sub initial {
    S_h[setting,disease] <- (setting == 0 && disease == 0 ? p_initial_susceptible : 1) * p_N_h[setting]
    E_h[setting,disease,delta_erlang_h] <- 0
    I_h[setting,disease] <- 0
    R_h[setting,disease] <- 0
    Z_h[setting,disease] <- 0
    E_m[setting,disease,delta_erlang_m] <- 0
    S_m[setting,disease] <- 1
    I_m[setting,disease] <- 0
    next_obs[setting,disease] <- 0
    started[setting,disease] <- 0
  }

  sub transition {

    inline r_death_m = 1 / p_d_life_m
    inline r_births_m = 1 / p_d_life_m

    Z_h[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? 0 : Z_h[setting,disease])
    next_obs[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? t_next_obs : next_obs[setting,disease])

    // n_transmission[setting,disease] ~ gamma(shape = pow(10, p_vol_transmission[setting]), scale = 1 / pow(10, p_vol_transmission[setting]))
    n_transmission[setting,disease] ~ gaussian(mean = 0, std = p_vol_transmission[setting])

    R0_m_h[setting,disease] <- min(0, p_tau[setting] * p_b_m[disease] * pow(p_d_life_m, 2) / (p_d_life_m + p_d_inc_m[disease]) + n_transmission[setting,disease])
    R0_h_m[setting,disease] <- min(0, p_tau[setting] * p_b_h[disease] * pow(10, p_lm[setting]) * p_d_inf_h[disease] + n_transmission[setting,disease])

    ode {
      dS_h[setting,disease]/dt =
      - R0_h_m[setting,disease] / p_d_inf_h[disease] * I_m[setting,disease] * S_h[setting,disease]

      dE_h[setting,disease,delta_erlang_h]/dt =
      + (delta_erlang_h == 0 ? R0_h_m[setting,disease] / p_d_inf_h[disease] * I_m[setting,disease] * S_h[setting,disease] : e_delta_h * (1 / p_d_inc_h[disease]) * E_h[setting,disease,delta_erlang_h - 1])
      - e_delta_h * (1 / p_d_inc_h[disease]) * E_h[setting,disease,delta_erlang_h]

      dI_h[setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[setting,disease,e_delta_h - 1]
      - (1 / p_d_inf_h[disease]) * I_h[setting,disease]

      dR_h[setting,disease]/dt =
      + p_p_asymptomatic[disease] * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[setting,disease,e_delta_h - 1]
      + (1 / p_d_inf_h[disease]) * I_h[setting,disease]

      dZ_h[setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[setting,disease,e_delta_h - 1]

      dS_m[setting,disease]/dt =
      + r_births_m
      - R0_m_h[setting,disease] * (p_d_life_m + p_d_inc_m[disease]) / pow(p_d_life_m, 2) * I_h[setting,disease] / p_N_h[setting] * S_m[setting,disease]
      - r_death_m * S_m[setting,disease]

      dE_m[setting,disease,delta_erlang_m]/dt =
      + (delta_erlang_m == 0 ? R0_m_h[setting,disease] * (p_d_life_m + p_d_inc_m[disease]) / pow(p_d_life_m, 2) * I_h[setting,disease] / p_N_h[setting] * S_m[setting,disease] : e_delta_m * (1 / p_d_inc_m[disease]) * E_m[setting,disease,delta_erlang_m - 1])
      - e_delta_m * (1 / p_d_inc_m[disease]) * E_m[setting,disease,delta_erlang_m]
      - r_death_m * E_m[setting,disease,delta_erlang_m]

      dI_m[setting,disease]/dt =
      + e_delta_m * (1 / p_d_inc_m[disease]) * E_m[setting,disease,e_delta_m - 1]
      - r_death_m * I_m[setting,disease]
    }

    I_h[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : I_h[setting,disease])
    S_h[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? S_h[setting,disease] - 1 : S_h[setting,disease])
    Z_h[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : Z_h[setting,disease])

    started[setting,disease] <- ((t_now + 1) >= p_t_start[setting,disease] ? 1 : 0 * Z_h[setting,disease]) // put Z_h in so libbi doesn't rearrange

  }

  sub observation {
    Cases[obs_id] ~ truncated_gaussian(mean = p_rep[obs_id / 2] * Z_h[obs_id % 2,obs_id / 2], std = sqrt((p_rep[obs_id / 2] * (1 - p_rep[obs_id / 2]) * Z_h[obs_id % 2,obs_id / 2] + 1) / p_phi_mult[obs_id / 2]), lower = 0)
    Sero[obs_id] ~ gaussian(mean = R_h[obs_id % 2,obs_id / 2] / p_N_h[obs_id % 2], std = 0.09 / 3.98)
  }

}
