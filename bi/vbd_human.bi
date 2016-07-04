model vbd_human {

  const e_delta_h = 1

  const e_setting = 2 // 0 = yap, 1 = fais
  const e_disease = 2 // 0 = dengue, 1 = zika
  const e_obs_id = 3 // 0 = yap/dengue, 1 = fais/dengue, 2 = yap/zika
                     // setting = obs_id % 2, disease = obs_id / 2
  const e_patch = 2 // 2 patches

  dim delta_h_erlang(e_delta_h)
  dim setting(e_setting)
  dim disease(e_disease)
  dim obs_id(e_obs_id)

  param p_d_inc[disease]
  param p_d_inf[disease]

  param p_vol_transmission[setting]

  param p_p_asymptomatic[disease] // proportion of infections that are asymptomatic (and not transmitting)
  param p_rep[disease] // reporting rate
  param p_R0[setting,disease] // number of female vectors per human (log base 10)
  param p_N[setting]
  param p_initial_susceptible // proportion initially susceptible for dengue in Yap

  param p_t_start[setting,disease]

  // humans
  state S[setting,disease](has_output = 0) // susceptible
  state E[setting,disease,delta_h_erlang](has_output = 0) // incubating
  state I[setting,disease](has_output = 0) // infectious
  state R[setting,disease] // recovered
  state Z[setting,disease] // incidence

  state next_obs[setting,disease](has_output = 0) // time of next observation
  state started[setting,disease](has_output = 0) // outbreak start switch

  noise n_transmission[setting,disease](has_output = 0)

  obs Cases[obs_id]
  obs Sero[obs_id]

  sub parameter {
    p_d_inc[disease] ~ truncated_gaussian(mean = 15.7/7, std = 1.45/7, lower = 0)

    p_d_inf[disease] ~ truncated_gaussian(mean = 4.5/7, std = 1.78/7, lower = 0)

    p_p_asymptomatic[disease] ~ uniform(lower = 0, upper = 1)
    p_rep[disease] ~ uniform(lower = 0, upper = 1)

    p_R0[setting,disease] ~ uniform(lower = 0, upper = 50)
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 9)

    p_initial_susceptible ~ uniform(lower = 0, upper = 1)

    p_vol_transmission[setting] ~ uniform(lower = 0, upper = 3)
  }

  sub initial {
    S[setting,disease] <- (setting == 0 && disease == 0 ? p_initial_susceptible : 1) * p_N[setting]
    E[setting,disease,delta_h_erlang] <- 0
    I[setting,disease] <- 0
    R[setting,disease] <- 0
    Z[setting,disease] <- 0
    next_obs[setting,disease] <- 0
    started[setting,disease] <- 0
  }

  sub transition {

    Z[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? 0 : Z[setting,disease])
    next_obs[setting,disease] <- (t_next_obs > next_obs[setting,disease] ? t_next_obs : next_obs[setting,disease])

    n_transmission[setting,disease] ~ gamma(shape = pow(10, p_vol_transmission[setting]), scale = 1 / pow(10, p_vol_transmission[setting]))

    ode {
      dS[setting,disease]/dt =
      - (p_R0[setting,disease] / p_d_inf[disease]) * I[setting,disease] * S[setting,disease] / p_N[setting] * n_transmission[setting,disease]

      dE[setting,disease,delta_h_erlang]/dt =
      + (delta_h_erlang == 0 ? p_R0[setting,disease] / p_d_inf[disease] * I[setting,disease] * S[setting,disease] / p_N[setting] * n_transmission[setting,disease] : e_delta_h * (1 / p_d_inc[disease]) * E[setting,disease,delta_h_erlang - 1])
      - e_delta_h * (1 / p_d_inc[disease]) * E[setting,disease,delta_h_erlang]

      dI[setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc[disease]) * E[setting,disease,e_delta_h - 1]
      - (1 / p_d_inf[disease]) * I[setting,disease]

      dR[setting,disease]/dt =
      + p_p_asymptomatic[disease] * e_delta_h * (1 / p_d_inc[disease]) * E[setting,disease,e_delta_h - 1]
      + (1 / p_d_inf[disease]) * I[setting,disease]

      dZ[setting,disease]/dt =
      + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc[disease]) * E[setting,disease,e_delta_h - 1]

    }

    I[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : I[setting,disease])
    S[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? S[setting,disease] - 1 : S[setting,disease])
    Z[setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : Z[setting,disease])

    started[setting,disease] <- ((t_now + 1) >= p_t_start[setting,disease] ? 1 : 0 * Z[setting,disease]) // put Z in so libbi doesn't rearrange

  }

  sub observation {
    Cases[obs_id] ~ truncated_gaussian(mean = p_rep[obs_id / 2] * Z[obs_id % 2,obs_id / 2], std = max(sqrt(p_rep[obs_id / 2] * Z[obs_id % 2,obs_id / 2]), 1), lower = 0)
    Sero[obs_id] ~ gaussian(mean = R[obs_id % 2,obs_id / 2] / p_N[obs_id % 2], std = 0.09 / 3.98)
  }

}
