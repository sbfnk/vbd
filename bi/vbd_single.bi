model vbd {

  const e_delta_h = 1
  const e_delta_m = 1

  dim delta_erlang_h(e_delta_h)
  dim delta_erlang_m(e_delta_m)

  param p_d_inc_h
  param p_d_inc_m
  param p_d_inf_h
  param p_d_life_m

  param p_vol_transmission

  param p_p_asymptomatic // proportion of infections that are asymptomatic
  param p_N_h
  param p_initial_susceptible // proportion initially susceptible for dengue in Yap

  param p_rep // reporting rate

  param p_t_start

  // humans
  state S_h(has_output = 0) // susceptible
  state E_h[delta_erlang_h](has_output = 0) // incubating
  state I_h(has_output = 0) // infectious
  state R_h // recovered
  state Z_h // incidence

  // vectors
  state S_m(has_output = 0) // susceptible
  state E_m[delta_erlang_m](has_output = 0) // incubating
  state I_m(has_output = 0) // infectious

  state next_obs(has_output = 0) // time of next observation
  state started(has_output = 0) // outbreak start switch

  param p_R0_m_h // human-to-mosquito
  param p_R0_h_m // mosquito-to-human

  param p_lm

  noise R0_m_h
  noise R0_h_m

  obs Cases
  obs Sero

  sub parameter {
    p_d_inc_h ~ log_gaussian(mean = log(5.9/7), std = 0.07/7)
    p_d_inc_m ~ log_gaussian(mean = log(9.8/7), std = 0.36/7)

    p_d_life_m ~ uniform(lower = 2, upper = 4)
    p_d_inf_h ~ truncated_gaussian(mean = 4.5/7, std = 1.78/7, lower = 0)

    p_p_asymptomatic ~ uniform(lower = 0, upper = 1)

    p_R0_m_h ~ uniform(lower = 0, upper = 10)
    p_R0_h_m ~ uniform(lower = 0, upper = 10)

    p_rep ~ uniform(lower = 0, upper = 1)
    p_lm ~ uniform(lower = -2, upper = 2)

    p_t_start ~ uniform(lower = 0, upper = 9)

    // p_phi_mult ~ uniform(lower = 0, upper = 1)
    // p_phi_add ~ uniform(lower = 0, upper = 5)

    p_initial_susceptible ~ uniform(lower = 0, upper = 1)

    p_vol_transmission ~ uniform(lower = 0, upper = 3)
  }

  sub initial {
    S_h <- p_initial_susceptible * p_N_h
    E_h[delta_erlang_h] <- 0
    I_h <- 0
    R_h <- 0
    Z_h <- 0
    E_m[delta_erlang_m] <- 0
    S_m <- 1
    I_m <- 0
    next_obs <- 0
    started <- 0
  }

  sub transition {

    inline r_death_m = 1 / p_d_life_m
    inline r_births_m = 1 / p_d_life_m

    Z_h <- (t_next_obs > next_obs ? 0 : Z_h)
    next_obs <- (t_next_obs > next_obs ? t_next_obs : next_obs)

    R0_m_h ~ gamma(shape = pow(p_R0_m_h, 2) / pow(10, - 2 * p_vol_transmission), scale = pow(10, 2 * - p_vol_transmission) / p_R0_m_h)
    R0_h_m ~ gamma(shape = pow(p_R0_h_m, 2) / pow(10, - 2 * p_vol_transmission), scale = pow(10, 2 * - p_vol_transmission) / p_R0_h_m)

    ode {
      dS_h/dt =
      - R0_h_m * pow(10, p_lm) / p_d_inf_h * I_m * S_h

      dE_h[delta_erlang_h]/dt =
      + (delta_erlang_h == 0 ? R0_h_m * pow(10, p_lm) / p_d_inf_h * I_m * S_h : e_delta_h * (1 / p_d_inc_h) * E_h[delta_erlang_h - 1])
      - e_delta_h * (1 / p_d_inc_h) * E_h[delta_erlang_h]

      dI_h/dt =
      + (1 - p_p_asymptomatic) * e_delta_h * (1 / p_d_inc_h) * E_h[e_delta_h - 1]

      - (1 / p_d_inf_h) * I_h

      dR_h/dt =
      + p_p_asymptomatic * e_delta_h * (1 / p_d_inc_h) * E_h[e_delta_h - 1]

      + (1 / p_d_inf_h) * I_h

      dZ_h/dt =
      + (1 - p_p_asymptomatic) * e_delta_h * (1 / p_d_inc_h) * E_h[e_delta_h - 1]


      dS_m/dt =
      + r_births_m
      - R0_m_h / pow(10, p_lm) * (p_d_life_m + p_d_inc_m) / pow(p_d_life_m, 2) * I_h / p_N_h * S_m
      - r_death_m * S_m

      dE_m[delta_erlang_m]/dt =
      + (delta_erlang_m == 0 ? R0_m_h * pow(10, p_lm) * (p_d_life_m + p_d_inc_m) / pow(p_d_life_m, 2) * I_h / p_N_h * S_m : e_delta_m * (1 / p_d_inc_m) * E_m[delta_erlang_m - 1])
      - e_delta_m * (1 / p_d_inc_m) * E_m[delta_erlang_m]
      - r_death_m * E_m[delta_erlang_m]

      dI_m/dt =
      + e_delta_m * (1 / p_d_inc_m) * E_m[e_delta_m - 1]
      - r_death_m * I_m
    }

    I_h <- (started == 0 && (t_now + 1) >= p_t_start ? 1 : I_h)
    S_h <- (started == 0 && (t_now + 1) >= p_t_start ? S_h - 1 : S_h)
    Z_h <- (started == 0 && (t_now + 1) >= p_t_start ? 1 : Z_h)

    started <- ((t_now + 1) >= p_t_start ? 1 : 0 * Z_h) // put Z_h in so libbi doesn't rearrange

  }

  sub observation {
    Cases ~ truncated_gaussian(mean = p_rep * Z_h, std = max(sqrt(p_rep * Z_h), 1), lower = 0)
    Sero ~ gaussian(mean = R_h / p_N_h, std = 0.09 / 3.98)
  }

}
