model vbd {
  const e_delta_h = 1
  const e_delta_m = 1
  const e_setting = 2
  const e_disease = 2
  const e_obs_id = 3
  const e_patch = 2
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
  const p_tau = 7
  param p_vol_transmission[setting]
  param p_p_asymptomatic[disease]
  param p_lm[setting]
  param p_N_h[setting]
  param p_initial_susceptible
  param p_rep[disease]
  param p_b_h[disease]
  param p_b_m[disease]
  param p_lr_patch_yap
  param p_p_patch_yap
  param p_t_start[setting,disease]
  param p_phi_mult[disease]
  state S_h[patch,setting,disease](has_output = 0)
  state E_h[patch,setting,disease,delta_erlang_h](has_output = 0)
  state I_h[patch,setting,disease](has_output = 0)
  state R_h[patch,setting,disease]
  state Z_h[patch,setting,disease]
  state S_m[patch,setting,disease](has_output = 0)
  state E_m[patch,setting,disease,delta_erlang_m](has_output = 0)
  state I_m[patch,setting,disease](has_output = 0)
  state next_obs[setting,disease](has_output = 0)
  state started[setting,disease](has_output = 0)
  state S_h_move[patch,disease](has_output = 0)
  state E_h_move[patch,disease,delta_erlang_h](has_output = 0)
  state I_h_move[patch,disease](has_output = 0)
  state R_h_move[patch,disease](has_output = 0)
  const n_R_move = 0
  const n_I_move = 0
  const n_E_move = 0
  const n_S_move = 0
  noise n_transmission[setting,disease](has_output = 0)
  obs Cases[obs_id]
  obs Sero[obs_id]
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
    p_t_start[setting,disease] ~ uniform(lower = 0, upper = 9)
    p_phi_mult[disease] ~ uniform(lower = 0, upper = 1)
    p_initial_susceptible ~ uniform(lower = 0, upper = 1)
    p_lr_patch_yap ~ uniform(lower = -5, upper = -3)
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
    S_h_move[patch,disease] <- - max(0, min(floor(n_S_move + 0.5), S_h[patch,0,disease])) + max(0, min(floor(n_S_move + 0.5), S_h[1 - patch,0,disease]))
    E_h_move[patch,disease,delta_erlang_h] <- - max(0, min(floor(n_E_move + 0.5), E_h[patch,0,disease,delta_erlang_h])) + max(0, min(floor(n_E_move + 0.5), E_h[1 - patch,0,disease,delta_erlang_h]))
    I_h_move[patch,disease] <- - max(0, min(floor(n_I_move + 0.5), I_h[patch,0,disease])) + max(0, min(floor(n_I_move + 0.5), I_h[1 - patch,0,disease]))
    R_h_move[patch,disease] <- - max(0, min(floor(n_R_move + 0.5), R_h[patch,0,disease])) + max(0, min(floor(n_R_move + 0.5), R_h[1 - patch,0,disease]))
    n_transmission[setting,disease] ~ gamma(shape = pow(10, p_vol_transmission[setting]), scale = 1 / pow(10, p_vol_transmission[setting]))
    S_h[patch,0,disease] <- S_h[patch,0,disease] + S_h_move[patch,disease]
    E_h[patch,0,disease,delta_erlang_h] <- E_h[patch,0,disease,delta_erlang_h] + E_h_move[patch,disease,delta_erlang_h]
    I_h[patch,0,disease] <- I_h[patch,0,disease] + I_h_move[patch,disease]
    R_h[patch,0,disease] <- R_h[patch,0,disease] + R_h_move[patch,disease]
    ode {
      dS_h[patch,setting,disease]/dt = - (p_tau * p_b_h[disease] * pow(10, p_lm[setting])) * I_m[patch,setting,disease] * S_h[patch,setting,disease] * n_transmission[setting,disease]
      dE_h[patch,setting,disease,delta_erlang_h]/dt = + (delta_erlang_h == 0 ? (p_tau * p_b_h[disease] * pow(10, p_lm[setting])) * I_m[patch,setting,disease] * S_h[patch,setting,disease] * n_transmission[setting,disease] : e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h - 1])
      - e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,delta_erlang_h]
      dI_h[patch,setting,disease]/dt = + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      - (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]
      dR_h[patch,setting,disease]/dt = + p_p_asymptomatic[disease] * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      + (1 / p_d_inf_h[disease]) * I_h[patch,setting,disease]
      dZ_h[patch,setting,disease]/dt = + (1 - p_p_asymptomatic[disease]) * e_delta_h * (1 / p_d_inc_h[disease]) * E_h[patch,setting,disease,e_delta_h - 1]
      dS_m[patch,setting,disease]/dt = + r_births_m
      - p_tau * p_b_m[disease] * I_h[patch,setting,disease] / p_N_h[setting] * S_m[patch,setting,disease]
      - r_death_m * S_m[patch,setting,disease]
      dE_m[patch,setting,disease,delta_erlang_m]/dt = + (delta_erlang_m == 0 ? p_tau * p_b_m[disease] * I_h[patch,setting,disease] / p_N_h[setting] * S_m[patch,setting,disease] : e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,delta_erlang_m - 1])
      - e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,delta_erlang_m]
      - r_death_m * E_m[patch,setting,disease,delta_erlang_m]
      dI_m[patch,setting,disease]/dt = + e_delta_m * (1 / p_d_inc_m[disease]) * E_m[patch,setting,disease,e_delta_m - 1]
      - r_death_m * I_m[patch,setting,disease]
    }
    I_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : I_h[0,setting,disease])
    S_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? S_h[0,setting,disease] - 1 : S_h[0,setting,disease])
    Z_h[0,setting,disease] <- (started[setting,disease] == 0 && (t_now + 1) >= p_t_start[setting,disease] ? 1 : Z_h[0,setting,disease])
    started[setting,disease] <- ((t_now + 1) >= p_t_start[setting,disease] ? 1 : 0 * Z_h[0,setting,disease])
  }
  sub observation {
    Cases[obs_id] ~ truncated_gaussian(mean = p_rep[obs_id / 2] * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]), std = sqrt((p_rep[obs_id / 2] * (1 - p_rep[obs_id / 2]) * (Z_h[0,obs_id % 2,obs_id / 2] + Z_h[1,obs_id % 2,obs_id / 2]) + 1) / p_phi_mult[obs_id / 2]), lower = 0)
    Sero[obs_id] ~ gaussian(mean = (R_h[0,obs_id % 2,obs_id / 2] + R_h[1,obs_id % 2,obs_id / 2]) / p_N_h[obs_id % 2], std = 0.09 / 3.98)
  }
  sub proposal_initial {
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
  sub proposal_parameter {
    p_d_inc_h[disease] ~ gaussian(mean = p_d_inc_h[disease], std = 0.03125 * 0.0058368843836619)
    p_d_inc_m[disease] ~ gaussian(mean = p_d_inc_m[disease], std = 0.03125 * 0.00532655394470054)
    p_d_inf_h[disease] ~ truncated_gaussian(mean = p_d_inf_h[disease], std = 0.03125 * 0.133324329530883, lower = 0)
    p_d_life_m ~ truncated_gaussian(mean = p_d_life_m, std = 0.03125 * 0.000194467525422842, lower = 2, upper = 4)
    p_vol_transmission[setting] ~ truncated_gaussian(mean = p_vol_transmission[setting], std = 0.03125 * 0.350590450278125, lower = 0, upper = 3)
    p_p_asymptomatic[disease] ~ truncated_gaussian(mean = p_p_asymptomatic[disease], std = 0.03125 * 0.0301073099161244, lower = 0, upper = 1)
    p_lm[setting] ~ truncated_gaussian(mean = p_lm[setting], std = 0.03125 * 0.547434227711879, lower = -1, upper = 2)
    p_initial_susceptible ~ truncated_gaussian(mean = p_initial_susceptible, std = 0.03125 * 0.000186743501059629, lower = 0, upper = 1)
    p_rep[disease] ~ truncated_gaussian(mean = p_rep[disease], std = 0.03125 * 0.248630453642588, lower = 0, upper = 1)
    p_b_h[disease] ~ truncated_gaussian(mean = p_b_h[disease], std = 0.03125 * 0.0583241371758106, lower = 0, upper = 1)
    p_b_m[disease] ~ truncated_gaussian(mean = p_b_m[disease], std = 0.03125 * 0.00429256542421774, lower = 0, upper = 1)
    p_lr_patch_yap ~ truncated_gaussian(mean = p_lr_patch_yap, std = 0.03125 * 0.000327178122246328, lower = -5, upper = -3)
    p_p_patch_yap ~ truncated_gaussian(mean = p_p_patch_yap, std = 0.03125 * 0.000227220978178083, lower = 0, upper = 1)
    p_t_start[setting,disease] ~ truncated_gaussian(mean = p_t_start[setting,disease], std = 0.03125 * 1.94559060380429, lower = 0, upper = 9)
    p_phi_mult[disease] ~ truncated_gaussian(mean = p_phi_mult[disease], std = 0.03125 * 0.10255385208764, lower = 0, upper = 1)
  }
}
