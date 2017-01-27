model vbd {

  const N = 0 // to be set in R script

  // ** parameters
  param p_d_gen // generation interval
  param p_d_inf_h // infectious period of mosquitoes

  param p_p_risk // proportion of the population at risk
  param p_p_immune // proportion of the population at risk

  param p_R0 // human-to-human basic reproduction number

  param p_p_rep // proportion of cases reported
  param p_rep_over // overdispersion of case reporting

  param p_s_peak // seasonal peak week
  param p_s_amp // strength of seasonal forcing

  param initI // initial number of infectious

  // ** states
  // humans
  state S (has_output = 0) // susceptible
  state E (has_output = 0) // incubating
  state I (has_output = 0) // infectious
  state R (has_output = 0) // recovered
  state Z // incidence accumulator
  state beta

  // auxiliary variable
  state next_obs (has_output = 0) // time of next observation (recorded for incidence calculation)

  input serology_sample

  obs Incidence
  obs Serology

  sub parameter {
    // 95% approximately 2 * std away from the mean
    p_d_gen ~ truncated_gaussian(mean = 20/7, std = 2/7, lower = 0)
    // 95% approximately 2 * std away from the mean, 1.5d shift see Ferguson et al., science
    p_d_inf_h ~ truncated_gaussian(mean = (5.9 - 1.5)/7, std = 0.25/7, lower = 0)

    p_p_immune ~ gamma(shape = 1, scale = 0.06)
    p_p_risk ~ uniform(lower = 0.1, upper = 1)
    p_R0 ~ uniform(lower = 0, upper = 25)

    p_s_amp ~ uniform(lower = 0, upper = 1)
    p_s_peak ~ gaussian(mean = 20, std = 2)

    // uninformed prior
    p_p_rep ~ uniform(lower = 0, upper = 1)
    // regularising prior
    p_rep_over ~ gamma(shape = 1, scale = 0.01)

    // weak prior on I
    initI ~ gamma(shape = 1, scale = 10)
  }

  sub initial {
    S <- max(N * (1 - p_p_immune) * p_p_risk - initI, 0)
    E <- 0
    I <- initI
    R <- N * p_p_immune * p_p_risk
    Z <- 0
    next_obs <- 0
    beta <- p_R0/p_d_inf_h * (1 + p_s_amp*cos(6.283*(-p_s_peak)/52))
  }

  sub transition {

    inline incubation_rate = 1/(p_d_gen-p_d_inf_h)
    inline recovery_rate = 1/p_d_inf_h
    inline infection_rate = p_R0/p_d_inf_h

    // reset accumulator if t_next_obs > next_obs
    Z <- (t_next_obs > next_obs ? 0 : Z)
    // set next_obs to the time of the next observation
    next_obs <- (t_next_obs > next_obs ? t_next_obs : next_obs)

    beta <- infection_rate*(1+p_s_amp*cos(6.283*(t_now-p_s_peak)/52))

    ode (h=0.1,atoler=1e-4,rtoler=1e-4) {
      dS/dt = -beta*S*I/(N*p_p_risk)
      dE/dt = +beta*S*I/(N*p_p_risk)-incubation_rate*E
      dI/dt = +incubation_rate*E-recovery_rate*I
      dR/dt = +recovery_rate*I
      dZ/dt = +incubation_rate*E
    }
  }

  sub observation {
    // cases: (approximately) binomial
    Incidence ~ gaussian(mean = p_p_rep * Z, std = sqrt(p_p_rep * Z + p_p_rep**2 * Z**2 * p_rep_over**2))
    // serology: (approximately) binomial
    Serology ~ gaussian(mean = serology_sample * R / N, std = sqrt(serology_sample * R / N * (1 - R / N)))
  }
}
