model vbd {

  const N = 0 // to be set in R script

  param p_d_inc_h // dealy infection by a human -> infectiousness of next human case
  param p_d_inf_h // infectious period of humans

  param p_p_risk // proportion of the population at risk
  param p_p_immune // proportion of the population at risk

  param p_R0 // human-to-human basic reproduction number

  param p_p_rep // proportion of cases reported
  param p_p_over // multiplicative overdispersion of reporting

  param p_s_peak // seasonal peak week
  param p_s_amp // strength of seasonal forcing

  param initI // initial number of infectious

  state S (has_output = 0) // susceptible
  state E (has_output = 0) // incubating
  state I (has_output = 0) // infectious
  state R (has_output = 0) // recovered
  state Z (has_output = 0) // incidence accumulator
  state beta_track
  state Reff

  input serology_sample

  obs Incidence
  obs Serology

  sub parameter {
    // from Ferguson et al., Science 95% approximately 2 * std away from the mean
    // incubation = IIP + mosquito-to-human GT
    // error on p_d_inc_h via error propgation rule with mean distance from expectation
    p_d_inc_h ~ truncated_gaussian(mean = 17.8/7, std = 2.3/7, lower = 0)
    p_d_inf_h ~ truncated_gaussian(mean = 4.7/7, std = 1.2/7, lower = 0)

    p_p_immune ~ gamma(shape = 1, scale = 0.06)
    p_p_risk ~ uniform(lower = 0.1, upper = 1)
    p_R0 ~ uniform(lower = 0, upper = 25)

    p_s_amp ~ uniform(lower = 0, upper = 1)
    p_s_peak ~ gaussian(mean = 20, std = 2)

    // uninformed prior
    p_p_rep ~ uniform(lower = 0, upper = 1)

    // weak prior on I
    initI ~ gamma(shape = 1, scale = 10)

    // regularising prior
    p_p_over ~ beta(1, 10)
  }

  sub initial {
    S <- max(N * (1 - p_p_immune) * p_p_risk - initI, 0)
    E <- 0
    I <- initI
    R <- N * p_p_immune * p_p_risk
    Z <- 0
    beta_track <- p_R0/p_d_inf_h * (1 + p_s_amp*cos(6.283*(-p_s_peak)/52))
    Reff <- p_R0 * S / (N * p_p_risk) * (1 + p_s_amp*cos(6.283*(-p_s_peak)/52))
  }

  sub transition {

    inline incubation_rate = 1/p_d_inc_h
    inline recovery_rate = 1/p_d_inf_h
    inline infection_rate = p_R0/p_d_inf_h
    inline transmission_rate = infection_rate*(1+p_s_amp*cos(6.283*(t_now-p_s_peak)/52))

    beta_track <- transmission_rate

    Reff <- transmission_rate/recovery_rate * S/(N * p_p_risk)
    // reset accumulator 
    Z <- 0

    ode {
      dS/dt = -transmission_rate*S*I/(N*p_p_risk)
      dE/dt = +transmission_rate*S*I/(N*p_p_risk)-incubation_rate*E
      dI/dt = +incubation_rate*E-recovery_rate*I
      dR/dt = +recovery_rate*I
      dZ/dt = +incubation_rate*E
    }
  }

  sub observation {
    // cases: (approximately) Poisson
    // use max(0, Z) to mitigate numerical inaccuracies
    Incidence ~ truncated_gaussian(mean = p_p_rep * max(0, Z), std = sqrt(p_p_rep * max(0, Z) / (1 - p_p_over)), lower=0)
    // serology: binomial
    Serology ~ binomial(size = serology_sample, prob = R / (N * p_p_risk))
  }
}
