# Set parameter values (all rates defined per year) and initial states
param_values <-
  c(mu = 1/68,              # Background mortality rate per capita
    muVL = 365/189,         # Excess mortality rate due to untreated VL
    cohort_model = 0,       # Flag for running a cohort model without births (0 = no, other value = yes)
    
    rhoL1 = 365/140,        # 1 / Average duration of latent infection 1 (low titre)
    rhoL2 = 365/21,         # 1 / Average duration of latent infection 2 (low titre)
    rhoH = 365/70,          # 1 / Average duration of high-titre latent infection (120 days means ~95% exits within 1 year)
    
    rhoI1 = 365/32,         # 1 / Average duration until detection and treatment of VL in group I1
    rhoI2 = 365/243,        # 1 / Average duration until detection and treatment of VL in group I2
    rhoD = 12/21,           # 1 / Average duration of the dormant stage
    rhoP = 1/5,             # 1 / Average duration of PKDL
    rhoR = 1/5,             # 1 / Average duration until spontaneous seroconversion of recovered individuals
    
    fh = 0.60,              # Proportion of infected humans that develop high antibody titres
    fs = 0.125,             # Proportion of latent infections with high titres that progress to VL
    fd = 0,                 # Proportion of VL cases that are in the high-detection group (I1); rest go to I2.
    fp = 0.05,              # Proportion of VL cases that develop PKDL
    fR = 0.5,               # Proportion of recovered cases who, when reinfected, may eventually develop clinical infection (related to immunity) 
    
    beta_bite = 365/4,      # Sandfly biting rate
    betaP = 1.0,            # Infectivity of PKDL relative to VL
    betaH = 0,              # Infectivity of high-titre latent infections relative to VL
    foi_ext = 0,            # External force of infection that does not depend on state of the system
    muF = 365/14,           # 1 / Average sandfly lifespan
    muIRS = 0,              # Excess sandfly mortality rate due to IRS
    fIRS = 0,               # Reduction in sandfly abundance due to IRS
    rhoEF = 365/5,          # 1 / Average time until sandflies because infectious
    pH = 1.0,               # Transmission probability from fly to human
    
    season_t = 0,           # The particular month of 12 months of the year 
    season_amp = 0,         # Relative amplitude of seasonality 
    
    NF = 0                  # Sandfly-to-human ratio
  )


# Initial distribution of human population over compartments
  init_values <- c(S_0 = .99, L1_0 = 0, L2_0 = 0, H_0 = 0,
                   I11_0 = 0, I12_0 = 0, I13_0 = 0,
                   I21_0 = 0, I22_0 = 0, I23_0 = 0,
                   D_0 = 0, P_0 = .01, R_0 = 0,
                   SF_0 = 0, EF_0 = 0, IF_0 = 0)

  paramnames <- c(names(param_values),
                  names(init_values))
  
### END OF CODE ### 
