# Source code for deterministic and stochastic models of transmission,
# detection, and control of visceral leishmaniasis, implemented in the pomp
# package (partially observed Markov processes).

# The model backbone is based on a previously developed model (Coffeng et al
# 2020 J Infect Dis) that focused on competing risks of mortality due to VL
# (Erlang distributed) and detection delays (exponentially distributed) and the
# effect of increasing detection effort in terms of coverage and intensity.

# The model here is an extension that considers that:
# - IRS may have an effect on transmission probability via excess sandfly
#   mortality
# - previously infected and recovered people may spontaneously become 
#   seropositive again


# Key assumptions:
#  - fixed population size
#  - no age structure
#  - stochastic human population
#  - deterministic fly population

# Author: Ananthu James (stochsatic model) and Luc Coffeng (deterministic model) 
# Created: 13 July 2023 
  
# Deterministic model ----
statenames_det <- c("S", "L1", "L2", "H",
                    "I11", "I12", "I13",
                    "I21", "I22", "I23",
                    "D", "P", "R",
                    "SF", "EF", "IF")

vl_code_det <-
"
  // Declare local variables
  double N;         // Total human population size
  double lambdaHF;  // Force of infection from humans to flies
  double lambdaFH;  // Force of infection from flies to humans

  N = S + L1 + L2 + H + I11 + I12 + I13 + I21 + I22 + I23 + D + P + R;
  lambdaHF = beta_bite * (betaH * H +
                          I11 + I12 + I13 + I21 + I22 + I23 +
                          betaP * P) / N;
  lambdaFH = beta_bite * pH * IF + foi_ext;
  

  // System of ODEs
  DS = - (lambdaFH + mu) * S ;
  if (cohort_model == 0) {
    DS += mu * N + 3 * muVL * (I13 + I23);
  }
  
  DL1 = lambdaFH * (S + fR * R) - (rhoL1 + mu) * L1;
  DL2 = (1 - fh) * rhoL1 * L1 + ((1 - fR) * lambdaFH + rhoR) * R + (1 - fs) * rhoH * H - (rhoL2 + mu) * L2;
  DH = fh * rhoL1 * L1 - (rhoH + mu) * H;
  DI11 = fd * fs * rhoH * H - (rhoI1 + mu + 3 * muVL) * I11;
  DI12 = 3 * muVL * I11 - (rhoI1 + mu + 3 * muVL) * I12;
  DI13 = 3 * muVL * I12 - (rhoI1 + mu + 3 * muVL) * I13;
  DI21 = (1 - fd) * fs * rhoH * H - (rhoI2 + mu + 3 * muVL) * I21;
  DI22 = 3 * muVL * I21 - (rhoI2 + mu + 3 * muVL) * I22;
  DI23 = 3 * muVL * I22 - (rhoI2 + mu + 3 * muVL) * I23;
  DD = rhoI1 * (I11 + I12 + I13) + rhoI2 * (I21 + I22 + I23) - (rhoD + mu) * D;
  DP = fp * rhoD * D - (rhoP + mu) * P;
  DR = rhoL2 * L2 + (1 - fp) * rhoD * D + rhoP * P
     - (rhoR + mu + lambdaFH) * R;
       
  DSF = (muF + muIRS) * NF * (1 + (season_amp * cos((t - season_t) * M_2PI))) * (1 - fIRS) - (lambdaHF + muF + muIRS) * SF; 
  DEF = lambdaHF * SF - (rhoEF + muF + muIRS) * EF;
  DIF = rhoEF * EF - (muF + muIRS) * IF;
"
  

vl_det_init <-
"
  double N_0;  // local variable for sum of initial human compartments sizes

  N_0 = S_0 + L1_0 + L2_0 + H_0
      + I11_0 + I12_0 + I13_0
      + I21_0 + I22_0 + I23_0
      + D_0 + P_0 + R_0;

  S = S_0/N_0;
  L1 = L1_0/N_0;
  L2 = L2_0/N_0;
  H = H_0/N_0;
  I11 = I11_0/N_0;
  I12 = I12_0/N_0;
  I13 = I13_0/N_0;
  I21 = I21_0/N_0;
  I22 = I22_0/N_0;
  I23 = I23_0/N_0;
  D = D_0/N_0;
  P = P_0/N_0;
  R = R_0/N_0;
  
  SF = SF_0;
  EF = EF_0;
  IF = IF_0;
"


# Stochastic compartmental model ----
# N.B.: this is old code that still needs to be updated to reflect the new model
#       structure (i.e., the deterministic model above).
counternames_sto <- c("VL_symp",
                      "VL_treat",
                      "VL_death",
                      "person_years")
statenames_sto <- c(statenames_det, counternames_sto)

vl_code_sto <- 
" 
// Declaring local variables 

  double N; // total human population size 
  double lambdaHF; // Force of infection from humans to flies 
  double lambdaFH; // Force of infection from flies to humans 
  double dSF; // Susceptible flies 
  double dEF; // Infected flies 
  double dIF; // Infectious flies 

  double e_S[2];  // number of individuals exiting from Susceptibles (S) compartment 
  double e_Latent1[3];  // number of exits from latent infection with low serological titre compartment 1 
  double e_Latent2[2];  // latent infection with low serological titre compartment 2 
  double e_H[4];    // latent infection with high serological titre 
  double e_IHE1[3]; // symptomatic infection compartment 1 with high probability of detection 
  double e_IHE2[3]; // symptomatic infection compartment 2 with high probability of detection 
  double e_IHE3[3]; // symptomatic infection compartment 3 with high probability of detection 
  double e_ILE1[3]; // symptomatic infection compartment 1 with low probability of detection 
  double e_ILE2[3]; // symptomatic infection compartment 2 with low probability of detection 
  double e_ILE3[3]; // symptomatic infection compartment 3 with low probability of detection 
  double e_Dormant[3];   
  double e_PKDL[2];   
  double e_R[4]; // Recovered 
  
  double rate_S[2];    
  double rate_Latent1[3];    
  double rate_Latent2[2];     
  double rate_H[4];    
  double rate_IHE1[3];     
  double rate_IHE2[3];     
  double rate_IHE3[3];     
  double rate_ILE1[3];     
  double rate_ILE2[3];     
  double rate_ILE3[3];   
  double rate_Dormant[3];    
  double rate_PKDL[2];    
  double rate_R[4];    
  
// Defining local variables 

  N = S + L1 + L2 + H + I11 + I12 + I13 + I21 + I22 + I23 + D + P + R;
  lambdaHF = beta_bite * (betaH * H + I11 + I12 + I13 + I21 + I22 + I23 + betaP * P) / N;  
  lambdaFH = beta_bite * pH * IF + foi_ext;   

  dSF = dt * ( (muF + muIRS) * NF * (1 + (season_amp * cos((t - season_t) * M_2PI))) * (1 - fIRS) - (lambdaHF + muF + muIRS) * SF ) ; 
  dEF = dt * ( lambdaHF * SF - (rhoEF + muF + muIRS) * EF ) ;  
  dIF = dt * ( rhoEF * EF - (muF + muIRS) * IF ) ;    
  
// Transitions (of discrete numbers) of humans exiting each compartment

  // S (susceptibles)
  rate_S[0] = lambdaFH;
  rate_S[1] = mu;
  reulermultinom(2, S, &rate_S[0], dt, &e_S[0]);

  // L1 (latent infection with low serological titre compartment 1) 
  rate_Latent1[0] = fh * rhoL1;            
  rate_Latent1[1] = (1 - fh) * rhoL1; 
  rate_Latent1[2] = mu;             
  reulermultinom(3, L1, &rate_Latent1[0], dt, &e_Latent1[0]);

  // L2 (latent infection with low serological titre compartment 2) 
  rate_Latent2[0] = rhoL2; 
  rate_Latent2[1] = mu;     
  reulermultinom(2, L2, &rate_Latent2[0], dt, &e_Latent2[0]);

  // H (latent infection with high serological titre)  
  rate_H[0] = fd * fs * rhoH;         // progression to symptomatic phase 1    
  rate_H[1] = (1-fd) * fs * rhoH;      
  rate_H[2] = (1-fs) * rhoH;    
  rate_H[3] = mu;     
  reulermultinom(4, H, &rate_H[0], dt, &e_H[0]); 

  // IHE1 (symptomatic infection compartment 1 with high probability of detection)  
  rate_IHE1[0] = rhoI1;                    // progression to  dormant stage 
  rate_IHE1[1] = mu;                       // death due to other causes
  rate_IHE1[2] = 3 * muVL; 
  reulermultinom(3, I11, &rate_IHE1[0], dt, &e_IHE1[0]);

  // IHE2 (symptomatic infection compartment 2 with high probability of detection)  
  rate_IHE2[0] = rhoI1;                    // progression to  dormant stage  
  rate_IHE2[1] = mu;                      // death due to other causes
  rate_IHE2[2] = 3 * muVL;                // death due to VL
  reulermultinom(3, I12, &rate_IHE2[0], dt, &e_IHE2[0]);

  // IHE3 (symptomatic infection compartment 3 with high probability of detection)  
  rate_IHE3[0] = rhoI1;                   // progression to  dormant stage
  rate_IHE3[1] = mu;                      // death due to other causes
  rate_IHE3[2] = 3 * muVL;                // death due to VL 
  reulermultinom(3, I13, &rate_IHE3[0], dt, &e_IHE3[0]);

  // ILE1 (symptomatic infection compartment 1 with low probability of detection)
  rate_ILE1[0] = rhoI2;                   // progression to  dormant stage 
  rate_ILE1[1] = mu;                      // death due to other causes
  rate_ILE1[2] = 3 * muVL;                // death due to VL
  reulermultinom(3, I21, &rate_ILE1[0], dt, &e_ILE1[0]);

  // ILE2 (symptomatic infection compartment 2 with low probability of detection)
  rate_ILE2[0] = rhoI2;                  // progression to  dormant stage
  rate_ILE2[1] = mu;                     // death due to other causes
  rate_ILE2[2] = 3 * muVL;               // death due to VL
  reulermultinom(3, I22, &rate_ILE2[0], dt, &e_ILE2[0]);

  // ILE3 (symptomatic infection compartment 3 with low probability of detection)
  rate_ILE3[0] = rhoI2;                 // progression to  dormant stage
  rate_ILE3[1] = mu;                    // death due to other causes
  rate_ILE3[2] = 3 * muVL;              // death due to VL
  reulermultinom(3, I23, &rate_ILE3[0], dt, &e_ILE3[0]);  

  // D (dormant)
  rate_Dormant[0] = (1 - fp) * rhoD;       // proportion with full recovery 
  rate_Dormant[1] = fp * rhoD;             // proportion developing PKDL
  rate_Dormant[2] = mu;                    // death due to other causes  
  reulermultinom(3, D, &rate_Dormant[0], dt, &e_Dormant[0]); 

  // PKDL  
  rate_PKDL[0] = rhoP;                  // full recovery   
  rate_PKDL[1] = mu;                    // death due to other causes
  reulermultinom(2, P, &rate_PKDL[0], dt, &e_PKDL[0]);  

  // R (recovered)  
  rate_R[0] = fR * lambdaFH;                 
  rate_R[1] = (1 - fR) * lambdaFH;                
  rate_R[2] = rhoR;                  // seroreversion rate
  rate_R[3] = mu;                    // death due to other causes 
  reulermultinom(4, R, &rate_R[0], dt, &e_R[0]);  // System of ODEs

  // System of ODEs
  S += - e_S[0] - e_S[1]; 
  if (cohort_model == 0) {
    S += e_Latent1[2] + e_Latent2[1] + e_H[3] + e_IHE1[1] + e_IHE2[1] + e_IHE3[1] + e_ILE1[1] + e_ILE2[1] + e_ILE3[1] + e_Dormant[2] + e_PKDL[1] + e_R[3] + e_IHE3[2] + e_ILE3[2] + e_S[1]; 
  }
  L1 += e_S[0] + e_R[0] - e_Latent1[0] - e_Latent1[1] - e_Latent1[2]; 
  L2 += e_Latent1[1] + e_R[1] + e_R[2] + e_H[2] - e_Latent2[0] - e_Latent2[1];  
  H += e_Latent1[0] - e_H[0] - e_H[1] - e_H[2] - e_H[3]; 
  I11 += e_H[0] - e_IHE1[0] - e_IHE1[1] - e_IHE1[2];   
  I12 += e_IHE1[2] - e_IHE2[0] - e_IHE2[1] - e_IHE2[2];   
  I13 += e_IHE2[2] - e_IHE3[0] - e_IHE3[1] - e_IHE3[2]; 
  I21 += e_H[1] - e_ILE1[0] - e_ILE1[1] - e_ILE1[2];  
  I22 += e_ILE1[2] - e_ILE2[0] - e_ILE2[1] - e_ILE2[2];
  I23 += e_ILE2[2] - e_ILE3[0] - e_ILE3[1] - e_ILE3[2]; 
  D += e_IHE1[0] + e_IHE2[0] + e_IHE3[0] + e_ILE1[0] + e_ILE2[0] + e_ILE3[0] - e_Dormant[0] - e_Dormant[1] - e_Dormant[2];
  P += e_Dormant[1] - e_PKDL[0] - e_PKDL[1];   
  R += e_Latent2[0] + e_Dormant[0] + e_PKDL[0] - e_R[0] - e_R[1] - e_R[2] - e_R[3]; 

  // Sandfly compartments 
  SF += dSF ; 
  EF += dEF ;  
  IF += dIF ;    
  
  // Accumulators to record numbers of incident cases per time window
  VL_symp += e_H[0] + e_H[1];     // True incidence of visceral leishmaniasis 
  VL_treat += e_IHE1[0] + e_IHE2[0] + e_IHE3[0] + e_ILE1[0] + e_ILE2[0] + e_ILE3[0];   
  VL_death += e_IHE3[2] + e_ILE3[2];   // Deaths due to untreated visceral leishmaniasis
  person_years += dt * N;
" 

vl_sto_init <-
"
  // Compartments
  S = S_0;
  L1 = L1_0;
  L2 = L2_0;
  H = H_0;
  I11 = I11_0;
  I12 = I12_0;
  I13 = I13_0;
  I21 = I21_0;
  I22 = I22_0;
  I23 = I23_0;
  D = D_0;
  P = P_0;
  R = R_0;
  SF = SF_0;
  EF = EF_0;
  IF = IF_0;
  
  // Accumulators (numerators for incidences)
  VL_symp = 0;       // True incidence of visceral leishmaniasis
  VL_death = 0;      // Deaths due to untreated visceral leishmaniasis
  VL_treat = 0;      // Detected and treated cases
  person_years = 0;
"

### END OF CODE ### ----

