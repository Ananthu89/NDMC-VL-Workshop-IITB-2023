# Script to run and compare the deterministic and stochastic versions of the
# model for VL transmission with time-varying parameters

# Author: Ananthu James, Luc Coffeng
# Date created: 10 August 2022, updated 11 Aug 2023

# Key assumptions:
#  - fixed population size
#  - no age structure
#  - determinstic or stochastic human population
#  - deterministic fly population


### Prep session ----
rm(list = ls())

library(ggplot2)
library(rstudioapi)  # requires RStudio 
library(data.table)
library(gridExtra)
library(pomp) 

theme_set(theme_bw() + theme(panel.spacing = unit(2, "lines")))
current_dir <- dirname(getActiveDocumentContext()$path)
base_dir <- file.path(current_dir, "..")

setwd(file.path(base_dir, "00_Source"))
source("model.r")
source("param_values.r")


### Comparison of deterministic and stochastic model with sandfly transmission present ----
## Model setup ----
# General

# Dropping the variable(s) for time-dependent parameters        
param_values_no_fIRS <- param_values[!(names(param_values) %in% c("fIRS"))]   

# Defining "param_names" with the updated values 
param_names <- c(names(param_values_no_fIRS), names(init_values))

# Changing parameter values
param_values_no_fIRS["NF"] <- 0.3           # Sandfly to human ratio (allowed values: >=0) 
param_values_no_fIRS["season_amp"] <- 1      # Relative amplitude of seasonality (allowed values: 0->1) 
param_values_no_fIRS["season_t"] <- 8 / 12   # Timing of peak sandfly-abundance 

start_time <- 0
sim_duration <- 100 # years

# Set time-varying parameters representing changes 
policy_fIRS <- covariate_table(fIRS = c(0, 0.5 , 0.5),
                             order = "constant",
                             times = c(start_time, 50, sim_duration + 1))  

# Deterministic model
dt_output_det <- 1 / 365  # interval for output generation
time_seq_det <- seq(start_time, start_time + sim_duration, by = dt_output_det)
pars_det <- c(param_values_no_fIRS, init_values)

# Stochastic model
timestep_sto <- 1 / 365     # step size for simulation
dt_output_sto <- 1          # interval for output generation (in years)
time_seq_sto <- seq(start_time, start_time + sim_duration, by = dt_output_sto)

nsim <- 1 
seed <- 1914679908L

pop_size <- 1e6
init_values_sto <- pop_size * init_values 
pars_sto <- c(param_values_no_fIRS, init_values_sto) 

# Compile models
model_det <- pomp(data = data.frame(time = time_seq_det),
                  times = "time",
                  covar = policy_fIRS,
                  t0 = start_time,
                  skeleton = vectorfield(Csnippet(vl_code_det)),
                  rinit = Csnippet(vl_det_init),
                  statenames = statenames_det,
                  paramnames = param_names) 

model_sto <- pomp(data = data.frame(time = time_seq_sto), 
                  times = "time",    
                  t0 = start_time,
                  covar = policy_fIRS,
                  rinit = Csnippet(vl_sto_init),
                  rprocess = euler(step.fun = Csnippet(vl_code_sto), 
                                   delta.t = timestep_sto),
                  statenames = statenames_sto, 
                  accumvars = counternames_sto,
                  paramnames = param_names)


## Simulate trends with deterministic model ----
output_det <- trajectory(model_det, params = pars_det, format = "data.frame")
output_det <- data.table(output_det)

# Calculate prevalences
col_names_output_det <- names(output_det)[1:13]
output_det[, N := Reduce("+", mget(col_names_output_det))]
output_det[, I := I11 + I12 + I13 + I21 + I22 + I23]

col_names_output_det <- c(col_names_output_det, "I")
output_det[, (paste0(col_names_output_det, "_prev")) :=
             lapply(.SD, function(x) x / N * 100),
           .SDcols = col_names_output_det]

# Calculate incidences
with(as.list(pars_det), {
  
  output_det[, I11_detected := rhoI1 * I11]
  output_det[, I12_detected := rhoI1 * I12]
  output_det[, I13_detected := rhoI1 * I13]
  output_det[, I21_detected := rhoI2 * I21]
  output_det[, I22_detected := rhoI2 * I22]
  output_det[, I23_detected := rhoI2 * I23]
  
  output_det[, VL_symp_inc := (fs * rhoH * H) / N]
  
  output_det[, VL_death_1 := 3 * muVL * I13]
  output_det[, VL_death_2 := 3 * muVL * I23]
  
})

output_det[, VL_treat_inc := Reduce("+",
                                    mget(grep("_detected",
                                              names(output_det),
                                              value = TRUE))) / N]
output_det[, (grep("_detected", names(output_det), value = TRUE)) := NULL]

output_det[, VL_death_inc := (VL_death_1 + VL_death_2) / N]
output_det[, c("VL_death_1", "VL_death_2") := NULL]

# Aggregate incidences by same interval such that numbers represent events in
# the time interval that has just past, like in the stochastic model output
# (this makes the impact of the time-varying parameter more visible by filtering
# out the seasonal pattern).
col_names_inc_det <- c("VL_symp_inc", "VL_treat_inc", "VL_death_inc")
output_det[, time_ceiling := ceiling(time / dt_output_sto) * dt_output_sto]
output_det <- output_det[, (col_names_inc_det) := lapply(.SD, mean),
                         .SDcols = col_names_inc_det,
                         by = time_ceiling]
output_det <- output_det[, .SD[.N], by = time_ceiling]
output_det[, time_ceiling := NULL]
# Note: prevalences are not averaged over intervals, but still represent the
# state of the population at the exact time recorded in the output object


## Simulate trends with population model ----
output_sto <- simulate(object = model_sto,
                       params = pars_sto,
                       nsim = nsim,
                       seed = seed,
                       format = "data.frame")
output_sto <- data.table(output_sto)

# Calculate prevalences
col_names_output_sto <- names(output_sto)[3:15]
output_sto[, N := Reduce("+", mget(col_names_output_sto))]
output_sto[, I := I11 + I12 + I13 + I21 + I22 + I23]

col_names_output_sto <- c(col_names_output_sto, "I")
output_sto[, (paste0(col_names_output_sto, "_prev")) :=
             lapply(.SD, function(x) x / N * 100),
           .SDcols = col_names_output_sto]
output_sto[, N := N / max(N), by = .id]  # for comparability to deterministic model 

# Calculate incidences
col_names_inc_sto <- counternames_sto[counternames_sto != "person_years"]
output_sto[, (paste0(col_names_inc_sto, "_inc")) :=
             lapply(.SD, function(x) x / person_years),
           .SDcols = col_names_inc_sto]


## Function for plotting ----
# Function to plot two data.tables together in the same panel
plot_two_vars <- function(df_sto,
                          df_det,
                          x_variable,
                          y_variable,
                          x_label,
                          y_label,
                          guide = "none") {
  df_sto_subset <- df_sto[, mget(c(x_variable, y_variable))]
  df_det_subset <- df_det[, mget(c(x_variable, y_variable))]
  df_sto_subset$dataset <- "Stochastic"
  df_det_subset$dataset <- "Deterministic"
  combined_df <- rbindlist(list(df_sto_subset, df_det_subset),
                           use.names = TRUE, fill = TRUE)
  
  plot <- ggplot(data = combined_df,
         mapping = aes_string(x = x_variable,
                              y = y_variable,
                              linetype = "dataset",
                              color = "dataset")) +
    geom_line() + 
    expand_limits(y = 0) +  
    labs(x = x_label,
         y = y_label,
         linetype = "Model",
         color = "Model") +
    expand_limits(y = 0)
  
  if (guide == "none") {
    plot <- plot +
      scale_linetype_discrete(guide = "none") + 
      scale_color_discrete(guide = "none")
  }
  
  plot
  
}


## Create plots comparing deterministic and stochastic model ----
setwd(file.path(base_dir, "Figures_population_model"))

# Prevalences
plot_N <- plot_two_vars(output_sto, output_det,
                        "time", "N",
                        "Time (years)", "Population size (relative)")

plot_S <- plot_two_vars(output_sto, output_det,
                        "time", "S_prev",
                        "Time (years)", "Susceptible (%)")

plot_H <- plot_two_vars(output_sto, output_det,
                        "time", "H_prev",
                        "Time (years)", "High titre (%)")

plot_I <- plot_two_vars(output_sto, output_det,
                        "time", "I_prev",
                        "Time (years)", "Untreated VL (%)")

plot_L1 <- plot_two_vars(output_sto, output_det,
                         "time", "L1_prev",
                         "Time (years)", "Low titre 1 (%)")

plot_L2 <- plot_two_vars(output_sto, output_det,
                        "time", "L2_prev",
                        "Time (years)", "Low titre 2 (%)")

plot_P <- plot_two_vars(output_sto, output_det,
                        "time", "P_prev",
                        "Time (years)", "PKDL (%)")

plot_D <- plot_two_vars(output_sto, output_det,
                        "time", "D_prev",
                        "Time (years)", "Dormant (%)")

plot_R <- plot_two_vars(output_sto, output_det,
                        "time", "R_prev",
                        "Time (years)", "Recovered (%)") 


# Incidences
plot_VL_symp <- plot_two_vars(output_sto, output_det,
                              "time", "VL_symp_inc",
                              "Time (years)", "Symptomatic VL (inc.)")

plot_VL_treat <- plot_two_vars(output_sto, output_det,
                              "time", "VL_treat_inc",
                              "Time (years)", "Treated VL (inc.)")

plot_VL_treat_zoom <- plot_VL_treat + coord_cartesian(x = c(45, 55))

plot_VL_death <- plot_two_vars(output_sto, output_det,
                               "time", "VL_death_inc",
                               "Time (years)", "VL-deaths (inc.)")

plot_VL_death_zoom <- plot_VL_death + coord_cartesian(x = c(45, 55))


# Print to disk
ggsave(filename = sprintf("NF%.2f-betaH%.2f-foi_ext%.3f_det_vs_stoch_over_time_fIRS.pdf",
                           param_values_no_fIRS["NF"],
                           param_values_no_fIRS["betaH"],
                           param_values_no_fIRS["foi_ext"]),
              plot = grid.arrange(plot_N, plot_H,  
                           plot_L1, plot_L2, plot_I, plot_S, plot_R,   
                           plot_D, plot_P,
                           plot_VL_symp, plot_VL_treat, plot_VL_treat_zoom, 
                           plot_VL_death, plot_VL_death_zoom, 
                           ncol = 2),
       width = 8, height = 14)


### END OF CODE ### ----
