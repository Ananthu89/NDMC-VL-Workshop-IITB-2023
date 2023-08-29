# Script to run and compare the deterministic and stochastic versions of the
# model for VL transmission (fixed force of infection only, no sandfly dynamics)

# Author: Ananthu James, Luc Coffeng
# Date created: 13 July 2022, updated 4 Aug 2023

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
base_dir <- dirname(getActiveDocumentContext()$path)

for (folder_name in c("Figures_population_model", "Figures_cohort_model")) {
  setwd(base_dir)
  if (!file.exists(folder_name)) {
    dir.create(folder_name)
    cat(paste("Folder '", folder_name, "' created successfully.", sep = ""))
  } else {
    cat(paste("Folder '", folder_name, "' already exists.", sep = ""))
  }
}

setwd(file.path(base_dir, "00_Source"))
source("model.r")
source("param_values.r")


### Comparison of deterministic and stochastic model with external force of infection ----
## Model setup ----
# General

# Changing parameter values
param_values["cohort_model"] <- 1   # Toggle: 0 = population model; 1 = cohort model
param_values["pH"] <- 0             # Transmission probability from fly to human  
param_values["foi_ext"] <- 0.025    # External force of infection 

start_time <- 0
sim_duration <- 50  # years

# Deterministic model
dt_output_det <- 1 / 365  # interval for output generation
time_seq_det <- seq(start_time, start_time + sim_duration, by = dt_output_det)
pars_det <- c(param_values, init_values)

# Stochastic model
timestep_sto <- 1 / 365  # step size for simulation
dt_output_sto <- 1       # interval for output generation (in years)
time_seq_sto <- seq(start_time, start_time + sim_duration, by = dt_output_sto)

nsim <- 1 
seed <- 1914679908L

pop_size <- 1e6
init_values_sto <- pop_size * init_values
pars_sto <- c(param_values, init_values_sto)

# Compile models
model_det <- pomp(data = data.frame(time = time_seq_det),
                  times = "time",
                  t0 = start_time,
                  skeleton = vectorfield(Csnippet(vl_code_det)),
                  rinit = Csnippet(vl_det_init),
                  statenames = statenames_det,
                  paramnames = paramnames) 

model_sto <- pomp(data = data.frame(time = time_seq_sto), 
                  times = "time",    
                  t0 = start_time,
                  rinit = Csnippet(vl_sto_init),
                  rprocess = euler(step.fun = Csnippet(vl_code_sto), 
                                   delta.t = timestep_sto),
                  statenames = statenames_sto, 
                  accumvars = counternames_sto,
                  paramnames = paramnames
)


## Simulate with deterministic model ----
output_det <- trajectory(model_det,
                         params = pars_det,
                         format = "data.frame")
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


## Simulate with stochastic model ----
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
stopifnot(param_values["cohort_model"] %in% c(0, 1))
if (param_values["cohort_model"] == 0) {
  setwd(file.path(base_dir, "Figures_population_model"))
} else {
  setwd(file.path(base_dir, "Figures_cohort_model"))
}

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


# Print to disk
ggsave(filename = sprintf("NF%.2f-foi_ext%.3f_det_vs_stoch_over_time.pdf",
                          param_values["NF"],
                          param_values["foi_ext"]),
       plot = grid.arrange(plot_N, plot_H, plot_S, plot_H, plot_I, plot_D,
                           plot_R,
                           plot_VL_symp, plot_VL_treat,
                           ncol = 2),
       width = 8, height = 10)


### END OF CODE ### ----
