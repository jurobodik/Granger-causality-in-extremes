# Comparison of our method with other state-of-the-art methods on VAR and GARCH models.
# This script creates the final comparative plot.
# It uses the results output for each alternative method by the corresponding scripts.

# library(ExtremeGranger)
source("./R/utils.R")

library(EnvStats)
library(ExtremeRisks)
library(igraph)

library(tidyverse)



number_of_simulations <- 100

sequence_of_p <- c(2, 4, 7, 10, 20)
sequence_of_n <- c(500, 5000)
sequence_of_p_GPDC = c(2, 4, 7) # GPDC is too time-consuming so we only have results for p=2,4,7
dat_structures <- c('VAR', 'GARCH')

tables_folder <- './Results/Simulations_comparison/tables/' # folder where results are stored
plots_folder <- './Results/Simulations_comparison/plots/' # folder where plots will be stored


set.seed(0)

# ================ MAIN FIGURE ================

data <-data.frame(
  method = character(),
  structure = character(),
  heavy_tailed = logical(),
  n = integer(),
  p = integer(),
  replication = integer(),
  distance = numeric()
)

for (met in c("ExtremeCausal", "PCMCI_Cor", "PCMCI_GPDC")) {
  for (dat_struct in dat_structures) {
    for (ht in c("heavy", "light")) {
      table_name <- paste0('results_', met, '_', dat_struct, '_', ht, '.csv')
      # assign(paste0('results_', met, '_', dat_struct, '_', ht), read.csv(paste0(tables_folder, table_name)))
      temp_df <- read.csv(paste0(tables_folder, table_name))
      # temp_df$method <- met
      # temp_df$structure <- dat_struct
      # temp_df$heavy_tailed <- ifelse(ht == "heavy", TRUE, FALSE)
      data <- rbind(data, temp_df)
    }
  }
}
# replace "True" and "TRUE" in heavy_tailed column by TRUE and "False" and "FALSE" by FALSE
# data$heavy_tailed <- ifelse(data$heavy_tailed %in% c("True", "TRUE"), TRUE, ifelse(data$heavy_tailed %in% c("False", "FALSE"), FALSE, data$heavy_tailed))
data$heavy_tailed <- as.logical(data$heavy_tailed)


summarize_and_standardize_data <- function(df) {
  # group df by method, structure, heavy_tailed, n and p and compute mean distance and standard error for each group
  df_summary <- df %>%
    group_by(method, structure, heavy_tailed, n, p) %>%
    summarise(
      mean_distance = mean(distance),
      std_error = sd(distance) / sqrt(n()),
      # std_error = sd(distance) / sqrt(number_of_simulations),
      std_dev = sd(distance),
      q_up = quantile(distance, 0.9),
      q_lo = quantile(distance, 0.1),
      nb_replications = n(),
      .groups = 'drop' # to ungroup after summarise
    )

  # Here are the resulting numbers that you should get by running the previous codes
  # We divide everything by p(p-1) because there are p*(p-1) arrows to estimate
  standardisation <- sequence_of_p * (sequence_of_p - 1)
  standardisation_GPDC <- sequence_of_p_GPDC * (sequence_of_p_GPDC - 1)

  # standardise the mean_distance and std_error in df_summary by p(p-1)
  # use p=standardisation_GPDC when method is PCMCI_GPDC and p=standardisation otherwise
  df_summary <- df_summary %>%
    mutate(
      mean_distance = ifelse(method == "PCMCI_GPDC", 
                            mean_distance / standardisation_GPDC[match(p, sequence_of_p_GPDC)], 
                            mean_distance / standardisation[match(p, sequence_of_p)]),
      std_error = ifelse(method == "PCMCI_GPDC", 
                        std_error / standardisation_GPDC[match(p, sequence_of_p_GPDC)], 
                        std_error / standardisation[match(p, sequence_of_p)]),
      std_dev = ifelse(method == "PCMCI_GPDC", 
                        std_dev / standardisation_GPDC[match(p, sequence_of_p_GPDC)], 
                        std_dev / standardisation[match(p, sequence_of_p)]),
      q_up = ifelse(method == "PCMCI_GPDC", 
                    q_up / standardisation_GPDC[match(p, sequence_of_p_GPDC)],
                    q_up / standardisation[match(p, sequence_of_p)]),
      q_lo = ifelse(method == "PCMCI_GPDC",
                    q_lo / standardisation_GPDC[match(p, sequence_of_p_GPDC)],
                    q_lo / standardisation[match(p, sequence_of_p)])
    )
  
  # df_summary <- df_summary %>%
  #   mutate(
  #     mean_distance = mean_distance / standardisation[match(p, sequence_of_p)],
  #     std_error = std_error / standardisation[match(p, sequence_of_p)]
  #   )
  
  return(df_summary)
}

data_summary <- summarize_and_standardize_data(data)


# add the Random models to data_summary
Random_model_df <- expand.grid(method = "Random", structure = dat_structures, heavy_tailed = c(TRUE, FALSE), n = sequence_of_n, p = sequence_of_p) %>%
  mutate(
    mean_distance = 0.5, # we set the mean distance for random model to 0.5 (this is an approximation based on the results we got for random model)
    std_error = 0, # we set the standard error for random model to 0
    std_dev = 0, # we set the standard deviation for random model to 0
    # std_error = 0.025 / standardisation[match(p, sequence_of_p)] # we set the standard error for random model to 0.025 (this is an approximation based on the results we got for random model) and we also standardise it by p(p-1)
    q_up = 0.5,
    q_lo = 0.5,
    nb_replications = 0 # 0 number_of_simulations
  ) %>%
  arrange(method, structure, heavy_tailed, n, p)

data_summary <- bind_rows(data_summary, Random_model_df) %>%
  arrange(method, structure, heavy_tailed, n, p)

# Create method_label, setting_label and n_label columns in data_summary
data_summary <- data_summary %>%
  mutate(
    method_label = case_when(
      (method == "ExtremeCausal") ~ "Our method",
      (method == "PCMCI_Cor") ~ "PCMCI cor",
      (method == "PCMCI_GPDC") ~ "PCMCI gpdc",
      (method == "Random") ~ "Random"
    ),
    setting_label = case_when(
      (structure == "VAR" & heavy_tailed) ~ "VAR heavy-tailed",
      (structure == "GARCH" & heavy_tailed) ~ "GARCH heavy-tailed",
      (structure == "VAR" & !heavy_tailed) ~ "VAR Gaussian",
      (structure == "GARCH" & !heavy_tailed) ~ "GARCH Gaussian"
    ),
    n_label = case_when(
      (n == 500) ~ "n = 500",
      (n == 5000) ~ "n = 5000"
    )
  )

data_summary <- data_summary %>%
  mutate(
    method_label = factor(method_label, levels = c("Our method", "PCMCI cor", "PCMCI gpdc", "Random")),
    setting_label = factor(setting_label, levels = c("VAR heavy-tailed", "GARCH heavy-tailed", "VAR Gaussian", "GARCH Gaussian")),
    n_label = factor(n_label, levels = c("n = 500", "n = 5000"))
  )



# Plotting
compar_plot <- ggplot(data_summary, aes(x = p, y = mean_distance, color = method_label, fill = method_label)) +
  # geom_ribbon(aes(ymin = mean_distance - std_error, ymax = mean_distance + std_error), alpha = 0.2, color = NA) +
  # geom_ribbon(aes(ymin = mean_distance - std_dev, ymax = mean_distance + std_dev), alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin = q_lo, ymax = q_up), alpha = 0.2, color = NA) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Average error") +
  xlab("Number of variables") +
  facet_grid(n_label ~ setting_label) + # Transposed facets
  scale_colour_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
  scale_fill_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 12), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(linewidth=0, linetype = "blank", fill="grey90"),# element_blank(),
    
    legend.position = "bottom", # Move legend to the bottom
    # panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.major = element_line(size = 0.4, color = "grey80"),
    panel.grid.minor = element_blank(),
    legend.box.spacing = margin(0.6),
    # plot.margin = margin(0,0,0,0, unit="pt"),
    aspect.ratio = 1
  ) + # Ensure each facet is squared
  scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) + # Regular breaks and limit for y-axis
  # limits without cutting the error bars
  coord_cartesian(ylim = c(0, 0.5))


ggsave(paste0(plots_folder, "Simulations_comparison_plot_main.pdf"), plot = compar_plot, width = 8.5, height = 7)
knitr::plot_crop(paste0(plots_folder, "Simulations_comparison_plot_main.pdf"))
plot(compar_plot)



# ================ PARETO(2) SUPPLEMENTARY FIGURE ================

data2 <-data.frame(
  method = character(),
  structure = character(),
  heavy_tailed = logical(),
  n = integer(),
  p = integer(),
  replication = integer(),
  distance = numeric()
)

for (met in c("ExtremeCausal", "PCMCI_Cor", "PCMCI_GPDC")) {
  for (dat_struct in c("VAR")) {
    for (ht in c("Pa2")) {
      table_name <- paste0('results_', met, '_', dat_struct, '_', ht, '.csv')
      # assign(paste0('results_', met, '_', dat_struct, '_', ht), read.csv(paste0(tables_folder, table_name)))
      temp_df <- read.csv(paste0(tables_folder, table_name))
      # temp_df$method <- met
      # temp_df$structure <- dat_struct
      # temp_df$heavy_tailed <- ifelse(ht == "heavy", TRUE, FALSE)
      data2 <- rbind(data2, temp_df)
    }
  }
}
# replace "True" and "TRUE" in heavy_tailed column by TRUE and "False" and "FALSE" by FALSE
# data2$heavy_tailed <- ifelse(data2$heavy_tailed %in% c("True", "TRUE"), TRUE, ifelse(data2$heavy_tailed %in% c("False", "FALSE"), FALSE, data2$heavy_tailed))
data2$heavy_tailed <- as.logical(data2$heavy_tailed)

data_summary2 <- summarize_and_standardize_data(data2)


# add the Random models to data_summary2
Random_model_df2 <- expand.grid(method = "Random", structure = c("VAR"), heavy_tailed = c(TRUE), n = sequence_of_n, p = sequence_of_p) %>%
  mutate(
    mean_distance = 0.5, # we set the mean distance for random model to 0.5 (this is an approximation based on the results we got for random model)
    std_error = 0, # we set the standard error for random model to 0
    std_dev = 0, # we set the standard deviation for random model to 0
    # std_error = 0.025 / standardisation[match(p, sequence_of_p)] # we set the standard error for random model to 0.025 (this is an approximation based on the results we got for random model) and we also standardise it by p(p-1)
    q_up = 0.5,
    q_lo = 0.5,
    nb_replications = 0 # 0 number_of_simulations
  ) %>%
  arrange(method, structure, heavy_tailed, n, p)

data_summary2 <- bind_rows(data_summary2, Random_model_df2) %>%
  arrange(method, structure, heavy_tailed, n, p)

# Create method_label, setting_label and n_label columns in data_summary2
data_summary2 <- data_summary2 %>%
  mutate(
    method_label = case_when(
      (method == "ExtremeCausal") ~ "Our method",
      (method == "PCMCI_Cor") ~ "PCMCI cor",
      (method == "PCMCI_GPDC") ~ "PCMCI gpdc",
      (method == "Random") ~ "Random"
    ),
    setting_label = case_when(
      (structure == "VAR" & heavy_tailed) ~ "VAR Pareto(2)"
    ),
    n_label = case_when(
      (n == 500) ~ "n = 500",
      (n == 5000) ~ "n = 5000"
    )
  )

data_summary2 <- data_summary2 %>%
  mutate(
    method_label = factor(method_label, levels = c("Our method", "PCMCI cor", "PCMCI gpdc", "Random")),
    setting_label = factor(setting_label, levels = c("VAR Pareto(2)")),
    # setting_label = factor(setting_label, levels = c("VAR Pareto(2)", "GARCH Pareto(2)", "VAR Gaussian", "GARCH Gaussian")),
    n_label = factor(n_label, levels = c("n = 500", "n = 5000"))
  )



# Plotting
compar_plot <- ggplot(data_summary2, aes(x = p, y = mean_distance, color = method_label, fill = method_label)) +
  # geom_ribbon(aes(ymin = mean_distance - std_error, ymax = mean_distance + std_error), alpha = 0.2, color = NA) +
  # geom_ribbon(aes(ymin = mean_distance - std_dev, ymax = mean_distance + std_dev), alpha = 0.2, color = NA) +
  geom_ribbon(aes(ymin = q_lo, ymax = q_up), alpha = 0.2, color = NA) +
  geom_line(size = 0.7, linetype = "solid") +
  geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
  ylab("Average error") +
  xlab("Number of variables") +
  # facet_grid(n_label ~ setting_label) + # Transposed facets
  facet_grid(setting_label ~ n_label) + # Transposed facets
  scale_colour_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
  scale_fill_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.text = element_text(size = 10), # Smaller font size for facet labels
    axis.text = element_text(size = 12), # Smaller font size for axis text
    axis.title = element_text(size = 12), # Smaller font size for axis titles
    legend.text = element_text(size = 10), # Smaller font size for legend text
    legend.title = element_blank(), # Hide legend title
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(linewidth=0, linetype = "blank", fill="grey90"),# element_blank(),
    
    legend.position = "bottom", # Move legend to the bottom
    # panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
    panel.grid.major = element_line(size = 0.4, color = "grey80"),
    panel.grid.minor = element_blank(),
    legend.box.spacing = margin(0.6),
    # plot.margin = margin(0,0,0,0, unit="pt"),
    aspect.ratio = 1
  ) + # Ensure each facet is squared
  scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) + # Regular breaks and limit for y-axis
  # limits without cutting the error bars
  coord_cartesian(ylim = c(0, 0.5))


ggsave(paste0(plots_folder, "Simulations_comparison_plot_Pa2.pdf"), plot = compar_plot, width = 4.5, height = 5)
knitr::plot_crop(paste0(plots_folder, "Simulations_comparison_plot_Pa2.pdf"))
plot(compar_plot)




