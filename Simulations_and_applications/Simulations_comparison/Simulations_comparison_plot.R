# Comparison of our method with other state-of-the-art methods on VAR and GARCH models.
# This script creates the final comparative plot.
# It uses the results output for each alternative method by the corresponding scripts.

source("./R/Main_functions.R")
source("./R/utils.R")
library(EnvStats)
library(ExtremeRisks)
library(igraph)

library(tidyverse)
# library(ggplot2)
# library(dplyr)
# library(tidyr)



number_of_simulations <- 100

sequence_of_p <- c(2, 4, 7, 10, 20)
sequence_of_n <- c(500, 5000)
sequence_of_p_GPDC = c(2, 4, 7) # GPDC is too time-consuming so we only have results for p=2,4,7
dat_structures <- c('VAR', 'GARCH')

tables_folder <- './Results/Simulations_comparison/tables/' # folder where results are stored
plots_folder <- './Results/Simulations_comparison/plots/' # folder where plots will be stored


set.seed(0)

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


# group data by method, structure, heavy_tailed, n and p and compute mean distance and standard error for each group
data_summary <- data %>%
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

# standardise the mean_distance and std_error in data_summary by p(p-1)
# use p=standardisation_GPDC when method is PCMCI_GPDC and p=standardisation otherwise
data_summary <- data_summary %>%
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

# data_summary <- data_summary %>%
#   mutate(
#     mean_distance = mean_distance / standardisation[match(p, sequence_of_p)],
#     std_error = std_error / standardisation[match(p, sequence_of_p)]
#   )


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
      # (structure == "VAR" & heavy_tailed == TRUE) ~ "VAR heavy-tailed",
      # (structure == "GARCH" & heavy_tailed == TRUE) ~ "GARCH heavy-tailed",
      # (structure == "VAR" & heavy_tailed == FALSE) ~ "VAR Gaussian",
      # (structure == "GARCH" & heavy_tailed == FALSE) ~ "GARCH Gaussian"
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

# plot(compar_plot)

ggsave(paste0(plots_folder, "Simulations_comparison_plot_new.pdf"), plot = compar_plot, width = 8.5, height = 7)
knitr::plot_crop(paste0(plots_folder, "Simulations_comparison_plot_new.pdf"))
plot(compar_plot)


# plot.background = element_blank(),
# panel.background = element_blank(),
# legend.background = element_blank(),
# strip.background = element_blank(),
# strip.placement = "outside",
# strip.text = element_text(size = font_size_axes),
# strip.switch.pad.grid = unit(0, "pt"),
# strip.switch.pad.wrap = unit(0, "pt"),
# plot.caption=element_text(size=font_size_captions, hjust=0, 
#                           margin=margin(t=15)),
# text = element_text(size = font_size),
# axis.ticks = element_blank(),
# axis.text = element_text(size = font_size_axes),#
# panel.grid.major = element_line(size = 0.25),
# legend.position = "bottom",
# legend.box.spacing = margin(0.5),
# plot.margin=margin(0,0,0,0, unit="pt")

# compar_plot <- ggplot(data_summary, aes(x = p, y = mean_distance, color = method_label, fill = method_label)) +
#   geom_ribbon(aes(ymin = mean_distance - std_error, ymax = mean_distance + std_error), alpha = 0.2, color = NA) +
#   geom_line(size = 0.7, linetype = "solid") +
#   geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
#   ylab("Average error") +
#   xlab("Number of variables") +
#   facet_grid(setting_label ~ n_label) + # Transposed facets
#   scale_colour_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
#   scale_fill_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
#   theme(
#     strip.text = element_text(size = 10), # Smaller font size for facet labels
#     axis.text = element_text(size = 12), # Smaller font size for axis text
#     axis.title = element_text(size = 12), # Smaller font size for axis titles
#     legend.text = element_text(size = 10), # Smaller font size for legend text
#     legend.title = element_blank(), # Hide legend title
#     legend.position = "right", # Move legend to the right
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
#     panel.grid.minor = element_blank(),
#     aspect.ratio = 1
#   ) + # Ensure each facet is squared
#   scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
#   scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))



# compar_plot <- ggplot(data, aes(x = p, y = mean_distance, color = method_label)) +
#   geom_line(size = 0.7, linetype = "solid") +
#   geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
#   ylab("Average error") +
#   xlab("Number of variables") +
#   facet_grid(setting_label ~ n_label) + # Transposed facets
#   scale_colour_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
#   theme(
#     strip.text = element_text(size = 10), # Smaller font size for facet labels
#     axis.text = element_text(size = 12), # Smaller font size for axis text
#     axis.title = element_text(size = 12), # Smaller font size for axis titles
#     legend.text = element_text(size = 10), # Smaller font size for legend text
#     legend.title = element_blank(), # Hide legend title
#     legend.position = "right", # Move legend to the right
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
#     panel.grid.minor = element_blank(),
#     aspect.ratio = 1
#   ) + # Ensure each facet is squared
#   scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
#   scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
# plot(compar_plot)


# OLD dirty version with copy-pasted table from xl:

# # Define the data vectors
# data_vectors <- list(
#   # OUR Method (VAR heavy), n=500 and n=5000 respectively
#   c(0.0198, 0.1485, 0.732, 1.5346, 6.5754257) / standardisation,
#   c(0.000, 0.0693, 0.1287, 0.4554, 1.485148) / standardisation,
#   # OUR Method (GARCH heavy), n=500 and n=5000 respectively
#   c(0.079, 0.316, 1.336, 2.940, 13.019) / standardisation,
#   c(0.009, 0.069, 0.376, 1.128, 5.546) / standardisation,
#   # OUR Method (VAR non-heavy), n=500 and n=5000 respectively
#   c(0.69306, 1.93069, 4.5247, 7.5643, 23.831) / standardisation,
#   c(0.38, 1.31, 2.95, 4.32, 9.81) / standardisation,
#   # OUR Method (GARCH non-heavy), n=500 and n=5000 respectively
#   c(0.603, 1.831, 4.287, 7.079, 20.336) / standardisation,
#   c(0.366, 1.207, 2.742, 4.316, 8.266) / standardisation,

#   # PCMCI COR (VAR heavy), n=500 and n=5000 respectively
#   c(0.15, 1.61, 5.84, 12.26, 54.5) / standardisation,
#   c(0.18, 2.63, 8.7, 17.71, 66.75) / standardisation,
#   # PCMCI COR (GARCH heavy), n=500 and n=5000 respectively
#   c(0.74, 2.9, 9.28, 18.7, 73.60) / standardisation,
#   c(0.66, 3.21, 9.54, 18.62, 72) / standardisation,
#   # PCMCI COR (VAR non-heavy), n=500 and n=5000 respectively
#   c(0.12, 1.26, 5.42, 12.05, 52.87) / standardisation,
#   c(0.12, 1.25, 5.05, 11.86, 48) / standardisation,
#   # PCMCI COR (GARCH non-heavy), n=500 and n=5000 respectively
#   c(0.84, 3.52, 9.93, 19.89, 74.83) / standardisation,
#   c(0.8, 3.44, 9.92, 20.27, 76.2) / standardisation,

#   # PCMCI GPDC (VAR heavy), n=500 and n=5000 respectively
#   c(0.39, 2.91, 9.43, NA, NA) / (sequence_of_p_GPDC * (sequence_of_p_GPDC - 1)),
#   c(NA, NA, NA, NA, NA),
#   # PCMCI GPDC (GARCH heavy), n=500 and n=5000 respectively
#   c(0.57, 2.74, 7.87, NA, NA) / (sequence_of_p_GPDC * (sequence_of_p_GPDC - 1)),
#   c(NA, NA, NA, NA, NA),
#   # PCMCI GPDC (VAR non-heavy), n=500 and n=5000 respectively
#   c(0.11, 0.77, 3.09, NA, NA) / (sequence_of_p_GPDC * (sequence_of_p_GPDC - 1)),
#   c(NA, NA, NA, NA, NA),
#   # PCMCI GPDC (GARCH non-heavy), n=500 and n=5000 respectively
#   c(0.58, 2.74, 7.87, NA, NA) / (sequence_of_p_GPDC * (sequence_of_p_GPDC - 1)),
#   c(NA, NA, NA, NA, NA),
#   # Random data
#   rnorm(5, 0.5, 0.025), # 0.025 is approximatelly the variance of the random graph error.
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025),
#   rnorm(5, 0.5, 0.025)
# )

# # Fill the data frame
# data <- data.frame()
# for (i in seq_along(data_vectors)) {
#   # Determine method, setting_label, and n_label
#   method_index <- (i - 1) %/% 8 + 1
#   setting_index <- ((i - 1) %/% 2) %% 4 + 1
#   n_label_index <- ifelse(i %% 2 == 1, 1, 2)

#   # Create a temporary data frame for this vector
#   temp_df <- expand.grid(
#     p = sequence_of_p,
#     method = c("Our method", "PCMCI cor", "PCMCI gpdc", "Random")[method_index],
#     n_label = c("n = 500", "n = 5000")[n_label_index],
#     setting_label = c("VAR heavy-tailed", "GARCH heavy-tailed", "VAR Gaussian", "GARCH Gaussian")[setting_index]
#   )

#   # Assign the values from data_vectors to the temporary data frame
#   temp_df$structInterv_dist <- data_vectors[[i]]

#   # Append to the main data frame
#   data <- rbind(data, temp_df)
# }

# # Plotting
# ggplot(data, aes(x = p, y = structInterv_dist, color = method)) +
#   geom_line(size = 0.7, linetype = "solid") +
#   geom_point(size = 2, shape = 16) + # Use shape = 16 for solid circles
#   ylab("Average error") +
#   xlab("Number of variables") +
#   facet_grid(setting_label ~ n_label) + # Transposed facets
#   scale_colour_manual(values = c("black", "#EE6677", "#228833", "#1f77b4")) +
#   theme(
#     strip.text = element_text(size = 10), # Smaller font size for facet labels
#     axis.text = element_text(size = 12), # Smaller font size for axis text
#     axis.title = element_text(size = 12), # Smaller font size for axis titles
#     legend.text = element_text(size = 10), # Smaller font size for legend text
#     legend.title = element_blank(), # Hide legend title
#     legend.position = "right", # Move legend to the right
#     panel.grid.major = element_line(size = 0.5, linetype = "dashed", color = "grey80"),
#     panel.grid.minor = element_blank(),
#     aspect.ratio = 1
#   ) + # Ensure each facet is squared
#   scale_x_continuous(breaks = seq(0, 20, by = 5)) + # Regular breaks
#   scale_y_continuous(breaks = seq(0, 0.5, by = 0.1))
# # 7x10 in export pdf



