require(tidyverse)
require(data.table)
require(tidyr)
require(dplyr)
require(stringr)
require(reshape2)
require(lubridate)
require(ggplot2)
require(gridExtra)
require(rjags)
require(coda)
require(bayesplot)
require(posterior)
require(ggridges)


gc()
rm(list = ls())

nday <- 300
n_iter <- 3000

#-------------------------------------------------------------------------------

# Bring in model output
out_dat <- read.csv('/Users/theothompson/Library/CloudStorage/OneDrive-Personal/1_Masters/2_Salvins_nest_survival/R_nest_survival/Sensitivity_analysis/Outputs/Global_Uniform.csv')

n_iter <- 10000

# Find number of iterations per chain
out <- out_dat %>%
  group_by(.chain) %>%
  filter(.iteration > 50000) %>%
  ungroup()

#-------------------------------------------------------------------------------

# Create a date index for model_day
years <- 2019:2025

date_list <- map(years, ~ seq(as.Date(paste0(.x-1, "-08-01")), as.Date(paste0(.x, "-04-29")), by = "day"))

date_long <- map_dfr(seq_along(years), ~ tibble(
  year = as.character(years[.x]),
  date = date_list[[.x]],
  model_day = seq_along(date_list[[.x]])
))


#-------------------------------------------------------------------------------
out$breeding_suc = (out$phi.1.^73)*(out$phi.2.^out$global_change.2.)*(out$phi.3.^out$global_change.3.)
out$inc_suc = (out$phi.1.^73)
out$guard_suc = (out$phi.2.^out$global_change.2.)
out$pg_suc = (out$phi.3.^out$global_change.3.)
out$chick_suc = (out$phi.2.^out$global_change.2.)*(out$phi.3.^out$global_change.3.)
out$laying = (out$global_change.1.-73)
out$chick_dur = (out$global_change.2.+out$global_change.3.)

#-------------------------------------------------------------------------------


# Create a data frame with mean, lower CI, and upper CI for each parameter
posterior_dates_long <- tibble(
  parameter = c(
    "posterior_h_date",
    "posterior_pg_date",
    "posterior_f_date",
    "estimated_lay_date",
    "guard_duration",
    "pg_duration",
    "egg_detection",
    "guard_detection",
    "pg_detection",
    "dsr_inc",
    "dsr_guard",
    "dsr_pg",
    "breeding_suc",
    "inc_suc",
    "guard_suc",
    "pg_suc",
    "chick_suc",
    "chick_dur"
  ),
  mean = c(
    mean(out$global_change.1.),
    mean(out$global_change.1. + out$global_change.2.),
    mean(out$global_change.1. + out$global_change.2. + out$global_change.3.),
    mean(out$laying),
    mean(out$global_change.2.),
    mean(out$global_change.3.),
    mean(out$pi.1.),
    mean(out$pi.2.),
    mean(out$pi.3.),
    mean(out$phi.1.),
    mean(out$phi.2.),
    mean(out$phi.3.),
    mean(out$breeding_suc),
    mean(out$inc_suc),
    mean(out$guard_suc),
    mean(out$pg_suc),
    mean(out$chick_suc),
    mean(out$chick_dur)
  ),
  lower_ci = c(
    quantile(out$global_change.1., 0.025),
    quantile(out$global_change.1. + out$global_change.2., 0.025),
    quantile(out$global_change.1. + out$global_change.2. + out$global_change.3., 0.025),
    quantile(out$laying, 0.025),
    quantile(out$global_change.2., 0.025),
    quantile(out$global_change.3., 0.025),
    quantile(out$pi.1., 0.025),
    quantile(out$pi.2., 0.025),
    quantile(out$pi.3., 0.025),
    quantile(out$phi.1., 0.025),
    quantile(out$phi.2., 0.025),
    quantile(out$phi.3., 0.025),
    quantile(out$breeding_suc, 0.025),
    quantile(out$inc_suc, 0.025),
    quantile(out$guard_suc, 0.025),
    quantile(out$pg_suc, 0.025),
    quantile(out$chick_suc, 0.025),
    quantile(out$chick_dur, 0.025)
  ),
  upper_ci = c(
    quantile(out$global_change.1., 0.975),
    quantile(out$global_change.1. + out$global_change.2., 0.975),
    quantile(out$global_change.1. + out$global_change.2. + out$global_change.3., 0.975),
    quantile(out$laying, 0.975),
    quantile(out$global_change.2., 0.975),
    quantile(out$global_change.3., 0.975),
    quantile(out$pi.1., 0.975),
    quantile(out$pi.2., 0.975),
    quantile(out$pi.3., 0.975),
    quantile(out$phi.1., 0.975),
    quantile(out$phi.2., 0.975),
    quantile(out$phi.3., 0.975),
    quantile(out$breeding_suc, 0.975),
    quantile(out$inc_suc, 0.975),
    quantile(out$guard_suc, 0.975),
    quantile(out$pg_suc, 0.975),
    quantile(out$chick_suc, 0.975),
    quantile(out$chick_dur, 0.975)
  )
)

# Now round only numeric columns
posterior_dates_long <- posterior_dates_long %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))


posterior_dates <- posterior_dates_long[1:4,]
#-------------------------------------------------------------------------------

# Set up matrix
plot <- matrix(NA, nrow = nday, ncol = n_iter + 2)
colnames(plot) <- c("model_day","mean_surv", paste0("iter", 1:n_iter))
plot[,"model_day"] <- 1:nday


#--------------------------------- Function ------------------------------------
cumulative_survival <- function(dsr_egg, dsr_chick, dsr_pg, posterior_dates, nday = nday) {
  # Set up matrix
  temp <- c("model_day", "dsr", "prod")
  plot <- matrix(1, nrow = nday, ncol = 3, dimnames = list(NULL, temp))
  plot[, "model_day"] <- 1:nday
  
  # Extract the relevant dates from the posterior_dates table
  estimated_lay_date <- posterior_dates %>% filter(parameter == "estimated_lay_date") %>% pull(mean)
  posterior_h_date <- posterior_dates %>% filter(parameter == "posterior_h_date") %>% pull(mean)
  posterior_pg_date <- posterior_dates %>% filter(parameter == "posterior_pg_date") %>% pull(mean)
  posterior_f_date <- posterior_dates %>% filter(parameter == "posterior_f_date") %>% pull(mean)
  
  # Fill in the DSR column based on the date ranges
  for (i in 1:nday) {
    if (i >= estimated_lay_date && i < posterior_h_date) {
      plot[i, "dsr"] <- dsr_egg
    } else if (i >= posterior_h_date && i < posterior_pg_date) {
      plot[i, "dsr"] <- dsr_chick
    } else if (i >= posterior_pg_date && i <= posterior_f_date) {
      plot[i, "dsr"] <- dsr_pg
    }
  }
  
  # Compute cumulative survival
  for (i in 2:nday) {
    plot[i, "prod"] <- plot[i - 1, "prod"] * plot[i, "dsr"]
  }
  
  # Set values to NA before estimated_lay_date and after posterior_f_date
  for (i in 1:nday) {
    if (plot[i, "model_day"] < estimated_lay_date || plot[i, "model_day"] > posterior_f_date) {
      plot[i, "prod"] <- NA
    }
  }
  
  # Convert to dataframe for easier manipulation
  plot_df <- as.data.frame(plot)
  plot_df$prod <- as.numeric(plot_df$prod)
  
  return(data.frame(prod = plot_df$prod))
}
#-------------------------------------------------------------------------------


# Mean cumulative survival 
plot[,"mean_surv"] <- cumulative_survival(
  mean(out$phi.1.),  
  mean(out$phi.2.),
  mean(out$phi.3.),
  posterior_dates,
  nday = nday
)$prod

# Loop over iterations
set.seed(123)  # For reproducibility

sample_indices <- sample(nrow(out), n_iter, replace = FALSE)

for (i in 1:n_iter) {
  
  sampled_rows <- out[sample_indices[i], ]
  
  # Compute DSR for this sample
  dsr_egg <- sampled_rows$phi.1.
  dsr_chick <- sampled_rows$phi.2.
  dsr_pg <- sampled_rows$phi.3.
  
  # Compute cumulative survival
  plot[, paste0("iter", i)] <- cumulative_survival(dsr_egg, dsr_chick, dsr_pg, posterior_dates, nday = nday)$prod
}

# Convert to data frame for easier analysis
plot_df <- as.data.frame(plot)

#-------------------------------------------------------------------------------


# Pivot the iteration columns
plot_long <- plot_df %>%
  pivot_longer(
    cols = starts_with("iter"),  # Only pivot the iter columns
    names_to = "iteration",      # Create an iteration column
    values_to = "prod"           # Values go into prod column
  )

#------------------- Clean up for ggplot -------------------------

#Clean up plot_long
iters <- plot_long %>%
  select(-mean_surv) %>%  # Remove the 'mean_surv' column
  filter(!is.na(prod))    

#Clean up mean_surv
mean_surv <- plot_long %>%
  select(model_day, mean_surv) %>%
  filter(!is.na(mean_surv))

#Get CI bands
ci_bands <- plot_long %>%
  group_by(model_day) %>%
  summarise(
    lower_95 = quantile(prod, 0.025, na.rm = TRUE),
    upper_95 = quantile(prod, 0.975, na.rm = TRUE)
  )

#-------------------------------------------------------------------------------
# Create stage-specific table for plotting
surv_tbl <- posterior_dates_long %>%
  filter(parameter %in% c(
    "dsr_inc", "dsr_guard", "dsr_pg",
    "egg_detection", "guard_detection", "pg_detection"
  )) %>%
  mutate(
    stage = case_when(
      parameter == "dsr_inc"       ~ "Egg",
      parameter == "dsr_guard"     ~ "Guard",
      parameter == "dsr_pg"        ~ "Post-guard",
      parameter == "egg_detection" ~ "Egg",
      parameter == "guard_detection" ~ "Guard",
      parameter == "pg_detection"  ~ "Post-guard"
    ),
    
    metric = case_when(
      grepl("^dsr", parameter) ~ "dsr",
      grepl("detection", parameter) ~ "detect"
    )
  ) %>%
  select(stage, metric, mean, lower_ci, upper_ci) %>%
  tidyr::pivot_wider(
    names_from = metric,
    values_from = c(mean, lower_ci, upper_ci)
  ) %>%
  rename(
    dsr = mean_dsr,
    dsr_lcl = lower_ci_dsr,
    dsr_ucl = upper_ci_dsr,
    detect = mean_detect,
    detect_lcl = lower_ci_detect,
    detect_ucl = upper_ci_detect
  ) %>%
  mutate(
    stage = factor(stage, levels = c("Egg", "Guard", "Post-guard"))
  )

#-------------------------------------------------------------------------------

# Base theme
base_theme <- theme(
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 10),
  axis.title.y = element_text(margin = margin(r = 10)),
  axis.ticks.length = unit(0.1, "cm"),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  legend.position = "bottom",
  legend.background = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(color = "black", size = 1, fill = NA),
  panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.4),
  panel.grid.major.x = element_blank(),
  plot.title = element_text(hjust = 0, size = 12),
  plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
)

#-------------------------------------------------------------------------------

p_dsr <- ggplot(surv_tbl, aes(x = stage, y = dsr)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = dsr_lcl, ymax = dsr_ucl),
    width = 0,
    linewidth = 0.8
  ) +
  labs(
    title = "Daily Survival Rate (DSR)",
    x = "",
    y = "DSR"
  ) +
  scale_y_continuous(limits = c(0.980, 1)) +
  base_theme

p_dsr


p_detect <- ggplot(surv_tbl, aes(x = stage, y = detect)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = detect_lcl, ymax = detect_ucl),
    width = 0,
    linewidth = 0.8
  ) +
  labs(
    title = "Daily Detection Probability",
    x = "",
    y = "Detection probability"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  base_theme

p_detect


#-------------------------------------------------------------------------------
Cumulative_survival <- ggplot() +  
  # Iteration lines (grey, transparent but not in legend)
  geom_line(data = iters, aes(x = model_day, y = prod, group = iteration), 
            color = "lightgrey", alpha = 0.01, size = 1, show.legend = FALSE) +  
  
  # Add a single, visible grey line for "Model Iterations" in legend
  geom_line(data = data.frame(model_day = 400, prod = 0.5), 
            aes(x = model_day, y = prod, color = "Model Iterations"), 
            size = 2) +  
  
  # Mean survival line (black, included in legend)
  geom_line(data = mean_surv, aes(x = model_day, y = mean_surv, color = "Mean Survival"), 
            size = .8) +  
  
  # Vertical reference lines
  geom_vline(data = posterior_dates, aes(xintercept = mean), linetype = "dashed", size = 0.3) +  
  
  geom_line(data = ci_bands, aes(x = model_day, y = lower_95, color = "95% CI"),linetype = "dotted", size = 0.5) +
  geom_line(data = ci_bands, aes(x = model_day, y = upper_95, color = "95% CI"),linetype = "dotted", size = 0.5) +
  
  # Labels
  labs(
    x = NULL,
    y = "Survival Probability",
    title = "Cumulative Survival of Salvin's Albatross Eggs and Chicks",
    color = ""  
  ) +  
  
  # X-axis with bottom and top annotations
  scale_x_continuous(
    limits = c(1, 300),
    breaks = c(1, 50, 100, 150, 200, 250, 300),
    labels = c(
      expression(Aug~1^st), 
      expression(Sep~19^th), 
      expression(Nov~8^th), 
      expression(Dec~28^th), 
      expression(Feb~16^th), 
      expression(Apr~7^th), 
      expression(May~27^th)
    ),
    sec.axis = sec_axis(
      trans = ~., 
      breaks = posterior_dates$mean, 
      labels = c("Hatching       ","           Post-Guard","Fledging","Laying"))
  ) + 
  
  # Y-axis settings
  scale_y_continuous(
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.125)
  ) +  
  
  # Legend
  scale_color_manual(
    values = c("Mean Survival" = "black", "Model Iterations" = "grey", "95% CI" = "black"), 
    breaks = c("Mean Survival", "Model Iterations", "95% CI")
  ) +  
  base_theme

print(Cumulative_survival)


#-------------------------------------------------------------------------------

ggsave(filename = "/Users/theothompson/Library/CloudStorage/OneDrive-Personal/1_Masters/Journal_Submission/Cumulative_survival_6.5x4.5.jpg", 
       plot = Cumulative_survival,
       device = "jpeg",      # Saves as a JPEG file
       width = 6.5,            # Width in inches
       height = 4.5,           # Height in inches
       dpi = 300)            

#-------------------------------------------------------------------------------

DSR_Detect <- grid.arrange(
  p_dsr,
  p_detect,
  ncol = 2
)


ggsave(filename = "/Users/theothompson/Library/CloudStorage/OneDrive-Personal/1_Masters/Journal_Submission/DSR_Detect_6.5x3.jpg", 
       plot = DSR_Detect,
       device = "jpeg",      # Saves as a JPEG file
       width = 6.5,            # Width in inches
       height = 3,           # Height in inches
       dpi = 300) 

