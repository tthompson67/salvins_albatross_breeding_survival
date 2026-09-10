
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



gc()
rm(list = ls())


# Set to range required 1:193
nest_index <- 1:193

nind <- length(nest_index)
nday <- 300


# ------------------------------ Bring in data ---------------------------------

daily_obs <- read.csv('/Users/theothompson/Library/CloudStorage/OneDrive-Personal/1_Masters/2_Salvins_nest_survival/R_nest_survival/Data/all_data.csv')


#----------------- Lining up the different years of data -----------------------

# Set date formats
daily_obs <- daily_obs %>%
  mutate(
    date = as.Date(date, format = "%Y:%m:%d"),
    camera_stop = as.Date(camera_stop, format = "%Y-%m-%d"),
    camera_start = as.Date(camera_start, format = "%Y-%m-%d")
  )

# Create a date index for model_day
years <- 2019:2025

date_list <- map(years, ~ seq(as.Date(paste0(.x-1, "-08-01")), as.Date(paste0(.x, "-04-29")), by = "day"))

date_long <- map_dfr(seq_along(years), ~ tibble(
  year = as.character(years[.x]),
  date = date_list[[.x]],
  model_day = seq_along(date_list[[.x]])
))


# Add model_day columns to our data, relating to the dates of First_Check, Last_Check and Date
daily_obs <- daily_obs %>%
  left_join(date_long %>% select(date, obs_start = model_day), by = c("camera_start" = "date")) %>% 
  left_join(date_long %>% select(date, obs_stop = model_day), by = c("camera_stop" = "date")) %>% 
  left_join(date_long %>% select(date, model_day = model_day), by = c("date" = "date"))



#--------------------Re-shaping the data into a matrix--------------------------


# Find a list of UniqueID's
list_unique_id <- unique(daily_obs$unique_id)

# Create a blank df with the correct number or rows columns. 
m <- data.frame(
  nest_id = list_unique_id,
  obs_start = NA,
  obs_stop = NA, 
  matrix(NA, nrow = length(list_unique_id), ncol = nday)
)

# Set the row names
colnames(m)[4:ncol(m)] <- seq(4:ncol(m))


# Populate obs_start and obs_stop from daily_obs 
m <- m %>% 
  rowwise() %>% 
  mutate(
    obs_start = daily_obs %>%
      filter(unique_id == nest_id) %>%
      pull(obs_start) %>%
      first(),
    obs_stop = daily_obs %>%
      filter(unique_id == nest_id) %>%
      pull(obs_stop) %>%
      last()
  ) %>%
  ungroup()



# Loop over rows in list-uniqueID and populate the m matrix with the value from classification where the model_day matches the column header.
for (i in 1:nrow(m)) {
  
  # Subset daily_obs by unique_id
  sub <- daily_obs[daily_obs$unique_id == m$nest_id[i], ]
  
  # Check if the subset has rows
  if (nrow(sub) > 0) {
    
    # Iterate over each row in the subset dataframe
    for (j in 1:nrow(sub)) {
      
      # Get the ModelDay value for the current row
      model_day <- as.numeric(sub$model_day[j])
      
      # Get the column index corresponding to the model_day
      col_index <- model_day
      
      # Get the classification value for the current row
      classification_value <- sub$classification[j]
      
      # Fill in the matrix at the corresponding row and column for classification
      m[i, col_index + 3] <- classification_value  # Note the '+ 3' to account for the first two columns being obs_start and obs_stop
    }
  }
}

# Turn m into a matrix
m <- as.matrix(m)


#------------------------Correcting state observations--------------------------

# Look for "C" observations that occur after the first "PG" observation and correct them to "PG"
for (i in 1:nrow(m)) {
  # Find the first occurrence of "PG" in the row
  first_pg <- which(m[i, ] == "PG")[1]
  
  # If there is a "PG", check for "C" after that point
  if (!is.na(first_pg)) {
    # Replace "C" with "PG" after the first "PG"
    m[i, (first_pg+1):ncol(m)][m[i, (first_pg+1):ncol(m)] == "C"] <- "PG"
  }
}

# Look for "E" observations that occur after the first "C" observation and correct them to "C"
for (i in 1:nrow(m)) {
  # Find the first occurrence of "PG" in the row
  first_c <- which(m[i, ] == "C")[1]
  
  # If there is a "PG", check for "C" after that point
  if (!is.na(first_c)) {
    # Replace "C" with "PG" after the first "PG"
    m[i, (first_c+1):ncol(m)][m[i, (first_c+1):ncol(m)] == "E"] <- "C"
  }
}


#-------------------------Remove and empty rows---------------------------------


# Make a list of rows that contain only NA or NC (Not Checked) values, excluding the FirstCheck and LastCheck columns.
to_remove <- apply(m[, -c(1, 2, 3)], 1, function(x) all(is.na(x) | x == "NC"))


# Remove rows with only NA or "NC" values
m <- m[!to_remove, ]


#---------------------------------Year------------------------------------------

# Year vector for use in model
year <- as.numeric(factor(substr(m[nest_index, "nest_id"], 1, 4))) ##### Range1 and Range2 ######

year <- as.factor(year)

#Number of unique years
nyear <- length(unique(year))

#-------------------------------Obs matrix-------------------------------------- 

## A matrix of occasions where were we observed and alive egg/chick in the nest. 

obs <- m[nest_index,]

obs[obs == "E"] <- 1
obs[obs == "C"] <- 2
obs[obs == "PG"] <- 3

obs[obs %in% c("D", "NC", "F")] <- 0

obs[is.na(obs)] <- 0

# Convert to matrix, remove the first three columns and set values to numeric 
obs <- as.matrix(apply(obs[, !colnames(obs) %in% c("nest_id", "obs_start", "obs_stop")], 2, as.numeric))


#------------------------- Dead matrix -------------------------------------

## A dead recovery matrix.

dead <- m[nest_index,]

# Loop through each row
for (i in 1:nrow(dead)) {
  last_obs <- NA # Initialize the last observation
  
  for (j in 1:ncol(dead)) {
    if (!is.na(dead[i, j])) {
      if (dead[i, j] %in% c("E", "C", "PG")) {
        # Update the last observation
        last_obs <- switch(dead[i, j],
                           "E" = 1,
                           "C" = 2,
                           "PG" = 3)
      } else if (dead[i, j] == "D") {
        # Replace "D" with the corresponding number
        dead[i, j] <- last_obs
      }
    }
  }
}

dead[!dead %in% c("1", "2", "3")] <- "0"

dead[is.na(dead)] <- 0

# Convert to matrix, remove the first three columns and set values to numeric 
dead <- as.matrix(apply(dead[, !colnames(dead) %in% c("nest_id", "obs_start", "obs_stop")], 2, as.numeric))


#---------------------------- Effort Matrix ------------------------------------

# Copy m matrix to e
effort <- m[nest_index,]

# Loop through each row of the matrix to handle LastCheck and FirstCheck value observations. This represents when the cameras were active, outside of these parameters is is impossible for us to have checked. 
for (i in 1:nrow(effort)) {
  last_check <- as.numeric(effort[i, "obs_stop"]) + 3  # Get the LastCheck value as numeric and adjust for three columns
  
  # Set everything after last_check to "0" (excluding the first two columns)
  if (last_check < ncol(effort)) {
    effort[i, (last_check + 1):ncol(effort)] <- 0
  }
  
  # Use the value in the FirstCheck column for the first observation
  first_check <- as.numeric(effort[i, "obs_start"]) + 3  # Get the FirstCheck value as numeric and adjust for three columns
  
  # If a non-NA value is found, set everything before it to "0"
  if (!is.na(first_check) && first_check > 4) {
    effort[i, 4:(first_check - 1)] <- 0
  }
}

# Replace all occurrences of "Not Checked" with "0". Because we didn't look.
effort[effort == "NC"] <- 0

# Replace all observations with "1" because we had to have looked to see these classifications.
effort[effort %in% c("E", "C", "PG", "D", "F")] <- 1

# Replace all NA values with "1". Because the remaining NA values represent occasions where we did look, we just didn't see anything.
effort[is.na(effort)] <- 1

# Convert to matrix, remove the first three columns and set values to numeric 
effort <- as.matrix(apply(effort[, !colnames(effort) %in% c("nest_id", "obs_start", "obs_stop")], 2, as.numeric))

#effort.df <- as.data.frame(e.dat)
#write_csv(effort.df,"/Users/theodrewwardlethompson/Library/CloudStorage/OneDrive-Personal/1_Masters/R_Working/data/effort.csv")



#------------------------- Stage matrix ----------------------------------------

#------ Sampling

# Initialize matrices for transition periods and vectors for first and last capture
first_capture <- rep(NA, nind)
last_capture <- rep(NA, nind)
egg_chick <- matrix(NA, nind, 2)
chick_postguard <- matrix(NA, nind, 2)

# Set the earliest possible day for observations
start <- 1

# Loop over each individual to find transitions and last capture date
for (i in 1:nind) {
  # Combine observations and dead recoveries for each state
  egg <- which(obs[i, ] == 1 | dead[i, ] == 1)
  chick <- which(obs[i, ] == 2 | dead[i, ] == 2)
  post_guard <- which(obs[i, ] == 3 | dead[i, ] == 3)
  
  # Determine the end of the egg stage and start of the chick stage
  if (length(egg) > 0) {
    egg_chick[i, 1] <- max(egg)  # End of egg stage (last day observed/recovered in egg stage)
  } else {
    egg_chick[i, 1] <- start  # Default start if egg stage is missing
  }
  
  if (length(chick) > 0) {
    egg_chick[i, 2] <- min(chick)       # Start of chick stage
    chick_postguard[i, 1] <- max(chick) # End of chick stage
  } else {
    # If chick observations are missing, infer the transition based on post-guard data or end of egg stage
    egg_chick[i, 2] <- ifelse(length(post_guard) > 0, min(post_guard), nday)
    chick_postguard[i, 1] <- egg_chick[i, 1]  # Align with end of egg stage if chick data is missing
  }
  
  # Determine the start of the post-guard stage and the end of chick stage
  if (length(post_guard) > 0) {
    chick_postguard[i, 2] <- min(post_guard)  # Start of post-guard stage
  } else {
    chick_postguard[i, 2] <- nday  # Default end if post-guard is missing
  }
  
  # Find the last_capture
  last_capture[i] <- max(c(egg, chick, post_guard), na.rm = TRUE)
  # Find the first_capture
  first_capture[i]<- min(c(egg, chick, post_guard),na.rm = TRUE)
}


# Find state at the first capture (1=Egg, 2=Guard, 3=Post-Guard) 
first_stage <- obs[cbind(1:nrow(obs), first_capture)]


# Set all dead recoveries to 1 
dead = (dead>0)*1

#----------------------- Censoring matrix forthe model -------------------------
cens = matrix(1, nind, 3)

#--------------------------------- Model ---------------------------------------

model <- 'model {

    # Priors for global change points 
      global_change[1] ~ dunif(74, 134)    # Hatching day (mean = 104, variance = 1/0.04)
      global_change[2] ~ dunif(7, 47)      # Guard duration (mean = 27, variance = 1/0.062)
      global_change[3] ~ dunif(92, 152)    # Post-guard duration (mean = 122, variance = 1/0.062)

    # Precision parameters for individual-level variability 
      sd ~ dt(0, 0.01, 3)T(0,)
      tau <- 1 / sd^2

    # Priors for survival (phi), detection (pi), and dead recovery (lambda)
      for (h in 1:3) {
        phi[h] ~ dbeta(1, 1)          # Survival probability  
        pi[h] ~ dbeta(1, 1)           # Detection probability  
        lambda[h] ~ dbeta(1, 1)       # Dead recovery probability  
      }

    # Setting probabilities to 0 after fledging, as no further observations are expected.
      pi[4] <- 0                   # No detection after fledging
      lambda[4] <- 0               # No dead recovery after fledging
      phi[4] <- 0.7                # Survival probability after fledging - this shouldnt matter to the model.

    # Individual-specific variation around change points
      for (i in 1:nind) {
        change[i,1] ~ dnorm(global_change[1], tau)
        change[i,2] ~ dnorm(change[i,1] + global_change[2], tau)
        change[i,3] ~ dnorm(change[i,2] + global_change[3], tau)
      }

    # Interval censoring for observed transitions 
     for (i in 1:nind) {
        cens[i, 1] ~ dinterval(change[i, 1], egg_chick[i, 1:2])       # Censoring for hatching day
        cens[i, 2] ~ dinterval(change[i, 2], chick_postguard[i, 1:2]) # Censoring for guard duration
        cens[i, 3] ~ dinterval(change[i, 3], last_capture[i])         # Censoring for post-guard duration

    # Define latent state at first capture 
        z[i, first_capture[i]] <- 1                      # Individual is alive at first capture
        stage[i, first_capture[i]] <- first_stage[i]    # First observed stage

    # Modeling state transitions and observations
        for (j in (first_capture[i]+1):nday) {
            # Determine current stage based on change points.
            stage[i, j] <- dinterval(j, change[i, ]) + 1  # Stage based on change points

            # Latent state transition (live survival).
            z[i, j] ~ dbern(phi[stage[i, j - 1]] * z[i, j - 1])

            # Observation model for live detection.
            obs[i, j] ~ dbern(pi[stage[i, j]] * z[i, j] * effort[i, j])

            # Dead recovery model.
            dead[i, j] ~ dbern((z[i, j - 1] - z[i, j]) * lambda[stage[i, j - 1]])
        }
    }
}
'

#-----------------------------Initial values------------------------------------

init_func <- function(){ 
  min_hatch <- 50
  max_hatch <- 150
  min_pg <- 100
  max_pg <- 170
  min_flg <- 200
  max_flg <- 290  
  
  # Initialize matrices for transition ranges
  init_hatch_range <- matrix(NA, nind, 2)
  init_pg_range <- matrix(NA, nind, 2)
  init_fledge_range <- matrix(NA, nind, 2)
  
  # Hatch date
  init_hatch_range[, 1] <- pmax(egg_chick[, 1], min_hatch)
  init_hatch_range[, 2] <- pmin(egg_chick[, 2], chick_postguard[, 2], max_hatch)
  init_hatch_date <- runif(nind, init_hatch_range[, 1], init_hatch_range[, 2])
  
  # Post-guard date
  init_pg_range[, 1] <- pmax(chick_postguard[, 1], init_hatch_date, min_pg)
  init_pg_range[, 2] <- pmin(chick_postguard[, 2], max_pg)
  init_pg_date <- runif(nind, init_pg_range[, 1], init_pg_range[, 2])
  
  # Fledge date
  init_fledge_range[, 1] <- pmax(last_capture, init_pg_date, min_flg)
  init_fledge_range[, 2] <- max_flg
  init_fledge_date <- runif(nind, init_fledge_range[, 1], init_fledge_range[, 2])
  
  # Lengths of stages
  init_g_length <- init_pg_date - init_hatch_date
  init_pg_length <- init_fledge_date - init_pg_date
  
  # Cumulative change points
  change <- matrix(NA, nind, 3)
  change[, 1] <- init_hatch_date                            
  change[, 2] <- change[, 1] + init_g_length               
  change[, 3] <- change[, 2] + init_pg_length  
  
  # Latent state matrix (z)
  z <- matrix(0, nind, nday)
  for (i in 1:nind) {
    dead_flag <- FALSE
    
    for (j in 1:nday) {
      if (dead_flag) {
        z[i, j] <- 0
      } else if (j <= last_capture[i]) {
        z[i, j] <- 1
      } else {
        z[i, j] <- rbinom(1, 1, 0.95 * z[i, j - 1])
      }
      
      if (dead[i, j] == 1) {
        dead_flag <- TRUE
        z[i, j] <- 0
      }
    }
    
    z[i, 1:first_capture[i]] <- NA
  }
  
  # Initialize pi, phi, and lambda
  pi <- c(runif(3), NA)
  phi <- c(runif(3), NA)
  lambda <- c(runif(3), NA)
  
  return(list(
    change = change,
    pi = pi,
    phi = phi,
    lambda = lambda,
    z = z
  ))
}

#init_func()

#--------------------------------Compile----------------------------------------
m1 = jags.model(textConnection(model), 
                data = list(nind = nind, 
                            obs = (obs>0)*1, 
                            effort = effort,
                            first_capture = first_capture, 
                            nday = nday, 
                            egg_chick = egg_chick, 
                            chick_postguard = chick_postguard, 
                            last_capture = last_capture,
                            first_stage = first_stage,
                            dead = dead,
                            cens = cens),
                inits = init_func, 
                n.chains = 3, 
                n.adapt = 1000)


#-------------------------------Simulate----------------------------------------

# Runs 100,000 iterations as burn-in
# update(m1, n.iter = 100000)

# Runs 100,000 iterations and stores monitored parameters
out = coda.samples(m1, 
                   c('global_change', 
                     'pi', 
                     'phi', 
                     'lambda'),
                   n.iter = 100000)

#-------------------------------------------------------------------------------

# Summary 
summary(out)

# Store chains in draws_df
draws_df <- as_draws_df(out)

mcmc_trace(draws_df, regex_pars = c('phi'))
mcmc_trace(draws_df, regex_pars = c("global_change"))


# Save the draws_df for plotting! 
write_csv(draws_df, '/home/theo/R/out_model_v3.csv')
