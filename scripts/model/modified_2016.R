## Desc
# Refactored version run file

## Setup
#rm(list = ls())
options(mc.cores = 6)
n_chains <- 6
n_cores <- 6
n_sampling <- 500
n_warmup <- 500
n_refresh <- n_sampling*0.1

## Libraries
{
  library(tidyverse, quietly = TRUE)
  library(rstan, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(gridExtra, quietly = TRUE)
  library(pbapply, quietly = TRUE)
  library(parallel, quietly = TRUE)
  library(boot, quietly = TRUE)
  library(lqmm, quietly = TRUE) 
  library(gridExtra, quietly = TRUE)
  library(ggrepel, quietly = TRUE)
  library(dplyr)
}

## Functions
cov_matrix <- function(n, sigma2, rho){
    m <- matrix(nrow = n, ncol = n)
    m[upper.tri(m)] <- rho
    m[lower.tri(m)] <- rho
    diag(m) <- 1
    (sigma2^.5 * diag(n))  %*% m %*% (sigma2^.5 * diag(n))
}

mean_low_high <- function(draws, states, id){
  tmp <- draws
  draws_df <- data.frame(mean = inv.logit(apply(tmp, MARGIN = 2, mean)),
                                  high = inv.logit(apply(tmp, MARGIN = 2, mean) + 1.96 * apply(tmp, MARGIN = 2, sd)), 
                                  low  = inv.logit(apply(tmp, MARGIN = 2, mean) - 1.96 * apply(tmp, MARGIN = 2, sd)),
                               state = states, 
                                type = id)
  return(draws_df) 
}

check_cov_matrix <- function(mat,wt=state_weights){
  # get diagnoals
  s_diag <- sqrt(diag(mat))
  # output correlation matrix
  cor_equiv <- cov2cor(mat)
  diag(cor_equiv) <- NA
  # output equivalent national standard deviation
  nat_product <- sqrt(t(wt) %*% mat %*% wt) / 4
  
  # print & output
  hist(as_vector(cor_equiv),breaks = 10)
  
  hist(s_diag,breaks=10)
  
  
  print(sprintf('national sd of %s',round(nat_product,4)))
}

## Master variables
RUN_DATE <- ymd("2016-11-08")
#RUN_DATE <- ymd("2016-10-19")

election_day <- ymd("2016-11-08")
start_date <- as.Date("2016-03-01") # Keeping all polls after March 1, 2016


# wrangle polls -----------------------------------------------------------
all_polls <- read.csv("data/all_polls.csv", stringsAsFactors = FALSE, header = TRUE)

# select relevant columns from HufFPost polls
all_polls <- all_polls %>%
  dplyr::select(state, pollster, number.of.observations, population, mode, 
                start.date, 
                end.date,
                clinton, trump, undecided, other, johnson, mcmullin) %>%
  filter(ymd(end.date) <= RUN_DATE)

# basic mutations
df <- all_polls %>% 
  tbl_df %>%
  rename(n = number.of.observations) %>%
  mutate(begin = ymd(start.date),
         end   = ymd(end.date),
         t = end - (1 + as.numeric(end-begin)) %/% 2) %>%
  filter(t >= start_date & !is.na(t)
         & (population == "Likely Voters" | 
              population == "Registered Voters" | 
              population == "Adults") # get rid of disaggregated polls
         & n > 1) 

# pollster mutations
df <- df %>%
  mutate(pollster = str_extract(pollster, pattern = "[A-z0-9 ]+") %>% sub("\\s+$", "", .),
         pollster = replace(pollster, pollster == "Fox News", "FOX"), # Fixing inconsistencies in pollster names
         pollster = replace(pollster, pollster == "WashPost", "Washington Post"),
         pollster = replace(pollster, pollster == "ABC News", "ABC"),
         pollster = replace(pollster, pollster == "DHM Research", "DHM"),
         pollster = replace(pollster, pollster == "Public Opinion Strategies", "POS"),
         undecided = ifelse(is.na(undecided), 0, undecided),
         other = ifelse(is.na(other), 0, other) + 
           ifelse(is.na(johnson), 0, johnson) + 
           ifelse(is.na(mcmullin), 0, mcmullin))

# mode mutations
df <- df %>% 
  mutate(mode = case_when(mode == 'Internet' ~ 'Online poll',
                          grepl("live phone",tolower(mode)) ~ 'Live phone component',
                          TRUE ~ 'Other'))

# vote shares etc
df <- df %>%
  mutate(two_party_sum = clinton + trump,
         polltype = population,
         n_respondents = round(n),
         # clinton
         n_clinton = round(n * clinton/100),
         pct_clinton = clinton/two_party_sum,
         n_trump = round(n * trump/100),
         pct_trump = trump/two_party_sum)


## --- numerical indices
state_abb_list <- read.csv("data/potus_results_76_16.csv") %>%
  pull(state) %>% unique()

df <- df %>% 
  mutate(poll_day = t - min(t) + 1,
         # Factors are alphabetically sorted: 1 = --, 2 = AL, 3 = AK, 4 = AZ...
         index_s = as.numeric(factor(as.character(state),
                                     levels = c('--',state_abb_list))),
         index_s = ifelse(index_s == 1, 52, index_s - 1),
         index_t = 1 + as.numeric(t) - min(as.numeric(t)),
         index_p = as.numeric(as.factor(as.character(pollster))),
         index_m = as.numeric(as.factor(as.character(mode))),
         index_pop = as.numeric(as.factor(as.character(polltype)))) %>%
  # selections
  arrange(state, t, polltype, two_party_sum) %>% 
  distinct(state, t, pollster, .keep_all = TRUE) %>%
  select(
    # poll information
    state, t, begin, end, pollster, polltype, method = mode, n_respondents, 
    # vote shares
    pct_clinton, n_clinton, 
    pct_trump, n_trump, 
    poll_day, index_s, index_p, index_m, index_pop, index_t)

# useful vectors
all_polled_states <- df$state %>% unique %>% sort

# day indices
first_day <- min(df$begin)
ndays <- max(df$t) - min(df$t)
all_t <- min(df$t) + days(0:(ndays))
all_t_until_election <- min(all_t) + days(0:(election_day - min(all_t)))
# pollster indices
all_pollsters <- levels(as.factor(as.character(df$pollster)))


# getting state contextual information from 2012 -------------------------
states2012 <- read.csv("data/2012.csv", 
                       header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(score = obama_count / (obama_count + romney_count),
         national_score = sum(obama_count)/sum(obama_count + romney_count),
         delta = score - national_score,
         share_national_vote = (total_count*(1+adult_pop_growth_2011_15))
         /sum(total_count*(1+adult_pop_growth_2011_15))) %>%
  arrange(state) 

state_abb <- states2012$state
rownames(states2012) <- state_abb

# get state indices
all_states <- states2012$state
state_name <- states2012$state_name
names(state_name) <- state_abb

# set prior differences
prior_diff_score <- states2012$delta
names(prior_diff_score) <- state_abb

# set state weights
state_weights <- c(states2012$share_national_vote / sum(states2012$share_national_vote))
names(state_weights) <- state_abb

# electoral votes, by state:
ev_state <- states2012$ev
names(ev_state) <- state_abb


# create covariance matrices ----------------------------------------------
# start by reading in data
state_data <- read.csv("data/potus_results_76_16.csv")
state_data <- state_data %>% 
  select(year, state, dem) %>%
  group_by(state) %>%
  mutate(dem = dem ) %>% #mutate(dem = dem - lag(dem)) %>%
  select(state,variable=year,value=dem)  %>%
  ungroup() %>%
  na.omit() %>%
  filter(variable == 2016)

census <- read.csv('data/acs_2013_variables.csv')
census <- census %>%
  filter(!is.na(state)) %>% 
  select(-c(state_fips,pop_total,pop_density)) %>%
  group_by(state) %>%
  gather(variable,value,
         1:(ncol(.)-1))

state_data <- state_data %>%
  mutate(variable = as.character(variable)) %>%
  bind_rows(census)

# add urbanicity
urbanicity <- read.csv('data/urbanicity_index.csv') %>%
  dplyr::select(state,pop_density = average_log_pop_within_5_miles) %>%
  gather(variable,value,
         2:(ncol(.)))

state_data <- state_data %>%
  bind_rows(urbanicity)

# add pct white evangelical
white_evangel_pct <- read_csv('data/white_evangel_pct.csv') %>%
  gather(variable,value,
         2:(ncol(.)))

state_data <- state_data %>%
  bind_rows(white_evangel_pct)

# add region, as a dummy for each region
regions <- read_csv('data/state_region_crosswalk.csv') %>%
  select(state = state_abb, variable=region) %>%
  mutate(value = 1) %>%
  spread(variable,value)

regions[is.na(regions)] <- 0

regions <- regions %>%
  gather(variable,value,2:ncol(.))

#state_data <- state_data %>%
#  bind_rows(regions)

# scale and spread
state_data_long <- state_data %>%
  group_by(variable) %>%
  # scale all varaibles
  mutate(value = (value - min(value, na.rm=T)) / 
           (max(value, na.rm=T) - min(value, na.rm=T))) %>%
  #mutate(value = (value - mean(value)) / sd(value)) %>%
  # now spread
  spread(state, value) %>% 
  na.omit() %>%
  ungroup() %>%
  select(-variable)

# compute the correlation matrix
# formula is
# a*(lambda*C + (1-lambda)*C_1)
# where C is our correlation matrix with min 0
# and C_1 is a sq matrix with all 1's
# lambda=0 is 100% correlation, lambda=1 is our corr matrix

# save correlation 
C <- cor(state_data_long)  

# increase the baseline correlation of the matrix to correspond to national-level error
C[C < 0] <- 0 # baseline cor for national poll error


tmp_C <- C
diag(tmp_C) <- NA
mean(tmp_C,na.rm=T)

# mixing with matrix of 0.5s
lambda <- 0.75
C_1 <- matrix(data=1,nrow = 51,ncol=51)
a <- 1
new_C <- (lambda*C + (1-lambda)*C_1) %>% make.positive.definite()

tmp <- new_C
diag(tmp) <- NA
mean(tmp,na.rm=T)

state_correlation_polling <- new_C

# make pos definite
state_correlation_polling <- make.positive.definite(state_correlation_polling)

# covariance matrix for polling error
state_covariance_polling_bias <- cov_matrix(51, 0.078^2, 0.9) # 3.4% on elec day
state_covariance_polling_bias <- state_covariance_polling_bias * state_correlation_polling

(sqrt(t(state_weights) %*% state_covariance_polling_bias %*% state_weights) / 4)
mean(apply(MASS::mvrnorm(100,rep(0,51),state_covariance_polling_bias),2,sd) /4) 

# covariance for prior e-day prediction
state_covariance_mu_b_T <- cov_matrix(n = 51, sigma2 = 0.18^2, rho = 0.9) # 6% on elec day
state_covariance_mu_b_T <- state_covariance_mu_b_T * state_correlation_polling

(sqrt(t(state_weights) %*% state_covariance_mu_b_T %*% state_weights) / 4)
mean(apply(MASS::mvrnorm(100,rep(0,51),state_covariance_mu_b_T),2,sd) /4) 

# covariance matrix for random walks
state_covariance_mu_b_walk <- cov_matrix(51, (0.017)^2, 0.9)
state_covariance_mu_b_walk <- state_covariance_mu_b_walk * state_correlation_polling # we want the demo correlations for filling in gaps in the polls

(sqrt(t(state_weights) %*% state_covariance_mu_b_walk %*% state_weights) / 4) * sqrt(300)
mean(apply(MASS::mvrnorm(100,rep(0,51),state_covariance_mu_b_walk),2,sd) /4) * sqrt(300)

## MAKE DEFAULT COV MATRICES
# we're going to make TWO covariance matrix here and pass it
# and 3 scaling values to stan, where the 3 values are 
# (1) the national sd on the polls, (2) the national sd
# on the prior and (3) the national sd of the random walk
# make initial covariance matrix (using specified correlation)
state_covariance_0 <- cov_matrix(51, 0.07^2, 0.9)
state_covariance_0 <- state_covariance_0 * state_correlation_polling # we want the demo correlations for filling in gaps in the polls

# national error of:
sqrt(t(state_weights) %*% state_covariance_0 %*% state_weights) / 4

# save the inital scaling factor
national_cov_matrix_error_sd <- sqrt(t(state_weights) %*% state_covariance_0 %*% state_weights) %>% as.numeric()
national_cov_matrix_error_sd

# save the other scales for later
fit_rmse_day_x <- function(x){0.03 +  (10^-6.6)*(x)^2} # fit to error from external script
fit_rmse_day_x(0:300)
days_til_election <- as.numeric(difftime(election_day,RUN_DATE))
expected_national_mu_b_T_error <- fit_rmse_day_x(days_til_election)

polling_bias_scale <- 0.013 # on the probability scale -- we convert later down
mu_b_T_scale <- expected_national_mu_b_T_error # on the probability scale -- we convert later down
random_walk_scale <- 0.05/sqrt(300) # on the probability scale -- we convert later down

# gen fake matrices, check the math (this is recreated in stan)
national_cov_matrix_error_sd <- sqrt(t(state_weights) %*% state_covariance_0 %*% state_weights) %>% as.numeric()

ss_cov_poll_bias = state_covariance_0 * (polling_bias_scale/national_cov_matrix_error_sd*4)^2
ss_cov_mu_b_T = state_covariance_0 * (mu_b_T_scale/national_cov_matrix_error_sd*4)^2
ss_cov_mu_b_walk = state_covariance_0 * (random_walk_scale/national_cov_matrix_error_sd*4)^2

sqrt(t(state_weights) %*% ss_cov_poll_bias %*% state_weights) / 4 
sqrt(t(state_weights) %*% ss_cov_mu_b_T %*% state_weights) / 4 
sqrt(t(state_weights) %*% ss_cov_mu_b_walk %*% state_weights) / 4 * sqrt(300)



# checking parameters -----------------------------------------------------
par(mfrow=c(3,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
check_cov_matrix(state_covariance_polling_bias)
check_cov_matrix(state_covariance_mu_b_T)
check_cov_matrix(state_covariance_mu_b_walk)

# poll bias should be:
err <- c(0.5, 1.9, 0.8, 7.2, 1.0, 1.4, 0.1, 3.3, 3.4, 0.9, 0.3, 2.7, 1.0) / 2 
sqrt(mean(err^2))

# implied national posterior on e-day
1 / sqrt( 1/((sqrt(t(state_weights) %*% state_covariance_polling_bias %*% state_weights) / 4)^2) + 
            1/((sqrt(t(state_weights) %*% state_covariance_mu_b_T %*% state_weights) / 4)^2) )

# random walks
replicate(1000,cumsum(rnorm(300)*0.5)) %>% 
  as.data.frame() %>%
  mutate(innovation = row_number()) %>%
  gather(walk_number,cumsum,1:100) %>%
  ggplot(.,aes(x=innovation,y=cumsum,group=walk_number)) +
  geom_line() +
  labs(subtitle='normal(0,0.5) random walk')

replicate(1000,cumsum(rt(300,7)*0.5)) %>% 
  as.data.frame() %>%
  mutate(innovation = row_number()) %>%
  gather(walk_number,cumsum,1:100) %>%
  ggplot(.,aes(x=innovation,y=cumsum,group=walk_number)) +
  geom_line() +
  labs(subtitle='student_t(7,0,0.5) random walk')

# create priors -----------------------------------------------------------
# read in abramowitz data
abramowitz <- read.csv('data/abramowitz_data.csv') %>% filter(year < 2016)
prior_model <- lm(incvote ~  juneapp + q2gdp, data = abramowitz)

# make predictions
national_mu_prior <- predict(prior_model,newdata = tibble(q2gdp = 1.1, juneapp = 4))
# on correct scale
national_mu_prior <- national_mu_prior / 100
# Mean of the mu_b_prior
mu_b_prior <- logit(national_mu_prior + prior_diff_score)
# or read in priors if generated already
prior_in <- read_csv("data/state_priors_08_12_16.csv") %>%
  filter(date <= RUN_DATE) %>%
  group_by(state) %>%
  arrange(date) %>%
  filter(date == max(date)) %>%
  select(state,pred) %>%
  ungroup() %>%
  arrange(state)

mu_b_prior <- logit(prior_in$pred + 0.0)
names(mu_b_prior) <- prior_in$state
names(mu_b_prior) == names(prior_diff_score) # correct order?
national_mu_prior <- weighted.mean(inv.logit(mu_b_prior), state_weights)
cat(sprintf('Prior Clinton two-party vote is %s\nWith a national sd of %s\n', 
            round(national_mu_prior,3),round(mu_b_T_scale,3)))

# alpha for disconnect in state v national polls 
score_among_polled <- sum(states2012[all_polled_states[-1],]$obama_count)/
  sum(states2012[all_polled_states[-1],]$obama_count + 
        states2012[all_polled_states[-1],]$romney_count)
alpha_prior <- log(states2012$national_score[1]/score_among_polled)

# firms that adjust for party (after 2016)
adjusters <- c(
  "ABC",
  "Washington Post",
  "Ipsos",
  "Pew",
  "YouGov",
  "NBC"
)

df %>% filter((pollster %in% adjusters)) %>% pull(pollster) %>% unique()


# passing data to Stan ----------------------------------------------------
N_state_polls <- nrow(df %>% filter(index_s != 52))
N_national_polls <- nrow(df %>% filter(index_s == 52))
T <- as.integer(round(difftime(election_day, first_day)))
current_T <- max(df$poll_day)
S <- 51
P <- length(unique(df$pollster))
M <- length(unique(df$method))
Pop <- length(unique(df$polltype))

state <- df %>% filter(index_s != 52) %>% pull(index_s)
day_national <- df %>% filter(index_s == 52) %>% pull(poll_day) 
day_state <- df %>% filter(index_s != 52) %>% pull(poll_day) 
poll_national <- df %>% filter(index_s == 52) %>% pull(index_p) 
poll_state <- df %>% filter(index_s != 52) %>% pull(index_p) 
poll_mode_national <- df %>% filter(index_s == 52) %>% pull(index_m) 
poll_mode_state <- df %>% filter(index_s != 52) %>% pull(index_m) 
poll_pop_national <- df %>% filter(index_s == 52) %>% pull(index_pop) 
poll_pop_state <- df %>% filter(index_s != 52) %>% pull(index_pop) 

n_democrat_national <- df %>% filter(index_s == 52) %>% pull(n_clinton)
n_democrat_state <- df %>% filter(index_s != 52) %>% pull(n_clinton)
n_two_share_national <- df %>% filter(index_s == 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
n_two_share_state <- df %>% filter(index_s != 52) %>% transmute(n_two_share = n_trump + n_clinton) %>% pull(n_two_share)
unadjusted_national <- df %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s == 52) %>% pull(unadjusted)
unadjusted_state <- df %>% mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% filter(index_s != 52) %>% pull(unadjusted)
                        
# priors (on the logit scale)
sigma_measure_noise_national <- 0.04
sigma_measure_noise_state <- 0.04
sigma_c <- 0.06
sigma_m <- 0.04
sigma_pop <- 0.04
sigma_e_bias <- 0.02

polling_bias_scale <- as.numeric(polling_bias_scale) * 4
mu_b_T_scale <- as.numeric(mu_b_T_scale) * 4
random_walk_scale <- as.numeric(random_walk_scale) * 4

# put the data in a list to export to Stan
data <- list(
  N_national_polls = N_national_polls,
  N_state_polls = N_state_polls,
  T = T,
  S = S,
  P = P,
  M = M,
  Pop = Pop,
  state = state,
  state_weights = state_weights,
  day_state = as.integer(day_state),
  day_national = as.integer(day_national),
  poll_state = poll_state,
  poll_national = poll_national,
  poll_mode_national = poll_mode_national, 
  poll_mode_state = poll_mode_state,
  poll_pop_national = poll_pop_national, 
  poll_pop_state = poll_pop_state,
  unadjusted_national = unadjusted_national,
  unadjusted_state = unadjusted_state,
  n_democrat_national = n_democrat_national,
  n_democrat_state = n_democrat_state,
  n_two_share_national = n_two_share_national,
  n_two_share_state = n_two_share_state,
  sigma_measure_noise_national = sigma_measure_noise_national,
  sigma_measure_noise_state = sigma_measure_noise_state,
  mu_b_prior = mu_b_prior,
  sigma_c = sigma_c,
  sigma_m = sigma_m,
  sigma_pop = sigma_pop,
  sigma_e_bias = sigma_e_bias,
  # covariance matrices
  # ss_cov_mu_b_walk = state_covariance_mu_b_walk,
  # ss_cov_mu_b_T = state_covariance_mu_b_T,
  # ss_cov_poll_bias = state_covariance_polling_bias
  state_covariance_0 = state_covariance_0,
  polling_bias_scale = polling_bias_scale,
  mu_b_T_scale = mu_b_T_scale,
  random_walk_scale = random_walk_scale
)


# run the Stan model ------------------------------------------------------
message("Running model...")
# model options
#("scripts/model/poll_model_2020_no_partisan_correction.stan")
#("scripts/model/poll_model_2020_no_mode_adjustment.stan")
#("scripts/model/poll_model_2020.stan")

# if rstan, uncomment these lines:
model <- rstan::stan_model("scripts/model/poll_model_2020.stan")
out <- rstan::sampling(model, data = data,
                       refresh = n_refresh,
                       chains  = n_chains, iter = 500, warmup = 250
)

# else if cmdstan, uncomment these
# model <- cmdstanr::cmdstan_model("scripts/model/poll_model_2020.stan",compile=TRUE,force=TRUE)
# fit <- model$sample(
#   data = data,
#   seed = 1843,
#   parallel_chains = n_cores,
#   chains = n_chains,
#   iter_warmup = n_warmup,
#   iter_sampling = n_sampling,
#   refresh = n_refresh
# )

# out <- rstan::read_stan_csv(fit$output_files())
# rm(fit)
# gc()

# save model for today
write_rds(out, sprintf('models/stan_model_%s.rds',RUN_DATE),compress = 'gz')

# extracting results ----
out <- read_rds(sprintf('models/stan_model_%s.rds',RUN_DATE))
# write_csv(out, file = "OUTCHECK.csv")
posterior_samples <- extract(out, permuted = TRUE)  # as list
posterior_samples_df <- as.data.frame(posterior_samples)  # as data frame
write_csv(posterior_samples_df, file = "OUTCHECK.csv")

## --- priors
## mu_b_T
y <- MASS::mvrnorm(1000, mu_b_prior, Sigma = state_covariance_mu_b_T)
mu_b_T_posterior_draw <- rstan::extract(out, pars = "mu_b")[[1]][,,254]
mu_b_T_prior_draws     <- mean_low_high(y, states = colnames(y), id = "prior")
mu_b_T_posterior_draws <- mean_low_high(mu_b_T_posterior_draw, states = colnames(y), id = "posterior")
mu_b_T <- rbind(mu_b_T_prior_draws, mu_b_T_posterior_draws)
mu_b_t_plt <- mu_b_T %>% arrange(mean) %>%
  ggplot(.) +
    geom_point(aes(y = mean, x = reorder(state, mean), color = type), position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = low, ymax = high, x = state, color = type), width = 0, position = position_dodge(width = 0.5)) +
    coord_flip() +
    theme_bw()
mu_b_t_plt
## mu_c
mu_c_posterior_draws <- rstan::extract(out, pars = "mu_c")[[1]] 
mu_c_posterior_draws <- data.frame(draws = as.vector(mu_c_posterior_draws),
                                   index_p = sort(rep(seq(1, P), dim(mu_c_posterior_draws)[1])), 
                                   type = "posterior")
mu_c_prior_draws <- data.frame(draws = rnorm(P * 1000, 0, sigma_c),
                               index_p = sort(rep(seq(1, P), 1000)), 
                               type = "prior")
mu_c_draws <- rbind(mu_c_posterior_draws, mu_c_prior_draws) 
pollster <- df %>% select(pollster, index_p) %>% distinct()
mu_c_draws <- merge(mu_c_draws, pollster, by = "index_p", all.x = TRUE)
mu_c_draws <- mu_c_draws %>%
  group_by(pollster, type) %>%
  summarize(mean = mean(draws), 
            low = mean(draws) - 1.96 * sd(draws),
            high = mean(draws) + 1.96 * sd(draws))
mu_c_plt <- mu_c_draws %>% 
  arrange(mean) %>% 
  filter(pollster %in% (df %>% group_by(pollster) %>% 
                          summarise(n=n()) %>% filter(n>=5) %>% pull(pollster))) %>%
  ggplot(.) +
    geom_point(aes(y = mean, x = reorder(pollster, mean), color = type), 
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = low, ymax = high, x = pollster, color = type), 
                  width = 0, position = position_dodge(width = 0.5)) +
    coord_flip() +
    theme_bw()
mu_c_plt
# write_csv(mu_c_draws,'output/mu_c_draws_2016.csv')
## mu_m
mu_m_posterior_draws <- rstan::extract(out, pars = "mu_m")[[1]] 
mu_m_posterior_draws <- data.frame(draws = as.vector(mu_m_posterior_draws),
                                   index_m = sort(rep(seq(1, M), dim(mu_m_posterior_draws)[1])), 
                                   type = "posterior")
mu_m_prior_draws <- data.frame(draws = rnorm(M * 1000, 0, sigma_m),
                               index_m = sort(rep(seq(1, M), 1000)), 
                               type = "prior")
mu_m_draws <- rbind(mu_m_posterior_draws, mu_m_prior_draws) 
method <- df %>% select(method, index_m) %>% distinct()
mu_m_draws <- merge(mu_m_draws, method, by = "index_m", all.x = TRUE)
mu_m_draws <- mu_m_draws %>%
  group_by(method, type) %>%
  summarize(mean = mean(draws), 
            low = mean(draws) - 1.96 * sd(draws),
            high = mean(draws) + 1.96 * sd(draws))
mu_m_plt <- mu_m_draws %>% arrange(mean) %>%
  ggplot(.) +
  geom_point(aes(y = mean, x = reorder(method, mean), color = type), 
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = low, ymax = high, x = method, color = type), 
                width = 0, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw()
mu_m_plt
## mu_pop
mu_pop_posterior_draws <- rstan::extract(out, pars = "mu_pop")[[1]] 
mu_pop_posterior_draws <- data.frame(draws = as.vector(mu_pop_posterior_draws),
                                   index_pop = sort(rep(seq(1, M), dim(mu_pop_posterior_draws)[1])), 
                                   type = "posterior")
mu_pop_prior_draws <- data.frame(draws = rnorm(Pop * 1000, 0, sigma_pop),
                               index_pop = sort(rep(seq(1, Pop), 1000)), 
                               type = "prior")
mu_pop_draws <- rbind(mu_pop_posterior_draws, mu_pop_prior_draws) 
method <- df %>% select(polltype, index_pop) %>% distinct()
mu_pop_draws <- merge(mu_pop_draws, method, by = "index_pop", all.x = TRUE)
mu_pop_draws <- mu_pop_draws %>%
  group_by(polltype, type) %>%
  summarize(mean = mean(draws), 
            low = mean(draws) - 1.96 * sd(draws),
            high = mean(draws) + 1.96 * sd(draws))
mu_pop_plt <- mu_pop_draws %>% arrange(mean) %>%
  ggplot(.) +
  geom_point(aes(y = mean, x = reorder(polltype, mean), color = type), 
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = low, ymax = high, x = polltype, color = type), 
                width = 0, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw()
mu_pop_plt
## state error terms
polling_bias_posterior <- rstan::extract(out, pars = "polling_bias")[[1]]
polling_bias_posterior %>% apply(.,2,sd) / 4
polling_bias_posterior_draws <- data.frame(draws = as.vector(polling_bias_posterior),
                                   index_s = sort(rep(seq(1, S), dim(polling_bias_posterior)[1])), 
                                   type = "posterior")
y <- MASS::mvrnorm(1000, rep(0, S), Sigma = state_covariance_polling_bias)
polling_bias_prior_draws <- data.frame(draws = as.vector(y),
                                   index_s = sort(rep(seq(1, S), dim(y)[1])), 
                                    type = "prior")
polling_bias_draws <- rbind(polling_bias_posterior_draws, polling_bias_prior_draws) 
states <- data.frame(index_s = 1:51, states = rownames(state_correlation_polling))
polling_bias_draws <- merge(polling_bias_draws, states, by = "index_s", all.x = TRUE)
polling_bias_draws <- polling_bias_draws %>%
  group_by(states, type) %>%
  summarize(mean = mean(draws), 
            low = mean(draws) - 1.96 * sd(draws),
            high = mean(draws) + 1.96 * sd(draws))
polling_bias_plt <- polling_bias_draws %>%
  ggplot(.) +
    geom_point(aes(y = mean, x = reorder(states, mean), color = type), 
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = low, ymax = high, x = states, color = type), 
                  width = 0, position = position_dodge(width = 0.5)) +
    coord_flip() +
    theme_bw() 
polling_bias_plt
## Posterior
# poll terms
poll_terms <- rstan::extract(out, pars = "mu_c")[[1]]
non_adjusters <- df %>% 
  mutate(unadjusted = ifelse(!(pollster %in% adjusters), 1, 0)) %>% 
  select(unadjusted, index_p) %>%
  distinct() %>%
  arrange(index_p)

e_bias <- rstan::extract(out, pars = "e_bias")[[1]]
plt_adjusted <- lapply(1:100,
       function(x){
         tibble(e_bias_draw = e_bias[x,]
                - mean(poll_terms[x, non_adjusters[non_adjusters$unadjusted == 0, 2]$index_p])
                + mean(poll_terms[x, non_adjusters[non_adjusters$unadjusted == 1, 2]$index_p]),
                trial = x) %>%
           mutate(date = min(df$end) + row_number())
       }) %>%
  do.call('bind_rows',.) %>%
  ggplot(.,aes(x=date,y=e_bias_draw,group=trial)) +
  geom_line(alpha=0.2)
plt_unadjusted <- lapply(1:100,
       function(x){
         tibble(e_bias_draw = e_bias[x,],
                trial = x) %>%
           mutate(date = min(df$end) + row_number())
       }) %>%
  do.call('bind_rows',.) %>%
  ggplot(.,aes(x=date,y=e_bias_draw,group=trial)) +
  geom_line(alpha=0.2)
grid.arrange(plt_adjusted, plt_unadjusted)

rstan::extract(out, pars = "e_bias")[[1]] %>% apply(.,2,mean) %>% plot

# states
predicted_score <- rstan::extract(out, pars = "predicted_score")[[1]]


# state correlation?
single_draw <- as.data.frame(predicted_score[,dim(predicted_score)[2],])
names(single_draw) <- colnames(state_correlation_polling)
single_draw %>% 
  select(AL,CA,FL,MN,NC,NM,RI,WI) %>%  #NV,FL,WI,MI,NH,OH,IA,NC,IN
  cor 


pct_clinton <- pblapply(1:dim(predicted_score)[3],
                    function(x){
                      # pred is mu_a + mu_b for the past, just mu_b for the future
                      temp <- predicted_score[,,x]
                      
                      # put in tibble
                      tibble(low = apply(temp,2,function(x){(quantile(x,0.025))}),
                             high = apply(temp,2,function(x){(quantile(x,0.975))}),
                             mean = apply(temp,2,function(x){(mean(x))}),
                             prob = apply(temp,2,function(x){(mean(x>0.5))}),
                             state = x) 
                      
                    }) %>% do.call('bind_rows',.)


pct_clinton$state = colnames(state_correlation_polling)[pct_clinton$state]

pct_clinton <- pct_clinton %>%
  group_by(state) %>%
  mutate(t = row_number() + min(df$begin)) %>%
  ungroup()

# national
pct_clinton_natl <- pblapply(1:dim(predicted_score)[1],
                         function(x){
                           # each row is a day for a particular draw
                           temp <- predicted_score[x,,] %>% as.data.frame()
                           names(temp) <- colnames(state_correlation_polling)
                           
                           # for each row, get weigted natl vote
                           tibble(natl_vote = apply(temp,MARGIN = 1,function(y){weighted.mean(y,state_weights)})) %>%
                             mutate(t = row_number() + min(df$begin)) %>%
                             mutate(draw = x)
                         }) %>% do.call('bind_rows',.)

pct_clinton_natl <- pct_clinton_natl %>%
  group_by(t) %>%
  summarise(low = quantile(natl_vote,0.025),
            high = quantile(natl_vote,0.975),
            mean = mean(natl_vote),
            prob = mean(natl_vote > 0.5)) %>%
  mutate(state = '--')

# bind state and national vote
pct_clinton <- pct_clinton %>%
  bind_rows(pct_clinton_natl) %>%
  arrange(desc(mean))


options(max.print = .Machine$integer.max)


election_results <- c(AL = 0, AK = 0, AZ = 0, AR = 0, CA = 1, CO = 1, CT = 1,
                      DC = 1, DE = 1, FL = 0, GA = 0, HI = 1, ID = 0, IL = 1,
                      IN = 0, IA = 0, KS = 0, KY = 0, LA = 0, ME = 1, MD = 1,
                      MA = 1, MI = 0, MN = 1, MS = 0, MO = 0, MT = 0, NE = 0,
                      NV = 1, NH = 1, NJ = 1, NM = 1, NY = 1, NC = 0, ND = 0,
                      OH = 0, OK = 0, OR = 1, PA = 0, RI = 1, SC = 0, SD = 0,
                      TN = 0, TX = 0, UT = 0, VT = 1, VA = 1, WA = 1, WV = 0,
                      WI = 0, WY = 0)


state_abbreviations <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", "GA", "HI",
                         "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN",
                         "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH",
                         "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA",
                         "WI", "WV", "WY")

specified_states <- c('AZ', 'CO', 'FL', 'GA', 'IA', 'ME', 'MI', 'MN', 'MO', 'MS', 'NC', 'NH',
                      'NM', 'NV', 'OH', 'PA', 'SC', 'TX', 'VA', 'WI')


"
*** Brier Score Computation (Method-1): ***
"


num_simulations <- dim(predicted_score)[1]
num_days <- dim(predicted_score)[2]


brier1_df <- data.frame(state = integer(),
                        mean_combined_prob = numeric())

for (i in 1:length(state_abbreviations)) {
  state_data <- predicted_score[, 254, i]

  probVal <- as.numeric(state_data > 0.5)
  probMean <- mean(probVal)  # Brier score calculation
  state_abbr <- state_abbreviations[i]
  actual_outcome <- election_results[state_abbr]
  brier_score <- mean((probMean - actual_outcome)^2)

  actual_outcome <- round(actual_outcome, 5)  # Rounding-off values
  probMean <- round(probMean, 5)
  brier_score <- round(brier_score, 5)

  brier1_df <- rbind(brier1_df, data.frame(state = state_abbr,
                                           actual = actual_outcome,
                                           prob = probMean,
                                           brier_score = brier_score))
}

write.csv(brier1_df, "brier-1(Unfiltered).csv", row.names = FALSE)
brier1_df_filtered <- brier1_df[brier1_df$state %in% specified_states,]
write.csv(brier1_df_filtered, "brier-1(Filtered).csv", row.names = FALSE)


"
*** Brier Score Computation (Method-2): ***
"

state_pairs <- combn(51, 2)  # Creating a matrix containing unique combinations of states

brier2_df <- data.frame(state1 = integer(),
                        state2 = integer(),
                        mean_combined_prob = numeric())

for (i in 1:ncol(state_pairs)) {
  state1 <- state_pairs[1, i]
  state2 <- state_pairs[2, i]

  state1_data <- predicted_score[, 254, state1]
  state2_data <- predicted_score[, 254, state2]

  state1_abbr <- state_abbreviations[state1]
  state2_abbr <- state_abbreviations[state2]

  CCval <- as.numeric((state1_data > 0.5) & (state2_data > 0.5))
  CCprob <- as.numeric((election_results[state1_abbr]>0.5) & (election_results[state2_abbr]>0.5))
  CTval <- as.numeric((state1_data > 0.5) & (state2_data < 0.5))
  CTprob <- as.numeric((election_results[state1_abbr]>0.5) & (election_results[state2_abbr]<0.5))
  TCval <- as.numeric((state1_data < 0.5) & (state2_data > 0.5))
  TCprob <- as.numeric((election_results[state1_abbr]<0.5) & (election_results[state2_abbr]>0.5))
  TTval <- as.numeric((state1_data < 0.5) & (state2_data < 0.5))
  TTprob <- as.numeric((election_results[state1_abbr]<0.5) & (election_results[state2_abbr]<0.5))

  CCmean <- round(mean(CCval),5)
  CCProbMean <- round(mean(CCprob),5)
  CTmean <- round(mean(CTval),5)
  CTProbMean <- round(mean(CTprob),5)
  TCmean <- round(mean(TCval),5)
  TCProbMean <- round(mean(TCprob),5)
  TTmean <- round(mean(TTval),5)
  TTProbMean <- round(mean(TTprob),5)

  CCbrier <- ((CCProbMean - CCmean)^2)
  CTbrier <- ((CTProbMean - CTmean)^2)
  TCbrier <- ((TCProbMean - TCmean)^2)
  TTbrier <- ((TTProbMean - TTmean)^2)
  AvgBrier <- ((CCbrier + CTbrier + TCbrier + TTbrier)/4)

  brier2_df <- rbind(brier2_df, data.frame(state1 = state1_abbr,
                                           state2 = state2_abbr,
                                           CC = CCbrier,
                                           CT = CTbrier,
                                           TC = TCbrier,
                                           TT = TTbrier,
                                           Avg = AvgBrier))
}

write.csv(brier2_df, "brier-2(Unfiltered).csv", row.names = FALSE)

brier2_df <- brier2_df[brier2_df$state1 %in% specified_states & brier2_df$state2 %in% specified_states,]

write.csv(brier2_df, "brier-2(Filtered).csv", row.names = FALSE)



"
*** Brier Score Computation (Method-3): ***
"

brier1_df_filtered <- brier1_df[brier1_df$state %in% specified_states, ]
brier2_df_filtered <- brier2_df[brier2_df$state1 %in% specified_states & brier2_df$state2 %in% specified_states,]

brier3_df <- merge(brier2_df_filtered, brier1_df_filtered, by.x = "state1", by.y = "state", all.x = TRUE)
brier3_df <- merge(brier3_df, brier1_df_filtered, by.x = "state2", by.y = "state", all.x = TRUE, suffixes = c("_state1", "_state2"))

brier3_df <- brier3_df[, c("state1", "brier_score_state1", "state2", "brier_score_state2", "CC", "CT", "TC", "TT", "Avg")]

brier3_df$CCcomb <- (brier3_df$brier_score_state1 + brier3_df$brier_score_state2 + brier3_df$CC)/3
brier3_df$CTcomb <- (brier3_df$brier_score_state1 + brier3_df$brier_score_state2 + brier3_df$CT)/3
brier3_df$TCcomb <- (brier3_df$brier_score_state1 + brier3_df$brier_score_state2 + brier3_df$TC)/3
brier3_df$TTcomb <- (brier3_df$brier_score_state1 + brier3_df$brier_score_state2 + brier3_df$TT)/3
brier3_df$AVGcomb <- (brier3_df$brier_score_state1 + brier3_df$brier_score_state2 + brier3_df$Avg)/3

write.csv(brier3_df, "brier-3.csv", row.names = FALSE)

"
-------------------------------------------
"


ex_states <- c('AZ', 'CO', 'FL', 'GA', 'IA', 'ME', 'MI', 'MN', 'MO', 'MS', 'NH', 'NM', 'NV', 'OH', 'PA', 'SC', 'TX', 'VA', 'WI')


pct_clinton %>% filter(t == RUN_DATE,state %in% c(ex_states,'--')) %>% mutate(se = (high - mean)/1.96) %>% dplyr::select(-t) %>% print
pct_clinton %>% filter(t == election_day,state %in% c(ex_states,'--')) %>% mutate(se = (high - mean)/1.96) %>% dplyr::select(-t) %>% print

pct_clinton %>% filter(t == election_day) %>% select(state, clinton=prob) %>% write_csv('today_2016.csv')
