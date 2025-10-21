//
// This Stan program defines a model where fish are released on a particular day
// and some number are recovered on subsequent days.
//
functions {
  
}
data {
  // Sizes
  int<lower=1> N_release; // Number of release days
  int<lower=1> N_pred; // Number of days to predict pcap and abundance
  int<lower=1> N_avail_obs; // Need extra rows for normalization
  int<lower=1> N_rec_obs; // Number of recapture observations
  int<lower=1> N_unm_obs; // Number of unmarked observations
  
  // Trap information
  // Currently using a vector here to avoid errors when mixing vectors and arrays
  vector<lower=0, upper=1>[N_pred] op;
  
  // Release information
  array[N_release] int<lower=1> num_released; // Number of individuals released
  array[N_release] int<lower=0> num_lost; // Number lost from each release
  
  // Availability information
  array[N_release + 1] int avail_idx; // ragged index for availability by release
  array[N_rec_obs] int<lower=1, upper=N_avail_obs> avail_rec_idx; // associate availability probs with recap days
  
  // Recapture information
  array[N_release + 1] int rec_idx; // ragged index of recaptures within each release
  array[N_rec_obs] int<lower=1, upper=N_pred> pred_rec_idx; // Associate pcap with recapture days
  array[N_rec_obs] int<lower=0> num_recaptured; // Number recaptured
  
  // Unmarked capture observations
  array[N_unm_obs] int<lower=1, upper=N_pred> unm_idx; // associate pcaps and unmarked counts, eliminating no-trap days.
  // Number of unmarked individuals captured
  array[N_unm_obs] int<lower=0> unm_count;
  
  // Availability model dimensions
  int<lower=0> N_avail_fixed;
  int<lower=0> N_avail_prior;
  int<lower=0> N_avail_random;
  int<lower=0> N_avail_hprior;
  // Availability model 
  matrix[N_avail_obs, N_avail_fixed] avail_fixed;
  vector<lower=0>[N_avail_prior] avail_prior_sd;
  vector[N_avail_prior] avail_prior_mean;
  array[N_avail_fixed] int<lower=0, upper=N_avail_fixed> avail_prior_idx;
  matrix[N_avail_obs, N_avail_random] avail_random;
  vector<lower=0>[N_avail_hprior] avail_hprior;
  array[N_avail_random] int<lower=0, upper=N_avail_random> avail_hprior_idx;
  
  // Probability of capture model dimensions
  int<lower=0> N_pcap_fixed;
  int<lower=0> N_pcap_prior;
  int<lower=0> N_pcap_random;
  int<lower=0> N_pcap_hprior;
  // Availability model 
  matrix[N_pred, N_pcap_fixed] pcap_fixed;
  vector<lower=0>[N_pcap_prior] pcap_prior_sd;
  vector[N_pcap_prior] pcap_prior_mean;
  array[N_pcap_fixed] int<lower=0, upper=N_pcap_fixed> pcap_prior_idx;
  matrix[N_pred, N_pcap_random] pcap_random;
  vector<lower=0>[N_pcap_hprior] pcap_hprior;
  array[N_pcap_random] int<lower=0, upper=N_pcap_random> pcap_hprior_idx;
  
  // Runsize model dimensions
  int<lower=0> N_rs_fixed;
  int<lower=0> N_rs_prior;
  int<lower=0> N_rs_random;
  int<lower=0> N_rs_hprior;
  // Runsize model 
  matrix[N_pred, N_rs_fixed] rs_fixed;
  vector<lower=0>[N_rs_prior] rs_prior_sd;
  vector[N_rs_prior] rs_prior_mean;
  array[N_rs_fixed] int<lower=0> rs_prior_idx;
  matrix[N_pred, N_rs_random] rs_random;
  vector<lower=0>[N_rs_hprior] rs_hprior;
  array[N_rs_random] int<lower=0> rs_hprior_idx;
  
  // Sample mixture posterior for mixture importance sampling LOO?
  int mixis;
}
transformed data {
  vector[N_avail_fixed] avail_offset = avail_prior_mean[avail_prior_idx];
  vector[N_avail_fixed] avail_mult = avail_prior_sd[avail_prior_idx];
  
  vector[N_pcap_fixed] pcap_offset = pcap_prior_mean[pcap_prior_idx];
  vector[N_pcap_fixed] pcap_mult = pcap_prior_sd[pcap_prior_idx];
  
  vector[N_rs_fixed] rs_offset = rs_prior_mean[rs_prior_idx];
  vector[N_rs_fixed] rs_mult = rs_prior_sd[rs_prior_idx];
}
parameters {
  // Cannot declare both multipliers and limits, so 
  // Travel time parameters
  vector<offset=avail_offset, multiplier=avail_mult>[N_avail_fixed] avail_coef_fixed;
  vector<lower=0>[N_avail_hprior] avail_rand_scale;
  vector<multiplier=avail_rand_scale[avail_hprior_idx]>[N_avail_random] avail_coef_random;
  
  // Capture probability parameters (uc = unconstrained)
  vector<offset=pcap_offset, multiplier=pcap_mult>[N_pcap_fixed] pcap_coef_fixed;
  vector<lower=0>[N_pcap_hprior] pcap_rand_scale;
  vector<multiplier=pcap_rand_scale[pcap_hprior_idx]>[N_pcap_random] pcap_coef_random;
  
  // Runsize parameters
  vector<offset=rs_offset, multiplier=rs_mult>[N_rs_fixed] rs_coef_fixed;
  vector<lower=0>[N_rs_hprior] rs_rand_scale;
  vector<multiplier=rs_rand_scale[rs_hprior_idx]>[N_rs_random] rs_coef_random;
}
transformed parameters {
  vector[N_avail_obs] avail_uc; // Unconstrained availability probabilities
  vector[N_avail_obs] avail; // Availability probabilities
  vector[N_rec_obs] avail_rec; // Availability reindexed to match recaptures
  
  vector[N_pred] pcap_uc; // Unconstrained recapture probabilities
  vector[N_pred] pcap; // probablities of capture
  vector[N_rec_obs] pcap_rec; // Probability of recapture 
  vector[N_release] p_lost;
  vector[N_rec_obs] p_recap;
  
  vector[N_pred] log_lambda;
  vector[N_pred] log_pcap;
  
  // Pre-compute the probability that an individual fish arrives each day within
  // the recapture window. Precalculate `_coef_` values so that they can be
  // tracked directly.
  avail_uc = avail_fixed * avail_coef_fixed
             + avail_random * avail_coef_random;
  // Calculate the logit probability of capture - multiply by `pcap_sigma`
  // to use non-centered parameterization, will be useful later when smoothing
  // penalty is estimated rather than fixed.
  pcap_uc = pcap_fixed * pcap_coef_fixed + pcap_random * pcap_coef_random;
  // Convert to probability scale
  pcap = inv_logit(pcap_uc);
  // Index into the predicted probability of capture to associate with each
  // recapture
  pcap_rec = pcap[pred_rec_idx];
  
  // Calculate the probability of recapture after each release
  for (rel in 1 : N_release) {
    // Index into the availability probability vector; need the full recapture
    // window so that they can be comparably normalized
    int ai = avail_idx[rel];
    int aj = avail_idx[rel + 1] - 1;
    // Calculate the probability of availability for each day in the release
    // window
    avail[ai : aj] = softmax(avail_uc[ai : aj]);
  }
  // Index into the availability probability vector to 
  avail_rec = avail[avail_rec_idx];
  
  for (rel in 1 : N_release) {
    int i = rec_idx[rel];
    int j = rec_idx[rel + 1] - 1;
    // Calculate the probability (proportion expected) for each day after
    // release. 
    p_recap[i : j] = pcap_rec[i : j] .* avail_rec[i : j];
    // Final entry is for all unrecaptured individuals
    p_lost[rel] = 1 - sum(p_recap[i : j]);
  }
  
  // Run size and predicted capture efficiency for all days
  log_lambda = rs_fixed * rs_coef_fixed + rs_random * rs_coef_random;
  // op is 0 or 1. If the trap is operational, log(1) = 0, and if not log(0) =
  // -Inf, indicating that probability of capture that day is zero.
  log_pcap = log_inv_logit(pcap_uc) + log(op);
}
model {
  // Vectos to save individual observation log-likelihoods
  vector[N_release] mr_log_lik;
  vector[N_unm_obs] um_log_lik;
  vector[N_release + N_unm_obs] log_lik;
  
  // Availability model coefficients
  avail_coef_fixed ~ normal(avail_offset, avail_mult);
  avail_rand_scale ~ normal(0, avail_hprior);
  avail_coef_random ~ normal(0, avail_rand_scale[avail_hprior_idx]);
  
  // Prior on probability of capture coefficients over time
  pcap_coef_fixed ~ normal(pcap_offset, pcap_mult);
  pcap_rand_scale ~ normal(0, pcap_hprior);
  pcap_coef_random ~ normal(0, pcap_rand_scale[pcap_hprior_idx]);
  
  // Noncentered prior on run size estimates coefficients
  rs_coef_fixed ~ normal(rs_offset, rs_mult);
  rs_rand_scale ~ normal(0, rs_hprior);
  rs_coef_random ~ normal(0, rs_rand_scale[rs_hprior_idx]);
  
  // Observation likelihood; multinomial including the lost group
  for (rel in 1 : N_release) {
    int i = rec_idx[rel];
    int j = rec_idx[rel + 1] - 1;
    // Need to add 2 here! j - i + 1 is the length of i:j
    array[j - i + 2] int obs = append_array(num_recaptured[i : j],
                                            num_lost[rel : rel]);
    vector[j - i + 2] p = append_row(p_recap[i : j], p_lost[rel]);
    // Number of recaptures and probability of recapture is stored as a single
    // long vector. This is indexed into so that we can treat e.g. each species
    // with a different recapture window, effectively giving us a ragged array.
    if (mixis) {
      mr_log_lik[rel] = multinomial_lpmf(obs | p);
    } else {
      obs ~ multinomial(p);
    }
  }
  
  // Likelihood of observing thinned Poisson process. Note that we need
  // to use the vector-vector addition here; we can't use elementwise `.+`
  // operator or we get an "invalid character found" error when lexing.
  // unm_count ~ poisson_log(log_lambda[unm_idx] + log_pcap[unm_idx]);
  if (mixis) {
    for (i in 1 : N_unm_obs) {
      um_log_lik[i] = poisson_log_lpmf(unm_count[i] | log_lambda[unm_idx][i]
                                                      + log_pcap[unm_idx][i]);
    }
    
    log_lik = append_row(mr_log_lik, um_log_lik);
    target += sum(log_lik);
    
    target += log_sum_exp(-log_lik);
  } else {
    unm_count ~ poisson_log(log_lambda[unm_idx] + log_pcap[unm_idx]);
  }
}
generated quantities {
  // Posterior predictive quantities
  array[N_rec_obs] int recap_rep;
  array[N_release] int lost_rep;
  array[N_unm_obs] int unm_rep;
  
  for (rel in 1 : N_release) {
    int i = rec_idx[rel];
    int j = rec_idx[rel + 1] - 1;
    
    // Need to add 2 here! j - i + 1 is the length of i:j
    array[j - i + 2] int obs;
    vector[j - i + 2] p = append_row(p_recap[i : j], p_lost[rel]);
    obs = multinomial_rng(p, num_released[rel]);
    
    recap_rep[i : j] = head(obs, j - i + 1);
    lost_rep[rel] = obs[j - i + 2];
  }
  
  unm_rep = poisson_log_rng(log_lambda[unm_idx] + log_pcap[unm_idx]);
  
  // Individual observation log-likelihoods
  // vector[sum(num_released)] mr_loglik;
  vector[N_release] mr_log_lik;
  vector[N_unm_obs] um_log_lik;
  vector[N_release + N_unm_obs] log_lik;
  // int mr_idx = 1;
  
  for (rel in 1 : N_release) {
    int i = rec_idx[rel];
    int j = rec_idx[rel + 1] - 1;
    // Need to add 2 here! j - i + 1 is the length of i:j
    array[j - i + 2] int obs = append_array(num_recaptured[i : j],
                                            num_lost[rel : rel]);
    vector[j - i + 2] p = append_row(p_recap[i : j], p_lost[rel]);
    mr_log_lik[rel] = multinomial_lpmf(obs | p);
  }
  
  for (i in 1 : N_unm_obs) {
    um_log_lik[i] = poisson_log_lpmf(unm_count[i] | log_lambda[unm_idx][i]
                                                    + log_pcap[unm_idx][i]);
  }
  
  log_lik = append_row(mr_log_lik, um_log_lik);
}
