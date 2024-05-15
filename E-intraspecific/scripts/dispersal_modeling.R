# Libaries
library(rethinking)

# Functions ---------------------------------------------------------------

# Find approximate alpha and beta parameters for a gamma distribution based on
# observing three quantiles
# Based on post below  - many thanks to "Sextus Empiricus"
#https://stats.stackexchange.com/questions/596388/fit-gamma-distribution-based-on-median-interquartile-range

# fit_gamma_pars() returns the fit of the supplied parameters, observed quantile values, and quantile probability thresholds
# using the Kolmogorov statistic
# large values are poor fits 
fit_gamma_pars = function(par=c(1,1), obs_quantiles, p_thresholds=c(0.025,0.5,0.975)) {
  p_quantiles = pgamma(obs_quantiles, shape = par[1], scale = par[2])
  statistic = max(abs(p_quantiles - p_thresholds))
  return(statistic)
} 

#approx_gamma_pars() returns the approximate shape and scale parameters based on Nelder and Mead (1965) approximation
#par is the starting gamma distribution parameter values
approx_gamma_pars<-function(par  = c(1,1),obs_quantiles, p_thresholds=c(0.025,0.5,0.975)) {
  params<-optim(par, 
                fit_gamma_pars, 
                par,
                obs_quantiles,
                p_thresholds)$par
  return(params)
}

# Looping through all taxa with significant slopes
for (taxa in c("AA1", "AA2", "AH1", "AH2", "AH3")) {
  print(taxa)
  # Construct priors for density -------------------------------------------
  
  # CENSUS density estimates per plot (units are De/m2)
  DNC <- read.csv(paste0("../data/",
                         taxa, ".census.txt"))
  
  #Estimating the frequencies as a log-normal distribution -> forces them to be positive
  freqs <- list(DNC_log = log(DNC$Freq))
  
  # intercept only model: mean has normal prior distribution ~N(1,1); error is exponentially distributed
  intercept_m <- ulam(
    alist( 
      DNC_log ~ dnorm(mu, error), 
      mu<-a,
      a ~ dnorm(1 , 1),
      error ~ dexp(1)
    ), data=freqs , chains=4 , log_lik=TRUE )
  
  precis(intercept_m) # = posterior for intercept and error on log scale
  exp(precis(intercept_m)$mean[1]) #= mean posterior freq on regular scale
  exp(precis(intercept_m )$sd[1])
  
  #sample DNC directly from posterior predicted distribution to use in your error propagation
  DNC_post<-exp(extract.samples(intercept_m)$a) #in D/m2
  
  # Make table of results
  result_table <- data.frame(param = numeric(0),
                             lower.10 = numeric(0),
                             mid.50 = numeric(0),
                             upper.90 = numeric(0))
  
  print("DNC")
  print(signif(quantile(DNC_post, probs = c(0.1, 0.5, 0.9)), 2))
  dnc.row <- c("DNC", round(quantile(DNC_post, probs = c(0.1, 0.5, 0.9)), 2))
  result_table <- rbind(result_table, dnc.row)
  
  # KINSHIP Ne density - assume gamma distribution b/c upper tail can be large
  Ne_estimates <- read.csv("../data/agaricia_effective_density.txt")
  Ne_estimates <- Ne_estimates[Ne_estimates$taxa == taxa,]
  
  # Use the 99% percentile as upper limit for NE estimates due to low confidence
  if (taxa == "AH1" | taxa == "AH3") {
    Ne_estimates$Ne_high <- round(quantile(DNC_post, probs = c(0.99))*Ne_estimates$Area, 0)
  }
  
  DNE<-list()
  
  for(r in 1:nrow(Ne_estimates)) {
    De<-Ne_estimates[r,"Ne"]/(Ne_estimates[r,"Area"])
    De_low<-Ne_estimates[r,"Ne_low"]/(Ne_estimates[r,"Area"])
    De_high<-Ne_estimates[r,"Ne_high"]/(Ne_estimates[r,"Area"])
    shape<-approx_gamma_pars(obs_quantiles=c(De_low, De, De_high))[1]
    scale<-approx_gamma_pars(obs_quantiles=c(De_low, De, De_high))[2]
    DNE[[r]] <-rgamma(2000, shape=shape, scale=scale)
  }
  
  DNE<-unlist(DNE)
  hist(log10(DNE))
  
  print("DNE")
  print(round(quantile(DNE, probs = c(0.1, 0.5, 0.9)), 2))
  dne.row <- c("DNE", round(quantile(DNE, probs = c(0.1, 0.5, 0.9)), 2))
  result_table <- rbind(result_table, dne.row)
  
  # Slope - estimating gamma function for prior -----------------------------
  
  #Different using Nb for Rousset and Loiselle
  if (taxa == "AA1") {
    nb_quantiles_r <- c(25, 33, 46)
    nb_quantiles_l <- c(35, 45, 60)
    
  } else if (taxa == "AA2") {
    nb_quantiles_r <- c(274, 582, Inf)
    nb_quantiles_l <- c(361, 501, 818)
    
  } else if (taxa == "AH1") {
    nb_quantiles_l <- c(57, 82, 142)
  } else if (taxa == "AH3") {
    nb_quantiles_l <- c(17, 24, 38)
  }
  
  # Loiselle for all taxa
  beta_quantiles_l <- rev(1/nb_quantiles_l)
  beta_shape_l<-approx_gamma_pars(obs_quantiles=c(beta_quantiles_l))[[1]]
  beta_scale_l<-approx_gamma_pars(obs_quantiles=c(beta_quantiles_l))[[2]]
  
  # Create generative simulation model to simulate sigma ----------------------
  # Simulate sigma calculations from isolation by distance slopes
  
  # Model: 
  # sigma ~ Exp(mean) #LaPlacian or half-Gaussian possible too 
  # sigma = sqrt(1/(4piDb))
  # DNC: use posterior observations 
  # DNE: use simulated values (already created above) with parameters estimated from empirical kinship data
  # beta: simulate values from empirical slope estimates 
  
  # Simulate beta (loiselle)
  sims<-2000
  beta_l<-rgamma(sims, beta_shape_l, scale=beta_scale_l)
  
  # Estimate census sigma (loiselle)
  sigma_DNC<-sqrt(1/(4*pi*DNC_post*beta_l))
  hist(log10(sigma_DNC))
  print("Sigma DNC")
  print(round(quantile(sigma_DNC, probs=c(0.1, 0.5, 0.9)),2))
  s.dnc.row <- c("sigma_DNC", round(quantile(sigma_DNC, probs = c(0.1, 0.5, 0.9)), 2))
  result_table <- rbind(result_table, s.dnc.row)
  
  # Estimate effective sigma (loiselle)
  sigma_DNE<-sqrt(1/(4*pi*DNE*beta_l))
  hist(log10(sigma_DNE))
  print("Sigma DNE")
  print(round(quantile(sigma_DNE, probs=c(0.1, 0.5, 0.9)),2))
  s.dne.row <- c("sigma_DNE", round(quantile(sigma_DNE, probs = c(0.1, 0.5, 0.9)), 2))
  result_table <- rbind(result_table, s.dne.row)
  
  # Estimate Nb - just from beta so same for NE and NC
  Neighborhood<-4*pi*DNC_post*sigma_DNC^2
  hist(log(Neighborhood))
  print("Neighbourhood")
  print(round(quantile(Neighborhood, probs=c(0.1, 0.5, 0.9)), 0))
  nb.row <- c("Nb", round(quantile(Neighborhood, probs = c(0.1, 0.5, 0.9)), 0))
  result_table <- rbind(result_table, nb.row)
  
  # save output
  write.csv(DNC_post, paste0("../results/distributions/", taxa,
                             ".DNC_copy.txt"), quote = FALSE, row.names = FALSE)
  write.csv(DNE, paste0("../results/distributions/", taxa,
                        ".DNE_copy.txt"), quote = FALSE, row.names = FALSE)
  
  # loiselle only
  write.csv(Neighborhood,
            paste0("../results/distributions/", taxa,
                   ".Neighborhood-loiselle_copy.txt"),
            quote = FALSE, row.names = FALSE)
  write.csv(sigma_DNC, paste0("../results/distributions/", taxa,
                              ".sigma_DNC-loiselle_copy.txt"), quote = FALSE, row.names = FALSE)
  write.csv(sigma_DNE, paste0("../results/distributions/", taxa,
                              ".sigma_DNE-loiselle_copy.txt"), quote = FALSE, row.names = FALSE)
  
  # Rousset for AA1 and AA2
  if (taxa == "AA1" | taxa == "AA2") {
    beta_quantiles_r <- rev(1/nb_quantiles_r)
    beta_shape_r<-approx_gamma_pars(obs_quantiles=c(beta_quantiles_r))[[1]]
    beta_scale_r<-approx_gamma_pars(obs_quantiles=c(beta_quantiles_r))[[2]]
    beta_r<-rgamma(sims, beta_shape_r, scale=beta_scale_r)
    
    # Estimate census sigma (rousset)
    sigma_DNC<-sqrt(1/(4*pi*DNC_post*beta_r))
    hist(log10(sigma_DNC))
    print("Sigma DNC")
    print(round(quantile(sigma_DNC, probs=c(0.1, 0.5, 0.9)),2))
    s.dnc.row <- c("sigma_DNC_r", round(quantile(sigma_DNC, probs = c(0.1, 0.5, 0.9)), 2))
    result_table <- rbind(result_table, s.dnc.row)
    
    # Estimate effective sigma (rousset)
    sigma_DNE<-sqrt(1/(4*pi*DNE*beta_r))
    hist(log10(sigma_DNE))
    print("Sigma DNE")
    print(round(quantile(sigma_DNE, probs=c(0.1, 0.5, 0.9)),2))
    s.dne.row <- c("sigma_DNE_r", round(quantile(sigma_DNE, probs = c(0.1, 0.5, 0.9)), 2))
    result_table <- rbind(result_table, s.dne.row)
    
    # Estimate Nb - just from beta so same for NE and NC
    Neighborhood<-4*pi*DNC_post*sigma_DNC^2
    hist(log(Neighborhood))
    print("Neighbourhood")
    print(round(quantile(Neighborhood, probs=c(0.1, 0.5, 0.9)), 0))
    nb.row <- c("Nb_r", round(quantile(Neighborhood, probs = c(0.1, 0.5, 0.9)), 0))
    result_table <- rbind(result_table, nb.row)
    
    # Rousset only
    write.csv(Neighborhood,
              paste0("../results/distributions/", taxa,
                     ".Neighborhood-rousset_copy.txt"),
              quote = FALSE, row.names = FALSE)
    write.csv(sigma_DNC, paste0("../results/distributions/", taxa,
                                ".sigma_DNC-rousset_copy.txt"), quote = FALSE, row.names = FALSE)
    write.csv(sigma_DNE, paste0("../results/distributions/", taxa,
                                ".sigma_DNE-rousset_copy.txt"), quote = FALSE, row.names = FALSE)
 
    
  }
  colnames(result_table) <- c("param", "low10", "mid50", "high90")
  write.csv(result_table, paste0("../results/distributions/", taxa,
                                  ".results-summary_copy.txt"), quote = FALSE, row.names = FALSE)
  }
  
  