#-------------------------------------------------------------

# File 03: power law  GOF

# call packages
library('igraph')
library('poweRlaw')

#-------------------------------------------------------------


# File contents:
# - A set of three functions to test for power law scaling vs Poisson scaling:
#    1. A function to calculate the degree distribution from real networks
#    2. A function to calculate ER random networks
#    3. A function to calculate LR tests regarding whether scaling ... 
#    ... varies from the scaling in random networks with the same G(n, m)
#    4. A function testing whether data meets requirements RE power law distribution ...
# ... see criticisms RE these requirements: https://www.barabasilab.com/post/love-is-all-you-need



# comparing degree distributions --------------------------------------------------

# For g, generate the degree distributions we observe --------------------
d = function(graph){
  
  # A degree distribution of the random graph
  d = igraph::degree(
    graph = graph, 
    mode = "total", 
    loops = FALSE, 
    normalized = FALSE
  )
  d <- d[ d != 0 ] # drop if degree k = 0
  return(d) # return the degree distribution as the result
  
  # End function
}

# For g, enerate Endos-Renyi (ER) random graphs by G(n,M) --------------------
r = function(graph){
  
  # set seed for replication purposes
  set.seed(20190812)
  
  # Get dimensions for G(n,M)
  n <- length(igraph::V(graph)) # 'n' vertices
  m <- length(igraph::V(graph)) # 'm' edges
  
  # calculate 10,000 samples/iterations
  samples <- 10000
  
  # 10,000 random networks G(n,m)
  r = c() # create a results vector
  for(i in 1 : samples) {
    # Calculate a random graph
    random = igraph::erdos.renyi.game(
      n = n, 
      p.or.m = m, 
      type = "gnm", 
      directed = FALSE, 
      loops = FALSE
    )
    # p.or.m = e for graphs of type G(n,M)... where n vertices and m edges ... 
    # ... m edges taken from the uniform random distribution from the set of all ...
    # ... possible edges 
    
    # A degree distribution of the random graph
    d = igraph::degree(
      graph = random, 
      mode = "total", 
      loops = FALSE, 
      normalized = FALSE
    )
    d <- d[d != 0] # drop if degree k = 0
    
    # put degree distribution into results vector
    r = c(r, d)  
  }
  
  # Average frequencies from the sampling distribution
  r = data.frame(r)
  r$n <-1
  r = aggregate(
    r$n, 
    by = list(r$r), 
    FUN = sum
    )
  r$x <- r$x/samples # averaging over the samples
  r$x <- round(x = r$x, digits = 0)
  r   <- rep(r$Group.1, r$x) # transform frequency table into a vector
  
  # return the results vector
  return(r)
  
  # End function
}

# compute degree distributions -------------------------
  
  # auto theft networks
  d_r_siren    = d(graph = r_siren)
  r_r_siren    = r(graph = r_siren)
  d_r_togo     = d(graph = r_togo)
  r_r_togo     = r(graph = r_togo)
  
  # drug trafficking networks
  d_d_caviar   = d(graph = d_caviar)
  r_d_caviar   = r(graph = d_caviar)
  d_d_cielnet  = d(graph = d_cielnet)
  r_d_cielnet  = r(graph = d_cielnet)
  d_d_cocaine  = d(graph = d_cocaine)
  r_d_cocaine  = r(graph = d_cocaine)
  d_d_heroin   = d(graph = d_heroin)
  r_d_heroin   = r(graph = d_heroin)
  
  # gang ties
  d_g_ity      = d(graph = g_ity)
  r_g_ity      = r(graph = g_ity)
  d_g_ldn      = d(graph = g_ldn)
  r_g_ldn      = r(graph = g_ldn)
  d_g_mtl      = d(graph = g_mtl)
  r_g_mtl      = r(graph = g_mtl)
  
  # mafia ties
  d_m_infinito = d(graph = m_infinito)
  r_m_infinito = r(graph = m_infinito)


# compare Vuong's likelihood ratio test for power law v. Poisson models -------------------------
vuong = function(d){
  
  # where dd = degree distribution from the above random graph generator
  
  # Fit power law distribution
  pl_fit = poweRlaw::displ(d)
  pl_fit$setXmin(1) # Evaluate @ degree k = 1
  pl_par = poweRlaw::estimate_pars(pl_fit)
  pl_fit$setPars(pl_par)
  
  # Fit Poisson distribution
  Pois_fit = poweRlaw::dispois(d)
  Pois_fit$setXmin(1) # Evaluate at degree k = 1
  Pois_par = poweRlaw::estimate_pars(Pois_fit)
  Pois_fit$setPars(Pois_par)
  
  # A LR test to compare nested models
  lrtest = poweRlaw::compare_distributions(pl_fit, Pois_fit)
  
  # Vuong's LR test statistic 
  lrstat = lrtest$test_statistic # because 'd' was input into the above funtion before ...
  # ... 'r', a negative test statistic indicates 'r' describes the degree distribution ...
  # ... and, assuming 'r' does give a better model fit, if we input 'd' and 'r' in ...
  # ... the reverse order, we'd get the same test statistic, but positive
  print( paste0( "Vuong's LR test statistic = ", lrstat ) )
  
  # p-value (two-tailed)
  p = lrtest$p_two_sided 
  # if p < 0.05, reject H0: neither power law or Poisson distribution provides good fit
  # if p > 0.05, fail to reject H1: degree distribution = power law or Poisson distribution
  print( paste0( "p-value = ", p ) )
  
  # Interpretation:
  print( paste0( "Interpretation:"))
  print( paste0( "A postive (+) > negative (-) LR test statistic suggests 'd' > 'r', or vice versa"))
  print( paste0( "> LR test statistics suggests better GOF 'd' vs. 'r'") )
  print( paste0( "If p < 0.05 ... reject H0, intepret by +/- LR test statistic") )
  
  # End function
}

# test results -------------------------
  
  # auto theft networks
  vuong(d = d_r_siren)
  vuong(d = r_r_siren)
  vuong(d = d_r_togo)
  vuong(d = r_r_togo)
  
  # drug trafficking networks
  vuong(d = d_d_caviar)
  vuong(d = r_d_caviar)
  vuong(d = d_d_cielnet)
  vuong(d = r_d_cielnet)
  vuong(d = d_d_cocaine)
  vuong(d = r_d_cocaine)
  vuong(d = d_d_heroin)
  vuong(d = r_d_heroin)
  
  # gangs 
  vuong(d = d_g_ity)
  vuong(d = r_g_ity)
  vuong(d = d_g_ldn)
  vuong(d = r_g_ldn)
  vuong(d = d_g_mtl)
  vuong(d = r_g_mtl)
  
  # mafia ties
  vuong(d = d_m_infinito)
  vuong(d = r_m_infinito)

  
  
# Function testing if data meets requirements for power law distribution ---------------

# As used in Clauset's papers ...
# see https://scholar.google.ca/scholar?hl=en&as_sdt=0%2C5&q=Goldstein+2004+power+law&btnG=
pk_fit = function(graph, k){
  
  # set initial conditions, where we test from the full distribution
  D = igraph::degree( # set degree distribution
    graph = graph,
    v = igraph::V(graph), 
    mode = "total", 
    loops = FALSE, 
    normalized = FALSE
  )
  
  # power law fit for Y for a given k cut point
  # for k / xmin, only values > k are used to fit the power law distribution
  Y = igraph::fit_power_law(
    x = D, 
    xmin = k, 
    implementation = 'plfit'
  )
  print( paste0("MLE for Y for a given k cut point:") )
  print( paste0("Y = ", Y$alpha) )
  print( paste0("log likelihood statistic = ", Y$logLik ) )
  print( paste0("Kolmogorov-Smirnov test statistic = ", Y$KS.stat ) )
  print( paste0("p-vale for Kolmogorov-Smirnov test statistic = ", Y$KS.p ) )
  
  # Interpretation
  print( paste0("Interpretation:") )
  print( paste0("Kolmogorov-Smirnov test compares the fitted distribution with the input vector --- smaller scores indicate better fit") )
  print( paste0("p-value < 0.05 for Kolmogorov-Smirnov test ---> reject H0: the original data could have been drawn from the fitted power-law distribution") )
  
  # End function -------------------------
}

# e.g., show tests for siren auto theft network
pk_fit( # i.e., for all k > 0, where k = degree
  graph = r_siren, 
  k = 1
  )
pk_fit( # i.e., for all k > 2
  graph = r_siren,
  k = 2
  )
pk_fit( # i.e., for all k > 4
  graph = r_siren, 
  k = 4
  )
# ... when 'k' for xmin isn't specified, the algorithm chooses ...
# ... the optimum value for k
pk_fit(
  graph = r_siren, 
  k = NULL
  )



# ... close .R script