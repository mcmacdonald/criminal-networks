#  -----------------------------------------------------------------------------------

# file 03: estimate the hierarchical network model

# 'mlergm' package
# https://cran.r-project.org/web/packages/mlergm/mlergm.pdf

# don't run
# install.package('mlergm')

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------



# transform igraph objects to network objects to estimate models
graph <- function(e, v){
  g <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = v)
  g <- igraph::simplify(g) # force the graph to not be of type 'multiple'
  g <- intergraph::asNetwork(g)
  return(g)
}
g_siren  <- graph(siren, v = v_siren)
g_togo   <- graph(togo, v = v_togo)
g_caviar <- graph(caviar, v = v_caviar)
g_cielnet  <- graph(cielnet, v = v_cielnet)
g_cocaine  <- graph(cocaine, v = v_cocaine)
g_heroin   <- graph(heroin, v = v_heroin)
g_oversize <- graph(oversize, v = v_oversize)
g_montagna <- graph(montagna, v = v_montagna)



# plot complimentary cumulative degree distribution -----------------------------
cdf_degree <- function(g){
  par(mfrow=c(1,1)) # plot dimensions
  # normalized degree centrality
  d <- sna::degree(g, gmode = "graph") # degree centrality for each node
  d <- d/max(d); d <- round(x = d, digits = 2) # scores to range 0-1 and round to two decimal places
  d <- d[order(d, decreasing = TRUE)] # sort from largest to smallest
  # calculate complimentary cumulative distribution function (ccdf)
  cdf <- stats::ecdf(d) # function to caclculate the empirical cdf
  ccdf <- 1 - cdf(d) # ccdf
  # join ccdf and normalized degree centrality
  data <- cbind(ccdf, d); data <- as.data.frame(data)
  colnames(data) <- c("ccdf", "ndegree") # column names 
  # plot the ccdf
  plot(log(data$ndegree), log(data$ccdf))
}
cdf_degree(g_siren)
cdf_degree(g_togo)
cdf_degree(g_caviar)
cdf_degree(g_cielnet)
cdf_degree(g_cocaine)
cdf_degree(g_heroin)
cdf_degree(g_oversize)
cdf_degree(g_montagna)
cdf_degree(g_super)



# posterior parameter estimation ---------------------------------------------
bayes <- function(y, x){
  i <- nrow(y) * 2  # burn in iterations to begin the MCMC run 
  k <-    250  # sample iterations
  h <-    5 * 2   # chains in the posterior distribution ... approximately twice the number of model parameters 
  n <-    h * k   # per Caimo & Friel (2011), auxiliary chain = # chains (h) * sample iterations (k)
  # load 'bergm' package
  require('Bergm')
  set.seed(20110210) # Halle's birthday
  bayes <- Bergm::bergm(
    x, # equation
    gamma = 0.1 # empirically, gamma ranges 0.1 - 1.5, where smaller gamma leads to better acceptance in big graphs
    )
  return(bayes)
}

# goodness-of-fit --------------------------------------------------------------
GOF <- function(bayes){
  set.seed(20110211) # Halle's birthday
  n <- 100 # graph simulations
  i <- 15 # degree distribution range
  j <- 10 # geodisstance range 
  k <- 15 # edgewise-shared partners range
  fit <- Bergm::bgof(
    bayes,
    n.deg  = i,
    n.dist = j,
    n.esp  = k,
    directed = F,    # symmetric graph
    sample.size = n # random graph realizations
  )
  return(fit)
}

# estimate the model -----------------------------------------------------------
bayes_01.siren <- bayes(
  y = g_siren,
  x = g_siren ~ edges + 
    gwdegree(decay = 3.0, fixed = TRUE) +
    gwdsp(decay = 3.0, fixed = TRUE) +
    gwesp(decay = 3.0, fixed = TRUE) +
    degcor
    )
summary(bayes_01.siren)

# goodness-of-fit
gof_01.siren <- GOF(bayes_01.siren)


# estimate the model -----------------------------------------------------------
bayes_02.togo <- bayes(
  y = g_togo,
  x = g_togo ~ edges + 
    gwdegree(decay = 0.5, fixed = TRUE) +
    gwdsp(decay = 3.0, fixed = TRUE) +
    gwesp(decay = 3.0, fixed = TRUE) +
    degcor
    )
summary(bayes_02.togo)

# goodness-of-fit
gof_02.togo <- GOF(bayes_02.togo)


# estimate the model -----------------------------------------------------------
bayes_03.caviar <- bayes(
  y = g_caviar,
  x = g_caviar ~ edges + 
    gwdegree(decay = 2.0, fixed = TRUE) +
    gwdsp(decay = 1.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_03.caviar)

# goodness-of-fit
gof_03.caviar <- GOF(bayes_03.caviar)


# estimate the model -----------------------------------------------------------
bayes_04.cielnet <- bayes(
  y = g_cielnet,
  x = g_cielnet ~ edges + 
    gwdegree(decay = 1.0, fixed = TRUE) +
    gwdsp(decay = 0.2, fixed = TRUE) +
    gwesp(decay = 0.2, fixed = TRUE) +
    degcor
    )
summary(bayes_04.cielnet)

# goodness-of-fit
gof_04.cielnet <- GOF(bayes_04.cielnet)


# estimate the model -----------------------------------------------------------
bayes_05.cocaine <- bayes(
  y = g_cocaine,
  x = g_cocaine ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) +
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_05.cocaine)

# goodness-of-fit
gof_05.cocaine <- GOF(bayes_05.cocaine)


# estimate the model -----------------------------------------------------------
bayes_06.heroin <- bayes(
  y = g_heroin,
  x = g_heroin ~ edges + 
    gwdegree(decay = 2.5, fixed = TRUE) +
    gwdsp(decay = 1.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_06.heroin)

# goodness-of-fit
gof_06.heroin <- GOF(bayes_06.heroin)


# estimate the model -----------------------------------------------------------
bayes_07.oversize <- bayes(
  y = g_oversize,
  x = g_oversize ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) +
    gwdsp(decay = 1.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_07.oversize)

# goodness-of-fit
gof_07.oversize <- GOF(bayes_07.oversize)


# estimate the model -----------------------------------------------------------
bayes_08.montagna <- bayes(
  y = g_montagna,
  x = g_montagna ~ edges + 
    gwdegree(decay = 2.0, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_08.montagna)

# goodness-of-fit
gof_08.montagna <- GOF(bayes_08.montagna)







# close .r script






# write equations --------------------------------------------------------------
x_01 = g_mob ~
  edges                            + 
  gwdegree(decay = 3, fixed = T)   +
  gwdsp(decay = 1, fixed = T)      +
  gwesp(decay = 1, fixed = T)      +
  degcor                           +
  nodefactor('family')             +
  nodefactor('rank', levels=-6)    +
  dyadcov(g_Family)                +
  dyadcov(g_upper)                 +
  dyadcov(g_lower)                 +
  dyadcov(d_upper)                 +
  dyadcov(d_lower)                 +
  # hierarchical relations by subgraph -----------------------------------------
S(~edges + 
    gwdegree(decay = 1, fixed = T) + 
    gwdsp(decay = 1, fixed = T)    + 
    gwesp(decay = 1, fixed = T)    ,
  ~(family == c("1_bonanno" )
  )
)



# hierarchical relations by subgraph -----------------------------------------
S(~nodemix('rank_c'), ~(family_c == c("1_bonanno" )))     +
  # S(~nodemix('rank_c'), ~(family_c == c("2_decavalcante"))) +
  S(~nodemix('rank_c'), ~(family_c == c("3_gambino" )))     +
  S(~nodemix('rank_c'), ~(family_c == c("4_genovese")))     +
  S(~nodemix('rank_c'), ~(family_c == c("5_lucchese")))     +
  S(~nodemix('rank_c'), ~(family_c == c("6_profaci" )))     +
  
  
  
  



#################################################################################


# construct the supergraph as a 'mlergm' object
g_super <- mlergm::mlnet(
  network = g_super, 
  node_memb = network::get.vertex.attribute(g_super, "group")
  )
# check that graph is of type multi-level network
mlergm::is.mlnet(g_super)
# plot the super network
plot(g_super, arrow.size = 2.5, arrow.gap = 0.025)



# estimate the model
model <- mlergm::mlergm(
  g_super ~ edges + 
    gwdegree(decay = 1.0, cutoff = 100) + 
    gwdsp(decay = 1.0) + 
    gwesp(decay = 1.0),
  parameterization = 'offset',
  options = set_options(
    burnin = 10000,
    interval = 1000,
    sample_size = 100000,
    NR_tol = 1e-04,
    NR_max_iter = 50,
    MCMLE_max_iter = 10,
    do_parallel = TRUE,
    # NR_step_len = 10
    adaptive_step_len = TRUE
    ),
  verbose = 2, # = 2 prints the full output
  seed = 123
)
summary(model_est)

# goodness-of-fit plots
model_gof <- mlergm::gof(model)
plot(
  gof_res, 
  cutoff = 15, 
  pretty_x = TRUE
  )





# function to estimate exponential random graph models --------------------------------------------------------------------
ERGM <- function(g, x){
  n <- nrow(g)
  k <- n * 100000 # 100,000 iterations per node
  b <- ergm::ergm(
    x, # mdoel specifciation
    estimate = 'MPLE',
    control = ergm::control.ergm(
      main.method = 'MCMLE', # see Snijders & van Duijn (2002) for info on 'Robbins-Monro' method
      MCMC.burnin = k,   # 'more is better' ... = v x 100,000
      MCMC.interval = k, # 'more is better' ... = v x 100,000
      MCMC.prop.weights = 'TNT', # faster convergence when paired with 'MPLE'
      seed = 20110210 # to replicate ERGM
    ),
    verbose = TRUE    # ... for networks with overall low density
  )
  print(summary(b))  # print results
  print(confint(b, level = 0.95)) # 95% confidence intervals
  return(b) # return ergm object
}
# pass model specifciations through function to estimate models 
model_01.siren <- ERGM(
  g = g_siren,
  x = g_siren ~ 
    edges + 
    gwdegree(decay = 3, fixed = TRUE) +
    gwdsp(decay = 1, fixed = TRUE) +
    gwesp(decay = 1, fixed = TRUE) +
    degcor
)
model_02.togo <- ERGM(
  g = g_siren,
  x = g_siren ~ 
    edges + 
    gwdegree(decay = 3, fixed = TRUE) +
    gwdsp(decay = 1, fixed = TRUE) +
    gwesp(decay = 1, fixed = TRUE) +
    degcor
)



# Bayesian exponential random graphs -------------------------------------------

# first, compute Bayesian priors ---------------------------------------------

# function to compute variance-covariance structure
sigma <- function(b){
  n <- length(unlist(b))
  s <- matrix(0:0, nrow = n, ncol = n) # where 'input' = number ERGM parameters 
  diag(s) <- 1 # set diagonals = 1
  return(s) # return matrix
}

# don't run
# compute vector of the multivariate normal priors outside of function
# mu_05 <- model_05$coefficients # model parameters
# sd_05 <- sigma(mu = mu_05) # variance-covariance structure/matrix
# diag(sd_01) <- coef(summary(rep_01))[, "Std. Error"]

# function to compute vector of the multivariate normal priors
vcov_ergm <- function(model){
  b <- model$coefficients
  vcov <- sigma(b) # embed prior function
  diag(vcov) <- coef(summary(model))[, "Std. Error"]
  return(vcov)
}
vcov_01.siren <- vcov_ergm(model_01.siren)






# posterior parameter estimation ---------------------------------------------
bayes <- function(g, x, npar, p, s){
  i <- nrow(g) * 2  # burn in iterations to begin the MCMC run 
  k <-    250  # sample iterations
  h <-    ergm::nparam(npar) * 2   # chains in the posterior distribution ... approximately twice the ERGM parameters 
  n <-    h * k   # per Caimo & Friel (2011), auxiliary chain = # chains (h) * sample iterations (k)
  # load 'bergm' package
  require('Bergm')
  set.seed(20110210) # Halle's birthday
  bayes <- Bergm::bergm(
    x,              # equation
    burn.in = i,     # burn ins
    main.iters = k,  # sample iterations
    aux.iters = n,   # auxiliary chains 
    nchains = h,     # essentially the # chains in the posterior distribution
    prior.mean = p,  # prior means
    prior.sigma = s, # prior variance/covariance structure
    gamma = 0.1      # empirically, gamma ranges 0.1 - 1.5, where smaller gamma leads to better acceptance in big graphs
  )
  return(bayes)
}