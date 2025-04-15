#  -----------------------------------------------------------------------------------

# file 06: estimate the hierarchical network model

# 'mlergm' package
# https://cran.r-project.org/web/packages/mlergm/mlergm.pdf

# don't run
# install.package('mlergm')

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------



# posterior parameter estimation -----------------------------------------------
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

# function to calculate the decay parameter for the geometrically weighted degree statistic
gwd_decay <- function(g, decay_range = seq(0.1, 5, by = 0.1), method = c("sse", "cor")) {
  
  # method, which can take one of two options
  method <- match.arg(method)
  
  # if not an network object, transform to one
  if (network::is.network(g)==TRUE){
    g <- g
  } else {
    g <- intergraph::asIgraph(g)
  }
  
  # degree distribution
  d_stat <- sna::degree(g, gmode = "graph")
  d_table <- table(d_stat)
  d_freq  <- as.numeric(d_table)
  d_value <- as.numeric(names(d_table))
  
  # normalize degree distribution
  d_norm <- d_freq / sum(d_freq)
  
  # calculate different decay values
  scores <- sapply(decay_range, function(decay) {
    
    weights <- (1 - exp(-decay))^d_value
    gwdegree <- sum(weights * d_freq)
    
    # Compare modeled vs actual
    if (method == "sse") {
      # Sum of squared error between weights and actual distribution
      sum((weights - d_norm)^2)
    } else if (method == "cor") {
      # Negative correlation (since we want to maximize it)
      -cor(weights, d_norm)
    }
  })
  # don't run
  # print(scores)
  
  # find the best decay value
  decay <- decay_range[which.min(scores)]
  return(decay)
}
# don't run 
# gwd_decay(g_siren, decay_range = seq(0.1, 3, by = 0.1), method = "cor"); gwd_decay(g_siren, decay_range = seq(0.1, 3, by = 0.1), method = "sse")



# compute Bayes' factors -------------------------------------------------------
# https://cran.r-project.org/web/packages/BFpack/vignettes/vignette_BFpack.html
BF <- function(model, priors){ 
  
  # for replication
  set.seed(20240517) # Bugsy's birthday
  
  # don't run
  # this function assumes equal prior probabilities of p = 1/3 = 0.33
  # bf <- BFpack::BF(model)
  
  # function to calculate Bayes' factor
  # https://github.com/jomulder/BFpack
  # https://bfpack.info
  # https://cran.r-project.org/web/packages/BFpack/vignettes/vignette_BFpack.html
  bf <- BFpack::BF(
    model,
    prior.hyp = priors
  )
  
  # hypothesis test
  bf <- bf$PHP_exploratory
  bf <- round(bf, digits = 3) # round posterior probabilities to three decimals
  
  # label coefficients
  row.names(bf) <- c("gwdegree", "gwdsp", "gwesp", "degcor")
  
  # interpretation of the hypothesis tests
  # see Mulder et al., p.48: https://www.sciencedirect.com/science/article/pii/S0378873323000801?via%3Dihub
  cat(message("Bayesian hypothesis test -- posterior probabilities of the hypotheses for each model parameter:"))
  cat("\n")
  print(bf)
  cat("\n")
  cat(message("The columns indicate the probability that the coefficient = 0, < 0, and > 0.")) # e.g., P(B > 0 | Y)
  cat("\n")
  cat(message("Note: Bayesian hypothesis tests not provided for the edges coefficient because an improper flat prior is used for this parameter."))
}



# goodness-of-fit --------------------------------------------------------------
GOF <- function(bayes){
  set.seed(20110211) # Halle's birthday
  n <- 10000 # graph simulations
  i <- 15 # degree distribution range
  j <- 15 # geodisstance range 
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
    gwdegree(decay = gwd_decay(g_siren), fixed = TRUE) + # gwdegree(decay = 3.0, fixed = TRUE) 
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_01.siren)

# hypothesis test
BF(bayes_01.siren, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_01.siren <- GOF(bayes_01.siren)


# estimate the model -----------------------------------------------------------
bayes_02.togo <- bayes(
  y = g_togo,
  x = g_togo ~ edges + 
    gwdegree(decay = gwd_decay(g_togo), fixed = TRUE) + # gwdegree(decay = 0.5, fixed = TRUE)
    gwdsp(decay = 1.0, fixed = TRUE) +
    gwesp(decay = 1.0, fixed = TRUE) +
    degcor
    )
summary(bayes_02.togo)

# hypothesis test
BF(bayes_02.togo, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_02.togo <- GOF(bayes_02.togo)


# estimate the model -----------------------------------------------------------
bayes_03.caviar <- bayes(
  y = g_caviar,
  x = g_caviar ~ edges + 
    gwdegree(decay = 1.5, fixed = TRUE) + # gwdegree(decay = 2.0, fixed = TRUE)
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_03.caviar)

# hypothesis test
BF(bayes_03.caviar, priors = c(1/3, 1/3, 1/3))

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

# hypothesis test
BF(bayes_04.cielnet, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_04.cielnet <- GOF(bayes_04.cielnet)


# estimate the model -----------------------------------------------------------
bayes_05.cocaine <- bayes(
  y = g_cocaine,
  x = g_cocaine ~ edges + 
    gwdegree(decay = 1.2, fixed = TRUE) + # gwdegree(decay = 1.5, fixed = TRUE)
    gwdsp(decay = 0.5, fixed = TRUE) +
    gwesp(decay = 0.5, fixed = TRUE) +
    degcor
    )
summary(bayes_05.cocaine)

# hypothesis test
BF(bayes_05.cocaine, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_05.cocaine <- GOF(bayes_05.cocaine)


# estimate the model -----------------------------------------------------------
bayes_06.heroin <- bayes(
  y = g_heroin,
  x = g_heroin ~ edges + 
    gwdegree(decay = gwd_decay(g_heroin), fixed = TRUE) +  # gwdegree(decay = 2.5, fixed = TRUE)
    gwdsp(decay = 0.3, fixed = TRUE) +
    gwesp(decay = 0.3, fixed = TRUE) +
    degcor
   )
summary(bayes_06.heroin)

# hypothesis test
BF(bayes_06.heroin, priors = c(1/3, 1/3, 1/3))

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

# hypothesis test
BF(bayes_07.oversize, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_07.oversize <- GOF(bayes_07.oversize)


# estimate the model -----------------------------------------------------------
bayes_08.montagna <- bayes(
  y = g_montagna,
  x = g_montagna ~ edges + 
    gwdegree(decay = 2.0, fixed = TRUE) + # gwdegree(decay = 2.0, fixed = TRUE)
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
summary(bayes_08.montagna)

# hypothesis test
BF(bayes_08.montagna, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_08.montagna <- GOF(bayes_08.montagna)



# function to compute the posterior estimates ----------------------------------
posteriors <- function(model, g, name){
  
  # posterior estimates of the model
  theta <- as.data.frame(model$Theta)
  colnames(theta) <- c("edges", "gwdegree", "gwdsp", "gwesp", "degcor")

  # transform the degree correlation term
  theta <- dplyr::mutate(theta, degcor = degcor/100)
  
  # means of the posterior estimates
  mu <- colMeans(theta)
  
  # standard deviation of the posterior estimates
  sd <- matrixStats::colSds(as.matrix(theta))
  
  # naive standard errors
  se1 <- apply(theta, 2, function(x) sd(x) / sqrt(nrow(theta)))
  
  # function to calculate the standard errors because of the auto-correlation in the iterations
  tsseFUN <- function(samples){
    mcmc <- coda::mcmc(samples) # MCMC object for auto-correlation analysis
    ess <- coda::effectiveSize(mcmc) # effective sample size
    var <- stats::var(samples) # variance of the posterior estimates
    stderr <- sqrt(var / ess) # time-series standard error
    return(stderr)
  }
  se2 <- apply(theta, 2, tsseFUN)
  
  # join posterior means and standard errors
  data <- cbind(mu, sd, se1, se2)
  
  # names of the parameters
  rownames(data) <- c("DYADS", "GWDEGREE", "GWDSP", "GWESP", "ASSORTATIVITY")
  parameters <- as.data.frame(rownames(data))
  parameters <- dplyr::mutate(parameters, model = name)
  
  # join
  data <- cbind(parameters, data)
  colnames(data) <- c("parameters", "model", "mean", "sd", "stderr", "stderr.ts")
  rownames(data) <- NULL
  return(data)
}
posteriors_bergm <- rbind(
  posteriors(model = bayes_01.siren, g = g_siren, name = "siren"),
  posteriors(model = bayes_02.togo, g = g_togo, name = "togo"),
  posteriors(model = bayes_03.caviar, g = g_caviar, name = "caviar"),
  posteriors(model = bayes_04.cielnet, g = g_cielnet, name = "cielnet"),
  posteriors(model = bayes_05.cocaine, g = g_cocaine, name = "cocaine"),
  posteriors(model = bayes_06.heroin, g = g_heroin, name = "heroin"),
  posteriors(model = bayes_07.oversize, g = g_oversize, name = "oversize"),
  posteriors(model = bayes_08.montagna, g = g_montagna, name = "montagna")
  )



# network meta-analysis that computes the weighted averages for each of the model parameters i.e., the hierarchical model
metanet <- function(data, parameter){
  
  # required packages
  require("metafor")
  
  # first, subset data by the parameter
  data <- subset(data, parameters == parameter)
  
  # estimate the hierarchical model i.e., random effects model
  model <- metafor::rma.mv(
    yi = mean, # parameter estimate
    V = stderr.ts^2, # squared standard error i.e., level 1 of the model
    random = ~1 | model, # parameter from each model for all the criminal networks i.e., level 2 of the model
    data = data
    )
  return(model)
}


# meta-analysis of each of the parameters
effect <- function(data, parameter){
  
  # required packages
  require("metafor")
  
  # first, subset data by the parameter
  data <- subset(data, parameters == parameter)
  
  # calculate the effect sizes
  outcome <- metafor::escalc(
    # see pg. 89 for options: https://cran.r-project.org/web/packages/metafor/metafor.pdf
    # see pg. 95: use "GEN" for generic outcome measure that is not further specified
    measure = "SMN", # standardized mean
    yi = mean,
    sei = stderr.ts,
    slab = model,
    data = data
    )
  return(outcome)
}
effect(data = posteriors_bergm, parameter = "GWDEGREE")




# network meta-analysis that computes the weighted averages for each of the model parameters i.e., the hierarchical Bayesian model
metanet_bayes <- function(data, parameter){
  
  # set seed for replication
  set.seed(20110210) # Halle's birthday
  
  # required packages
  require("bayesmeta"); require("metafor")
  
  # first, subset data by the parameter
  data <- subset(data, parameters == parameter)
  
  # vector of effect sizes
  mean <- data$mean
  
  # vector of standard errors
  sigma <- if (parameter == "ASSORTATIVITY") {
    data$stderr.ts
    } else {
      data$stderr
  }
  
  # vector of labels 
  model <- data$model
  
  # sample weights
  weights <- 1/(sigma^2)
  
  # prior mean
  prior.mu <- sum(weights * mean)/sum(weights) 
  
  # prior standard deviation
  w_var <- sum(weights * (mean - prior.mu)^2)/(sum(weights) - sum(weights^2)/sum(weights)) # weighted sample variance
  prior.sd <- sqrt(w_var) #  weighted standard deviation is the square root of the weighted variance
  
  # tau prior
  prior.tau <- function(t){bayesmeta::dhalfcauchy(t, scale = 1)}
  # prior.TAU <- function(tau) {
    # ifelse(tau >= 0.001 & tau <= 3, bayesmeta::dhalfnormal(tau, scale = 1), 0)
  # }
  
  # estimate the hierarchical Bayesian model
  model <- bayesmeta::bayesmeta(
    y = mean,
    sigma = sigma,
    mu.prior.mean = prior.mu,
    mu.prior.sd = prior.sd,
    tau.prior = prior.tau,
    label = model
    )
  return(model)
}
meta_01 <- metanet_bayes(data = posteriors_bergm, parameter = "DYADS")
meta_02 <- metanet_bayes(data = posteriors_bergm, parameter = "GWDEGREE")
meta_03 <- metanet_bayes(data = posteriors_bergm, parameter = "GWDSP")
meta_04 <- metanet_bayes(data = posteriors_bergm, parameter = "GWESP")
meta_05 <- metanet_bayes(data = posteriors_bergm, parameter = "ASSORTATIVITY")



# compile results of the meta analysis into one dataset to plot the results
reshape_metabayes <- function(model, effect){
  
  # model parameter
  effect <- effect

  # reconfigure average effects into dataframe
  data <- as.data.frame(model$summary)
  rownames(data) <- NULL
  
  # transpose
  data <- as.data.frame(t(data))
  
  # columns
  colnames(data) <- c("mode", "median", "mean", "sd", "ci.lo", "ci.hi")
  
  # columns
  data <- dplyr::select(data, mean, sd, ci.lo, ci.hi)
  
  # compute quantile range for plot
  Q1 <- data$mean + (stats::qnorm(0.25) * data$sd)
  Q3 <- data$mean + (stats::qnorm(0.75) * data$sd)
  
  # join
  data <- cbind(data, Q1, Q3)
  
  # join
  data <- cbind(effect, data)
  data <- as.data.frame(data[2,])
  rownames(data) <- NULL
  
  # return the dataset for plotting
  return(data)
}
meta_results2 <- rbind(
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "DYADS"), effect = "DYADS"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWDEGREE"), effect = "GWDEGREE"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWDSP"), effect = "GWDSP"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWESP"), effect = "GWESP"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "ASSORTATIVITY"), effect = "ASSORTATIVITY")
  )



# compile results of the meta analysis into one dataset to plot the results
reshape_metafor <- function(model, effect){
  
  # model parameter
  effect <- effect
  
  # multilevel model parameters
  fetch <- function(thing){
    thing <- as.numeric(thing)
    return(thing)
  }
  beta <- fetch(model$beta)
  stderr <- fetch(model$se)
  pvalue <- fetch(model$pval)
  ci.lo <- fetch(model$ci.lb)
  ci.hi <- fetch(model$ci.ub)
  
  # compute quantile range for plot
  Q1 <- beta + (stats::qnorm(0.25) * stderr); Q3 <- beta + (stats::qnorm(0.75) * stderr)
  
  # join
  data <- as.data.frame(cbind(beta, stderr, pvalue, ci.lo, ci.hi, Q1, Q3))
  data <- cbind(effect, data)
  
  # function to label p-values with asterisks
  label <- function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("n.s.") # uniciode characters for n.s.: \u207F\u00B7\u02E2\u00B7
  }
  data$p.signif <- label(data$pvalue)
  
  # return the dataset for plotting
  return(data)
}
meta_results <- rbind(
  reshape_metafor(metanet(data = posteriors_bergm, parameter = "DYADS"), effect = "DYADS"),
  reshape_metafor(metanet(data = posteriors_bergm, parameter = "GWDEGREE"), effect = "GWDEGREE"),
  reshape_metafor(metanet(data = posteriors_bergm, parameter = "GWDSP"), effect = "GWDSP"),
  reshape_metafor(metanet(data = posteriors_bergm, parameter = "GWESP"), effect = "GWESP"),
  reshape_metafor(metanet(data = posteriors_bergm, parameter = "ASSORTATIVITY"), effect = "ASSORTATIVITY")
  )



# boxplot of the weighted average model parameters across the models

# set seed for replication
set.seed(20190816) # Maeve's birthday

# function to set the labels for the y-axis
scaleFUN <- function(x) sprintf("%.2f", x)

# level factors so that the function plots them in the correct order
meta_results2$effect <- factor(meta_results2$effect, levels = c("DYADS", "GWDEGREE", "GWDSP", "GWESP", "ASSORTATIVITY"))

# construct the boxplot
fig_meta <- ggplot2::ggplot(meta_results2, aes(x = effect, y = mean)) +
  # horizontal line to marker where y = 0
  ggplot2::geom_hline(yintercept = 0.00, color = "red", linetype = "dashed", size = 1, alpha = 0.8) +
  # error bars to visualize the confidence intervals
  ggplot2::geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black") +
  # data points
  ggplot2::geom_jitter(
    data = posteriors_bergm, 
    mapping = ggplot2::aes(x = parameters, y = mean), 
    width = 0.2, 
    size = 2, 
    shape = 21, # circles
    color = "black",
    fill = "white",
    alpha = 1.0
    ) +
  # boxplot
  ggplot2::geom_boxplot(
    ggplot2::aes(lower = Q1, upper = Q3, middle = mean, ymin = ci.lo, ymax = ci.hi), 
    stat = "identity",
    fill = "lightgrey", color = "black", width = 0.5, alpha = 1.0
    ) +
  # don't run
  # circle points for the mean
  # ggplot2::geom_point(size = 3, color = "black") +
  ggplot2::scale_y_continuous(limits = c(-6, 2), breaks = seq(-6, 2, 1), labels = scaleFUN) +
  ggplot2::labs(
    # title = "NETWORK META-ANALYSIS",
    x = "MODEL PARAMETERS",
    y = "WEIGHTED POSTERIOR ESTIMATES"
    ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 8),
    axis.title.x = ggplot2::element_text(size = 10, face = "plain"),
    axis.text.y = ggplot2::element_text(angle =  0, vjust = 0.5, hjust = 0.5, size = 8),
    axis.title.y = ggplot2::element_text(size = 10, face = "plain")
    )

# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop", 
    width = 5, 
    height = 5, 
    device = 'png', 
    dpi = 700
    )
  }
output(plot = fig_meta, filename = "fig_meta.png")






# construct the boxplot
fig_meta <- ggplot2::ggplot(meta_results2, aes(x = effect, y = mean)) +
  # horizontal line to marker where y = 0
  ggplot2::geom_hline(yintercept = 0.00, color = "red", linetype = "dashed", size = 1, alpha = 0.8) +
  # error bars to visualize the confidence intervals
  ggplot2::geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black") +
  # data points
  ggplot2::geom_jitter(
    data = posteriors_bergm, 
    mapping = ggplot2::aes(x = parameters, y = mean), 
    width = 0.2, 
    size = 2, 
    shape = 21, # circles
    color = "black",
    fill = "white",
    alpha = 1.0
  ) +
  # boxplot
  ggplot2::geom_boxplot(
    ggplot2::aes(lower = Q1, upper = Q3, middle = mean, ymin = ci.lo, ymax = ci.hi), 
    stat = "identity",
    fill = "lightgrey", color = "black", width = 0.5, alpha = 1.0
  ) +
  # don't run
  # circle points for the mean
  # ggplot2::geom_point(size = 3, color = "black") +
  ggplot2::scale_y_continuous(limits = c(-6, 2), breaks = seq(-6, 2, 1), labels = scaleFUN) +
  ggplot2::labs(
    # title = "NETWORK META-ANALYSIS",
    x = "MODEL PARAMETERS",
    y = "WEIGHTED POSTERIOR ESTIMATES"
  ) +
  ggthemes::theme_base() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 8),
    axis.title.x = ggplot2::element_text(size = 10, face = "plain"),
    axis.text.y = ggplot2::element_text(angle =  0, vjust = 0.5, hjust = 0.5, size = 8),
    axis.title.y = ggplot2::element_text(size = 10, face = "plain")
  ) +
  # markers that indicate statistical significance
  ggplot2::geom_text(
    data = meta_results, 
    ggplot2::aes(x = effect, y = 1.9, label = p.signif), 
    size = 3, 
    color = "black", 
    vjust = 0
  )



# boxplot of the average posterior estimates across the models
ggplot2::ggplot(posteriors_bergm, ggplot2::aes(x = parameters, y = mean, fill = "white")) +
  ggplot2::geom_boxplot(color = "black", alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.2, size = 2, color = "black", fill = "white", alpha = 0.9) +
  ggplot2::scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1)) +
  ggplot2::labs(
    # title = "TITLE",
    x = NULL,
    y = "POSTERIOR ESTIMATES"
    ) +
  ggthemes::theme_base() +
  ggplot2::theme(legend.position = "none")





# estimate the model -----------------------------------------------------------
g_super <- intergraph::asNetwork(igraph::simplify(intergraph::asIgraph(g_super)))
bayes_09.super <- bayes(
  y = g_super,
  x = g_super ~ edges + 
    gwdegree(decay = 1.0, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor + 
    nodematch("group")
    )
summary(bayes_09.super)

# hypothesis test
BF(bayes_09.super, priors = c(1/3, 1/3, 1/3))

# goodness-of-fit
gof_09.super <- GOF(bayes_09.super)


# estimate the model -----------------------------------------------------------
# note: can't estimate degree assortativity i.e., 'degcor' in the subgraph estimates
bayes_10.super <- bayes(
  y = g_super,
  x = g_super ~ edges + 
    gwdegree(decay = 0.1, fixed = TRUE) +
    gwdsp(decay = 0.1, fixed = TRUE) +
    gwesp(decay = 1.0, fixed = TRUE) +
    degcor +
  # subgraph estimates for different 'groups' of the super network -------------
  # siren auto theft ring
  S(~edges + 
      gwdegree(decay = 3.0, fixed = TRUE) + 
      gwdsp(decay = 3.0, fixed = TRUE)    + 
      gwesp(decay = 3.0, fixed = TRUE)    ,
    ~(group == c("siren")
    )
  ) +
  # togo auto theft ring
  S(~edges + 
      gwdegree(decay = 0.5, fixed = TRUE) + 
      gwdsp(decay = 3.0, fixed = TRUE)    + 
      gwesp(decay = 3.0, fixed = TRUE)    ,
    ~(group == c("togo")
    )
  # caviar drug trafficking organization
  ) +
  S(~edges + 
      gwdegree(decay = 2.0, fixed = TRUE) + 
      gwdsp(decay = 1.0, fixed = TRUE)    + 
      gwesp(decay = 2.0, fixed = TRUE)    ,
    ~(group == c("caviar")
    )
  ) +
  # cielnet drug trafficking organization
  S(~edges + 
      gwdegree(decay = 1.0, fixed = TRUE) + 
      gwdsp(decay = 0.2, fixed = TRUE)    + 
      gwesp(decay = 0.2, fixed = TRUE)    ,
    ~(group == c("cielnet")
    )
  ) +
  # La Cosa Nostra cocaine trafficking outfit
  S(~edges + 
      gwdegree(decay = 1.5, fixed = TRUE) + 
      gwdsp(decay = 0.5, fixed = TRUE)    + 
      gwesp(decay = 0.5, fixed = TRUE)    ,
    ~(group == c("cocaine")
     )
    ) +
  # New York City heroin traffickers
  S(~edges + 
      gwdegree(decay = 2.5, fixed = TRUE) + 
      gwdsp(decay = 1.5, fixed = TRUE)    + 
      gwesp(decay = 0.5, fixed = TRUE)    ,
    ~(group == c("heroin")
      )
    ) +
  # 'Ndrangheta wiretap records -- operation oversize
  S(~edges + 
      gwdegree(decay = 1.5, fixed = TRUE) + 
      gwdsp(decay = 1.5, fixed = TRUE)    + 
      gwesp(decay = 0.5, fixed = TRUE)    ,
    ~(group == c("oversize")
      )
    ) +
  # Cosa Nostra wiretap records -- operation montagna 
  S(~edges + 
      gwdegree(decay = 2.0, fixed = TRUE) + 
      gwdsp(decay = 2.0, fixed = TRUE)    + 
      gwesp(decay = 2.0, fixed = TRUE)    ,
    ~(group == c("montagna")
    )
  ) +
  # the french connection - federal bureau of narcotics
  S(~edges + 
      gwdegree(decay = 2.0, fixed = TRUE) + 
      gwdsp(decay = 2.0, fixed = TRUE)    + 
      gwesp(decay = 2.0, fixed = TRUE)    ,
    ~(group == c("tfc")
    )
  )
)
summary(bayes_10.super)

# goodness-of-fit
gof_10.super <- GOF(bayes_10.super)



# construct results table from the posterior estimates
thetas <- bayes_10.super$Theta # posterior estimates
colnames(thetas) <- c( # parameter labels
  "edges",
  "gwdegree 0.1",
  "gwdsp 0.1",
  "gwesp 1.0",
  "degcor",
  "edges - siren",
  "gwdegree 3.0 - siren",
  "gwdsp 3.0 - siren",
  "gwesp 3.0 - siren",
  "edges - togo",
  "gwdegree 0.5 - togo",
  "gwdsp 3.0 - togo",
  "gwesp 3.0 - togo",
  "edges - caviar",
  "gwdegree 2.0 - caviar",
  "gwdsp 1.0 - caviar",
  "gwesp 2.0 - caviar",
  "edges - cielnet",
  "gwdegree 1.0 - cielnet",
  "gwdsp 0.2 - cielnet",
  "gwesp 0.2 - cielnet",
  "edges - cocaine",
  "gwdegree 1.5 - cocaine",
  "gwdsp 0.5 - cocaine",
  "gwesp 0.5 - cocaine",
  "edges - heroin",
  "gwdegree 2.5 - heroin",
  "gwdsp 1.5 - heroin",
  "gwesp 0.5 - heroin",
  "edges - oversize",
  "gwdegree 1.5 - oversize",
  "gwdsp 1.5 - oversize",
  "gwesp 0.5 - oversize",
  "edges - montagna",
  "gwdegree 2.0 - montagna",
  "gwdsp 2.0 - montagna",
  "gwesp 2.0 - montagna"
  )

# compute means and standard deviations of the postieor estimates

# don't run
# mean(thetas[, 1]) # estimate for the edges parameter
# sd(thetas[, 1])

# mean
thetas_mu <- colMeans(thetas)

# standard deviation
thetas_sd <- matrixStats::colSds(thetas) # install.packages("matrixStats")

# 
thetas_table <- cbind(thetas_mu, thetas_sd)
thetas_table <- as.data.frame(thetas_table)
thetas_names <- rownames(thetas_table)
thetas_table <- cbind(thetas_names, thetas_table)
rownames(thetas_table) <- c()
colnames(thetas_table) <- c("parameter", "mean", "SD")

# close .r script








  
  
  
  

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
model_02.caviar <- ERGM(
  g = g_caviar,
  x = g_caviar ~ 
    edges + 
    gwdegree(decay = 2.0, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )
# pass model specifciations through function to estimate models 
model_09.tfc <- ERGM(
  g = g_tfc,
  x = g_tfc ~ 
    edges + 
    gwdegree(decay = 2.0, fixed = TRUE) +
    gwdsp(decay = 2.0, fixed = TRUE) +
    gwesp(decay = 2.0, fixed = TRUE) +
    degcor
    )

# spectral goodness-of-fit measure to critique the model fit to the observed data
# installation guidelines
# https://github.com/blubin/SpectralGOF/blob/master/instructions/spectralGOFwalkthrough.pdf
# install.packages("~/Desktop/spectralGOF", repos = NULL, type = "source")
spectralGOF::SGOF(model_09.tfc, nsim = 100)



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


# hierarchical relations by subgraph -----------------------------------------
S(~nodemix('rank_c'), ~(family_c == c("1_bonanno" )))     +
  # S(~nodemix('rank_c'), ~(family_c == c("2_decavalcante"))) +
  S(~nodemix('rank_c'), ~(family_c == c("3_gambino" )))     +
  S(~nodemix('rank_c'), ~(family_c == c("4_genovese")))     +
  S(~nodemix('rank_c'), ~(family_c == c("5_lucchese")))     +
  S(~nodemix('rank_c'), ~(family_c == c("6_profaci" )))     +
  
  
  
  
  
  # function to average the posterior estimates
  mod_avg <- function(posteriors){
    
    # required packages
    require(c("dplyr", "tidyr", "magrittr"))
    `%>%` <- magrittr::`%>%`
    
    # calculate weights from the inverse squared variance
    posteriors$weight <- 1/(posteriors$stderr.ts^2)
    
    # calculate means of the posterior estimates across models
    means <- posteriors %>% 
      dplyr::group_by(parameters) %>%
      dplyr::summarise(
        w_mean = sum(weight * mean)/sum(weight),
        w_var = 1/sum(weight),
        w_stderr = sqrt(w_var)
      )
    print(means)
    
    # names of the parameters
    parameters <- rownames(posteriors)
    parameters <- c("DYADS", "GWDEGREE", "GWDSP", "GWESP", " ASSORTATIVITY")
    
    # join means to dataset
    posteriors <- cbind(posteriors, means)
    
    # join model parameters to dataset
    posteriors <- cbind(parameters, posteriors)
    
    # name columns
    colnames(posteriors) <- c("parameters", "siren", "togo", "caviar", "cielnet", "cocaine", "heroin", "oversize", "montagna", "means")
    
    # reshape wide to long
    posteriors <- as.data.frame(posteriors)
    posteriors <- tidyr::pivot_longer(
      posteriors,
      cols = -parameters,
      names_to = "model",
      values_to = "posterior"
    )
    posteriors <- dplyr::filter(posteriors, model != "means")
    posteriors <- dplyr::mutate(posteriors, posterior = as.numeric(posterior))
    # return the posteriors
    return(posteriors)
  }
posteriors <- mod_avg(posteriors_bergm)

  