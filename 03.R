#  -----------------------------------------------------------------------------------

# file 03: estimate the Bayesian hierarchical model

# last updated: 23/04/2025

# ------------------------------------------------------------------------------------



# function to compute the posterior estimates
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
  
  # return
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
  # Gelman recommends the half-Caucy priors (0, 5), where scale parameter = 5
  # https://projecteuclid.org/journals/bayesian-analysis/volume-1/issue-3/Prior-distributions-for-variance-parameters-in-hierarchical-models-comment-on/10.1214/06-BA117A.full
  prior.tau <- function(t){bayesmeta::dhalfcauchy(t, scale = 1)}
  # don't run
  # let the tau prior vary within defined range
  # prior.TAU <- function(tau) {
  # ifelse(tau >= 0.001 & tau <= 3, bayesmeta::dhalfnormal(tau, scale = 1), 0)
  # }
  
  # estimate the hierarchical Bayesian model
  # sources: 
  # https://cran.r-project.org/web/packages/bayesmeta/vignettes/bayesmeta.html
  # https://cran.r-project.org/web/packages/bayesmeta/bayesmeta.pdf
  # https://rshiny.gwdg.de/apps/bayesmeta/
  model <- bayesmeta::bayesmeta(
    y = mean,
    sigma = sigma,
    mu.prior.mean = prior.mu,
    mu.prior.sd = prior.sd,
    tau.prior = prior.tau,
    label = model
    )
  
  # return
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
meta_results <- rbind(
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "DYADS"), effect = "DYADS"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWDEGREE"), effect = "GWDEGREE"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWDSP"), effect = "GWDSP"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "GWESP"), effect = "GWESP"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "ASSORTATIVITY"), effect = "ASSORTATIVITY")
  )



# boxplot of the weighted average model parameters across the models

  # set seed for replication
  set.seed(20190816) # Maeve's birthday
  
  # function to set the labels for the y-axis
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # level factors so that the function plots them in the correct order
  meta_results$effect <- factor(meta_results$effect, levels = c("DYADS", "GWDEGREE", "GWDSP", "GWESP", "ASSORTATIVITY"))

# construct the boxplot
fig_meta <- ggplot2::ggplot(meta_results, ggplot2::aes(x = effect, y = mean)) +
  # horizontal line to marker where y = 0
  ggplot2::geom_hline(yintercept = 0.00, color = "red", linetype = "dashed", size = 1.0, alpha = 0.8) +
  # error bars to visualize the confidence intervals
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black") +
  # data points
  ggplot2::geom_jitter(
    data = posteriors_bergm, 
    mapping = ggplot2::aes(x = parameters, y = mean), 
    width = 0.2, 
    size = 2.0, 
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
  ggthemes::theme_few() +
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
    path = "~/Desktop/Super", 
    width = 5, 
    height = 5, 
    device = 'pdf', 
    dpi = 700
    )
}
output(plot = fig_meta, filename = "fig2.pdf")



# simulate networks from the effect sizes of the hierarchical model

  # this code simulates 120,000 networks in total
  # 10,000 simulations x four networks of different sizes x three densities

  # posterior distribution for the effect sizes from the hierarchical model
  posterior_effects <- data.frame(
    dyads = meta_01$rposterior(n = 10000)[, "mu"],
    gwdegree = meta_02$rposterior(n = 10000)[, "mu"],
    gwdsp = meta_03$rposterior(n = 10000)[, "mu"],
    gwesp = meta_04$rposterior(n = 10000)[, "mu"],
    degcor = meta_05$rposterior(n = 10000)[, "mu"]
    )
  
  # priors
  prior <- brms::prior(normal(0, 2), class = "b") + brms::prior(normal(0, 5), class = "Intercept")
  posterior_scaled <- as.data.frame(scale(posterior_effects[, c("dyads", "gwdegree", "gwdsp", "gwesp", "degcor")]))
  
  # regression model
  fit <- brms::brm(
    formula = dyads ~ gwdegree + gwdsp + gwesp + degcor,
    data = posterior_scaled,
    prior = prior,
    chains = 4,
    cores = 4,
    iter = 5000,
    seed = 20240517 # Malone's birthday
    )
  brms::bayes_R2(fit)
  plot(fit)
  
  summary(brms::VarCorr(fit, summary = TRUE))

  
# join data
thetas <- rbind(
  cbind(as.data.frame(bayes_01.siren$Theta), network = "siren"),
  cbind(as.data.frame(bayes_02.togo$Theta), network = "togo"),
  cbind(as.data.frame(bayes_03.caviar$Theta), network = "caviar"),
  cbind(as.data.frame(bayes_04.cielnet$Theta), network = "cielnet"),
  cbind(as.data.frame(bayes_05.cocaine$Theta), network = "cocaine"),
  cbind(as.data.frame(bayes_06.heroin$Theta), network = "heroin"),
  cbind(as.data.frame(bayes_07.oversize$Theta), network = "oversize"),
  cbind(as.data.frame(bayes_08.montagna$Theta), network = "montagna")
  )  
colnames(thetas) <- c("dyads", "gwdegree", "gwdsp", "gwesp", "degcor", "model")   



# bayes model
correlate <- function(thetas){
  
  # posterior estimates
  thetas <- dplyr::mutate(thetas, degcor = degcor/100)
  
  # model names
  model <- thetas$model
  
  # scale
  thetas <- as.data.frame(scale(thetas[, c("dyads", "gwdegree", "gwdsp", "gwesp", "degcor")]))
  
  # join
  thetas <- cbind(thetas, model = model)
  
  # regression model
  fit <- brms::brm(
    formula = gwdegree ~ gwdsp + gwesp + degcor + (1 | model),
    data = thetas,
    # prior = prior,
    chains = 4,
    cores = 4,
    iter = 5000,
    seed = 20240517 # Malone's birthday
    )
  print(summary(fit))
  
  # posterior fitted values (posterior means)
  thetas$fitted <- fitted(fit)[, "Estimate"]
  
  # Bayesian R²
  r2 <- brms::bayes_R2(fit)[1]
  
  # decompose the variance across model components
  vc <- brms::VarCorr(fit, summary = TRUE)
  print(vc)  
  
  # plot actual vs fitted with line of best fit
  best <- ggplot2::ggplot(thetas, aes(x = fitted, y = gwdegree)) +
    ggplot2::geom_point(
      shape = 21, 
      size = 1, 
      color = "black", 
      fill = "white", 
      alpha = 0.5
      ) +
    ggplot2::scale_y_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "firebrick1", alpha = 0.80) +
    ggplot2::labs(
      title = "PREDICTED LINE OF BEST FIT",
      subtitle = paste0("Bayesian R² = ", round(r2, digits = 2)),
      x = "PREDICTED GWDEGREE",
      y = "GWDEGREE STATISTICS FROM THE MODEL"
      ) +
    ggthemes::theme_few()
  return(best)
}
correlate(thetas)






# bayes model
correlate <- function(theta, title){
  
  # posterior estimates
  theta <- as.data.frame(theta)
  colnames(theta) <- c("dyads", "gwdegree", "gwdsp", "gwesp", "degcor")
  theta <- dplyr::mutate(theta, degcor = degcor/100)
  
  # scale
  theta <- as.data.frame(scale(theta[, c("dyads", "gwdegree", "gwdsp", "gwesp", "degcor")]))

  # regression model
  fit <- brms::brm(
    formula = gwdegree ~ gwdsp + gwesp + degcor,
    data = theta,
    prior = prior,
    chains = 4,
    cores = 4,
    iter = 5000,
    seed = 20240517 # Malone's birthday
    )
    print(summary(fit))
    
    # posterior fitted values (posterior means)
    theta$fitted <- fitted(fit)[, "Estimate"]
    theta$lower <- fitted(fit)[, "Q2.5"]
    theta$upper <- fitted(fit)[, "Q97.5"]
    
    # Bayesian R²
    r2 <- brms::bayes_R2(fit)[1]
    
    # variance decomposition
    summary(brms::VarCorr(fit, summary = TRUE))
    
    # plot actual vs fitted with line of best fit
    best <- ggplot2::ggplot(theta, aes(x = fitted, y = gwdegree)) +
      ggplot2::geom_point(
        shape = 21, 
        size = 1, 
        color = "black", 
        fill = "white", 
        alpha = 0.5
        ) +
      ggplot2::geom_line(ggplot2::aes(y = fitted), color = "firebrick1", alpha = 1.00, linewidth = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "firebrick1", alpha = 0.50) +
      ggplot2::scale_y_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
      ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
      ggplot2::annotate(
        "text", x = -2.75, y = 2.75, 
        label = bquote(italic(R)^2 == .(round(r2, digits = 2))), 
        hjust = 1, 
        size = 3
        ) +
      ggplot2::labs(
        title = title,
        x = "FITTED VALUES",
        y = "GWDEGREE STATISTICS"
        ) +
      ggthemes::theme_few()
    return(best)
}
correlate(bayes_01.siren$Theta, title = "(A) SIREN AUTO THEFT NETWORK")
correlate(bayes_02.togo$Theta)
correlate(bayes_03.caviar$Theta)
correlate(bayes_04.cielnet$Theta)
correlate(bayes_05.cocaine$Theta)
correlate(bayes_06.heroin$Theta)
correlate(bayes_07.oversize$Theta)
correlate(bayes_08.montagna$Theta)







  # number of nodes and density of networks to simulate 
  n <- c(30, 50, 75, 100)
  p <- c(0.02, 0.05, 0.10)
  
  # dyad counts 
  dyads <- expand.grid(n = n, p = p)
  dyads$dyads <- round(dyads$n * (dyads$n - 1) / 2 * dyads$p)
  
  # simulate
  simulator <- function(nodes, theta){
    
    # initialize network of n number of nodes
    g <- network::network.initialize(nodes, directed = FALSE)
    
    # formula
    f <- g ~ edges + gwesp(0.5, fixed = TRUE) + gwdsp(0.5, fixed = TRUE) + gwdegree(0.5, fixed = TRUE) + degcor()
    
    # simulate
    sim <- tryCatch({
      # wrapper function for igraph local average transitivity
      stats::simulate(
        f, # model specification
        coef = unname(as.numeric(theta)), # posterior estimates
        seed = 20240517, # Malone's birthday
        nsim = 1, # number of simulations
        output = "network" # thing to simulate
        )
    }, 
    error = function(e) {
      message("simulation failed.")
      return(NULL)
    })
    
    # return
    return(sim)
  }
  
  results <- list()
  counter <- 1
  
  for (i in 1:nrow(posterior_effects)) {
    theta <- posterior_effects[i, ]
    
    for (j in 1:nrow(dyads)) {
      
      nodes <- dyads$n[j]
      density <- dyads$p[j]
      links <- dyads$dyads[j]
      
      g_sim <- simulator(nodes, theta)
      
      results[[counter]] <- list(
        draw = i,
        nodes = nodes,
        density = density,
        edges = links,
        net = g_sim
        )
      counter <- counter + 1
    }
  }
  
  
  
  





# close .r script





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



# construct the boxplot
fig_meta <- ggplot2::ggplot(meta_results2, ggplot2::aes(x = effect, y = mean)) +
  # horizontal line to marker where y = 0
  ggplot2::geom_hline(yintercept = 0.00, color = "red", linetype = "dashed", size = 1, alpha = 0.8) +
  # error bars to visualize the confidence intervals
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci.lo, ymax = ci.hi), width = 0.2, color = "black") +
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



