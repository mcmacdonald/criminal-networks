#  -----------------------------------------------------------------------------------

# file 03: estimate the Bayesian hierarchical model

# ------------------------------------------------------------------------------------



# function to compute the posterior estimates
posteriors <- function(model, g, name){
  
  # posterior estimates of the model
  theta <- as.data.frame(model$Theta)
  colnames(theta) <- c("edges", "gwdegree", "gwdsp", "gwesp", "degcor")
  
  # rescale the degree correlation term
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
  rownames(data) <- c("EDGES", "CENTRALIZATION", "COMMON FRIENDS", "TRIADIC CLOSURE", "ASSORTATIVITY")
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
  
  # vector of standard deviations
  sigma <- data$sd
  
  # vector of labels 
  model <- data$model
  
  # tau prior
  # see Gelman's paper
  # https://projecteuclid.org/journals/bayesian-analysis/volume-1/issue-3/Prior-distributions-for-variance-parameters-in-hierarchical-models-comment-on/10.1214/06-BA117A.full
  prior.tau <- function(t){bayesmeta::dhalfcauchy(t, scale = 5)}
  
  # estimate the hierarchical Bayesian model
  # sources: 
  # https://cran.r-project.org/web/packages/bayesmeta/vignettes/bayesmeta.html
  # https://cran.r-project.org/web/packages/bayesmeta/bayesmeta.pdf
  # https://rshiny.gwdg.de/apps/bayesmeta/
  model <- bayesmeta::bayesmeta(
    y = mean,
    sigma = sigma,
    tau.prior = prior.tau,
    label = model
    )
  
  # return
  return(model)
}
meta_01 <- metanet_bayes(data = posteriors_bergm, parameter = "EDGES")
meta_02 <- metanet_bayes(data = posteriors_bergm, parameter = "COMMON FRIENDS")
meta_03 <- metanet_bayes(data = posteriors_bergm, parameter = "TRIADIC CLOSURE")
meta_04 <- metanet_bayes(data = posteriors_bergm, parameter = "ASSORTATIVITY")
meta_05 <- metanet_bayes(data = posteriors_bergm, parameter = "CENTRALIZATION")



# compile results of the meta analysis into one dataset to plot the results
reshape_metabayes <- function(model, effect){
  
  # required packages
  require('dplyr'); require('stats')
  
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
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "EDGES"), effect = "EDGES"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "COMMON FRIENDS"), effect = "COMMON FRIENDS"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "TRIADIC CLOSURE"), effect = "TRIADIC CLOSURE"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "ASSORTATIVITY"), effect = "ASSORTATIVITY"),
  reshape_metabayes(model = metanet_bayes(data = posteriors_bergm, parameter = "CENTRALIZATION"), effect = "CENTRALIZATION")
  )



# boxplot of the weighted average model parameters across the models

  # set seed for replication
  set.seed(20190816) # Maeve's birthday
  
  # function to set the labels for the y-axis
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  # level factors so that the function plots them in the correct order
  meta_results$effect <- factor(
    meta_results$effect, 
    levels = c("EDGES", "COMMON FRIENDS", "TRIADIC CLOSURE", "ASSORTATIVITY", "CENTRALIZATION"
              )
    )
  meta_results <- dplyr::filter(meta_results, effect != "EDGES")
  posteriors_bergm <- dplyr::filter(posteriors_bergm, parameters != "EDGES")

# construct the boxplot
fig2 <- ggplot2::ggplot(meta_results, ggplot2::aes(x = effect, y = mean)) +
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
  ggplot2::scale_y_continuous(limits = c(-8, 2.00), breaks = seq(-8, 2.00, 2), labels = scaleFUN) +
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
    device = 'png', 
    dpi = 700
    )
}
output(plot = fig2, filename = "fig2.png")



# close .R script

