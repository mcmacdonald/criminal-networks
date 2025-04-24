#  -----------------------------------------------------------------------------------

# file 03: estimate the shape of the degree distributions of each graph

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------



# plot complimentary cumulative degree distribution ----------------------------------
cdf_degree <- function(g, cmode){
  par(mfrow=c(1,1)) # plot dimensions
  # normalized degree centrality
  d <- sna::degree(g, gmode = "graph", cmode = "freeman") # degree centrality for each node
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



# first, estimate the shape of the degree distribution and test it against different statistical distributions
# start from the assumption that the degree is log-normal... Broido and Clauset 2019 find that most social networks produce log-normal degree distributions
# 1) compare to the power law distribution, another heavy tailed distribution
# 2) compare to the exponential distribution, a statistical distribution that does not have a heavy tail
# 3) compare to the Poisson distribution, a statistical distribution that does not have a heavy tail



# Vuong's likelihood ratio test for goodness-of-fit for the log normal distribution ----------------------------------------------------------------
vuong = function(g, model){
  
  # required packages
  require("poweRlaw"); require("sna")
  
  # interpretation of the goodness-of-fit test
  cat("\n") # space
  message("Vuong's likelihood ratio test is a sign test:")
  message("* A likelihood-ratio test statistic > 0 suggests the degree distribution more closely resembles the shape of the log normal distribution.")
  message("* A likelihood-ratio test statistic < 0 suggests the degree distribution more closely resembles the shape of the comparison distribution.")
  message("* A larger test statistic, in either dirtection, suggests better goodness-of-fit.")
  cat("\n") # space
  message("Hypothesis statements for Vuong's likelihood ratio test:")
  message("    H0: The degree distribution does not resemble the shape of the log normal or comparison distribution.")
  message("    H1: The degree distribution resembles the shape of the log normal or comparison distribution.")
  cat("\n") # space
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[d!=0] # drop nodes that do not have sender or receiver ties
  d <- d[order(d, decreasing = FALSE)]
  
  # fit log-normal distribution
  y = poweRlaw::dislnorm(d)
  y$setXmin(1) # for all degree k
  y$setPars(poweRlaw::estimate_pars(y))
  
  # fit comparison distribution
  x = model(d)
  x$setXmin(1) # for all degree k
  x$setPars(poweRlaw::estimate_pars(x))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(y, x)
  # notes on the interpretation of Vuong's likelihood-ratio test statistic: 
  # the sign of the test statistic (i.e., +/-) has meaning for interpretation (Vuong's formula is a sign test)
  # because 'y' is the first input into the poweRlaw::compare_distributions() function and 'x' is the second:
  # ... a positive (+) test statistic suggests the degree distribution more so resembles the log normal distribution
  # ... a negative (-) test statistic suggests the degree distribution more so resembles the comparison distribution
  # ... reversing the input order in poweRlaw::compare_distributions() (i.e., 'x' before 'y') computes the same test statistic, but in the opposite direction
  
  # test results
  message("Results of Vuong's likelihood ratio test:")
  cat("\n") # space
  
  # Vuong's likelihood-ratio test statistic 
  lrstat = lrtest$test_statistic
  cat("Vuong's likelihood-ratio test statistic = "); cat(lrstat ); cat("\n"); cat("\n")
  
  # p-value (one-tailed)
  p = lrtest$p_one_sided # one sided p-values because Broido & Clauset (2019) find that the degree distributions of most social networks resemble log-normal distributions
  # if p < 0.05, reject H0: the degree distribution neither resembles the power law distribution or Poisson distribution
  # if p > 0.05, fail to reject H1: degree distribution resemebles power law distribution or Poisson distribution (see notes on Vuong's likelihood-ratio test)
  cat( "p-value = "); cat(p); cat("\n"); cat("\n")
  
  # interpretation:
  if(lrstat > 0){
    message("The degree distribution more closely resembles the log normal distribution.")
    if(p < 0.10){
      cat("\n")
      message("Reject the null hypothesis that the degree distribution does not resemble the log normal distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the log normal distribution."}
  } else {
    message("The degree distribution more closely resembles the comparison distribution.")
    if(p < 0.10){
      cat("\n") 
      message("Reject the null hypothesis that the degree distribution does not resemble the comparison distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the comparison distribution."}
  }
  
  # return results 
  p <- as.data.frame(p)
  colnames(p) <- c("p")
  return(p)
}
# compare log normal distribution to the power law distribution
vuong_siren <- vuong(g_siren, model = poweRlaw::displ)
vuong_togo <- vuong(g_togo, model = poweRlaw::displ)
vuong_caviar <- vuong(g_caviar, model = poweRlaw::displ)
vuong_cielnet <- vuong(g_cielnet, model = poweRlaw::displ)
vuong_cocaine <- vuong(g_cocaine, model = poweRlaw::displ)
vuong_heroin <- vuong(g_heroin, model = poweRlaw::displ)
vuong_oversize <- vuong(g_oversize, model = poweRlaw::displ)
vuong_montagna <- vuong(g_montagna, model = poweRlaw::displ)

# compare log normal distribution to the exponential distribution
vuong(g_siren, model = poweRlaw::disexp)
vuong(g_togo, model = poweRlaw::disexp)
vuong(g_caviar, model = poweRlaw::disexp)
vuong(g_cielnet, model = poweRlaw::disexp)
vuong(g_cocaine, model = poweRlaw::disexp)
vuong(g_heroin, model = poweRlaw::disexp)
vuong(g_oversize, model = poweRlaw::disexp)
vuong(g_montagna, model = poweRlaw::disexp)

# compare log normal distribution to the Poisson distribution
vuong(g_siren, model = poweRlaw::dispois)
vuong(g_togo, model = poweRlaw::dispois)
vuong(g_caviar, model = poweRlaw::dispois)
vuong(g_cielnet, model = poweRlaw::dispois)
vuong(g_cocaine, model = poweRlaw::dispois)
vuong(g_heroin, model = poweRlaw::dispois)
vuong(g_oversize, model = poweRlaw::dispois)
vuong(g_montagna, model = poweRlaw::dispois)



# Vuong's likelihood ratio test for goodness-of-fit for the log normal distribution ----------------------------------------------------------------
vuong_simulator = function(simulations, model1, model2){
  
  # required packages
  require("poweRlaw"); require("sna")
  
  # number of simulations
  samples <- length(simulations)
  
  # results
  results <- c()
  
  for (i in 1:samples){
  
  # degree distribution for the criminal networks
  d <- sna::degree(simulations[[i]], gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[d!=0] # drop nodes that do not have sender or receiver ties
  d <- d[order(d, decreasing = FALSE)]
  
  # fit distribution
  y = model1(d)
  y$setXmin(min(d)) # for all degree k
  y$setPars(poweRlaw::estimate_pars(y))
  
  # fit comparison distribution
  x = model2(d)
  x$setXmin(min(d)) # for all degree k
  x$setPars(poweRlaw::estimate_pars(x))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(y, x)
  
  # Vuong's likelihood-ratio test statistic 
  lrstat = lrtest$test_statistic
 
  # p-value (one-tailed)
  p = lrtest$p_one_sided
  
  # join p-values to results vector
  results <- c(results, p)
  }
  
  # return results
  results <- as.data.frame(results)
  colnames(results) <- c("p")
  return(results)
}
# compare log normal distribution to the power law distribution
vuong_siren.sim <- vuong_simulator(simulations = g_siren.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_togo.sim <- vuong_simulator(simulations = g_togo.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_caviar.sim <- vuong_simulator(simulations = g_caviar.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_cielnet.sim <- vuong_simulator(simulations = g_cielnet.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_cocaine.sim <- vuong_simulator(simulations = g_cocaine.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_heroin.sim <- vuong_simulator(simulations = g_heroin.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_oversize.sim <- vuong_simulator(simulations = g_oversize.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)
vuong_montagna.sim <- vuong_simulator(simulations = g_montagna.sim, model1 = poweRlaw::dislnorm, model2 = poweRlaw::displ)

# compare the power law distribution to the log normal distribution
vuong_siren.sim <- vuong_simulator(simulations = g_siren.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_togo.sim <- vuong_simulator(simulations = g_togo.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_caviar.sim <- vuong_simulator(simulations = g_caviar.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_cielnet.sim <- vuong_simulator(simulations = g_cielnet.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_cocaine.sim <- vuong_simulator(simulations = g_cocaine.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_heroin.sim <- vuong_simulator(simulations = g_heroin.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_oversize.sim <- vuong_simulator(simulations = g_oversize.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)
vuong_montagna.sim <- vuong_simulator(simulations = g_montagna.sim, model1 = poweRlaw::displ, model2 = poweRlaw::dislnorm)




# compare log normal distribution to the exponential distribution
vuong_siren.sim <- vuong_simulator(simulations = g_siren.sim, model = poweRlaw::disexp)
vuong_togo.sim <- vuong_simulator(simulations = g_togo.sim, model = poweRlaw::disexp)
vuong_caviar.sim <- vuong_simulator(simulations = g_caviar.sim, model = poweRlaw::disexp)
vuong_cielnet.sim <- vuong_simulator(simulations = g_cielnet.sim, model = poweRlaw::disexp)
vuong_cocaine.sim <- vuong_simulator(simulations = g_cocaine.sim, model = poweRlaw::disexp)
vuong_heroin.sim <- vuong_simulator(simulations = g_heroin.sim, model = poweRlaw::disexp)
vuong_oversize.sim <- vuong_simulator(simulations = g_oversize.sim, model = poweRlaw::disexp)
vuong_montagna.sim <- vuong_simulator(simulations = g_montagna.sim, model = poweRlaw::disexp)

# compare log normal distribution to the Poisson distribution
vuong_siren.sim <- vuong_simulator(simulations = g_siren.sim, model = poweRlaw::dispois)
vuong_togo.sim <- vuong_simulator(simulations = g_togo.sim, model = poweRlaw::dispois)
vuong_caviar.sim <- vuong_simulator(simulations = g_caviar.sim, model = poweRlaw::dispois)
vuong_cielnet.sim <- vuong_simulator(simulations = g_cielnet.sim, model = poweRlaw::dispois)
vuong_cocaine.sim <- vuong_simulator(simulations = g_cocaine.sim, model = poweRlaw::dispois)
vuong_heroin.sim <- vuong_simulator(simulations = g_heroin.sim, model = poweRlaw::dispois)
vuong_oversize.sim <- vuong_simulator(simulations = g_oversize.sim, model = poweRlaw::dispois)
vuong_montagna.sim <- vuong_simulator(simulations = g_montagna.sim, model = poweRlaw::dispois)





# plot the distribution of p-values from the Vuong likelihood ratio test
vuong_plot <- function(sim.data, obs.data, title){
  
  # required packages
  require('ggplot2'); require('scales'); library('ggplot2')
  
  # label to annotate the p-value for the Vuong test of the actual criminal networks
  sig <- obs.data
  sig <- sig$p
  label <- round(sig, digits = 2)
  
  # flag statistical significance
  sim.data$tail_flag <- sim.data$p < 0.10
  
  # p-values for the Vuong test from the networks simulated from the ERGMs
  histogram <- ggplot2::ggplot(sim.data, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(count)/sum(ggplot2::after_stat(count)), fill = tail_flag),
      binwidth = 0.05,
      color = "black",
      alpha = 0.80
      ) +
    # shaded tail for p < 0.10
    ggplot2::scale_fill_manual(values = c("TRUE" = "skyblue1", "FALSE" = "white"), guide = "none") +
    # transform y-axis to percentage scale
    ggplot2::scale_y_continuous(
      name = "PROBABILITY DENSITY FUNCTION (PDF)", 
      labels = scales::percent_format(accuracy = 1L), # 2L to round to one decimal place, 3L to round to two decimal places, etc.
      limits = c(0.00, 1.00), 
      breaks = seq(0.00, 1.00, by = 0.20)
      ) +
    ggplot2::scale_x_continuous(
      name = expression(italic(P)*"-VALUES"),
      labels = scales::label_number(accuracy = 0.01),
      breaks = seq(0.00, 1.00, by = 0.1)
      ) +
    ggplot2::coord_cartesian(xlim = c(0.00, 1.00)) +
    # p-value for the Vuong test from the actual criminal networks
    ggplot2::geom_vline(
      data = obs.data, 
      ggplot2::aes(xintercept = p), 
      color = "firebrick3", 
      linetype = "dashed", 
      linewidth = 1,
      alpha = 0.80
      ) +
    # annotation for the p-value for the actual networks
    ggplot2::annotate(geom = "label", x = label, y = 0.50, label = sprintf('%0.2f', label), size = 2) +
    ggthemes::theme_few() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 5, face = "plain"),
      axis.text.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain"),
      axis.text.y = ggplot2::element_text(color = "black", size = 5, hjust = 1.0, vjust = 0.0, face = "plain"),  
      axis.title.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.0, face = "plain"),
      axis.title.y = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain")
      )
  plot(histogram)
  return(histogram)
}

# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 2.5, 
    height = 2, 
    device = 'pdf', 
    dpi = 700
    )
}
output(vuong_plot(
  sim.data = vuong_siren.sim, 
  obs.data = vuong_siren, 
  title = "(A) SIREN AUTO THEFT RING"
  ),
  filename = "figA7a.pdf"
  )
output(vuong_plot(
  sim.data = vuong_togo.sim, 
  obs.data = vuong_togo, 
  title = "(B) TOGO AUTO THEFT RING"
  ),
  filename = "figA7b.pdf"
  )
output(vuong_plot(
  sim.data = vuong_caviar.sim, 
  obs.data = vuong_caviar, 
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
  ),
  filename = "figA7c.pdf"
  )
output(vuong_plot(
  sim.data = vuong_cielnet.sim, 
  obs.data = vuong_cielnet, 
  title = "(D) CIEL DRUG TRAFFICKING ORGANIZATION"
  ),
  filename = "figA7d.pdf"
  )
output(vuong_plot(
  sim.data = vuong_cocaine.sim, 
  obs.data = vuong_cocaine, 
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS"
  ),
  filename = "figA7e.pdf"
  )
output(vuong_plot(
  sim.data = vuong_heroin.sim, 
  obs.data = vuong_heroin, 
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT"
  ),
  filename = "figA7f.pdf"
  )
output(vuong_plot(
  sim.data = vuong_oversize.sim, 
  obs.data = vuong_oversize, 
  title = "(G) 'NDRANGHETA DRUG TRAFFICKING OPERATION" 
  ),
  filename = "figA7g.pdf"
  )
output(vuong_plot(
  sim.data = vuong_montagna.sim, 
  obs.data = vuong_montagna, 
  title = "(H) COSA NOSTRA BID-RIGGING CONSPIRACY"
  ),
  filename = "figA7h.pdf"
  )





# estimate the shape of the degree distribution from the model -------------------------------------------------------------------
mle <- function(g, model, n){
  
  # required packages
  require("poweRlaw"); require("sna")
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = FALSE)]
  
  # fit model to estimate the shape of the cumulative distribution function
  ccdf = model(d)
  
  # if-else statement to set the cut-point for the distribution
  if (any(class(model) %in% c("dislnorm", "disexp", "dispois"))) {
    ccdf$setXmin(min(d)) # estimate the model for all degree k > 0
  } else{
    # let the model set the cut-off point for degree k for the power law distribution
    k <- poweRlaw::estimate_xmin(ccdf)
    ccdf$setXmin(k)
  }
  
  # model parameters
  ccdf$setPars(poweRlaw::estimate_pars(ccdf))
  
  # if-else statement to estimates bootstrapped confidence intervals for the shape of the distribution
  if (any(class(model) %in% c("dislnorm", "disexp", "dispois"))) {
    ci <- poweRlaw::bootstrap(
      ccdf,
      xmins = min(d),
      no_of_sims = n,
      threads = 4,
      seed = 16052024
    )
  } else {
    ci <- poweRlaw::bootstrap(
      ccdf,
      # xmins = seq(min(d), max(d), 1),
      xmins = ccdf$xmin,
      no_of_sims = n,
      threads = 4,
      seed = 16052024
    )
  }
  
  # if-else statement to calculate the p-values for the hypothesis test
  if (any(class(model) %in% c("dislnorm", "disexp", "dispois"))) {
    p <- poweRlaw::bootstrap_p(
      ccdf,
      xmins = min(d),
      no_of_sims = n, 
      threads = 4, 
      seed = 16052024
    )
  } else {
    p <- poweRlaw::bootstrap_p(
      ccdf,
      # xmins = seq(min(d), max(d), 1),
      xmins = ccdf$xmin,
      no_of_sims = n, 
      threads = 4, 
      seed = 16052024
    )
  }
  
  # return model fit
  fit <- list(ccdf, ci, p)
  return(fit)
}
# log normal curve
mle_siren.mu <- mle(g_siren, model = poweRlaw::dislnorm, n = 10000)
mle_togo.mu  <- mle(g_togo, model = poweRlaw::dislnorm, n = 10000)
mle_caviar.mu <- mle(g_caviar, model = poweRlaw::dislnorm, n = 10000)
mle_cielnet.mu <- mle(g_cielnet, model = poweRlaw::dislnorm, n = 10000)
mle_cocaine.mu <- mle(g_cocaine, model = poweRlaw::dislnorm, n = 10000)
mle_heroin.mu <- mle(g_heroin, model = poweRlaw::dislnorm, n = 10000)
mle_oversize.mu <- mle(g_oversize, model = poweRlaw::dislnorm, n = 10000)
mle_montagna.mu <- mle(g_montagna, model = poweRlaw::dislnorm, n = 10000)

# power law distribution
mle_siren.pl <- mle(g_siren, model = poweRlaw::displ, n = 10000)
mle_togo.pl  <- mle(g_togo, model = poweRlaw::displ, n = 10000)
mle_caviar.pl <- mle(g_caviar, model = poweRlaw::displ, n = 10000)
mle_cielnet.pl <- mle(g_cielnet, model = poweRlaw::displ, n = 10000)
mle_cocaine.pl <- mle(g_cocaine, model = poweRlaw::displ, n = 10000)
mle_heroin.pl <- mle(g_heroin, model = poweRlaw::displ, n = 10000)
mle_oversize.pl <- mle(g_oversize, model = poweRlaw::displ, n = 10000)
mle_montagna.pl <- mle(g_montagna, model = poweRlaw::displ, n = 10000)

# exponential curve
mle_siren.lambda <- mle(g_siren, model = poweRlaw::disexp, n = 10000)
mle_togo.lambda  <- mle(g_togo, model = poweRlaw::disexp, n = 10000)
mle_caviar.lambda <- mle(g_caviar, model = poweRlaw::disexp, n = 10000)
mle_cielnet.lambda <- mle(g_cielnet, model = poweRlaw::disexp, n = 10000)
mle_cocaine.lambda <- mle(g_cocaine, model = poweRlaw::disexp, n = 10000)
mle_heroin.lambda <- mle(g_heroin, model = poweRlaw::disexp, n = 10000)
mle_oversize.lambda <- mle(g_oversize, model = poweRlaw::disexp, n = 10000)
mle_montagna.lambda <- mle(g_montagna, model = poweRlaw::disexp, n = 10000)



# estimate the shape of the simulated degree distribution from the log-normal model

# loop through networks and extract fitted mean log and sd log for each
fit_mle <- function(simulations){
  
  fit <- lapply(simulations, function(g) {
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[d > 0]  # optional: remove zeros to avoid issues
  
  # fit the log-normal model
  model <- poweRlaw::dislnorm$new(d)
  model$setXmin(min(d))  # set cut-off point manually for all degree k
  result <- poweRlaw::estimate_pars(model)
  
  # return vector of parameters
  return(result$pars)
  })
  # data frame
  results <- as.data.frame(do.call(rbind, fit))
  names(results) <- c("meanlog", "sdlog")
  return(results)
}
line.siren.sim <- fit_mle(g_siren.sim)
line.togo.sim <- fit_mle(g_togo.sim)
line.caviar.sim <- fit_mle(g_caviar.sim)
line_cielnet.sim <- fit_mle(g_cielnet.sim)
line_cocaine.sim <- fit_mle(g_cocaine.sim)
line_heroin.sim <- fit_mle(g_heroin.sim)
line_oversize.sim <- fit_mle(g_oversize.sim)
line_montagna.sim <- fit_mle(g_montagna.sim)



# calculate ccdf for each fitted log-normal
ccdf_fit <- function(g, simulations){
  
  # get the max for the degree distribution
  d <- sna::degree(g, gmode = "graph", cmode = "freeman")
  k <- 1:max(d)
    
  # cumulative distribution function
  ccdf <- apply(simulations, 1, function(stats) {
    1 - stats::plnorm(k, meanlog = stats[1], sdlog = stats[2])
    })
  # return results
  return(list(ccdf = ccdf, k = k))
}
line.siren.sim <- ccdf_fit(g = g_siren, simulations = line.siren.sim)
line.togo.sim <- ccdf_fit(g = g_togo, simulations = line.togo.sim)
line.caviar.sim <- ccdf_fit(g = g_caviar, simulations = line.caviar.sim)
line.cielnet.sim <- ccdf_fit(g = g_cielnet, simulations = line.cielnet.sim)
line.cocaine.sim <- ccdf_fit(g = g_cocaine, simulations = line.cocaine.sim)
line.heroin.sim <- ccdf_fit(g = g_heroin, simulations = line.heroin.sim)
line.oversize.sim <- ccdf_fit(g = g_oversize, simulations = line.oversize.sim)
line.montagna.sim <- ccdf_fit(g = g_montagna, simulations = line.montagna.sim)




# data points
summarize <- function(simulations){
  ccdf <- simulations$ccdf
  k <- simulations$k
  data <- data.frame(
    k = k,
    median = apply(ccdf, 1, stats::quantile, probs = 0.500, na.rm = TRUE),
    lower  = apply(ccdf, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
    upper  = apply(ccdf, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
    )
  return(data)
}
line.siren.sim <- summarize(line.siren.sim)





# plot the degree distribution for each of the criminal networks
mle_plot <- function(model){
  plot <- plot(model, draw = F)
}
plot_siren <- mle_plot(mle_siren.mu[[1]])
plot_togo <- mle_plot(mle_togo.mu[[1]])
plot_caviar <- mle_plot(mle_caviar.mu[[1]])
plot_cielnet <- mle_plot(mle_cielnet.mu[[1]])
plot_cocaine <- mle_plot(mle_cocaine.mu[[1]])
plot_heroin <- mle_plot(mle_heroin.mu[[1]])
plot_oversize <- mle_plot(mle_oversize.mu[[1]])
plot_montagna <- mle_plot(mle_montagna.mu[[1]])



# plot the maximum likelihood estimates for each of the different models
mle_line <- function(model){
  line <- lines(model, draw = F)
}
# log normal curve
line_siren.mu <- mle_line(mle_siren.mu[[1]])
line_togo.mu <- mle_line(mle_togo.mu[[1]])
line_cavair.mu <- mle_line(mle_caviar.mu[[1]])
line_cielnet.mu <- mle_line(mle_cielnet.mu[[1]])
line_cocaine.mu <- mle_line(mle_cocaine.mu[[1]])
line_heroin.mu <- mle_line(mle_heroin.mu[[1]])
line_oversize.mu <- mle_line(mle_oversize.mu[[1]])
line_montagna.mu <- mle_line(mle_montagna.mu[[1]])



# power law distribution
line_siren.pl <- mle_line(mle_siren.pl[[1]])
line_togo.pl <- mle_line(mle_togo.pl[[1]])
line_cavair.pl <- mle_line(mle_caviar.pl[[1]])
line_cielnet.pl <- mle_line(mle_cielnet.pl[[1]])
line_cocaine.pl <- mle_line(mle_cocaine.pl[[1]])
line_heroin.pl <- mle_line(mle_heroin.pl[[1]])
line_oversize.pl <- mle_line(mle_oversize.pl[[1]])
line_montagna.pl <- mle_line(mle_montagna.pl[[1]])



# exponential curve
line_siren.lambda <- mle_line(mle_siren.lambda[[1]])
line_togo.lambda <- mle_line(mle_togo.lambda[[1]])
line_cavair.lambda <- mle_line(mle_caviar.lambda[[1]])
line_cielnet.lambda <- mle_line(mle_cielnet.lambda[[1]])
line_cocaine.lambda <- mle_line(mle_cocaine.lambda[[1]])
line_heroin.lambda <- mle_line(mle_heroin.lambda[[1]])
line_oversize.lambda <- mle_line(mle_oversize.lambda[[1]])
line_montagna.lambda <- mle_line(mle_montagna.lambda[[1]])



# plot the log normal curves for each of the criminal networks estimated from the model ------------------------------------------------------------------------------
ccdf_plot <- function(plot, line1, line2, sims, title){
  require("ggplot2"); require("scales"); require("ggthemes")
  # tutorial on how to plot poweRlaw() objects with ggplot2
  # https://rpubs.com/lgadar/power-law
  ccdf <- ggplot2::ggplot(plot) + 
  ggplot2::geom_point(
    ggplot2::aes(x = x, y = y), 
    shape = 21, 
    size = 3, 
    color = "black", 
    fill = "white"
    ) +
  # fitted power law
  ggplot2::geom_line(
    data = line2, 
    ggplot2::aes(x = x, y = y), 
    color = "grey80", 
    linetype = "solid", 
    linewidth = 1, 
    alpha = 0.50
    ) +
  # fitted log-normal distribution
  ggplot2::geom_line(
    data = line1, 
    ggplot2::aes(x = x, y = y), 
    color = "firebrick1", 
    linetype = "solid", 
    linewidth = 1, 
    alpha = 1.00
    ) +
  # fitted log-normal distribution (simulations)
  ggplot2::geom_ribbon(
      data = sims, 
      ggplot2::aes(x = k, ymin = lower, ymax = upper), 
      fill = "skyblue1", 
      alpha = 0.5
      ) +
  ggplot2::geom_line(
    data = sims,
    ggplot2::aes(x = k, y = median), 
    color = "skyblue1", 
    linewidth = 1,
    alpha = 1.00
    ) +
    # tutorial to add log scales and tick marks
    # https://www.datanovia.com/en/blog/ggplot-log-scale-transformation/
    # ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
    # ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
  ggplot2::scale_x_log10(limits = c(1, 30)) +
  ggplot2::scale_y_log10(limits = c(0.01, 1.00)) +
  ggplot2::annotation_logticks(sides = "trbl") +
  ggplot2::labs(
    title = title,
    y = "LOGGED CUMULATIVE DISTRIBUTION FUNCTION (CDF)", 
    x = "LOGGED DEGREE"
    ) +
  ggthemes::theme_few() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 5, hjust = 0.5, face = "plain", color = "black"),
    axis.title = ggplot2::element_text(color = "black", size = 5, face = "plain"),
    axis.text.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain"),
    axis.text.y = ggplot2::element_text(color = "black", size = 5, hjust = 1.0, vjust = 0.0, face = "plain"),  
    axis.title.x = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.0, face = "plain"),
    axis.title.y = ggplot2::element_text(color = "black", size = 5, hjust = 0.5, vjust = 0.5, face = "plain")
    )
  print(ccdf)
  return(ccdf)
}
ccdf_siren <- ccdf_plot(
  plot = plot_siren, 
  line1 = line_siren.mu,
  line2 = line_siren.pl,
  sims = test,
  title = "(A) SIREN AUTO THEFT RING"
  )
ccdf_togo <- ccdf_plot(
  plot = plot_togo, 
  line1 = line_togo.mu,
  line2 = line_togo.pl,
  title = "(B) TOGO AUTO THEFT RING"
  )
ccdf_caviar <- ccdf_plot(
  plot = plot_caviar, 
  line1 = line_cavair.mu,
  line2 = line_cavair.pl,
  title = "(C) CAVAIR DRUG TRAFFICKING ORGANIZATION"
  )
ccdf_cielnet <- ccdf_plot(
  plot = plot_cielnet, 
  line1 = line_cielnet.mu, 
  line2 = line_cielnet.pl,
  title = "(D) CIEL DRUG COCAINE TRAFFICKING ORGANIZATION"
  )
ccdf_cocaine <- ccdf_plot(
  plot = plot_cocaine, 
  line1 = line_cocaine.mu,
  line2 = line_cocaine.pl,
  title = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS"
  )
ccdf_heroin <- ccdf_plot(
  plot = plot_heroin, 
  line1 = line_heroin.mu,
  line2 = line_heroin.pl,
  title = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT"
  )
ccdf_oversize <- ccdf_plot(
  plot = plot_oversize, 
  line1 = line_oversize.mu,
  line2 = line_oversize.pl,
  title = "(G) 'NDRANGHETA DRUG TRAFFICKING OPERATION"
  )
ccdf_montagna <- ccdf_plot(
  plot = plot_montagna, 
  line1 = line_montagna.mu,
  line2 = line_montagna.pl,
  title = "(H) COSA NOSTRA BID-RIGGING CONSPIRACY"
  )

# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 2.5, 
    height = 2, 
    device = 'pdf', 
    dpi = 700
  )
}
output(plot = ccdf_siren, filename = "figA8a.pdf")
output(plot = ccdf_togo, filename = "figA8b.pdf")
output(plot = ccdf_caviar, filename = "figA8c.pdf")
output(plot = ccdf_cielnet, filename = "figA8d.pdf")
output(plot = ccdf_cocaine, filename = "figA8e.pdf")
output(plot = ccdf_heroin, filename = "figA8f.pdf")
output(plot = ccdf_oversize, filename = "figA8g.pdf")
output(plot = ccdf_montagna, filename = "figA8h.pdf")



# close .r script



