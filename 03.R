#  -----------------------------------------------------------------------------------

# file 03: estimate the shape of the degree distributions of each graph

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------



# plot complimentary cumulative degree distribution ----------------------------------
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


# first, estimate the shape of the degree distribution and test it against different statistical distributions
# start from the assumption that the degree is log-normal... Broido and Clauset 2019 find that most social networks produce log-normal degree distributions
# 1) compare to the power law distribution, another heavy tailed distribution
# 2) compare to the exponential distribution, a statistical distribution that does not have a heavy tail
# 3) compare to the Poisson distribution, a statistical distribution that does not have a heavy tail


# Vuong's likelihood ratio test for goodness-of-fit for the power law distribution ----------------------------------------------------------------
vuong_power.law = function(g){
  
  # required packages
  require("poweRlaw")
  
  # interpretation of the goodness-of-fit test
  cat("\n") # space
  message("Vuong's likelihood ratio test is a sign test:")
  message("* A likelihood-ratio test statistic > 0 suggests the degree distribution more closely resembles the shape of the power law distribution.")
  message("* A likelihood-ratio test statistic < 0 suggests the degree distribution more closely resembles the shape of the log normal distribution.")
  message("* A larger test statistic, in either dirtection, suggests better goodness-of-fit.")
  cat("\n") # space
  message("Hypothesis statements for Vuong's likelihood ratio test:")
  message("    H0: The degree distribution does not resemble the shape of the power law or log normal distribution.")
  message("    H1: The degree distribution resembles the shape of the power law or log normal distribution.")
  cat("\n") # space
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = TRUE)]
  
  # fit power law distribution
  alpha = poweRlaw::displ(d)
  alpha$setXmin(1) # for all degree k >= 1
  alpha$setPars(poweRlaw::estimate_pars(alpha))
  
  # fit log normal distribution
  mu = poweRlaw::dislnorm(d)
  mu$setXmin(1) # for all degree k => 1
  mu$setPars(poweRlaw::estimate_pars(mu))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(alpha, mu)
  # notes on the interpretation of Vuong's likelihood-ratio test statistic: 
  # the sign of the test statistic (i.e., +/-) has meaning for interpretation (Vuong's formula is a sign test)
  # because 'alpha' is the first input into the poweRlaw::compare_distributions() function and 'mu' is the second:
  # ... a positive (+) test statistic suggests the degree distribution more so resembles the power law distribution
  # ... a negative (-) test statistic suggests the degree distribution more so resembles the log normal distribution
  # ... reversing the input order in poweRlaw::compare_distributions() (i.e., 'mu' before 'alpha') computes the same test statistic, but in the opposite direction
  
  # test results
  message("Results of Vuong's likelihood ratio test:")
  cat("\n") # space
  
  # Vuong's likelihood-ratio test statistic 
  lrstat = lrtest$test_statistic
  cat("Vuong's likelihood-ratio test statistic = "); cat(lrstat ); cat("\n"); cat("\n")
  
  # p-value (two-tailed)
  p = lrtest$p_two_sided
  # if p < 0.05, reject H0: the degree distribution neither resembles the power law distribution or Poisson distribution
  # if p > 0.05, fail to reject H1: degree distribution resemebles power law distribution or Poisson distribution (see notes on Vuong's likelihood-ratio test)
  cat( "p-value = "); cat(p); cat("\n"); cat("\n")
  
  # interpretation:
  if(lrstat > 0){
  message("The degree distribution more closely resembles the power law distribution.")
    if(p < 0.05){
  cat("\n")
  message("Reject the null hypothesis that the degree distribution does not resemble the power law distribution.")
     } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the power law distribution."}
    } else {
  message("The degree distribution more closely resembles the log normal distribution.")
    if(p < 0.05){
    cat("\n") 
    message("Reject the null hypothesis that the degree distribution does not resemble the log normal distribution.")
     } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the log normal distribution."}
    }
}
# test results for the criminal networks
vuong_power.law(g_siren)
vuong_power.law(g_togo)
vuong_power.law(g_caviar)
vuong_power.law(g_cielnet)
vuong_power.law(g_cocaine)
vuong_power.law(g_heroin)
vuong_power.law(g_oversize)
vuong_power.law(g_montagna)
vuong_power.law(g_super)



# Vuong's likelihood ratio test for goodness-of-fit for the exponential distribution ----------------------------------------------------------------
vuong_exponential = function(g){
  
  # required packages
  require("poweRlaw")
  
  # interpretation of the goodness-of-fit test
  cat("\n") # space
  message("Vuong's likelihood ratio test is a sign test:")
  message("* A likelihood-ratio test statistic > 0 suggests the degree distribution more closely resembles the shape of the exponential distribution.")
  message("* A likelihood-ratio test statistic < 0 suggests the degree distribution more closely resembles the shape of the log normal distribution.")
  message("* A larger test statistic, in either dirtection, suggests better goodness-of-fit.")
  cat("\n") # space
  message("Hypothesis statements for Vuong's likelihood ratio test:")
  message("    H0: The degree distribution does not resemble the shape of the exponential or log normal distribution.")
  message("    H1: The degree distribution resembles the shape of the exponential or log normal distribution.")
  cat("\n") # space
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = TRUE)]
  
  # fit exponential distribution
  lambda = poweRlaw::disexp(d)
  lambda$setXmin(1) # for all degree k >= 1
  lambda$setPars(poweRlaw::estimate_pars(lambda))
  
  # fit log normal distribution
  mu = poweRlaw::dislnorm(d)
  mu$setXmin(1) # for all degree k => 1
  mu$setPars(poweRlaw::estimate_pars(mu))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(lambda, mu)
  # notes on the interpretation of Vuong's likelihood-ratio test statistic: 
  # the sign of the test statistic (i.e., +/-) has meaning for interpretation (Vuong's formula is a sign test)
  # because 'alpha' is the first input into the poweRlaw::compare_distributions() function and 'mu' is the second:
  # ... a positive (+) test statistic suggests the degree distribution more so resembles the power law distribution
  # ... a negative (-) test statistic suggests the degree distribution more so resembles the log normal distribution
  # ... reversing the input order in poweRlaw::compare_distributions() (i.e., 'mu' before 'alpha') computes the same test statistic, but in the opposite direction
  
  # test results
  message("Results of Vuong's likelihood ratio test:")
  cat("\n") # space
  
  # Vuong's likelihood-ratio test statistic 
  lrstat = lrtest$test_statistic
  cat("Vuong's likelihood-ratio test statistic = "); cat(lrstat ); cat("\n"); cat("\n")
  
  # p-value (two-tailed)
  p = lrtest$p_two_sided
  # if p < 0.05, reject H0: the degree distribution neither resembles the power law distribution or Poisson distribution
  # if p > 0.05, fail to reject H1: degree distribution resemebles power law distribution or Poisson distribution (see notes on Vuong's likelihood-ratio test)
  cat( "p-value = "); cat(p); cat("\n"); cat("\n")
  
  # interpretation:
  if(lrstat > 0){
    message("The degree distribution more closely resembles the exponential distribution.")
    if(p < 0.05){
      cat("\n")
      message("Reject the null hypothesis that the degree distribution does not resemble the exponential distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the exponential distribution."}
  } else {
    message("The degree distribution more closely resembles the log normal distribution.")
    if(p < 0.05){
      cat("\n") 
      message("Reject the null hypothesis that the degree distribution does not resemble the log normal distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the log normal distribution."}
  }
}
# test results for the criminal networks
vuong_exponential(g_siren)
vuong_exponential(g_togo)
vuong_exponential(g_caviar)
vuong_exponential(g_cielnet)
vuong_exponential(g_cocaine)
vuong_exponential(g_heroin)
vuong_exponential(g_oversize)
vuong_exponential(g_montagna)
vuong_exponential(g_super)



# Vuong's likelihood ratio test for goodness-of-fit for the Poisson law distribution ----------------------------------------------------------------
vuong_Poisson = function(g){
  
  # required packages
  require("poweRlaw")
  
  # interpretation of the goodness-of-fit test
  cat("\n") # space
  message("Vuong's likelihood ratio test is a sign test:")
  message("* A likelihood-ratio test statistic > 0 suggests the degree distribution more closely resembles the shape of the Poisson distribution.")
  message("* A likelihood-ratio test statistic < 0 suggests the degree distribution more closely resembles the shape of the log normal distribution.")
  message("* A larger test statistic, in either dirtection, suggests better goodness-of-fit.")
  cat("\n") # space
  message("Hypothesis statements for Vuong's likelihood ratio test:")
  message("    H0: The degree distribution does not resemble the shape of the Poisson or log normal distribution.")
  message("    H1: The degree distribution resembles the shape of the Poisson or log normal distribution.")
  cat("\n") # space
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = TRUE)]
  
  # fit exponential distribution
  lambda = poweRlaw::dispois(d)
  lambda$setXmin(1) # for all degree k >= 1
  lambda$setPars(poweRlaw::estimate_pars(lambda))
  
  # fit log normal distribution
  mu = poweRlaw::dislnorm(d)
  mu$setXmin(1) # for all degree k => 1
  mu$setPars(poweRlaw::estimate_pars(mu))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(lambda, mu)
  # notes on the interpretation of Vuong's likelihood-ratio test statistic: 
  # the sign of the test statistic (i.e., +/-) has meaning for interpretation (Vuong's formula is a sign test)
  # because 'alpha' is the first input into the poweRlaw::compare_distributions() function and 'mu' is the second:
  # ... a positive (+) test statistic suggests the degree distribution more so resembles the power law distribution
  # ... a negative (-) test statistic suggests the degree distribution more so resembles the log normal distribution
  # ... reversing the input order in poweRlaw::compare_distributions() (i.e., 'mu' before 'alpha') computes the same test statistic, but in the opposite direction
  
  # test results
  message("Results of Vuong's likelihood ratio test:")
  cat("\n") # space
  
  # Vuong's likelihood-ratio test statistic 
  lrstat = lrtest$test_statistic
  cat("Vuong's likelihood-ratio test statistic = "); cat(lrstat ); cat("\n"); cat("\n")
  
  # p-value (two-tailed)
  p = lrtest$p_two_sided
  # if p < 0.05, reject H0: the degree distribution neither resembles the power law distribution or Poisson distribution
  # if p > 0.05, fail to reject H1: degree distribution resemebles power law distribution or Poisson distribution (see notes on Vuong's likelihood-ratio test)
  cat( "p-value = "); cat(p); cat("\n"); cat("\n")
  
  # interpretation:
  if(lrstat > 0){
    message("The degree distribution more closely resembles the Poisson distribution.")
    if(p < 0.05){
      cat("\n")
      message("Reject the null hypothesis that the degree distribution does not resemble the Poisson distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the Poisson distribution."}
  } else {
    message("The degree distribution more closely resembles the log normal distribution.")
    if(p < 0.05){
      cat("\n") 
      message("Reject the null hypothesis that the degree distribution does not resemble the log normal distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the log normal distribution."}
  }
}
# test results for the criminal networks
vuong_Poisson(g_siren)
vuong_Poisson(g_togo)
vuong_Poisson(g_caviar)
vuong_Poisson(g_cielnet)
vuong_Poisson(g_cocaine)
vuong_Poisson(g_heroin)
vuong_Poisson(g_oversize)
vuong_Poisson(g_montagna)
vuong_Poisson(g_super)



# estimate the shape of the degree distribution from the log-normal model -------------------------------------------------------------------
mle <- function(g, sims){
  require("poweRlaw")
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = TRUE)]
  # d <- d/max(d)
  
  # fit log normal distribution to estimate the decay parameter
  mu = poweRlaw::dislnorm(d) #  mu = poweRlaw::conlnorm(d)
  mu$setXmin(min(d)) # for all degree k => 1
  mu$setPars(poweRlaw::estimate_pars(mu))
  
  # bootstrapped confidence intervals for the shape of the distribution
  ci <- poweRlaw::bootstrap(
    mu,
    xmins = min(d),
    # xmins = seq(1, k, 1),
    no_of_sims = sims,
    threads = 4,
    seed = 16052024
    )
  
  # p-value for the log normal distribution 
  p <- poweRlaw::bootstrap_p(
    mu,
    xmins = min(d),
    # xmins = seq(1, k, 1),
    no_of_sims = sims, 
    threads = 4, 
    seed = 16052024
    )
  
  # return model fit
  fit <- list(mu, ci, p)
  return(fit)
}
mle_siren <- mle(g_siren, sims = 10000)
mle_togo  <- mle(g_togo, sims = 10000)
mle_caviar <- mle(g_caviar, sims = 10000)
mle_cielnet <- mle(g_cielnet, sims = 10000)
mle_cocaine <- mle(g_cielnet, sims = 10000)
mle_heroin <- mle(g_heroin, sims = 10000)
mle_oversize <- mle(g_oversize, sims = 10000)
mle_montagna <- mle(g_montagna, sims = 10000)
mle_super <- mle(g_super, sims = 10000)


# plot the degree distribution for each of the criminal networks
mle_plot <- function(model){
  plot <- plot(model, draw = F)
}
plot_siren <- mle_plot(mle_siren[[1]])
plot_togo <- mle_plot(mle_togo[[1]])
plot_caviar <- mle_plot(mle_caviar[[1]])
plot_cielnet <- mle_plot(mle_cielnet[[1]])
plot_cocaine <- mle_plot(mle_cocaine[[1]])
plot_heroin <- mle_plot(mle_heroin[[1]])
plot_oversize <- mle_plot(mle_oversize[[1]])
plot_montagna <- mle_plot(mle_montagna[[1]])
plot_super <- mle_plot(mle_super[[1]])


# plot the maximum likelihood estimate for each of the criminal networks
mle_line <- function(model){
  line <- lines(model, draw = F)
}
line_siren <- mle_line(mle_siren[[1]])
line_togo <- mle_line(mle_togo[[1]])
line_cavair <- mle_line(mle_caviar[[1]])
line_cielnet <- mle_line(mle_cielnet[[1]])
line_cocaine <- mle_line(mle_cocaine[[1]])
line_heroin <- mle_line(mle_heroin[[1]])
line_oversize <- mle_line(mle_oversize[[1]])
line_montagna <- mle_line(mle_montagna[[1]])
line_super <- mle_line(mle_super[[1]])


# plot the log normal curves for each of the criminal networks estimated from the model ------------------------------------------------------------------------------
ccdf_plot <- function(plot, line, title){
  require(ggplot2); require(scales)
  # tutorial on how to plot poweRlaw() objects with ggplot2
  # https://rpubs.com/lgadar/power-law
  ccdf <- ggplot2::ggplot(plot) + 
  ggplot2::geom_point(ggplot2::aes(x = x, y = y), shape = 21, size = 5, color = "black", fill = "white") +
  ggplot2::geom_line(data = line, ggplot2::aes(x = x, y = y), color = "firebrick1", linewidth = 2) +
  ggplot2::ggtitle(title) +
  ggplot2::labs(y = "LOGGED CUMULATIVE DISTRIBUTION FUNCTION (CDF)", x = "LOGGED DEGREE") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain"),
    axis.text.y = ggplot2::element_text(color = "black", size = 8, hjust = 1.0, vjust = 0.0, face = "plain"),  
    axis.title.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.0, face = "plain"),
    axis.title.y = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain")
    ) +
  # tutorial to add log scales and tick marks
  # https://www.datanovia.com/en/blog/ggplot-log-scale-transformation/
  ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", math_format(10^.x))) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", math_format(10^.x))) +
  ggplot2::annotation_logticks(sides = "trbl") +
  ggplot2::theme_bw()
}
ccdf_siren <- ccdf_plot(plot = plot_siren, line = line_siren, title = "(A) SIREN AUTO THEFT RING")
ccdf_togo <- ccdf_plot(plot = plot_togo, line = line_togo, title = "(B) TOGO AUTO THEFT RING")
ccdf_caviar <- ccdf_plot(plot = plot_caviar, line = line_cavair, title = "(C) CAVAIR DRUG TRAFFICKING ORGANIZATION")
ccdf_cielnet <- ccdf_plot(plot = plot_cielnet, line = line_cielnet, title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION")
ccdf_cocaine <- ccdf_plot(plot = plot_cocaine, line = line_cocaine, title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT")
ccdf_heroin <- ccdf_plot(plot = plot_heroin, line = plot_heroin, title = "(F) NEW YORK CITY HEROIN TRAFFICKERS")
ccdf_oversize <- ccdf_plot(plot = plot_oversize, line = plot_oversize, title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE")
ccdf_montagna <- ccdf_plot(plot = plot_montagna, line = plot_montagna, title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA")
ccdf_super <- ccdf_plot(plot = plot_super, line = line_super, title = "(I) SUPER POPULATION OF CRIMINAL NETWORKS")



# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop", 
    width = 5, 
    height = 5, 
    device = 'pdf', 
    dpi = 700
    )
}
output(plot = ccdf_siren, filename = "fig2a.pdf")
output(plot = ccdf_togo, filename = "fig2b.pdf")
output(plot = ccdf_caviar, filename = "fig2c.pdf")
output(plot = ccdf_cielnet, filename = "fig2d.pdf")
output(plot = ccdf_cocaine, filename = "fig2e.pdf")
output(plot = ccdf_heroin, filename = "fig2f.pdf")
output(plot = ccdf_oversize, filename = "fig2g.pdf")
output(plot = ccdf_montagna, filename = "fig2h.pdf")
output(plot = ccdf_super, filename = "fig2i.pdf")



# close .r script


