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
cdf_degree(g_tfc)
cdf_degree(g_super)



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
  d <- d[order(d, decreasing = TRUE)]
  
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
  
  # p-value (two-tailed)
  p = lrtest$p_two_sided
  # if p < 0.05, reject H0: the degree distribution neither resembles the power law distribution or Poisson distribution
  # if p > 0.05, fail to reject H1: degree distribution resemebles power law distribution or Poisson distribution (see notes on Vuong's likelihood-ratio test)
  cat( "p-value = "); cat(p); cat("\n"); cat("\n")
  
  # interpretation:
  if(lrstat > 0){
    message("The degree distribution more closely resembles the log normal distribution.")
    if(p < 0.05){
      cat("\n")
      message("Reject the null hypothesis that the degree distribution does not resemble the log normal distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the log normal distribution."}
  } else {
    message("The degree distribution more closely resembles the comparison distribution.")
    if(p < 0.05){
      cat("\n") 
      message("Reject the null hypothesis that the degree distribution does not resemble the comparison distribution.")
    } else {"Fail to reject the null hypothesis that the degree distribution does not resemble the comparison distribution."}
  }
}
# compare log normal distribution to the power law distribution
vuong(g_siren, model = poweRlaw::displ)
vuong(g_togo, model = poweRlaw::displ)
vuong(g_caviar, model = poweRlaw::displ)
vuong(g_cielnet, model = poweRlaw::displ)
vuong(g_cocaine, model = poweRlaw::displ)
vuong(g_heroin, model = poweRlaw::displ)
vuong(g_oversize, model = poweRlaw::displ)
vuong(g_montagna, model = poweRlaw::displ)
vuong(g_tfc, model = poweRlaw::displ)
vuong(g_super, model = poweRlaw::displ)

# compare log normal distribution to the exponential distribution
vuong(g_siren, model = poweRlaw::disexp)
vuong(g_togo, model = poweRlaw::disexp)
vuong(g_caviar, model = poweRlaw::disexp)
vuong(g_cielnet, model = poweRlaw::disexp)
vuong(g_cocaine, model = poweRlaw::disexp)
vuong(g_heroin, model = poweRlaw::disexp)
vuong(g_oversize, model = poweRlaw::disexp)
vuong(g_montagna, model = poweRlaw::disexp)
vuong(g_tfc, model = poweRlaw::disexp)
vuong(g_super, model = poweRlaw::disexp)

# compare log normal distribution to the Poisson distribution
vuong(g_siren, model = poweRlaw::dispois)
vuong(g_togo, model = poweRlaw::dispois)
vuong(g_caviar, model = poweRlaw::dispois)
vuong(g_cielnet, model = poweRlaw::dispois)
vuong(g_cocaine, model = poweRlaw::dispois)
vuong(g_heroin, model = poweRlaw::dispois)
vuong(g_oversize, model = poweRlaw::dispois)
vuong(g_montagna, model = poweRlaw::dispois)
vuong(g_tfc, model = poweRlaw::dispois)
vuong(g_super, model = poweRlaw::dispois)



# estimate the shape of the degree distribution from the model -------------------------------------------------------------------
mle <- function(g, model, n){
  
  # required packages
  require("poweRlaw"); require("sna")
  
  # degree distribution for the criminal networks
  d <- sna::degree(g, gmode = "graph", cmode = "freeman", rescale = FALSE)
  d <- d[order(d, decreasing = TRUE)]
  
  # fit model to estimate the decay parameter
  mu = model(d)
  mu$setXmin(min(d)) # for all degree k
  mu$setPars(poweRlaw::estimate_pars(mu))
  
  # bootstrapped confidence intervals for the shape of the distribution
  ci <- poweRlaw::bootstrap(
    mu,
    xmins = min(d),
    # xmins = seq(1, k, 1),
    no_of_sims = n,
    threads = 4,
    seed = 16052024
    )
  
  # p-value for hypothesis test
  p <- poweRlaw::bootstrap_p(
    mu,
    xmins = min(d),
    # xmins = seq(1, k, 1),
    no_of_sims = n, 
    threads = 4, 
    seed = 16052024
    )
  
  # return model fit
  fit <- list(mu, ci, p)
  return(fit)
}
# log normal curve
mle_siren.mu <- mle(g_siren, model = poweRlaw::dislnorm, n = 10000)
mle_siren.mu <- mle(g_siren, model = poweRlaw::dislnorm, n = 10000)
mle_togo.mu  <- mle(g_togo, model = poweRlaw::dislnorm, n = 10000)
mle_caviar.mu <- mle(g_caviar, model = poweRlaw::dislnorm, n = 10000)
mle_cielnet.mu <- mle(g_cielnet, model = poweRlaw::dislnorm, n = 10000)
mle_cocaine.mu <- mle(g_cielnet, model = poweRlaw::dislnorm, n = 10000)
mle_heroin.mu <- mle(g_heroin, model = poweRlaw::dislnorm, n = 10000)
mle_oversize.mu <- mle(g_oversize, model = poweRlaw::dislnorm, n = 10000)
mle_montagna.mu <- mle(g_montagna, model = poweRlaw::dislnorm, n = 10000)
mle_tfc.mu <- mle(g_tfc, model = poweRlaw::dislnorm, n = 10000)
mle_super.mu <- mle(g_super, model = poweRlaw::dislnorm, n = 10000)

# power law distribution
mle_siren.pl <- mle(g_siren, model = poweRlaw::displ, n = 10000)
mle_siren.pl <- mle(g_siren, model = poweRlaw::displ, n = 10000)
mle_togo.pl  <- mle(g_togo, model = poweRlaw::displ, n = 10000)
mle_caviar.pl <- mle(g_caviar, model = poweRlaw::displ, n = 10000)
mle_cielnet.pl <- mle(g_cielnet, model = poweRlaw::displ, n = 10000)
mle_cocaine.pl <- mle(g_cielnet, model = poweRlaw::displ, n = 10000)
mle_heroin.pl <- mle(g_heroin, model = poweRlaw::displ, n = 10000)
mle_oversize.pl <- mle(g_oversize, model = poweRlaw::displ, n = 10000)
mle_montagna.pl <- mle(g_montagna, model = poweRlaw::displ, n = 10000)
mle_tfc.pl <- mle(g_tfc, model = poweRlaw::displ, n = 10000)
mle_super.pl <- mle(g_super, model = poweRlaw::displ, n = 10000)

# exponential curve
mle_siren.lambda <- mle(g_siren, model = poweRlaw::disexp, n = 10000)
mle_siren.lambda <- mle(g_siren, model = poweRlaw::disexp, n = 10000)
mle_togo.lambda  <- mle(g_togo, model = poweRlaw::disexp, n = 10000)
mle_caviar.lambda <- mle(g_caviar, model = poweRlaw::disexp, n = 10000)
mle_cielnet.lambda <- mle(g_cielnet, model = poweRlaw::disexp, n = 10000)
mle_cocaine.lambda <- mle(g_cielnet, model = poweRlaw::disexp, n = 10000)
mle_heroin.lambda <- mle(g_heroin, model = poweRlaw::disexp, n = 10000)
mle_oversize.lambda <- mle(g_oversize, model = poweRlaw::disexp, n = 10000)
mle_montagna.lambda <- mle(g_montagna, model = poweRlaw::disexp, n = 10000)
mle_tfc.lambda <- mle(g_tfc, model = poweRlaw::disexp, n = 10000)
mle_super.lambda <- mle(g_super, model = poweRlaw::disexp, n = 10000)



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
plot_tfc <- mle_plot(mle_tfc.mu[[1]])
plot_super <- mle_plot(mle_super.mu[[1]])

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
line_tfc.mu <- mle_line(mle_tfc.mu[[1]])
line_super.mu <- mle_line(mle_super.mu[[1]])

# power law distribution
line_siren.pl <- mle_line(mle_siren.pl[[1]])
line_togo.pl <- mle_line(mle_togo.pl[[1]])
line_cavair.pl <- mle_line(mle_caviar.pl[[1]])
line_cielnet.pl <- mle_line(mle_cielnet.pl[[1]])
line_cocaine.pl <- mle_line(mle_cocaine.pl[[1]])
line_heroin.pl <- mle_line(mle_heroin.pl[[1]])
line_oversize.pl <- mle_line(mle_oversize.pl[[1]])
line_montagna.pl <- mle_line(mle_montagna.pl[[1]])
line_tfc.pl <- mle_line(mle_tfc.pl[[1]])
line_super.pl <- mle_line(mle_super.pl[[1]])

# exponential curve
line_siren.lambda <- mle_line(mle_siren.lambda[[1]])
line_togo.lambda <- mle_line(mle_togo.lambda[[1]])
line_cavair.lambda <- mle_line(mle_caviar.lambda[[1]])
line_cielnet.lambda <- mle_line(mle_cielnet.lambda[[1]])
line_cocaine.lambda <- mle_line(mle_cocaine.lambda[[1]])
line_heroin.lambda <- mle_line(mle_heroin.lambda[[1]])
line_oversize.lambda <- mle_line(mle_oversize.lambda[[1]])
line_montagna.lambda <- mle_line(mle_montagna.lambda[[1]])
line_tfc.lambda <- mle_line(mle_tfc.lambda[[1]])
line_super.lambda <- mle_line(mle_super.lambda[[1]])



# plot the log normal curves for each of the criminal networks estimated from the model ------------------------------------------------------------------------------
ccdf_plot <- function(plot, line, title){
  require("ggplot2"); require("scales"); require("ggthemes")
  # tutorial on how to plot poweRlaw() objects with ggplot2
  # https://rpubs.com/lgadar/power-law
  ccdf <- ggplot2::ggplot(plot) + 
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), shape = 21, size = 5, color = "black", fill = "white") +
    ggplot2::geom_line(data = line, ggplot2::aes(x = x, y = y), color = "firebrick1", linetype = "solid", linewidth = 2) +
    ggplot2::ggtitle(title) +
    # tutorial to add log scales and tick marks
    # https://www.datanovia.com/en/blog/ggplot-log-scale-transformation/
    ggplot2::scale_x_log10(name = "LOGGED DEGREE", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
    ggplot2::scale_y_log10(name = "LOGGED CUMULATIVE DISTRIBUTION FUNCTION (CDF)", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
    ggplot2::annotation_logticks(sides = "trbl") +
    ggthemes::theme_base() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "plain"),
      axis.text.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain"),
      axis.text.y = ggplot2::element_text(color = "black", size = 8, hjust = 1.0, vjust = 0.0, face = "plain"),  
      axis.title.x = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.0, face = "plain"),
      axis.title.y = ggplot2::element_text(color = "black", size = 8, hjust = 0.5, vjust = 0.5, face = "plain")
      )
}
ccdf_siren <- ccdf_plot(
  plot = plot_siren, 
  line = line_siren.mu,
  title = "(A) SIREN AUTO THEFT RING"
  )
ccdf_togo <- ccdf_plot(
  plot = plot_togo, 
  line = line_togo.mu,
  title = "(B) TOGO AUTO THEFT RING"
  )
ccdf_caviar <- ccdf_plot(
  plot = plot_caviar, 
  line = line_cavair.mu,
  title = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION"
  )
ccdf_cielnet <- ccdf_plot(
  plot = plot_cielnet, 
  line = line_cielnet.mu, 
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
  )
ccdf_cocaine <- ccdf_plot(
  plot = plot_cocaine, 
  line = line_cocaine.mu, 
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
  )
ccdf_heroin <- ccdf_plot(
  plot = plot_heroin, 
  line = line_heroin.mu, 
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
  )
ccdf_oversize <- ccdf_plot(
  plot = plot_oversize, 
  line = line_oversize.mu, 
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
  )
ccdf_montagna <- ccdf_plot(
  plot = plot_montagna, 
  line = line_montagna.mu, 
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
  )
ccdf_tfc <- ccdf_plot(
  plot = plot_tfc,
  line = line_tfc.mu, 
  title = "(I) THE FRENCH CONNECTION - FEDERAL BUREAU OF NARCOTICS"
  )
ccdf_super <- ccdf_plot(
  plot = plot_super, 
  line = line_super.mu, 
  title = "(J) SUPER POPULATION OF CRIMINAL NETWORKS"
  )


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
output(plot = ccdf_tfc, filename = "fig2i.pdf")
output(plot = ccdf_super, filename = "fig2j.pdf")



# close .r script



# plot the log normal curves for each of the criminal networks estimated from the model ------------------------------------------------------------------------------
ccdf_plot <- function(plot, line1, line2, line3, title){
  require(ggplot2); require(scales)
  # tutorial on how to plot poweRlaw() objects with ggplot2
  # https://rpubs.com/lgadar/power-law
  ccdf <- ggplot2::ggplot(plot) + 
  ggplot2::geom_point(ggplot2::aes(x = x, y = y), shape = 21, size = 5, color = "black", fill = "white") +
  ggplot2::geom_line(data = line1, ggplot2::aes(x = x, y = y), color = "firebrick1", linetype = "solid", linewidth = 2) +
  ggplot2::geom_line(data = line2, ggplot2::aes(x = x, y = y), color = "skyblue1", linetype = "solid", linewidth = 2) +
  ggplot2::geom_line(data = line3, ggplot2::aes(x = x, y = y), color = "green3", linetype = "solid", linewidth = 2) +
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
  ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", math_format(10^.x))) +
  ggplot2::annotation_logticks(sides = "trbl") +
  ggthemes::theme_base()
}
ccdf_siren <- ccdf_plot(
  plot = plot_siren, 
  line1 = line_siren.mu,
  line2 = line_siren.lambda,
  line3 = line_siren.pl,
  title = "(A) SIREN AUTO THEFT RING"
  )
ccdf_togo <- ccdf_plot(
  plot = plot_togo, 
  line1 = line_togo.mu, 
  line2 = line_togo.lambda, 
  line3 = line_togo.pl,
  title = "(B) TOGO AUTO THEFT RING"
  )
ccdf_caviar <- ccdf_plot(
  plot = plot_caviar, 
  line1 = line_cavair.mu,
  line2 = line_cavair.lambda,
  line3 = line_cavair.pl,
  title = "(C) CAVAIR DRUG TRAFFICKING ORGANIZATION"
  )
ccdf_cielnet <- ccdf_plot(
  plot = plot_cielnet, 
  line1 = line_cielnet.mu, 
  line2 = line_cielnet.lambda,
  line3 = line_cielnet.pl,
  title = "(D) CIELNET DRUG TRAFFICKING ORGANIZATION"
  )
ccdf_cocaine <- ccdf_plot(
  plot = plot_cocaine, 
  line1 = line_cocaine.mu,
  line2 = line_cocaine.lambda,
  line3 = line_cocaine.pl,
  title = "(E) LA COSA NOSTRA COCAINE TRAFFICKING OUTFIT"
  )
ccdf_heroin <- ccdf_plot(
  plot = plot_heroin, 
  line1 = line_heroin.mu, 
  line2 = line_heroin.lambda,
  line3 = line_heroin.pl,
  title = "(F) NEW YORK CITY HEROIN TRAFFICKERS"
  )
ccdf_oversize <- ccdf_plot(
  plot = plot_oversize, 
  line1 = line_oversize.mu,
  line2 = line_oversize.lambda,
  line3 = line_oversize.pl,
  title = "(G) 'NDRANGHETA WIRETAPS - OPERATION OVERSIZE"
  )
ccdf_montagna <- ccdf_plot(
  plot = plot_montagna, 
  line1 = line_montagna.mu,
  line2 = line_montagna.lambda,
  line3 = line_montagna.pl,
  title = "(H) COSA NOSTRA WIRETAPS - OPERATION MONTAGNA"
  )
ccdf_tfc <- ccdf_plot(
  plot = plot_tfc,
  line1 = line_tfc.mu, 
  line2 = line_tfc.lambda,
  line3 = line_tfc.pl,
  title = "(I) THE FRENCH CONNECTION - FEDERAL BUREAU OF NARCOTICS"
  )
ccdf_super <- ccdf_plot(
  plot = plot_super, 
  line1 = line_super.mu,
  line2 = line_super.lambda,
  line3 = line_super.pl,
  title = "(J) SUPER POPULATION OF CRIMINAL NETWORKS"
  )


