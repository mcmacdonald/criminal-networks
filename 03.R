#  -----------------------------------------------------------------------------------

# file 03: estimate the shape of the degree distributions of each graph

# 'mlergm' package
# https://cran.r-project.org/web/packages/mlergm/mlergm.pdf

# don't run
# install.package('mlergm')

# last updated: 27/02/2025

# ------------------------------------------------------------------------------------




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



# compare Vuong's likelihood ratio test for goodness-of-fit of power law distribution
vuong = function(g){
  
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
  lambda = poweRlaw::dislnorm(d)
  lambda$setXmin(1) # for all degree k => 1
  lambda$setPars(poweRlaw::estimate_pars(lambda))
  
  # Vuong's likelihood-ratio test to compare models
  lrtest = poweRlaw::compare_distributions(alpha, lambda)
  # notes on the interpretation of Vuong's likelihood-ratio test statistic: 
  # the sign of the test statistic (i.e., +/-) has meaning for interpretation (Vuong's formula is a sign test)
  # because 'alpha' is the first input into the poweRlaw::compare_distributions() function and 'lambda' is the second:
  # ... a positive (+) test statistic suggests the degree distribution more so resembles the power law distribution
  # ... a negative (-) test statistic suggests the degree distribution more so resembles the Poisson distribution
  # ... reversing the input order in poweRlaw::compare_distributions() (i.e., 'lambda' before 'alpha') computes the same test statistic, but in the opposite direction
  
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
vuong(g_siren)
vuong(g_togo)
vuong(g_caviar)
vuong(g_cielnet)
vuong(g_cocaine)
vuong(g_heroin)
vuong(g_oversize)
vuong(g_montagna)
vuong(g_super)










