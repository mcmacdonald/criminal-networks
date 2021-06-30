#-------------------------------------------------------------

# File 02: B-A models to fit power law distributions

# call packages
library('igraph')

#-------------------------------------------------------------


# File contents:
# - A function to estimate degree scaling exponents in networks


# Function to fit the power law distribution ---------------------------------------------
pk = function(graph) {
  
  # calculate degree
  d = igraph::degree(
    graph = graph, 
    v = igraph::V(graph), 
    mode = "all", 
    loops = FALSE, 
    normalized = FALSE
    )
  dd = igraph::degree.distribution(
    graph = graph, 
    mode = "all", 
    cumulative = FALSE
    )
  degree = 1:max(d) # range of centralities from the network
  
  # calculate probability of the degree frequency distribution, where p = probability
  p = dd[-1] 
  pk1 = which(p != 0) # remove degree = 0
  p = p[pk1]
  
  # overwrite degree distribution when degree centrality > 0
  degree = degree[pk1]
  
  # scaling coefficients
  m = lm(log(p) ~ log(degree))
  b = coef(m)
  print(summary(m))
  print(confint.lm(m)) # 95% CIs
  print(paste("N =", nobs(m) 
              ) 
        )
  plfit = function(x) exp(b[[1]] + b[[2]] * log(x))
  alpha = -b[[2]]
  R2    = summary(m)$r.squared
  print(paste("Alpha =", round(x = alpha, digits = 3) 
              ) 
        )
  print(paste("R2 =", round(x = R2, digits = 3)
              )
        )
  
  # plot the linear regression
  k1 <- 1
  kn <- max(
    igraph::degree(
      graph = graph, 
      v = igraph::V(graph), 
      mode = "total", 
      loops = FALSE, 
      normalized = FALSE
      )
    ) 
  options(scipen = 999) # turn off scientific notation 
  plot(p ~ degree, 
       log = "xy",
       xlim = c(k1, kn), # scale x-axis 
       ylim = c(0.0001, 1),
       xlab = "degree k (log)", 
       ylab = "probability k (log)", 
       main = "SCALING IN DEGREE CENTRALITY",
       pch = 1, 
       cex = 2)
  curve(plfit, 
        col = "black", 
        lwd = 2, 
        add = TRUE, 
        n = length(d)
        )
  
  # End function
}

# plot results -------------------------
pk(graph = r_siren)
pk(graph = r_togo)
pk(graph = d_caviar)
pk(graph = d_cocaine)
pk(graph = d_heroin)
pk(graph = d_cielnet)
pk(graph = g_ity)
pk(graph = g_ldn)
pk(graph = g_mtl)
pk(graph = m_infinito)



# ... close .R script