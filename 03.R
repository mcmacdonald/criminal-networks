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




# compute Bayes Factors for simulated networks ---------------------------------

# function to simulate networks from Bayesian models
simulator <- function(thetas, model, nsim){
  
  # required packages
  require('stats')
  
  # model parameters
  theta <- colMeans(thetas)
  
  # simulate
  sims <- stats::simulate(
    model, # model specification
    coef = theta, # model parameters
    seed = 20240517, # Malone's birthday
    nsim = nsim, # number of simulations
    output = "network" # thing to simulate
    )
  
  # return 
  return(sims)
}
g_siren.sim <- simulator(
  thetas = bayes_01.siren$Theta,
  model = bayes_01.siren$formula,
  nsim = 100
  )
g_togo.sim <- simulator(
  thetas = bayes_02.togo$Theta,
  model = bayes_02.togo$formula,
  nsim = 100
  )
g_caviar.sim <- simulator(
  thetas = bayes_03.caviar$Theta,
  model = bayes_03.caviar$formula,
  nsim = 100
  )
g_cielnet.sim <- simulator(
  thetas = bayes_04.cielnet$Theta,
  model = bayes_04.cielnet$formula,
  nsim = 100
  )
g_cocaine.sim <- simulator(
  thetas = bayes_05.cocaine$Theta,
  model = bayes_05.cocaine$formula,
  nsim = 100
  )
g_heroin.sim <- simulator(
  thetas = bayes_06.heroin$Theta,
  model = bayes_06.heroin$formula,
  nsim = 100
  )
g_oversize.sim <- simulator(
  thetas = bayes_07.oversize$Theta,
  model = bayes_07.oversize$formula,
  nsim = 100
  )
g_montagna.sim <- simulator(
  thetas = bayes_08.montagna$Theta,
  model = bayes_08.montagna$formula,
  nsim = 100
  )



# estimate the model for the simulations
bayes_simulator <- function(y_list, alpha_deg, alpha_dsp, alpha_esp){
  
  # required packages
  require('Bergm')
  
  # set seed for replication
  set.seed(20190812) # Hayes's birthday
  
  # number of simulations
  n <- length(y_list)
  
  # list to store models
  models <- vector("list", n)
  
  # loop function to estimate models from the exponential random graph model
  for (i in seq_len(n)) { 
    
    # store inside the function
    y <- y_list[[i]]
    
    # formula
    formula <- stats::as.formula(paste0("y ~ edges + ","gwdegree(", alpha_deg, ", fixed=TRUE) + ", "gwdsp(", alpha_dsp, ", fixed=TRUE) + ", "gwesp(", alpha_esp, ", fixed=TRUE) + ", "degcor"))
    
    # model estimation  
    fit <- tryCatch(
      {
        bayes(y = y, formula = formula)
      },
      error = function(e) {
        message(paste("Error in model", i, ":", conditionMessage(e)))
        return(NULL)
      # },
      # warning = function(w) {
        # message(paste("Warning in model", i, ":", conditionMessage(w)))
        # return(invisible(NULL))
      }
    )
    
    # join ith model to list of the models
    models[[i]] <- fit
    
  }
  
  return(models)
}
bayes_01.siren.sim <- bayes_simulator(
  y_list = g_siren.sim,
  alpha_deg = 0.7,
  alpha_dsp = 2.0,
  alpha_esp = 2.0
  )
bayes_02.togo.sim <- bayes_simulator(
  y_list = g_togo.sim, 
   alpha_deg = 0.7, 
   alpha_dsp = 1.7, 
   alpha_esp = 1.0
   )
bayes_03.caviar.sim <- bayes_simulator(
  y_list = g_caviar.sim, 
   alpha_deg = 2.0, 
   alpha_dsp = 2.0, 
   alpha_esp = 2.0
   )
bayes_04.cielnet.sim <- bayes_simulator(
  y_list = g_cielnet.sim, 
   alpha_deg = 1.7, 
   alpha_dsp = 0.2, 
   alpha_esp = 0.2
   )
bayes_05.cocaine.sim <- bayes_simulator(
  y_list = g_cocaine.sim, 
   alpha_deg = 2.5, 
   alpha_dsp = 0.5, 
   alpha_esp = 0.5
   )
bayes_06.heroin.sim <- bayes_simulator(
  y_list = g_heroin.sim, 
   alpha_deg = 0.5, 
   alpha_dsp = 0.1, 
   alpha_esp = 1.0
   )
bayes_07.oversize.sim <- bayes_simulator(
  y_list = g_oversize.sim, 
   alpha_deg = 1.5, 
   alpha_dsp = 0.5, 
   alpha_esp = 1.7
   )
bayes_08.montagna.sim <- bayes_simulator(
  y_list = g_montagna.sim, 
   alpha_deg = 1.2, 
   alpha_dsp = 0.5, 
   alpha_esp = 0.5
   )



# estimate Bayes Factors for the 100 simulations
BF <- function(y_list, sim_list, alpha_deg, alpha_dsp, alpha_esp, hypothesis, effect, model) {
  
  # required packages
  require('stats'); require('BFpack'); require('Bergm')
  
  # set seed for replication
  set.seed(20190812) # Hayes's birthday
  
  # number of simulations
  n <- length(y_list)
  
  # list to store models
  bf_test <- vector("list", n)
  
  # loop function to estimate models from the exponential random graph model
  for (i in seq_len(n)) { 
    
    # store inside the function
    y <- y_list[[i]]
    
    # model
    fit <- sim_list[[i]]
    
    # ERGM formula
    # formula <- stats::as.formula(paste0("y ~ edges + ","gwdegree(", alpha_deg, ", fixed=TRUE) + ", "gwdsp(", alpha_dsp, ", fixed=TRUE) + ", "gwesp(", alpha_esp, ", fixed=TRUE) + ", "degcor"))
    
    # pass through prior function to estimate Bayesian exponential random graph model
    bf_test[i] <- tryCatch({
      
    # model estimation  
    # fit <- bayes(y = y, formula = formula)
    
    # check that theta samples exist
    if (is.null(fit$Theta)) stop("bayes() failed: Theta is NULL")
    
    # posteriors 
    posterior <- colMeans(fit$Theta)
    
    # variance
    Sigma <- stats::cov(fit$Theta)
    
    # make sure the call uses 'y' from global env
    # fit$call$formula[[2]] <- quote(y)
    
    # temporarily call to global environment so that BFpack:":BF() can call the simulated data
    assign("y", y, envir = .GlobalEnv)
    
    # hypothesis test
    test <- BFpack::BF(
      fit, 
      posterior = posterior, 
      Sigma = Sigma, 
      hypothesis = hypothesis
      )
    
    # remove to ensure there is no error in the temporarily assignment
    rm("y", envir = .GlobalEnv)
    
    # Bayes factor
    test <- test$BFtable_confirmatory[1,7]
    
    # return
    as.numeric(test)
    
    # in case of error, don't interrupt/stop
    }, error = function(e) {
      warning(sprintf("Model %d failed: %s", i, e$message))
      NA_real_
    })
  }
  # classify the strength of the evidence, per Rafferty (1995) ... see Table 6 (p. 139)
  # https://www.jstor.org/stable/271063?origin=crossref&seq=29
  evidence <- sapply(bf_test, function(bf) {
    if (is.na(bf)) { # error
      "NA" 
    } else if (bf < 1) { # evidence to support null hypothesis
      "NULL" 
    } else if (bf < 3) { # weak evidence for hypothesis
      "WEAK" 
    } else if (bf < 20) { # moderate evidence for hypothesis
      "MODERATE"
    } else if (bf < 150) { # strong evidence for hypothesis
      "STRONG"
    } else { # very strong evidence for hypothesis
      "VERY STRONG"
    }
  })
  # return tidy data frame
  results <- data.frame(
    BF = as.numeric(bf_test),
    evidence = factor(
      evidence, 
      levels = c(
        "NULL", 
        "WEAK", 
        "MODERATE", 
        "STRONG", 
        "VERY STRONG", 
        "NA"
        )
      ),
    effect = rep(effect, n), # repeat for each simulation
    model = rep(model, n), # repeat for each simulation
    stringsAsFactors = FALSE
    )

  # return the results
  return(results)
}

# centralization i.e., geometrically weighted degree
bf_centralization <- rbind(
  BF(y_list = g_siren.sim, 
     sim_list = bayes_01.siren.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "SIREN"
     ),
  BF(y_list = g_togo.sim, 
     sim_list = bayes_02.togo.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 1.7, 
     alpha_esp = 1.0,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "TOGO"
     ),
  BF(y_list = g_caviar.sim,
     sim_list = bayes_03.caviar.sim,
     alpha_deg = 2.0, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "CAVIAR"
     ),
  BF(y_list = g_cielnet.sim, 
     sim_list = bayes_04.cielnet.sim,
     alpha_deg = 1.7, 
     alpha_dsp = 0.2, 
     alpha_esp = 0.2,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "CIELNET"
     ),
  BF(y_list = g_cocaine.sim, 
     sim_list = bayes_05.cocaine.sim,
     alpha_deg = 2.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "COCAINE"
     ),
  BF(y_list = g_heroin.sim, 
     sim_list = bayes_06.heroin.sim,
     alpha_deg = 0.5, 
     alpha_dsp = 0.1, 
     alpha_esp = 1.0,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "HEROIN"
     ),
  BF(y_list = g_oversize.sim, 
     sim_list = bayes_07.oversize.sim,
     alpha_deg = 1.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 1.7,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "OVERSIZE"
     ),
  BF(y_list = g_montagna.sim, 
     sim_list = bayes_08.montagna.sim,
     alpha_deg = 1.2, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta2 < 0",
     effect = "CENTRALIZATION",
     model = "MONTAGNA"
     )
  )

# common co-conspirators i.e., geometrically weighted dyadwise shared partners
bf_dsp <- rbind(
  BF(y = g_siren.sim,
     sim_list = bayes_01.siren.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "SIREN"
     ),
  BF(y = g_togo.sim, 
     sim_list = bayes_02.togo.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 1.7, 
     alpha_esp = 1.0,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "TOGO"
     ),
  BF(y = g_caviar.sim, 
     sim_list = bayes_03.caviar.sim,
     alpha_deg = 2.0, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "CAVIAR"
     ),
  BF(y = g_cielnet.sim, 
    sim_list = bayes_04.cielnet.sim,
    alpha_deg = 1.7, 
    alpha_dsp = 0.2, 
    alpha_esp = 0.2,
    hypothesis = "theta3 > 0",
    effect = "COMMON FRIENDS",
    model = "CIELNET"
    ),
  BF(y = g_cocaine.sim,
     sim_list = bayes_05.cocaine.sim,
     alpha_deg = 2.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "COCAINE"
     ),
  BF(y = g_heroin.sim, 
     sim_list = bayes_06.heroin.sim,
     alpha_deg = 0.5, 
     alpha_dsp = 0.1, 
     alpha_esp = 1.0,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "HEROIN"
     ),
  BF(y = g_oversize.sim, 
     sim_list = bayes_07.oversize.sim,
     alpha_deg = 1.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 1.7,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "OVERSIZE"
     ),
  BF(y = g_montagna.sim, 
     sim_list = bayes_08.montagna.sim,
     alpha_deg = 1.2, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta3 > 0",
     effect = "COMMON FRIENDS",
     model = "MONTAGNA"
     )
  )

# triadic closure i.e., geometrically weighted edgewise shared partners
bf_esp <- rbind(
  BF(y = g_siren.sim,
     sim_list = bayes_01.siren.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "SIREN"
     ),
  BF(y = g_togo.sim, 
     sim_list = bayes_02.togo.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 1.7, 
     alpha_esp = 1.0,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "TOGO"
     ),
  BF(y = g_caviar.sim, 
     sim_list = bayes_03.caviar.sim,
     alpha_deg = 2.0, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "CAVIAR"
     ),
  BF(y = g_cielnet.sim, 
     sim_list = bayes_04.cielnet.sim,
     alpha_deg = 1.7, 
     alpha_dsp = 0.2, 
     alpha_esp = 0.2,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "CIELNET"
     ),
  BF(y = g_cocaine.sim, 
     sim_list = bayes_05.cocaine.sim,
     alpha_deg = 2.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "COCAINE"
     ),
  BF(y = g_heroin.sim, 
     sim_list = bayes_06.heroin.sim,
     alpha_deg = 0.5, 
     alpha_dsp = 0.1, 
     alpha_esp = 1.0,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "HEROIN"
     ),
  BF(y = g_oversize.sim, 
     sim_list = bayes_07.oversize.sim,
     alpha_deg = 1.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 1.7,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "OVERSIZE"
     ),
  BF(y = g_montagna.sim, 
     sim_list = bayes_08.montagna.sim,
     alpha_deg = 1.2, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta4 > 0",
     effect = "TRIADIC CLOSURE",
     model = "MONTAGNA"
     )
  )

# degree assortativity i.e., degree correlation
bf_degcor <- rbind(
  BF(y = g_siren.sim,
     sim_list = bayes_01.siren.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "SIREN"
     ),
  BF(y = g_togo.sim, 
     sim_list = bayes_02.togo.sim,
     alpha_deg = 0.7, 
     alpha_dsp = 1.7, 
     alpha_esp = 1.0,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "TOGO"
     ),
  BF(y = g_caviar.sim, 
     sim_list = bayes_03.caviar.sim,
     alpha_deg = 2.0, 
     alpha_dsp = 2.0, 
     alpha_esp = 2.0,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "CAVIAR"
     ),
  BF(y = g_cielnet.sim, 
     sim_list = bayes_04.cielnet.sim,
     alpha_deg = 1.7, 
     alpha_dsp = 0.2, 
     alpha_esp = 0.2,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "CIELNET"
     ),
  BF(y = g_cocaine.sim, 
     sim_list = bayes_05.cocaine.sim,
     alpha_deg = 2.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "COCAINE"
     ),
  BF(y = g_heroin.sim, 
     sim_list = bayes_06.heroin.sim,
     alpha_deg = 0.5, 
     alpha_dsp = 0.1, 
     alpha_esp = 1.0,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "HEROIN"
     ),
  BF(y = g_oversize.sim, 
     sim_list = bayes_07.oversize.sim,
     alpha_deg = 1.5, 
     alpha_dsp = 0.5, 
     alpha_esp = 1.7,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "OVERSIZE"
     ),
  BF(y = g_montagna.sim, 
     sim_list = bayes_08.montagna.sim,
     alpha_deg = 1.2, 
     alpha_dsp = 0.5, 
     alpha_esp = 0.5,
     hypothesis = "theta5 > 0",
     effect = "ASSORTATIVITY",
     model = "MONTAGNA"
     )
  )

# join the bayes factors for the simulated networks together
`%>%` <- magrittr::`%>%`
bf <- dplyr::bind_rows(
  bf_centralization,
  bf_dsp,
  bf_esp,
  bf_degcor
  ) %>%
  dplyr::filter(!is.na(BF))  # drop any failed simulations

# set the order of facets
bf$model <- factor(bf$model, levels = c("SIREN", "TOGO", "CAVIAR", "CIELNET", "COCAINE", "HEROIN", "OVERSIZE", "MONTAGNA"))

# rename the labels
labels <- c(
  SIREN = "(A) SIREN AUTO THEFT RING",
  TOGO = "(B) TOGO AUTO THEFT RING",
  CAVIAR = "(C) CAVIAR DRUG TRAFFICKING ORGANIZATION",
  CIELNET = "(D) CIEL DRUG TRAFFICKING ORGANIZATION",
  COCAINE = "(E) CARTEL SATELLITE COCAINE TRAFFICKERS",
  HEROIN = "(F) LA COSA NOSTRA HEROIN TRAFFICKING OUTFIT",
  OVERSIZE = "(G) OVERSIZE - 'NDRANGHETA RACKETEERING",
  MONTAGNA = "(H) MONTAGNA - COSA NOSTRA BID-RIGGING CONSPIRACY"
  )
bf$model <- factor(bf$model, levels = names(labels))

# set order of the model parameters
bf$effect <- factor(bf$effect, levels = c("COMMON FRIENDS", "TRIADIC CLOSURE", "ASSORTATIVITY", "CENTRALIZATION"))

# don't run
# save
# saveRDS(bf, file = "~/Desktop/super/BF.rds")

# summarize statistics for boxplot
bf_box <- bf %>%
  dplyr::group_by(effect) %>%
  dplyr::summarise(
    ymin = stats::quantile(BF, 0.025),
    lower = stats::quantile(BF, 0.25),
    middle = mean(BF),           # substitute for median if you want
    upper = stats::quantile(BF, 0.75),
    ymax = stats::quantile(BF, 0.975)
    )


# plot the box and whisker plot
fig3 <- ggplot2::ggplot(bf, ggplot2::aes(x = effect, y = BF, color = evidence)) +
  ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 1) +
  ggplot2::geom_boxplot(
    data = bf_box,
    ggplot2::aes(
      x = effect,
      ymin = ymin,
      lower = lower,
      middle = middle,
      upper = upper,
      ymax = ymax
      ),
    stat = "identity",
    fill = NA,
    color = "black",
    width = 0.5,
    inherit.aes = FALSE
    ) +
  ggplot2::geom_errorbar(
    data = bf_box,
    ggplot2::aes(x = effect, ymin = upper, ymax = ymax),
    width = 0.15,
    color = "black",
    inherit.aes = FALSE
    ) +
  ggplot2::geom_errorbar(
    data = bf_box,
    ggplot2::aes(x = effect, ymin = ymin, ymax = lower),
    width = 0.15,
    color = "black",
    inherit.aes = FALSE
    ) +
  ggplot2::scale_y_log10(
    limits = c(1, 10),
    breaks = c(1, 3, 10, 30, 100, 150),
    labels = c("1", "3", "10", "30", "100", "150")
    ) +
  ggplot2::scale_color_manual(values = c(
    "NULL" = "#f0f0f0", "WEAK" = "#bdbdbd", "MODERATE" = "#969696",
    "STRONG" = "#636363", "VERY STRONG" = "#252525"
     )) +
  ggplot2::labs(
    x = "MODEL PARAMETERS",
    y = "BAYES FACTORS",
    color = "Evidence Strength"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(
      # family = "serif",
      size = 6, 
      face = "plain"
    ),
    axis.title.x = ggplot2::element_text(color = "black", size = 8),
    axis.text.x = ggplot2::element_text(color = "black", size = 6, hjust = 0.5, vjust = 0.5, face = "plain", angle = 45),
    axis.title.y = ggplot2::element_text(color = "black", size = 8),
    axis.text.y = ggplot2::element_text(color = "black", size = 6, hjust = 1.0, vjust = 0.0, face = "plain", angle =  0),
    panel.grid.major.y = ggplot2::element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
    panel.grid.major.x = ggplot2::element_line(color = "grey80", linetype = "dotted")
    )


# plot the box and whisker plot
fig3 <- ggplot2::ggplot(bf, ggplot2::aes(x = effect, y = BF, color = evidence)) +
  ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 1) +  # dot plot
  ggplot2::geom_boxplot(ggplot2::aes(group = effect), fill = NA, color = "black", width = 0.5, outlier.shape = NA) +  # boxplot overlay
  # ggplot2::facet_wrap(~model, nrow = 4, ncol = 2, labeller = ggplot2::labeller(model = labels)) + 
  ggplot2::scale_y_log10(
    limits = c(1, 10),
    breaks = c(1, 3, 10, 30, 100, 150), 
    labels = c("1", "3", "10", "30", "100", "150")
    ) +
  # greyscale palette (dark for stronger evidence)
  ggplot2::scale_color_manual(values = c(
    "NULL" = "#f0f0f0",
    "WEAK" = "#bdbdbd",
    "MODERATE" = "#969696",
    "STRONG" = "#636363",
    "VERY STRONG" = "#252525"
     )) +
  ggplot2::labs(
    # title = "Bayes Factors by Network Effect (Log Scale)",
    x = "MODEL PARAMETERS",
    y = "BAYES FACTORS", # y = expression("BAYES FACTOR (log"[10]*")"),
    color = "Evidence Strength"
    ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(
      # family = "serif",
      size = 6, 
      face = "plain"
      ),
    axis.title.x = ggplot2::element_text(color = "black", size = 8),
    axis.text.x = ggplot2::element_text(color = "black", size = 6, hjust = 0.5, vjust = 0.5, face = "plain", angle = 45),
    axis.title.y = ggplot2::element_text(color = "black", size = 8),
    axis.text.y = ggplot2::element_text(color = "black", size = 6, hjust = 1.0, vjust = 0.0, face = "plain", angle =  0),
    panel.grid.major.y = ggplot2::element_line(color = "grey80", linetype = "dotted"),
    panel.grid.minor.y = ggplot2::element_line(color = "grey90", linetype = "dotted"),
    panel.grid.major.x = ggplot2::element_line(color = "grey80", linetype = "dotted")
    )

# output high resolution images
output <- function(plot, filename){
  ggplot2::ggsave(
    filename,
    plot,
    path = "~/Desktop/super/", 
    width = 5, 
    height = 5, 
    device = 'png', 
    dpi = 700
    )
}
output(plot = fig3, filename = "fig3.png")



# close . R script






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
    data$stderr.ts
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
  prior.tau <- function(t){bayesmeta::dhalfcauchy(t, scale = 0.5)}
  
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