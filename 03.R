#  -----------------------------------------------------------------------------------

# file 03: estimate the hierarchical network model

# 'mlergm' package
# https://cran.r-project.org/web/packages/mlergm/mlergm.pdf

# don't run
# install.package('mlergm')

# ------------------------------------------------------------------------------------



# construct the supergraph as a 'mlergm' object
g_super <- mlergm::mlnet(
  network = g_super, 
  node_memb = network::get.vertex.attribute(g_super, "group")
  )
# plot the supergraph
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



# close .r script


