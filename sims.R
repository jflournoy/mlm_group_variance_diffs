.libPaths('/R_LIB/')
library(lme4)
library(lmerTest)
library(parallel)

#generate data in which there is variance in intercepts and slopes and a
#difference in mean intercept and slope between groups. Optionally, create
#different variance of REs for the different groups.

sim <- function(delta_b0 = 0, delta_b1 = 0, sd_b0_g1 = 1, delta_sd_b0 = 1, sd_b1_g1 = .2, delta_sd_b1 = .4, J_grp = 250, N_per_j = 3){
  RE_coef <- function(mean, sd, N){
    r <- Map(function(m, s, n){
      rnorm(n = n, mean = m, sd = s)
    }, mean, sd, N)
    return(unlist(r))
  }
  
  b0 <- RE_coef(mean = c(0, delta_b0), sd = c(sd_b0_g1, sd_b0_g1 + delta_sd_b0), N = c(J_grp, J_grp))
  b1 <- RE_coef(mean = c(0, delta_b1), sd = c(sd_b1_g1, sd_b1_g1 + delta_sd_b1), N = c(J_grp, J_grp))
  N_total <- J_grp * 2 * N_per_j
  X <- split(rnorm(N_total), 1:length(b0))
  
  y <- Map(function(b0_j, b1_j, X_j){
    rnorm(n = N_per_j, mean = b0_j + b1_j * X_j, sd = 1)
  }, b0, b1, X)
  
  d <- data.frame(x = unlist(X), y = unlist(y), 
                  id = rep(1:length(b0), each = N_per_j),
                  group = rep(c('control', 'treat'), each = J_grp * N_per_j))
  
  m <- lmer(y ~ 1 + x * group + (1 + x | id), data = d)
  
  return(coef(summary(m))[, 'Pr(>|t|)'])
}

NCPU <- detectCores() - 1
cl <- makeForkCluster(NCPU)  
# get library support needed to run the code
nada <- clusterEvalQ(cl, {library(lme4); library(lmerTest)})
# put objects in place that might be needed for the code
clusterExport(cl, c('sim'))
# Set a different seed on each member of the cluster (just in case)
clusterSetRNGStream(cl)
asim <- parallel::parSapply(cl = cl, X = 1:5000*NCPU, FUN = function(i) sim())
stopCluster(cl)

apply(asim, 1, function(x) sum(x < .05) / length(x))
