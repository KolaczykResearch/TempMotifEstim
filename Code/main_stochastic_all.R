library(Rcpp) 
library('ggplot2')
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# or Additional modules to load field, type gcc/8.3.0
sourceCpp("functions.cpp") 

# Simulate from Poisson process uniform model 
# generate a sequence of times (event) from poisson(lamba), then sample uniformly an edge
# return a 3-column matrix, each row contains a directed edge with time stamp (i, j, t)
sim_PPunif = function(n, lambda, TT){
  # simulate arrival times from Poisson processes within time interval [0,TT]
  event_times = cumsum(rexp(max(1, TT*lambda), lambda))
  while (tail(event_times, n=1)<TT){
    more_event_times = cumsum(rexp(max(1, TT*lambda*0.2), lambda))
    event_times = c(event_times, tail(event_times, n=1)+more_event_times)
  }
  event_times = event_times[event_times<=TT]
  m = length(event_times)
  data = matrix(0, m, 3)
  # uniformly sample edges for each time
  edges_list = c(1:n^2)
  edges_sampled = sample(edges_list, m, replace=TRUE)
  row_ids = edges_sampled%%n
  col_ids = edges_sampled%/%n+1
  data[,1] = row_ids
  data[,2] = col_ids
  data[,3] = event_times
  return(data)
}

# Simulate from Poisson process stochastic block model 
# return a 3-column matrix, each row contains a directed edge with time stamp (i, j, t)
sim_homoPPBM = function(n, TT, pi, lambda){
  data = NULL
  # simulate n(n-1) poisson processes within [0, TT]
  for(i in 1:n){
    for(j in 1:n){
      if (i!=j){
        # simulate a poisson process on edge (i, j) within [0, TT] 
        event_times = cumsum(rexp(max(1, TT*lambda[pi[i],pi[j]]), lambda[pi[i],pi[j]]))
        while (tail(event_times, n=1)<TT){
          more_event_times = cumsum(rexp(max(1, TT*lambda[pi[i],pi[j]]*0.2), lambda[pi[i],pi[j]]))
          event_times = c(event_times, tail(event_times, n=1)+more_event_times)
        }
        event_times = event_times[event_times<=TT]
        data_singe_edge = matrix(0, length(event_times), 3)
        data_singe_edge[,1] = rep(i, length(event_times)) 
        data_singe_edge[,2] = rep(j, length(event_times)) 
        data_singe_edge[,3] = event_times
        data = rbind(data, data_singe_edge)
      }
    }
  }
  data = data[order(data[,3]),]
  return(data)
}

###############################################################
#### Simulation Experiments ###################################
###############################################################
# -------------------------------------------------------
# - Stochastic case 1: Poisson process uniform model
# ------------------------------------------------------
# one run generate 1 data set from sim_PPunif(), then compute C and C_hat
onerun = function(n, motif, delta, p, lambda, TT){
  seed = sample(1e9, 1)
  # n: number of individuals (nodes)
  # TT: length of time
  # Get temporal edges, simulated from Poisson process uniform model 
  data = sim_PPunif(n, lambda, TT)

  C_true = estimate_motif_counts_var(data, motif, delta, p=1, seed=seed)[1]
  output = estimate_motif_counts_var(data, motif, delta, p, seed=seed)
  C_hat = output[1]
  sig_2_hat = output[2]
  low = C_hat-sqrt(sig_2_hat)*qnorm(0.975)
  up = C_hat+sqrt(sig_2_hat)*qnorm(0.975)
  bool = (C_true<=up & C_true>=low)*1
  # return (c(bool, C_hat, sqrt(sig_2_hat), C_true, low, up))
  return (C_hat/C_true)
}

# ------------
# Consistency - stoc case 1
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
delta = 2
p=0.03
lambda = 30
n = 100
n_rep = 100
TT_set = c(200, 300, 400, 600, 800, 1600, 3200, 5000)

set.seed(60)
results = NULL
for (TT in TT_set){
  cat('--> TT = ', TT, '\n')
  ss = replicate(n_rep, onerun(n, motif, delta, p, lambda, TT))
  row = c(TT, mean(ss), sd(ss))
  cat(row, '\n')
  results = rbind(results, row)
}

# plot
results = results[!is.na(results[,2]),] # remove NA
df = data.frame(results)
colnames(df) = c('tau', 'mean', 'se')
pd = position_dodge(0) # move them .05 to the left and right
p1_st = ggplot(df, aes(x=tau, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=50, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=1.5) + theme_bw() + 
  ylab(expression(hat(C)(H, T(tau))/C(H, T(tau)))) + 
  xlab(expression(tau)) +
  geom_hline(yintercept=1)
p1_st

# ------------
# CLT - stoc case 1
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
delta = 2
p=0.03
lambda = 30
n = 100
TT = 9000

system.time({
  n_rep = 5000
  set.seed(60)
  ss = replicate(n_rep, onerun(n, motif, delta, p, lambda, TT))
  hist(ss)
  hist(ss, breaks=seq(min(ss)*0, max(ss)*1.05, 0.04), probability = TRUE, 
       main='', xlab=expression(hat(C)(H, T(tau))/C(H, T(tau))))
  curve(dnorm(x, mean=mean(ss), sd=sd(ss)),
        col="darkblue", lwd=2, add=TRUE, yaxt="n")
})

# -------------------------------------------------------
# - Stochastic case 2: Poisson process SBM
# ------------------------------------------------------
# one run generate 1 data set from sim_homoPPBM(), then compute C and C_hat
onerun = function(n, motif, delta, p, TT, inten_off_diag, inten_diag){
  seed = sample(1e9, 1)
  # n: number of individuals (nodes)
  # TT: length of time
  pi = rep(c(1,2), each=n/2) # memberships
  
  # intensity functions (homogeneous)
  lambda = matrix(0, 2, 2)
  # off-diagonal intensity function
  lambda[1,2] = lambda[2,1] = inten_off_diag
  # diagonal intensity function
  lambda[1,1] = lambda[2,2] = inten_diag
  
  # Get temporal edges, simulated from Poisson process SBM 
  data = sim_homoPPBM(n, TT, pi, lambda)
  # cat('nrow:', nrow(data))
  C_true = estimate_motif_counts_var(data, motif, delta, p=1, seed=seed)[1]
  output = estimate_motif_counts_var(data, motif, delta, p, seed=seed)
  C_hat = output[1]
  sig_2_hat = output[2]
  low = C_hat-sqrt(sig_2_hat)*qnorm(0.975)
  up = C_hat+sqrt(sig_2_hat)*qnorm(0.975)
  bool = (C_true<=up & C_true>=low)*1
  return (C_hat/C_true)
}

# ------------
# Consistency - stoc case 2
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
delta = 2
p=0.03
n = 100
n_rep = 100
TT_set = c(5, 10, 15, 20, 30, 40, 80, 160)

inten_off_diag = 0.06
inten_diag = 0.02

system.time({
set.seed(60)
results = NULL
for (TT in TT_set){
  cat('--> TT = ', TT, '\n')
  ss = replicate(n_rep, onerun(n, motif, delta, p, TT, inten_off_diag, inten_diag))
  row = c(TT, mean(ss), sd(ss))
  cat(row, '\n')
  results = rbind(results, row)
}
})
# plot
results = results[!is.na(results[,2]),] # remove NA
df = data.frame(results)
colnames(df) = c('tau', 'mean', 'se')
pd = position_dodge(0) # move them .05 to the left and right
p3_st = ggplot(df, aes(x=tau, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=1, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=1.5) + theme_bw() + 
  ylab(expression(hat(C)(H, T(tau))/C(H, T(tau)))) + 
  xlab(expression(tau)) +
  geom_hline(yintercept=1)
p3_st

# ------------
# CLT - stoc case 2
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
delta = 2
p=0.03
n = 100
TT = 10
inten_off_diag = 0.06
inten_diag = 0.02

system.time({
  n_rep = 5000
  set.seed(60)
  ss = replicate(n_rep, onerun(n, motif, delta, p, TT, inten_off_diag, inten_diag))
  hist(ss)
  hist(ss, breaks=seq(min(ss)*0, max(ss)*1.2, 0.04), probability = TRUE, 
       main='', xlab=expression(hat(C)(H, T(tau))/C(H, T(tau))))
  curve(dnorm(x, mean=mean(ss), sd=sd(ss)),
        col="darkblue", lwd=2, add=TRUE, yaxt="n")
})

