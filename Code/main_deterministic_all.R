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
# - Deterministic case 1: Poisson process uniform model
# ------------------------------------------------------
# Get data from sim_homoPPBM()
get_data = function(n, TT, m, lambda){
  # n: number of individuals (nodes)
  # TT: length of time
  # Get temporal edges, simulated from Poisson process uniform model
  data = sim_PPunif(n, lambda, TT)
  # cat('length:', nrow(data), '\n')
  data = data[1:m, ]
  return(data)
}

# Estimate motifs on the generated data
onerun = function(data, motif, delta, p){
  seed = sample(1e9, 1)
  C_hat = estimate_motif_counts_var(data, motif, delta, p, seed=seed)[1]
  C = estimate_motif_counts_var(data, motif, delta, p=1, seed=seed)[1]
  return (C_hat/C)
}

# data = sim_PPunif(n=100, lambda=30, TT=250)
# cat('length:', nrow(data), '\n')

# ------------
# Consistency - det case 1
# ------------
# Define motif: directed cyclic triangles 
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
n=100; m=7000; TT=280
delta = 2
p=0.03
n_rep = 100
lambda_set = c(25, 30, 35, 50, 70, 100, 150, 250)

set.seed(36)
results = NULL
for (lambda in lambda_set){
  ptm = proc.time()
  cat('--> lambda = ', lambda, '\n')
  data = get_data(n, TT, m, lambda)
  ss = replicate(n_rep, onerun(data, motif, delta, p))
  row = c(lambda, mean(ss), sd(ss))
  cat(row, '\n')
  results = rbind(results, row)
  time = proc.time() - ptm
  # cat('time:', time, '\n')
}

# Con Plot for deterministic case 1: PP Uniform
df = data.frame(results)
colnames(df) = c('lambda', 'mean', 'se')
pd = position_dodge(0) # move them .05 to the left and right
p1 = ggplot(df, aes(x=lambda, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=2.5, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=1.5) + theme_bw() + 
  ylab(expression(hat(C)(H, T[m])/C(H, T[m]))) + 
  xlab(expression(lambda)) + 
  geom_hline(yintercept=1)
p1

# ------------
# CLT - det case 1
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
# n=100; m=7000; TT=30
n=100; m=200000; TT=820
delta = 2
p=0.03
n_rep = 5000
lambda = 250

# CLT Plot for deterministic case 1: PP Uniform
system.time({
set.seed(5)
data = get_data(n, TT, m, lambda)
ss = replicate(n_rep, onerun(data, motif, delta, p))
})
hist(ss, breaks=seq(min(ss)*0.9, max(ss)*1.1, 0.01), probability = TRUE, 
     main='', xlab=expression(hat(C)(H, T[m])/C(H, T[m])), ylim=c(0,12))
curve(dnorm(x, mean=mean(ss), sd=sd(ss)),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

# -------------------------------------------------------
# - Deterministic case 2: Poisson process SBM
# ------------------------------------------------------
# Get data from sim_homoPPBM()
get_data = function(n, TT, m, inten_off_diag=1, inten_diag){
  # n = 10 # number of individuals (nodes)
  # TT = 50 # length of time
  pi = rep(c(1,2), each=n/2) # memberships
  
  # intensity functions (homogeneous)
  lambda = matrix(0, 2, 2)
  # off-diagonal intensity function
  lambda[1,2] = lambda[2,1] = inten_off_diag
  # diagonal intensity function
  lambda[1,1] = lambda[2,2] = inten_diag
  
  # Get temporal edges, simulated from Poisson process stochastic block model 
  data = sim_homoPPBM(n, TT, pi, lambda)
  cat('length:', nrow(data), '\n')
  data = data[1:m, ]
  return(data)
}

# Estimate motifs on the generated data
onerun = function(data, motif, delta, p){
  # cat('Data length:', nrow(data), '\n')
  seed = sample(1e9, 1)
  C_hat = estimate_motif_counts_var(data, motif, delta, p, seed=seed)[1]
  C = estimate_motif_counts_var(data, motif, delta, p=1, seed=seed)[1]
  return (C_hat/C)
}

# ------------
# Consistency - det case 2
# ------------
# Define motif: directed cyclic triangles 
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
n=100; TT=25; m=7000
delta = 2
p=0.03
n_rep = 100
inten_diag_set = c(0.001, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

set.seed(60)
results = NULL
for (inten_diag in inten_diag_set){
  ptm = proc.time()
  cat('--> inten_diag = ', inten_diag, '\n')
  data = get_data(n, TT, m, inten_off_diag=0.06, inten_diag=inten_diag)
  ss = replicate(n_rep, onerun(data, motif, delta, p))
  row = c(inten_diag, mean(ss), sd(ss))
  cat(row, '\n')
  results = rbind(results, row)
  time = proc.time() - ptm
  cat('time:', time, '\n')
}

# Con Plot for deterministic case 2: PP SBM
df = data.frame(results)
colnames(df) = c('diagonal_intensity', 'mean', 'se')
pd = position_dodge(0) # move them .05 to the left and right
p3 = ggplot(df, aes(x=diagonal_intensity, y=mean)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.002, position=pd) +
  # geom_line(position=pd) +
  geom_point(position=pd, size=1.5) + theme_bw() + 
  ylab(expression(hat(C)(H, T[m])/C(H, T[m]))) + 
  xlab('Diagonal intensity') + 
  geom_hline(yintercept=1)
p3

# ------------
# CLT - det case 2
# ------------
motif = matrix(c(0, 1, 1,
                 1, 2, 2,
                 2, 0, 3), 3, 3, byrow=TRUE)
n=100; TT=7; m=7000 
delta = 2
p=0.03
n_rep = 5000
inten_diag_set = 0.2

# CLT Plot for deterministic case 2: PP SBM
ptm = proc.time()
set.seed(5)
data = get_data(n, TT, m, inten_off_diag=0.06, inten_diag=0.2)
ss2 = replicate(n_rep, onerun(data, motif, delta, p))
time = proc.time() - ptm
cat('time:', time, '\n') # 302 2s

hist(ss2, breaks=seq(min(ss2)*0.7, max(ss2)*1.1, 0.02), probability = TRUE, 
     main='', xlab=expression(hat(C)(H, T[m])/C(H, T[m])), ylim=c(0,5))
curve(dnorm(x, mean=mean(ss2), sd=sd(ss2)),
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
