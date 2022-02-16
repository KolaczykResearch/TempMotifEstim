library(Rcpp) 
library('ggplot2')
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# or Additional modules to load field, type gcc/8.3.0
sourceCpp("functions.cpp") 

# Obtain data ------------------------------------------------
data = read.table('../Data/CollegeMsg.txt', header = FALSE)
data = as.matrix(data)

# Specify motif ---------------------------------------------
motif1 = matrix(c(0, 1, 1,
                  0, 1, 2,
                  1, 0, 3), 3, 3, byrow=TRUE)
motif2 = matrix(c(0, 1, 1,
                  1, 2, 2,
                  2, 0, 3), 3, 3, byrow=TRUE)
motif3 = matrix(c(0, 1, 1,
                  0, 2, 2,
                  3, 1, 3,
                  3, 2, 4), 4, 3, byrow=TRUE)
delta = 86400
motifs = list()
motifs[[1]] = motif1
motifs[[2]] = motif2
motifs[[3]] = motif3

# Run multiple replicates to assess CI ------------------------------------------
onerun = function(data, motif, delta, p, C_true){
  # cat('Data length:', nrow(data), '\n')
  seed = sample(1e9, 1)
  output = estimate_motif_counts_var(data, motif, delta, p, seed=seed)
  C_hat = output[1]
  sig_2_hat = output[2]
  low = C_hat-sqrt(sig_2_hat)*qnorm(0.975)
  up = C_hat+sqrt(sig_2_hat)*qnorm(0.975)
  bool = (C_true<=up & C_true>=low)*1
  return (bool)
}

#--------------
# RUN ALL
#--------------
# For motifs = c(type1, type2, type3); p=c(0.01, 0.03, 0.05, 0.1, 0.2)
# true counts: 381720, 9850, 271022
system.time({
set.seed(5)
n_rep = 5000
results = NULL
for (i in 1:length(motifs)){
  cat('---------------------\n')
  cat('For motif type', i, '\n')
  cat('---------------------\n')
  
  motif = motifs[[i]]
  motif_type = i
  C_true = estimate_motif_counts_var(data, motif, delta, p=1, seed=22)[1]
  for (p in c(0.01, 0.03, 0.05, 0.1, 0.2)){
    ptm = proc.time()
    cat('--> p = ', p, '\n')
    ss = replicate(n_rep, onerun(data, motif, delta, p, C_true))
    row = c(motif_type, p, mean(ss), sd(ss))
    cat(row, '\n')
    results = rbind(results, row)
    time = proc.time() - ptm
    cat('time:', time, '\n')
  }
  write.csv(results, file = "results.csv", row.names = FALSE)
}

})
results
# results = read.csv('results.csv', header=TRUE)

# Visualization ========================
colnames(results) = c('Motif', 'p', 'RF', 'se')
rownames(results) = NULL
results_plot = data.frame(results)
results_plot$Motif = as.factor(results_plot$Motif)

plt1 = ggplot(results_plot, aes(x=p, y=RF, group=Motif)) +
  geom_hline(yintercept=0.95, linetype=5, color = "grey") +
  geom_line(aes(color=Motif), linetype=2) +
  geom_point(aes(color=Motif, shape=Motif), size=3) + 
  labs(x='Sampling ratio p', y='Estimated coverage probability',
       color = "Motif H", shape = "Motif H") +
  theme_bw()
plt1


### --------------------
### Continue with more motifs
### ---------------------
# Obtain data ------------------------------------------------
data = read.table('../Data/CollegeMsg.txt', header = FALSE)
data = as.matrix(data)

# Specify motif ---------------------------------------------
motif1 = matrix(c(0, 1, 1,
                  1, 2, 2,
                  0, 2, 3), 3, 3, byrow=TRUE)
motif2 = matrix(c(0, 1, 1,
                  0, 1, 2,
                  2, 1, 3), 3, 3, byrow=TRUE)
motif3 = matrix(c(0, 1, 1,
                  2, 1, 2,
                  1, 0, 3), 3, 3, byrow=TRUE)
delta = 86400
motifs = list()
motifs[[1]] = motif1
motifs[[2]] = motif2
motifs[[3]] = motif3

# Run multiple replicates to assess CI ------------------------------------------
onerun = function(data, motif, delta, p, C_true){
  # cat('Data length:', nrow(data), '\n')
  seed = sample(1e9, 1)
  output = estimate_motif_counts_var(data, motif, delta, p, seed=seed)
  C_hat = output[1]
  sig_2_hat = output[2]
  low = C_hat-sqrt(sig_2_hat)*qnorm(0.975)
  up = C_hat+sqrt(sig_2_hat)*qnorm(0.975)
  bool = (C_true<=up & C_true>=low)*1
  return (bool)
}

#--------------
# RUN ALL
#--------------
# For motifs = c(type1, type2, type3); p=c(0.01, 0.03, 0.05, 0.1, 0.2)
# true counts: 381720, 9850, 271022, 
# 16064, 1201092, 295970
system.time({
  set.seed(50)
  n_rep = 5000
  results = NULL
  for (i in 1:length(motifs)){
    cat('---------------------\n')
    cat('For motif type', i, '\n')
    cat('---------------------\n')
    
    motif = motifs[[i]]
    motif_type = i
    C_true = estimate_motif_counts_var(data, motif, delta, p=1, seed=22)[1]
    for (p in c(0.01, 0.03, 0.05, 0.1, 0.2)){
      ptm = proc.time()
      cat('--> p = ', p, '\n')
      ss = replicate(n_rep, onerun(data, motif, delta, p, C_true))
      row = c(motif_type, p, mean(ss), sd(ss))
      cat(row, '\n')
      results = rbind(results, row)
      time = proc.time() - ptm
      cat('time:', time, '\n')
    }
    write.csv(results, file = "results2.csv", row.names = FALSE)
  }
  
})
results
# results = read.csv('results2.csv', header=TRUE)

# Visualization ========================
colnames(results) = c('Motif', 'p', 'RF', 'se')
rownames(results) = NULL
results_plot = data.frame(results)
results_plot$Motif = as.factor(results_plot$Motif)

plt1 = ggplot(results_plot, aes(x=p, y=RF, group=Motif)) +
  geom_hline(yintercept=0.95, linetype=5, color = "grey") +
  geom_line(aes(color=Motif), linetype=2) +
  geom_point(aes(color=Motif, shape=Motif), size=3) + 
  labs(x='Sampling ratio p', y='Estimated coverage probability',
       color = "Motif H", shape = "Motif H") +
  theme_bw()
plt1

###------------------------------
### Combine both results
###-----------------------------
results1 = read.csv('results.csv', header=TRUE)
results2 = read.csv('results2.csv', header=TRUE)
results1[,1]=rep(c(1,4,6), each=5)
results2[,1]=rep(c(5,2,3), each=5)
results = rbind(results1, results2)
results[,4] = sqrt(results[,3]*(1-results[,3]))/sqrt(5000)

colnames(results) = c('Motif', 'p', 'RF', 'se')
rownames(results) = NULL
results_plot = data.frame(results)
results_plot$Motif = as.factor(results_plot$Motif)

plt1 = ggplot(results_plot, aes(x=p, y=RF, group=Motif)) +
  geom_hline(yintercept=0.95, linetype=5, color = "grey") +
  # geom_errorbar(aes(ymin=RF-se, ymax=RF+se), width=.005) +
  geom_line(aes(color=Motif), linetype=2) +
  geom_point(aes(color=Motif, shape=Motif), size=3) + 
  labs(x='Sampling ratio p', y='Estimated coverage probability',
       color = "Motif H", shape = "Motif H") +
  theme_bw()
plt1
