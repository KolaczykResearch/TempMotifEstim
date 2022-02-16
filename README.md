# TempMotifEstim 
This repository contains code supporting the article ''Quantifying Uncertainty for Temporal Motif Estimation in Graph Streams under Sampling''. 

## Data 
We use CollegeMsg temporal network data set, which is available [here](https://snap.stanford.edu/data/CollegeMsg.html). 

## Code
The corresponding code for repoducing the numerical results in the paper is provided in `Code/`. The code is built upon [this implementation](https://github.com/jingjing-cs/Temporal-Motif-Counting) of the temporal motif sampling and counting algorithm. 

* `Code/main_deterministic_all.R`: code for Figure 1 and 2.
* `Code/main_stochastic_all.R/`: code for Figure 3 and 4. 
* `Code/main_app_CI_fig.R/`: code for Figure 6. 

For issues encountered in running the code, please contact Xiaojing Zhu at xiaojzhu@bu.edu
