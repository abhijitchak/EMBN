# ERG
Fortarn codes for Markov chain monte carlo estimation of exponential random graph models
These are Fortarn codes for Markov chain monte carlo estimation of exponential random graph (ERG) models as proposed by Snijders, T. A. in Journal of Social Structure, 3, 1-40 (2002).
These codes are for Bipartite network (Bank-firm network). 
It estimates the ERG parameter $theta$'s for (1)  Bernoulli ERG model and (2)  Two star ERG model.

Command to compile and execute: 

 (1) Bernouli model: 
 gfortran -C -o Bernoulli BN_MCMC_Bernoulli.f 
 
 ./Bernoulli 
 
 (2) Two star model:
 gfortran -C -o twostar BN_MCMC_2Star.f
 
 ./twostar
