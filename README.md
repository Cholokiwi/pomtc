# Model-based clustering for longitudinal ordinal data #


### Computer scripts to reproduce simulation in Costilla et al 2018 ###
This repository contains the R and C++ binary files to reproduce results presented in Table 2 (section 3.4 Model validation using simulated data). C++ source files are also included for convenience. Scripts run in Linux and Windows 64 bits versions (x86-64).

A brief description of these files can be found below. 

**R scripts**

+ [pomtc.sim.r](https://github.com/Cholokiwi/pomtc/blob/master/pomtc.sim.r)
Simulates data, estimates the model using a Metropolis-Hastings sampler, and relabels MCMC chains.

+ [functions.pomtc.r](https://github.com/Cholokiwi/pomtc/blob/master/functions.pomtc.r)
Functions to perform all above taks, including loading C++ binaries.


**C++ files**

*Binaries, x86-64 versions*

+ Linux: [theta.pomtc.so](https://github.com/Cholokiwi/pomtc/blob/master/theta.pomtc.so), [Zrc_tt.so](https://github.com/Cholokiwi/pomtc/blob/master/Zrc_tt.so)

+ Windows: [theta.pomtc.dll](https://github.com/Cholokiwi/pomtc/blob/master/theta.pomtc.dll), [Zrc_tt.dll](https://github.com/Cholokiwi/pomtc/blob/master/Zrc_tt.dll)


*Source files*

+ [theta.pomtc.cpp](https://github.com/Cholokiwi/pomtc/blob/master/theta.pomtc.cpp). Calculates cell probabilities given parameters mu, alpha, beta, and gamma.

+ [Zrc_tt.cpp](https://github.com/Cholokiwi/pomtc/blob/master/Zrc_tt.cpp). Numerator for z_ig in logs. Used to compute likelihood as a function of theta, pi, and the data.

  
  
### Running instructions ###

To run the simulation please:

1. [Download](https://github.com/Cholokiwi/pomtc/archive/master.zip) or clone the [git repository](https://github.com/Cholokiwi/pomtc). Uncompress it to your local computer.
2. Run "pomtc.sim.r" within R. Note that you need to access to several R libraries.
3. Table 2 contents will be saved as a csv file. 

For the simulation in the paper (n=1000, p=15, q=5, R=3) programs take about an hour to run using R 3.3.3 in a Xeon E5-2680 2.50GHz CPU. Depending on your computer specifications this time might vary. Running time includes simulating the data, estimating the model using 3 MCMC chains and relabeling. In addition to that, traceplots for the original and relabelled chains are also produced and the R session saved. The complete output is available [here](https://rawgit.com/Cholokiwi/pomtc/master/pomtc.sim.html) in R markdown format.


#### References ####
Costilla, Liu, Arnold, Fernandez (2018). A Bayesian model-based approach to estimate latent groups in longitudinal ordinal data (submitted).

### Comments/questions to ###
Roy Costilla

PhD Statistics

Institute for Molecular Biosciences. Univerity of Queensland

r.costilla@imb.uq.edu.au

https://www.researchgate.net/profile/Roy_Costilla

https://twitter.com/CmRoycostilla

