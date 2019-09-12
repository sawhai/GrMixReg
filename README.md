# GrMixReg
## Group Mixture of Regressions

This is an extension to the regular Finite Mixture of Regressions (FMR), which allows for known group structure among observations, in addition to their posibly unknown heterogeneity. The algorithm assigns groups (instead of individual observations) to clusters. Expectation Maximization (EM) is implemented and applied for this setup. In addition, a Maximum a Posterior (MAP) prediction density is developed. 

The code is to run a simulation study used in the paper.


plot_predict_comparison.R is the code for plotting the comparison between MMCL++ and GMR.

The .RData files are the results saved and used for plotting.
