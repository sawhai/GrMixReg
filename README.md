# GrMixReg
Group Mixture of Regressions

The source code to reproduce the simulation study for GMR.

GMR.K2.d2.R is the simulation code for the case k = 2 and dimension = 2.
GMR.K4.d4.R is the simulation code for the case k = 4 and dimension = 4.
GMR_data_gen.R is the code for generating the data (called by GMR.K2.d2.R and GMR.K4.d4.R)
NMI.AMI.calc.fn.R is the code for NMI and AMI calculation (called by GMR.K2.d2.R and GMR.K4.d4.R)
fit_GMR.R is the code to fit GMR (called by GMR.K2.d2.R and GMR.K4.d4.R)
plot_Gr_EM2.R is the code for plotting simulation result.
plot_predict_comparison.R is the code for plotting the comparison between MMCL++ and GMR.
The .RData files are the results saved and used for plotting.
