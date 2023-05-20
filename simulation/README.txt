This folder consists of the R codes for the algorithm used to perform the simulation studies in the manuscript. You will need the following R files to run the simulations. In this README file, I will explain each R file.

###########################################################

1. simple_DataGen.R: the R code used to generate data for the simulation studies presented in the manuscript.

2. initial_fit.R: the R code used to find an initial estimation for the parameters of the semi-parametric functional quantile regression (FQR) with one functional covariate. The output is used as the starting point of the algorithm for computing the estimations of the proposed CLoSE method.

3. CLoSE_fit.R: the R code used to perform the proposed CLoSE method for the FQR model with one functional covariate. The output includes the estimated B-spline coefficients for the unknown slope function, the estimated parameters for the parametric part of the semi-parametric FQR, the estimated null subregions and non-null subregions.

4. step2_fit: the R code used to perform the SQL method based on the output of CLoSE_fit.R. More specifically, this R code is used to perform SQL method on the processed data, for which the estimated null subregions of the functional covariates have been removed from the estimation procedure. 

5. sim.R: the R code used to conduct the simulation studies presented in the manuscript.
###########################################################