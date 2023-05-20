This folder consists of the R codes for the algorithm used to perform the real data analysis and the R files used to generate figures related to the real data analysis in the manuscript.

The cleaned soybean data are saved in three RDS files: "soy_x.RDS", "soy_y.RDS" and "soy_add.RDS". "soy_x.RDS" contains functional covariates observations: the daily minimum temperature and daily maximum temperature observations. "soy_y.RDS" contains the annual soybean yields of different counties of different years. "soy_add.RDS" contains the annual precipitations and proportions of irrigated area of each county in different years.

The estimations obtained from the proposed CLoSE method and the proposed wild bootstrap procedure are saved in the zipped folders "0.25-quantile.zip", "0.5-quantile.zip" and "0.75-quantile.zip" for different quantiles respectively.

In this README file, I will explain each R file.
###########################################################
1. initial_fit_2fvars.R: the R code used to find an initial estimation for the parameters of the semi-parametric functional quantile regression (FQR) with two functional covariates. The output is used as the starting point of the algorithm for computing the estimations of the proposed CLoSE method.

2. CLoSE_fit_2fvars.R: the R code used to perform the proposed CLoSE method for the semi-parametric FQR with two functional covariates. The output includes the estimated B-spline coefficients for the unknown slope functions, the estimated parameters for the parametric part of the semi-parametric FQR, the estimated null subregions and non-null subregions of each slope function.

3. step2_fit_2fvars: the R code used to perform the SQL method based on the output of CLoSE_fit_2fvars.R. More specifically, this R code is used to perform SQL method on the processed data, for which the estimated null subregions of the functional covariates have been removed from the estimation procedure. 

4. soybean_fit.R: the R code used to perform the CLoSE method on the soybean data set introduced in the manuscript.

5. soybean_summary.R: the R code used to generate Figure 3 and Figure 4 in the manuscript.
###########################################################