Results for generalized potentials of the form: 
U(r) = 4* eps * [(sig/r)**m - (sig/r)**n]. 

Python and shell script files to generate the results for arbitraty potentials are included. 

1. First calculate isotherms for a range of temperatures on both liquid and vapor sides: 

    This needs an estimate of the Tc (or eps_critical) - this is obtained from running a isochore calculation at rho = 0.5 from low to high epsilon (or high to low Temperature) and identifying where the HNC solver first fails to converge (i.e., we enter the two phase region). This is an estimate as the critical temperature as the eps spacing for this is fairly coarse (0.05) and the critical density is assumed to be 0.50. 

2. Then calculate the spinodal by running the second script file - extrapolation of S(k=0) by fitting to a power law divergence. 



