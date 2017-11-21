All files are numbered according to when they entered the paper. 

Sometimes, a,b,c, etc are used to indicate code chunks that refer to the same part. For example, the Monte Carlo was broken up into three chunks: one to set up the monte carlo for batch processing, one to actually batch it, one to  actually perform the calculations, one to aggregate the results, and one to report the results. 

0.parameters_setup - this file sets up the parameters (e.g., means, variances, iterations) used in the demo. I call these values both from the simulation and the paper, to avoid errors
1. demonstration_Only_CseIII - this file produces the output in Figure 1
2. monte_carlo_run_runif - this file performs a preliminary simulation where each parameter is varied using a random uniform distribution. The purpose is simply to determine which variables affect bias. 
3 (a-d) These files perform the Monte Carlo for the results section
4. report results caseiii - this file reports the results of the monte carlo