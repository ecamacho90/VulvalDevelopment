# VulvalDevelopment
Github repository containing the code to simulate and fit the model of vulval development in C. elegans.
The repository contains two folders, each corresponding to a version of the model: **Type I (Model_v1)** and **Type II (Model_v2)**. The folders contain code to fit the models to the Training Data by using ABC SMC. All functions in the two folders are equivalent except the functions containing the priors and constraints for the parameters.

Given a model type, we can divide the code into two parts, a group of functions which, given a parameter vector, simulate the model and gives the percentage of times that each cell took each fate, and a group of functions which run ABC SMC. In the following, we will explain each group.

## Functions that simulate the model
Given a parameter vector, the model is simulated in two steps:

### Step 1
Linear Noise Approximation to determine the starting distribution around the initial state (the attractor corresponding to tertiary fate). This is done with function **lna_v10.m** which also uses functions: 
  1. **solution_det_eqtn_v10.m**: solves function **cusp_and_saddlenode_model_singlecell_v10.m**, which contains the model for an isolated cell without receiving any signals, i.e. s=l=0.
  2. **covariancemat_direct_solution_v10.m**: solves function **covariance_matrix_ode_v10.m**, using LNA one can find the covariance matrix around the deterministic solution, by solving **model_jacoft_v10.m**.

### Step 2
Given a sample of initial conditions obtained from Step 1, the function **simulationEulerACablation_Competence_v10_vec** simulates the model taking advantage of Euler-Maruyama algorithm. It uses the following functions:
  1. **cusp_and_saddlenode_model_v10_vec.m** which contains the differential equations of the tristable binary switch landscape for P4-6.p.
  2. **computefates_PostCompetence.m** which, given the end point 

## Functions for ABC SMC





