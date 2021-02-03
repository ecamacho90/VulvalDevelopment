# ABC fitting of the binary flip with cusp landscape for vulval development
## Elena Camacho-Aguilar, Aryeh Warmflash, David A. Rand

This respository contains the code associated with the following manuscript:

> [Camacho-Aguilar E, Warmflash A, Rand DA. "Quantifying cell transitions in C. elegans with data-fitted landscape models." bioRxiv (2021).](https://www.biorxiv.org/content/10.1101/2021.01.22.426019v1#ref-5)

The repository contains two folders, each corresponding to a version of the model: **Type I (Model_v1)** and **Type II (Model_v2)**. The folders contain code to fit the models to the Training Data by using ABC SMC. All functions in the two folders are equivalent except the functions containing the priors and constraints for the parameters.

Given a model type, we can divide the code into two parts, a group of functions which, given a parameter vector, simulate the model and gives the percentage of times that each cell took each fate, and a group of functions which run ABC SMC. In the following, we will explain each group.

## Functions that simulate the model
Given a parameter vector, the model is simulated in two steps:

### Step 1
Linear Noise Approximation to determine the starting distribution around the initial state (the attractor corresponding to tertiary fate). This is done with function **lna_v10.m** which also uses functions: 
  1. **solution_det_eqtn_v10.m**: solves function **cusp_and_saddlenode_model_singlecell_v10.m**, which contains the model for an isolated cell without receiving any signals, i.e. s=l=0.
  2. **covariancemat_direct_solution_v10.m**: solves function **covariance_matrix_ode_v10.m**, using LNA one can find the covariance matrix around the deterministic solution, by solving **model_jacoft_v10.m**.

### Step 2
Given a sample of initial conditions obtained from Step 1, the function **simulationEulerACablation_Competence_v10_vec.m** simulates the model taking advantage of Euler-Maruyama algorithm. It uses the following functions:
  1. **cusp_and_saddlenode_model_v10_vec.m** which contains the differential equations of the tristable binary switch landscape for P4-6.p.
  2. **computefates_PostCompetence.m** which, given the end point of the trajectories of the cells, quantifies the basin of attraction at which they are located in order to score the fates. This function calls the following functions:
      1. **findfate_1_2_deg_PostCompetence.m**: In a landscape in which fates 1 or 2 are degenerate (the discriminant is very close to 0), given the coordinates of the end point of the trajectory of the cell, it looks for the attractor towards which the trajectory will converge in case we take that point as initial condition. It calls functions singlecell_cusp_and_saddlenode_model.m and myEventsFcnDegbmin.m to find the basin of attraction in which the end point of the trajectory is. 
      2. **findfate_1_2_nondeg_PostCompetence.m**: In a landscape in which fates 1 or 2 are not degenerate (the discriminant is not 0), given the coordinates of the end point of the trajectory of the cell, it looks for the attractor towards which the trajectory will converge in case we take that point as initial condition. It calls functions singlecell_cusp_and_saddlenode_model.m and myEventsFcnNonDeg.m to find the basin of attraction in which the end point of the trajectory is. 

## Functions for ABC SMC
Given a threshold e, the number of particles N, the data to fit and the parameters to be fitted, the function **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L.m** (local version) (or **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance.m** (parallel version) runs ABC SMC to find the N particles which simulations are e-close to the data. The details regarding the threshold, the number of particles, parameters etc is set with the function **Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL.m** (or **Call_Parallel_function_AbsDist_Modelv10_v1.m**); in this function, loopept controls the step of the algorithm, and EpT vector and names cell need to be updated accordingly every time one is running a new step of the algorithm.
The functions **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L.m** (or **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance.m**), in turn, call the following functions:
  1. **Vulval_Development_Modelv10_v1_Priors.m** that contains the priors for the parameters and samples from them.
  2. **Vulval_Development_Modelv10_v1_EvalPriors.m** that given a value for a parameter checks the probability of being sampled from its prior.
  2. **Vulval_Development_Modelv10_v1_Constraints.m** that checks constraints on the parameter values.
  3. **Vulval_Development_Modelv10_v1_Relations_Constraints1.m** that checks relationships between the parameters.
  4. **findm22m32.m** finds m22 and m33 as functions of the other matrix values.
  5. **Vulval_Development_Modelv10_v1_Relations_Constraints2.m** checks final constraints.
  6. **equilibria_fates_1_2.m** computes the values of the critical points on the x-axis.
  7. Calls **lna_v10.m** to find the initial conditions.
  8. **Vulval_Development_Modelv10_v1_AbsDistance_Fates.m** that simulates each mutant (using the function **simulationEulerACablation_Competence_v10_vec.m**) and computes the distance between the simulation and the data for each of the mutants being fitted. 
  9. **Vulval_Development_Modelv10_v1_AbsDistance_EvalKernels.m** evaluates the density function kernel of a perturbed particle to find the weight of the new particle if accepted.
  
In order to run ABC, one needs to set up the function **Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL.m** (or **Call_Parallel_function_AbsDist_Modelv10_v1.m**) which the settings for the particular step of the algorithm and then run **Parallel_function_AbsDist_Modelv10_v1_LOCAL.m**, which will compute the local covariance matrices for each particle (if the step of the algorithm is bigger than one) with the function **Vulval_Development_Modelv10_v1_CovarianceMatricesOLCM.m** and then will call **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L.m** (local version) (or **ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance.m** if parallelised). 

The data will be saved in __.mat__ format, one for each parallel job. If one wants to combine all the subfiles from the parallel jobs into one, then use function **SaveData.m**. 
  




