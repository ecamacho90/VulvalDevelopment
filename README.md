# VulvalDevelopment
Github repository containing the code to simulate and fit the model of vulval development in C. elegans.
The repository contains two folders, each corresponding to a version of the model: Type I (Model_v1) and Type II (Model_v2). The folders contain code to fit the models to the Training Data by using ABC SMC. All functions in the two folders are equivalent except the functions containing the priors and constraints for the parameters.

Given a model type, we can divide the code into two parts, a group of functions which, given a parameter vector, simulate the model and gives the percentage of times that each cell took each fate, and a group of functions which run ABC SMC. In the following, we will explain each group.

# Functions that simulate the model


# Functions for ABC SMC





