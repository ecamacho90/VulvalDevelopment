function kernelvalue = Vulval_Development_Modelv10_v1_EvalKernels(newparticle,particleprevstep,CovMat)

kernelvalue = mvnpdf(newparticle,particleprevstep,CovMat);