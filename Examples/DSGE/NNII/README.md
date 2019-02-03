# DSGE_Neural_Net_Estimation_Example

NN_Estimate.m
This code shows how the parameters of a small nonlinear DSGE model may be estimated with good accuracy, very, very quickly (e.g., in around one second).

The parameters \theta are estimated using the limited information posterior mean E(\theta|Z), where Z is a vector of statistics. 

The posterior mean function has previously been fit using a neural net, using the code at https://github.com/mcreel/NeuralNetsForIndirectInference.jl  Implementing this fit only requires prior information, so this part can be done before the actual sample is observed.

The model is solved using Dynare, with a 3rd order solution. The model is the same one described in  "Bayesian Indirect Inference and the ABC of GMM" by Creel, Gao, Hong and Kristensen, http://arxiv.org/abs/1512.07385 The Dynare .mod file is provided, see it for the details.

The present code requires Octave or Matlab, plus Dynare. I believe that it is completely self-contained, apart from these dependencies. The code has been checked using Octave. If there are problems running it with Matlab, please let me know.
