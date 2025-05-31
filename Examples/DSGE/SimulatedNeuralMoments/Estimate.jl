#= 

This is written to be used interactively, with VScode, to explain the
methods, step by step. For good performance, it is better to wrap 
everything into a function. See the example in the DSGE directory of
https://github.com/mcreel/SNM  for how to do that.

Start julia as "julia --proj -t auto" to use threads

=#


using Pkg
cd(@__DIR__)
Pkg.activate("../../..")
using SimulatedNeuralMoments, Flux, SolveDSGE, MCMCChains
using Distributions, StatsPlots, CSV, PrettyTables, DataFrames
using BSON:@save
using BSON:@load

# get the things to define the structure for the model
include("CKlib.jl") # contains the functions for the DSGE model
n = 160
lb,ub = PriorSupport()
model = SNMmodel("DSGE example", n, lb, ub, GoodData, InSupport, Prior, PriorDraw, auxstat)

# train the net (you probably don't want to do this!)
# nnmodel, nninfo, params, stats, transf_stats = MakeNeuralMoments(model, TrainTestSize=100000, Epochs=1000)
# @save "neuralmodel.bson" nnmodel nninfo # save the trained model

# load the trained model
@load "neuralmodel.bson" nnmodel nninfo # use this to load a trained net
# fill in the structure that defines the model
n = 160
lb, ub = PriorSupport()
model = SNMmodel("DSGE example", n, lb, ub, GoodData, InSupport, Prior, PriorDraw, auxstat)

## see how the NN estimator works with some random parameter draws
S = 100
errs = zeros(S,7)
Threads.@threads for i = 1:S
    # generate some date and define the neural moments using the data
    θtrue = PriorDraw()
    ok = false
    while !ok
        data = dgp(θtrue, dsge, 1, rand(1:Int64(1e10)))[1]
        ok = GoodData(auxstat(data))
        if ok
            θnn =  NeuralMoments(auxstat(data), model, nnmodel, nninfo)[:]
            errs[i,:] = θnn - θtrue
            pretty_table([θtrue θnn], header = ["θtrue", "θnn"])
        else
            println("bad draw, repeating")
        end    
    end    
end

## see how the NN estimator works at the "true parameters" for the DSGE example
# using the common Monte Carlo data sets
S = 1000
errs = zeros(S, size(lb,1))
θtrue = TrueParameters()
Threads.@threads for  i = 1:S
    # generate some date and define the neural moments using the data
    data = Matrix(CSV.read("../GenData/MCdata/mcdata-design-$i.csv", DataFrame))
    θnn =  NeuralMoments(auxstat(data), model, nnmodel, nninfo)[:]
    errs[i,:] = θnn - θtrue
end
b = mean(errs, dims=1)[:]
m = b .+ θtrue
s = std(errs, dims=1)[:]
r = sqrt.(b.^2 + s.^2)[:]
names = ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"] 
printstyled("Monte Carlo SNM neural net results for the DSGE model, 1000 reps\n", color=:green)
pretty_table([names round.([TrueParameters() m b s r],digits=4)],
 header=["parameter", "True", "mean", "bias", "st. dev.", "rmse"])


## Now, let's move on to Bayesian MSM using either the typical data set, or generate a new one
data = CSV.File("dsgedata.csv") |> CSV.Tables.matrix 
#θtrue = TrueParameters() # the parameters of the example data
#data = dgp(θtrue, dsge, 1, rand(1:Int64(1e10)))[1]



## get the NN estimate from the data, using the trained net
θnn = NeuralMoments(auxstat(data), model, nnmodel, nninfo)[:]
pretty_table(round.([TrueParameters() θnn], digits=4), header = (["True", "θnn"]))

## settings for MCMC
S = 50
covreps = 1000
tuninglength = 500
finallength = 1000
burnin = 100
verbosity = 100 # show results every X draws
tuning = 0.1 # pre-checked to work ok

## define the proposal and the log-likelihood
junk, Σp = mΣ(θnn, covreps, model, nnmodel, nninfo)
while !isposdef(Σp)
    for i = 1:size(Σp,1)
        Σp[i,i] += 1e-5
    end
end    
proposal(θ) = rand(MvNormal(θ, tuning*Σp))
lnL = θ -> snmobj(θ, θnn, S, model, nnmodel, nninfo)

## run a short chain to improve proposal
# tuning the chain and creating a good proposal may
# need care - this is just an example!
chain = SimulatedNeuralMoments.mcmc(θnn, tuninglength, lnL, model, nnmodel, nninfo, proposal, burnin, verbosity)
acceptance = mean(chain[:,end])
start = 0.
# update proposal until acceptance rate is good
while acceptance < 0.2 || acceptance > 0.3
    global tuning, chain, acceptance, start
    acceptance < 0.2 ? tuning *= 0.5 : nothing
    acceptance > 0.3 ? tuning *= 1.5 : nothing
    proposal(θ) = rand(MvNormal(θ, tuning*Σp))
    start = mean(chain[:,1:end-2], dims=1)[:]
    chain = SimulatedNeuralMoments.mcmc(start, tuninglength, lnL, model, nnmodel, nninfo, proposal, burnin, verbosity)
    acceptance = mean(chain[:,end])
end

## final long chain
start = mean(chain[:,1:end-2], dims=1)[:]
chain = SimulatedNeuralMoments.mcmc(start, finallength, lnL, model, nnmodel, nninfo, proposal, burnin, verbosity)

## visualize results
chn = Chains(chain[:,1:end-2], ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"])
plot(chn)
#savefig("chain.png") # saved one used 5000 draws, current settings are fewer
pretty_table([TrueParameters() θnn mean(chain[:,1:end-2],dims=1)[:]], header = (["θtrue", "θnn", "θmcmc"]))
display(chn)

