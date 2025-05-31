using Dynare, Serialization, Statistics
cd(@__DIR__)
context = @dynare "CK.mod" ; 

chain = deserialize("CK/output/mcmc_chain_1.jls")
run(`./cleanup`)
