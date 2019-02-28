# shows some data analysis using DataFrames
using CSV, DataFrames, GLM
dir = dirname(dirname(pathof(Econometrics)))
nerlove = CSV.read(dir*"/Examples/Data/nerlove.csv")
nerlove[:lnC] = log.(nerlove[:cost])
nerlove[:lnQ] = log.(nerlove[:output])
nerlove[:lnPL] = log.(nerlove[:labor])
nerlove[:lnPF] = log.(nerlove[:fuel])
nerlove[:lnPK] = log.(nerlove[:capital])
f = @formula(lnC ~ 1 + lnQ + lnPL + lnPF + lnPK)
lm(f, nerlove)
