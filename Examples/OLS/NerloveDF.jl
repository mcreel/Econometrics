# shows some data analysis using DataFrames
using CSV, DataFrames, GLM
nerlove = CSV.read("nerlove.csv")
nerlove[:lnC] = log.(nerlove[:cost])
nerlove[:lnQ] = log.(nerlove[:output])
nerlove[:lnPL] = log.(nerlove[:labor])
nerlove[:lnPF] = log.(nerlove[:fuel])
nerlove[:lnPK] = log.(nerlove[:capital])
f = @formula(lnC ~ 1 + lnQ + lnPL + lnPF + lnPK)
lm(f, nerlove)
