# Basic fixed effect panel data
# using the Arellano-Bond data (exported from GRETL)

## packages
using CSV, DataFrames, ShiftedArrays, Term, StatsPlots, FixedEffectModels

## read in the data
cd(@__DIR__)
ab = CSV.read("../Examples/Data/abdata.csv", DataFrame; missingstring="NA")
show(describe(ab))

## drop missings (use different name, to be able to go back)
ab_nm = dropmissing(ab)

## examine the data
show(describe(ab_nm))
@df ab_nm boxplot(:n, label="n")
@df ab_nm boxplot!(:w, label="w")
@df ab_nm boxplot!(:k, label="k")
@df ab_nm boxplot!(:ys, label="ys")
plot!(legend=:bottomright)

## simple FE model: time and firm fixed effects
println()
println(@green "simple FE model with year and firm effects")
display(reg(ab_nm, @formula(n ~ w + k + ys + fe(YEAR) + fe(unit)), Vcov.cluster(:unit)))

## now, prepare DPD estimation
# first, get first differences of the variables
transform!(groupby(ab, :unit), :n => (x -> x .- lag(x)) => :Δn)
transform!(groupby(ab, :unit), :w => (x -> x .- lag(x)) => :Δw)
transform!(groupby(ab, :unit), :k => (x -> x .- lag(x)) => :Δk)
transform!(groupby(ab, :unit), :ys => (x -> x .- lag(x)) => :Δys)
# lagged dep var, after differencing
transform!(groupby(ab, :unit), :Δn => (x -> lag(x)) => :Δnlag)
# instruments: lags of levels of dep var
transform!(groupby(ab, :unit), :n => (x -> lag(x,2)) => :nlag2)
transform!(groupby(ab, :unit), :n => (x -> lag(x,3)) => :nlag3)
transform!(groupby(ab, :unit), :n => (x -> lag(x,4)) => :nlag4)

## examine the new variables
show(describe(ab))

## drop missings, to run the regression on the available observations
dropmissing!(ab)

## DPD: time and firm fixed effects
println()
println(@green "DPD with year and firm effects")
# add more instruments if desired, e.g., (nlag2+nlag3) 
# or (nlag2+nlag3+nlag4)
#=
Note: these results don't coincide with GRETL, because GRETL replaces missings with 0.
This code drops observations with missings. Thus, the number of observations is 
considerably different. Also, the covariance estimator is not the same, as the
specialized DPD covariance estimator is not used, so reported standard deviations
are different.
=#
display(reg(ab, @formula(Δn ~ (Δnlag ~ (nlag2)) + Δw + Δk + Δys + fe(YEAR) + fe(unit)), Vcov.cluster(:unit)))


