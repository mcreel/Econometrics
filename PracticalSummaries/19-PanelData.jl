# Basic fixed effect panel data, using the Arellano-Bond
# data (exported from GRETL)


## read in the data and drop missings
using CSV, DataFrames, ShiftedArrays, Term
cd(@__DIR__)
ab = CSV.read("../Examples/Data/abdata.csv", DataFrame; missingstring="NA")
show(describe(ab))

## drop missings (use different name, to be able to go back)
ab_nm = dropmissing(ab)

## examine the data
show(describe(ab_nm))
using StatsPlots
@df ab_nm boxplot(:n, label="n")
@df ab_nm boxplot!(:w, label="w")
@df ab_nm boxplot!(:k, label="k")
@df ab_nm boxplot!(:ys, label="ys")
plot!(legend=:bottomright)

## run a model: time and firm fixed effects
println(@green "simple FE model with year and firm effects")
using FixedEffectModels
display(reg(ab_nm, @formula(n ~ w + k + ys + fe(YEAR) + fe(unit))))

## now, prepare DPD estimation
# first,  add first differences of the variables
println(@green "adding first differences")
transform!(groupby(ab, :unit), :n => (x -> x .- lag(x)) => :Δn)
transform!(groupby(ab, :unit), :w => (x -> x .- lag(x)) => :Δw)
transform!(groupby(ab, :unit), :k => (x -> x .- lag(x)) => :Δk)
transform!(groupby(ab, :unit), :ys => (x -> x .- lag(x)) => :Δys)
show(describe(ab))
printstyled(@green "note: there are still missings")


