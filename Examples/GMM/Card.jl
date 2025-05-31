## this script shows how to do IV estimation, in the GMM form, using
#  the well-known Card data on returns to schooling

##
using Econometrics, CSV, DataFrames, DataFramesMeta, Statistics, Term
cd(@__DIR__)
card = CSV.read("../Julia/cooked.csv", DataFrame)   # this is the prepared data
                                                    # from the script BasicDataAnalysis.jl
display(first(card,6))

## first, OLS
y = Matrix(@select(card, :lnwage))
x = Matrix(@select(card, :educ, :exper, :expsq, :black, :smsa, :south))
x = [ones(size(y,1)) x]
names = ["const", "educ", "exp", "expsq", "black", "smsa", "south"]
βols, Vols, junk = ols(y,x,names=names);

## define instruments
w = Matrix(@select(card, :nearc4, :age, :agesq, :black, :smsa, :south))
w, junk = stnorm(w) # this type of scaling the istruments is good for numeric accuracy
w = [ones(size(w,1)) w]

## define moments
moments = θ -> w.*(y-x*θ)

## CUE
βcue, junk, Vcue, junk =  gmmresults(moments, βols, "", "CUE GMM", names, true);

## two-step
βgmm1, objvalue, junk =  gmm(moments, βols, eye(7))
ms = moments(βgmm1)
W = inv(cov(ms))
βgmm, junk =  gmmresults(moments, βgmm1, W, "two step GMM", names, true);

## Note: this is all just identified, so the following should 
# all be the same. Without the scaling instruments, step1 will 
# be different, which means it didn't actually find the minimizer,
# or that the poor scaling causes lack of identification. How to
# figure out which it is?
prettyprint([βcue βgmm1 βgmm], ["cue", "step1", "step2"])

## Are the instruments ok? Look at first stage regressions
# we are regressing the 3 endogenous variables (educ, exper and exper^2)
# on the instruments, and testing that the outside instruments (nearc4, age, age^2)
# are contributing to the fit, on top of what the included exogenous variables
# account for. This is important: if they really made no contribution to fit, then the
# purged regressors (X hat) would not be linearly independent
R = [zeros(3) eye(3) zeros(3,3)]
r = zeros(3)
y = Matrix{Float64}(@select(card, :educ))
F1 = (TestStatistics(y, w, R, r; silent=true)[1]) # first return of TestStatistics is F
y = Matrix(@select(card, :exper))
F2 = (TestStatistics(y, w, R, r; silent = true)[1])
y = Matrix(@select(card, :expsq))
F3 = (TestStatistics(y, w, R, r; silent=true)[1])
println()
PrintDivider()
println(@green "F test of instrument strength for the 3 endog variables")
prettyprint([F1 F2 F3],["educ", "exper", "expsq"])
println("according to the simple rule that F should be ≥ 10")
println("there is concern that education is not well-instrumented")

## Standard Hausman test: compare CUE and OLS
# this is assuming that OLS is efficient, which
# is a bit of a stretch (possible HET, at least)
using LinearAlgebra, Distributions
e = βcue-βols
H = round(dot(e, inv(Vcue-Vols), e),digits=4)
df = rank(Vcue-Vols)
pval = round(1. - cdf(Chisq(df),H), digits=4)
PrintDivider()
println(@green "Hausman test")
println("statistic: $H, degrees of freedom: $df, p-value: $pval")

##
# to summarize, given the weakness of the nearc4 instrument,
# and the Hausman test, which gives some support to OLS,
# there's not too much reason to favor the IV results over
# the OLS results