using Econometrics, CSV, DataFrames, DataFramesMeta, Statistics
card = CSV.read("../Julia/cooked.csv", DataFrame)   # this is the prepared data
                                                    # from the script BasicDataAnalysis.jl
display(first(card,6))
# first, OLS
y = Matrix(@select(card, :lnwage))
x = Matrix(@select(card, :educ, :exper, :expsq, :black, :smsa, :south))
x = [ones(size(y,1)) x]
names = ["const", "educ", "exp", "expsq", "black", "smsa", "south"]
βols, junk = ols(y,x,names=names);

# define instruments
w = Matrix(@select(card, :nearc4, :age, :agesq, :black, :smsa, :south))
w = [ones(size(y,1)) w]
# define moments
moments = θ -> w.*(y-x*θ)
# CUE
βcue, junk =  gmmresults(moments, βols, "", "CUE GMM", names, true)
# two-step
βgmm1, junk =  gmm(moments, βols, eye(7))
ms = moments(βgmm1)
W = inv(cov(ms))
βgmm, junk =  gmmresults(moments, βols, W, "two step GMM", names, true)


