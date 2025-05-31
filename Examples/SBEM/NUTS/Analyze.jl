include("SV.jl")
posterior, chain = main()
# Extract the parameter and plot the posterior
σe_hat = [i[1][1] for i in posterior]
#σu_hat = 2.0*log.(σu_hat)
dstats(σe_hat)
p1 = npdensity(σe_hat)
# Extract the parameter and plot the posterior
ρhat = [i[2][1] for i in posterior]
dstats(ρhat)
p2 = npdensity(ρhat)
# Extract the parameter and plot the posterior
σu_hat = [i[3][1] for i in posterior]
dstats(σu_hat)
p3 = npdensity(σu_hat)
plot(p1, p2, p3)
plot!(legend=false)
gui()
#savefig("nuts.svg")
