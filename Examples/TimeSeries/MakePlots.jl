using Econometrics, Statistics, Distributions, StatsPlots, DelimitedFiles
function main()
# weekly close price of NSYE, data provided with GRETL
data = readdlm("nysewk.txt")
# compute weekly percentage growth
y = 100.0 * log.(data[2:end] ./ data[1:end-1])
p1 = plot(y, leg=false, show=true)
p2 = histogram(y, normed=true)
plot!(p2, Distributions.Normal(mean(y),std(y)), show=true, leg=false)
p3 = npdensity(y)
plot!(p3,legend=:left)
plot(p1,p2,p3, layout = (3,1))
gui()
savefig("nyse.png")
return
end
main()
