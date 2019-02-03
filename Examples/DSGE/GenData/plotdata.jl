using Plots
data = readdlm("dsgedata.txt")
plot(data, label = ["y", "c", "n", "r", "w"], show=true)
savefig("dsgedata.svg")
