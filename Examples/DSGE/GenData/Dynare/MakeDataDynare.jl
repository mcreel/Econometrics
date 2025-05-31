using Dynare, CSV, DataFrames


@views function main()
c = @dynare "GenData.mod";
d = c.results.model_results[1].simulations[1].data
data = [d.y d.c d.n d.r d.w]

burnin = 500
n = 160
for rep = 1:1000
    start = (rep-1)*(burnin+n)+burnin+1
    stop = rep*(burnin+n)
    temp = data[start:stop,:]
    df = DataFrame(temp, ["y", "c", "n","r","w"])
    CSV.write("./MCdata/mcdata-design-$rep.csv", df)
end

end
main()



