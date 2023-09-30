## computes returns for several data series, to estimate a CAPM model
using Glob, CSV, DataFrames, DataFramesMeta, StatsPlots, Dates, DelimitedFiles
cd(@__DIR__)

## read the data into a DataFrame
files = glob("*.csv", ".")
CAPM = DataFrame()
global i = 0
for f in files
    global i +=1
    temp = CSV.read(f, DataFrame)
    rename!(temp, :var"Adj Close" => :close)
    ## compute percentage returns
    name = f[3:end-4]
    rets = vcat(0.0,[100.0 *(log(temp.close[i]) - log(temp.close[i-1])) for i=2:size(temp,1)])
    i == 1 ? CAPM[!, "Date"] = temp.Date[2:end] : nothing
    CAPM[!,name*"rets"] = rets[2:end]
    CAPM[!,name*"close"] = temp.close[2:end]
end
CSV.write("CAPMdata.csv", CAPM) # write as CSV

