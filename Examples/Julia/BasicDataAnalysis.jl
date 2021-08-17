#=
This shows basic data analysis using the Card data
on returns to education

See the github pages of the packages for documentation

=#

using CSV, DataFrames, DataFramesMeta, DelimitedFiles, StatsPlots, Econometrics
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)#pyplot(reuse=false, show=true)
dir = dirname(Base.load_path()[1])
card = CSV.read(dir*"/Examples/Data/card.csv", DataFrame) # how to read CSV data into DataFrame
size(card)
names(card)
# see the first 6 lines
show(first(card,6))
# see some statistics using DataFrames functionality
show(describe(card))
# see some statistics using Econometrics functionality
println()
dstats(Matrix{Float64}(card))

# a simple way of looking at geographical and racial discrimination
# select north and white
nw = @subset(card, :south .== 0, :black .== 0)
@df nw scatter(:educ, :wage,label="North-White")
# select south and black
sb = @subset(card, :south .== 1, :black .== 1)
@df sb scatter!(:educ, :wage, label="South-Black",title="Wage vs. Educ",xlabel="Educ", ylabel="Wage")
gui()
wait_for_key("press enter to continue")
# the following is similar, but maybe better, as it plots all points
@df card scatter(:educ, :wage, zcolor=:south+:black)
gui()
wait_for_key("press enter to continue")

##
# density plots by groups
@df card density(:wage, group=(:black), fill=true, fillalpha=0.3,legend=:topright)
gui()
wait_for_key("press enter to continue")
##
# violin and box plots
@df card boxplot(:wage) # puts is beside
@df card violin!(:wage)
gui()
wait_for_key("press enter to continue")
# violin plot of wage, by education
@df card violin(:educ, :wage)
# add a boxplot on top
@df card boxplot!(:educ, :wage)
gui()
wait_for_key("press enter to continue")
##
# convert the dataframe to an ordinary array,
# as this is usually better for more advanced
# modeling
data = Matrix{Float64}(card)
writedlm("card.ascii", data) # write to plain ASCII file
