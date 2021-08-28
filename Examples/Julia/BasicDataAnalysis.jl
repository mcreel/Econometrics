#=
This shows basic data analysis using the Card data
on returns to education

If running from the REPL, it will pause after each
plot, and wait for an <enter>

If using VSCode, just run the blocks.
=#

using CSV, DataFrames, DataFramesMeta, DelimitedFiles, StatsPlots
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)#pyplot(reuse=false, show=true)
card = CSV.read("../Data/card.csv", DataFrame) # how to read CSV data into DataFrame
size(card)
names(card)
# see the first 6 lines
show(first(card,6))
# see some statistics using DataFrames functionality
show(describe(card))

# a simple way of looking at geographical and racial discrimination
# select north and white
nw = @subset(card, :south .== 0, :black .== 0)
@df nw scatter(:educ, :wage,label="North-White")
# select south and black
sb = @subset(card, :south .== 1, :black .== 1)
@df sb scatter!(:educ, :wage, label="South-Black",title="Wage vs. Educ",xlabel="Educ", ylabel="Wage")

# this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

# the following is similar, but maybe better, as it plots all points
@df card scatter(:educ, :wage, zcolor=:south+:black)

# this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

# density plots by groups
@df card density(:wage, group=(:black,:south), fill=true, fillalpha=0.3,legend=:topright, title="wage density by south and\n black categorical variables")

@df card marginalkde(:wage, :educ)
# this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

# violin and box plots
@df card boxplot(:wage) # puts is beside
@df card violin!(:wage)

# this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

# violin plot of wage, by education
@df card violin(:educ, :wage)
# add a boxplot on top
@df card boxplot!(:educ, :wage)

# this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

# convert the dataframe to an ordinary array,
# as this is usually better for more advanced
# modeling
data = Matrix{Float64}(card)
# writedlm("card.ascii", data) # write to plain ASCII file
