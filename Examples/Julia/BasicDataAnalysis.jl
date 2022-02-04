#=
This shows basic data analysis using the Card data
on returns to education

If running from the REPL, it will pause after each
plot, and wait for an <enter>

If using VSCode, just run the blocks.
=#

##
using CSV, DataFrames, DataFramesMeta, DelimitedFiles, StatsPlots
cd(@__DIR__)
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)#pyplot(reuse=false, show=true)
card = CSV.read("../Data/card.csv", DataFrame) # how to read CSV data into DataFrame
size(card)
names(card)
## see the first 6 lines
show(first(card,6))
# see some statistics using DataFrames functionality
show(describe(card))

## add some variables
card[!,:lnwage] = round.(log.(card[!,:wage]), digits=4) # round to save storage space
card[!,:expsq] = (card[!,:exper].^2)/100
card[!,:agesq] = card[!,:age].^2

## select the ones we want
cooked = @select(card, :lnwage, :educ, :exper, :expsq, :black, :south, :smsa, :age, :agesq, :nearc4 )
display(first(cooked,6))

## write as CSV
CSV.write("cooked.csv", cooked)

## write as plain ASCII, using limited precision to reduce file size
data = Matrix{Float32}(cooked)
writedlm("cooked.txt", data)

## a simple way of looking at geographical and racial discrimination
# select north and white
nw = @subset(card, :south .== 0, :black .== 0)
@df nw scatter(:educ, :wage,label="North-White")
# select south and black
sb = @subset(card, :south .== 1, :black .== 1)
@df sb scatter!(:educ, :wage, label="South-Black",title="Wage vs. Educ",xlabel="Educ", ylabel="Wage")

## this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

## the following is similar, but maybe better, as it plots all points
@df card scatter(:educ, :wage, zcolor=:south+:black)

## this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

## density plots by groups
@df card density(:wage, group=(:black,:south), fill=true, fillalpha=0.3,legend=:topright, title="wage density by south and\n black categorical variables")

@df card marginalkde(:wage, :educ)
## this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

## violin and box plots
@df card boxplot(:wage) # puts is beside
@df card violin!(:wage)
gui()
wait_for_key("press enter to continue")

## violin plot of wage, by education
@df card violin(:educ, :wage)
# add a boxplot on top
@df card boxplot!(:educ, :wage)

## this is here in case running from REPL
gui()
wait_for_key("press enter to continue")

