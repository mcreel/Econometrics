#=
This shows basic data analysis using the Card data
on returns to education

See the github pages of the packages for documentation

=#

using CSV, DataFrames, DataFramesMeta, DelimitedFiles, StatsPlots
pyplot(reuse=false, show=true)
card = CSV.read("../Data/card.csv") # how to read CSV data into DataFrame
size(card)
names(card)
show(first(card,6))
show(describe(card))

# select north and white
nw = @where(card, :south .== 0, :black .== 0)
@df nw scatter(:educ, :wage,label="North-White")
# select south and black
sb = @where(card, :south .== 1, :black .== 1)
@df sb scatter!(:educ, :wage, label="South-Black",title="Wage vs. Educ",xlabel="Educ", ylabel="Wage")
##
# density plots by groups
@df card density(:wage, group=(:black), fill=true, fillalpha=0.3,legend=:topright)

##
# violin and box plots
@df card boxplot(:wage) # puts is beside
@df card violin!(:wage)
# violin plot of wage, by education
@df card violin(:educ, :wage)
# add a boxplot on top
@df card boxplot!(:educ, :wage)
##
# convert the dataframe to an ordinary array,
# as this is usually better for more advanced
# modeling
data = convert(Array{Float64}, card)
writedlm("card.ascii", data) # write to plain ASCII file
