using Plots
plotlyjs()
include("FoggyMountainObj.jl")

n = 200
x = range(-20., stop=20., length=n)
y = x

##
p1 = surface(x,y,(x,y)->-FoggyMountainObj(x,y))
savefig("FoggySurface.png")
gui()

##
p2 = contour(x,y,(x,y)->-FoggyMountainObj(x,y))
gui()
savefig("FoggyContour.png")
