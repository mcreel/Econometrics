# this is the funtion that appears in the
# FoggySurface and FoggyContour graphics,
# It is the Octave sombrero, tilted a bit
function FoggyMountainObj(θ₁, θ₂)
   	r = sqrt(θ₁^2.0 + θ₂^2.0) + 1e-10
	z = sin(r) / (2*r)
	z = z + θ₁/40.0 - 2*(θ₁/40.0)^2.0 + θ₂/40.0 - 2*(θ₂/40.0)^2.0 
	z = -z/10.0 # switch to minimization
end

function FoggyMountainObj(θ)
	θ₁ = θ[1]
	θ₂ = θ[2]
    FoggyMountainObj(θ₁, θ₂)
end
