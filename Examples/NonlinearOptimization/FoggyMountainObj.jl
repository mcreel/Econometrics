# this is the funtion that appears in the
# FoggySurface and FoggyContour graphics,
# It is the Octave sombrero, tilted a bit
function FoggyMountainObj(theta)
	theta1 = theta[1]
	theta2 = theta[2]
   	r = sqrt(theta1^2.0 + theta2^2.0) + 1e-10
	z = sin(r) / r
	z = z + theta1/80.0 - 0.1*(theta1/16.0)^2.0
	z = -z/10.0 # switch to minimization
end
