# M. Creel's modification of sombrero, to illustrate problems in
# optimization. The function is tilted slightly, so that the
# concentric rings each have a unique local maximum.


function FoggyMountainPictures(n)

	if (nargin != 1)
	usage ("sombrero (n)");
	endif

	if (n > 1)
		x = y = linspace (-20, 20, n)';
		[xx, yy] = meshgrid (x, y);
		r = sqrt (xx .^ 2 + yy .^ 2) + eps;
		z = sin (r) ./ r;
		z = z + xx/80 - 0.1*(xx/16).^2;
		z = z/10;

		mesh(x, y, z);
		print("FoggySurface.png", "-dpng"); # let's not overwrite that file!
		figure;
		contour(x, y, z);
		grid("on");
		grid("minor");
		print("FoggyContour.png", "-dpng"); # let's not overwrite that file!

		else
		error ("FoggyMountainPictures: number of grid lines must be greater than 1");
	endif

endfunction
