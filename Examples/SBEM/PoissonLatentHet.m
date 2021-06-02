# this calculates the log likelihood of a
# Poisson model where lambda contains a latent variable
# in this example, the latent variable is normally distributed
function [logdensity, score, prediction] = PoissonLatentHet(theta, data, otherargs)
	k = otherargs{1};

	y = data(:,1);
	x = data(:,2:k+1);
	randdraws = data(:,k+2:columns(data));
	S = columns(randdraws);
	n = rows(x);

	beta = theta(1:k,:);
	sig = exp(theta(k+1,:));
	lambda = exp(repmat(x*beta,1,S) + sig*randdraws);
	Y = repmat(y,1,S);
	density = exp(-lambda) .* (lambda .^ Y) ./ factorial(Y); # compare to Poisson
	density = mean(density, 2);
	logdensity = slog(density);
	score = "na"; # we don't have the scores
	prediction = "na";
endfunction
