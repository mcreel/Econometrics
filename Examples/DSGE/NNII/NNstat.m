# this function take a raw 40 element vector 
# of stats for the DSGE model and outputs
# the fitted 300_40.jl neural net fittted 
# value
# Z enters as a row vector, output stat is also a row vector
function fit = NNstat(Z)
        # load the preprocess items and NN coefs
        load NNitems;
        # get the NN fit
        x = (Z - mZs) ./ sZs;
        preprocess = [ones(size(x,1),1) x]*beta;
        h1 = max(0,alpha1' + x*beta1);
        h2 = max(0,alpha2' + h1*beta2);
        fit = alpha3' + h2*beta3;
        fit = fit.*sErrors + preprocess;
endfunction        
