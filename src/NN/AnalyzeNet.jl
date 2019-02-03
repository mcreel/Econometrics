using PyPlot
# this script makes predictions, gets RMSEs, and plots the figure
function AnalyzeNet(savefile, epochs, data, trainsize, noutputs; title="", params="", doplot=false)
    Y = copy(data[trainsize+1:end,1:noutputs])
    data, m, s = stnorm(data)
    X = data[trainsize+1:end,noutputs+1:end]'
    model = mx.load_checkpoint(savefile, epochs, mx.FeedForward) # load trained model
    # obtain predictions
    provider = mx.ArrayDataProvider(:data => X)
    fit = mx.predict(model, provider)
    fit = fit'
    fit = fit.*s[:,1:noutputs]
    fit = fit.+m[:,1:noutputs]
    # compute RMSE
    error = fit - Y
    bias = mean(error,1)
    mse = mean(error.^2.0,1)
    rmse = sqrt.(mse)
    #rsq = 1.0 .- mse./mean((Y .- mean(Y,1)).^2.0,1)
    rsq = 1.0 .- var(error,1) ./ var(Y,1)
    names = ["bias","rmse","mse","R²"]
    println()
    println("__________________________________________________")
    println("Neural net indirect inference results")
    if title != "" println(title) end
    if params != ""
        prettyprint([bias' rmse' mse' rsq'], names, params)
    else
        prettyprint([bias' rmse' mse' rsq'], names)
    end     
    println("__________________________________________________")
    if doplot
        # get the first layer parameters for influence analysis
        model = mx.load_checkpoint(savefile, epochs) # load trained model
        beta = copy(model[2][:fullyconnected0_weight])
        z = maximum(abs.(beta),2);
        cax3 = matshow(z', interpolation="nearest")
        colorbar(cax3)
        ninputs = size(z,1)
        xlabels = [string(i) for i=1:ninputs]
        xticks(0:ninputs-1, xlabels)
        show()
        println("") # stop PyPlot screen spam
    end
    return fit
end

