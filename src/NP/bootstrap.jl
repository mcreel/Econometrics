# simple iid bootstrap
function bootstrap(data)
    n = size(data,1)
    resampled = similar(data)
    for i = 1:n
        j = rand(1:n)
        resampled[i,:] = data[j,:]
    end
    return resampled
end    


