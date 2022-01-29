# simple iid bootstrap
function bootstrap(data)
    n = size(data,1)
    data[rand(1:n,n),:]
end    


