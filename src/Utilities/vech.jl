using LinearAlgebra
function vech(x)
    k = size(x,1)
    a = zeros((k^2)/2 + k)
    m = 1
    for i = 1:k
        for j = 1:i
            a[m] = x[i,j]
            m += 1
        end
    end
end    
