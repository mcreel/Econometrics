using LinearAlgebra
function vech(x)
    k = size(x,1)
    a = zeros(Int((k^2-k)/2 + k))
    m = 1
    for i = 1:k
        for j = 1:i
            a[m] = x[i,j]
            m += 1
        end
    end
    a
end

# for testing vech()
function vech()
    a = eye(3)
    b = vech(a)
    c = [1.0; 0.0; 1.0; 0.0; 0.0; 1.0]
    b == c
end    
