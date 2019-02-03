using LinearAlgebra
function vech(x)
    a = triu(x)
    a = a[findall(x->x!=0,a)]
end    
