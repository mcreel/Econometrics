using LinearAlgebra: I

function eye(m::Int64)
    Matrix(1.0I, m, m)
end

function eye()
    a = eye(3)
    sum(sum(a)) == 3
end    
