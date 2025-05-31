using Ipopt, MathProgBase
"""
    qreg(y,x,tau=0.5)

Estimate a quantile regression, quantile tau

"""
function qreg(y,x,q=0.5)
    # thanks to Paul Soderlind for suggesting this.
    # TBD: estimated standard errors.
    (T,K) = (size(x,1),size(x,2))
    c   = vcat(fill(q,T),fill(1-q,T),zeros(K))    #c'[u;v;bet]
    A   = hcat(Matrix(1.0I,T,T),-Matrix(1.0I,T,T),x)  #u[t] - v[t] + x[t,:]*b  = y[t]
    b   = copy(vec(y))
    lb  = vcat(zeros(T),zeros(T),fill(-Inf,K))        #u>=0,v>=0
    ub  = fill(Inf,2*T+K)
    Sol = linprog(c,A,'=',b,lb,ub,IpoptSolver(print_level=1))
    if Sol.status == :Optimal
    par = Sol.sol[2*T+1:end]          #extract regression coeffs
    else
    par = NaN
    end
    return par
end
