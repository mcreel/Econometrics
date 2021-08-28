 using Optim
 junk=2.
 # shows use of obj. fun. as a closure
 function sse(x) # very simple quadratic objective function
 objvalue = junk + sum(x.*x)
 end
 k = 5
 x = rand(k,1)
 lb = -ones(k,1)
 ub = -lb
 xopt= (Optim.optimize(sse, lb, ub, x, SAMIN(verbosity=2),Optim.Options(iterations=10^6))).minimizer

