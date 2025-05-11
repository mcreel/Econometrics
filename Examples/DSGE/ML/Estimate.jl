using Dynare

cd(@__DIR__)

context = @dynare "ML.mod" ;
run(`./cleanup`);
