system("cp CGHK_ml.mod junk.mod");
system("cp CGHK_ml_steadystate.m junk_steadystate.m");
dynare junk.mod
system("./cleanup");
