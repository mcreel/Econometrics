system("cp CGHK_Bayes.mod junk.mod");
system("cp CGHK_Bayes_steadystate.m junk_steadystate.m");
dynare junk.mod
system("./cleanup");
