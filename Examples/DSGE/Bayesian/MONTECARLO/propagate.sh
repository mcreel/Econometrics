#!/bin/bash
for i in {1..1000}
do
	cp TEMPLATE.mod MonteCarlo.mod
    sed -i "s/REPLACE/$i/"g MonteCarlo.mod
    julia --proj --eval 'using Dynare, DelimitedFiles, Serialization, Statistics; c = @dynare "MonteCarlo.mod"; writedlm("posmean", mean((deserialize("./MonteCarlo/output/mcmc_chain_1.jls").value)[:,1:end-1],dims=1))'
    if [ $i -eq 1 ]; then mv posmean posmeans; else cat posmeans posmean >> posmeans; fi
done
rm posmean
rm -R -f MonteCarlo
rm MonteCarlo.log
