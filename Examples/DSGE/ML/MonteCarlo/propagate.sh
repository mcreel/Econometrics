#!/bin/bash
for i in {1..1000}
do
	cp TEMPLATE.mod MonteCarlo.mod
    sed -i "s/REPLACE/$i/"g MonteCarlo.mod
    julia --proj --eval 'using Dynare, DelimitedFiles; c = @dynare "MonteCarlo.mod"; writedlm("posmode", transpose(c.results.model_results[1].estimation.posterior_mode))'
    if [ $i -eq 1 ]; then mv posmode posmodes; else cat posmodes posmode >> posmodes; fi
done
rm posmode
rm -R -f MonteCarlo
rm MonteCarlo.log
