# run this to make a precompiled system image with the packages used
# in Econometrics pre-compiled, and with some commonly used functions
# also pre-compiled.
# 
# once it's done, start Julia from the Econometrics directory with
# julia --proj -JEconometricsSysimage.so"
#
# just start without the switch to return to normal. You will
# need to re-do this if you update your packages.

using Pkg, PackageCompiler
Pkg.activate(".")
PackageCompiler.create_sysimage(:Econometrics; sysimage_path="JuliaSysimage.so",
                                                precompile_execution_file="warmup.jl")

